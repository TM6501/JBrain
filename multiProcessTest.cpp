
#ifdef __linux__
// Only used by Linux for multiprocessing:
// #include <unistd.h>
#include <sys/wait.h>

pid_t startSingleChild(std::vector<int>& workVector)
{
    pid_t id = fork();
    if (id > 0)
    {
        workVector.erase(
            workVector.begin(), workVector.begin() +
            std::min(1, static_cast<int>(workVector.size())));
    }
    return id;
}

void doWork(const std::vector<int>& workVector, unsigned int workSize)
{
    std::this_thread::sleep_for(std::chrono::seconds(3));
    for (unsigned int i = 0; i < workSize && i < workVector.size(); ++i)
        std::cout << workVector[i] << std::endl;
    std::cout << "Child process ending." << std::endl;
    exit(0);
}

void multiProcessTest(unsigned int maxConcurrentChildren)
{
    unsigned int workSize = 1;
    unsigned int currentConcurrentChildren = 0;
    pid_t id = 1;

    std::vector<int> valsToOutput{ 1, 2, 3, 4, 5, 6, 7, 8, 9,
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 };

    // Assuming valsToOutput size > maxConcurrentChildren, get X children
    // started:
    while (currentConcurrentChildren < maxConcurrentChildren && id > 0)
    {
        id = startSingleChild(valsToOutput);
        if (id == 0)
        {
            doWork(valsToOutput, workSize);
        }
        currentConcurrentChildren += 1;
    }

    // Now we wait on each child, starting a new one when each finishes:
    while (valsToOutput.size() > 0)
    {
        wait(NULL);
        id = startSingleChild(valsToOutput);
        if (id == 0)
            doWork(valsToOutput, workSize);
    }

    std::cout << "Parent done creating children. Waiting on all children." << std::endl;
    while (currentConcurrentChildren > 0)
    {
        wait(NULL);
        currentConcurrentChildren -= 1;
    }
    std::cout << "Parent exiting." << std::endl;
}

int main()
{
    multiProcessTest(4);
    return 0;
}

#elif _WIN32
#include <windows.h>

bool startSingleChild(std::vector<STARTUPINFOA>& startInfos,
    std::vector<PROCESS_INFORMATION>& processInfos,
    std::vector<int>& workToDo)
{
    STARTUPINFOA tempSi;
    PROCESS_INFORMATION tempPi;

    ZeroMemory(&tempSi, sizeof(tempSi));
    tempSi.cb = sizeof(tempSi);
    ZeroMemory(&tempPi, sizeof(tempPi));
    startInfos.push_back(tempSi);
    processInfos.push_back(tempPi);
    char exeBuffer[500] = "Starter.exe doWork ";
    char numBuffer[5];

    // Get the bit of work to be done:
    _itoa_s(workToDo[0], numBuffer, 5, 10);
    strcat_s(exeBuffer, 500, numBuffer);
    workToDo.erase(workToDo.begin(), workToDo.begin() + 1);

    if (CreateProcessA(
        NULL, // No module name
        exeBuffer, // Command line to execute
        NULL, // Process handle not inheritable
        NULL, // Thread handle not inheritable
        FALSE, // Handle inheritance is false
        0, // No creation flags
        NULL, // Use parent's environment block
        NULL, // Use parent's starting directory
        &startInfos[startInfos.size() - 1],  // pointer to startup info
        &processInfos[processInfos.size() - 1] // pointer to process_information
    ))
    {
        return true;
    }
    else
    {
        startInfos.pop_back();
        processInfos.pop_back();
        return false;
    }
}

bool checkRunningAndRemove(std::vector<STARTUPINFOA>& startInfos,
    std::vector<PROCESS_INFORMATION>& processInfos,
    unsigned int processToCheck)
{
    DWORD exitCode;
    GetExitCodeProcess(processInfos[processToCheck].hProcess, &exitCode);

    // Process is still running:
    if (exitCode == STILL_ACTIVE)
    {
        return true;
    }
    // Not running. Clean it up and remove it:
    else
    {
        CloseHandle(processInfos[processToCheck].hProcess);
        CloseHandle(processInfos[processToCheck].hThread);
        startInfos.erase(startInfos.begin() + processToCheck);
        processInfos.erase(processInfos.begin() + processToCheck);
        return false;
    }
}

void multiProcessTest(unsigned int maxConcurrentChildren)
{
    unsigned int currentConcurrentChildren = 0;
    unsigned int loopWaitMilliseconds = 500;

    std::vector<int> valsToOutput{ 1, 2, 3, 4, 5, 6, 7, 8, 9,
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 };

    std::vector<STARTUPINFOA> startInfos;
    std::vector<PROCESS_INFORMATION> processInfos;


    // Start as many as we need in the beginning:
    for (currentConcurrentChildren = 0;
        currentConcurrentChildren < maxConcurrentChildren && valsToOutput.size() > 0;
        ++currentConcurrentChildren)
    {
        if (!startSingleChild(startInfos, processInfos, valsToOutput))
        {
            std::cout << "Failed to start child process. Exiting." << std::endl;
            exit(-1);
        }
    }

    // Keep creating new worker processes until the work is done:
    while (valsToOutput.size() > 0)
    {
        // Sleep to not spin too quickly on this loop:
        std::this_thread::sleep_for(std::chrono::milliseconds(loopWaitMilliseconds));

        // Start children if we need them:
        while (currentConcurrentChildren < maxConcurrentChildren && valsToOutput.size() > 0)
        {
            startSingleChild(startInfos, processInfos, valsToOutput);
            ++currentConcurrentChildren;
        }

        // Check if each process is running:
        for (unsigned int i = 0; i < processInfos.size(); ++i)
        {
            // Process is done? Remove it and decrement the number of
            // children currently running:
            if (!checkRunningAndRemove(startInfos, processInfos, i))
                --currentConcurrentChildren;
        }
    }

    std::cout << "All children created. Waiting for them to finish." << std::endl;
    while (currentConcurrentChildren > 0)
    {
        for (unsigned int i = 0; i < processInfos.size(); ++i)
        {
            // Process is done? Remove it and decrement the number of
            // children currently running:
            if (!checkRunningAndRemove(startInfos, processInfos, i))
                --currentConcurrentChildren;
        }
    }
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cout << "Usage: <AppName> <'main'|'doWork'> <numberOfConcurrentProcesses|workToDo>" << std::endl;
    }
    // Main runs all the other workers:
    if (strcmp(argv[1], "main") == 0)
    {
        std::cout << "Beginning work with " << argv[2] << " child processes." << std::endl;
        multiProcessTest(static_cast<unsigned int>(atoi(argv[2])));
        std::cout << "Parent exiting." << std::endl;
    }
    else if (strcmp(argv[1], "doWork") == 0)
    {
        std::this_thread::sleep_for(std::chrono::seconds(5));
        std::cout << "Did work: " << argv[2] << "..." << " Child exiting." << std::endl;
    }
    else
    {
        std::cout << "Unrecognized command: " << argv[1] << std::endl;
    }

    return 0;
}

#endif