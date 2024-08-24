# JBrain
Simulated brains learning in reinforcement learning environments

This code compiles and runs on my Windows 10 machine. I have an (outdated as of 2024-08-24) Makefile that compiled the project for Linux, but in the interest of getting things to work in as short of a time as possible, I haven't updated it in a long time.

Also, making C++ call into Python libraries is a massive pain and the exact location and version of files is critical.  Once this project reaches as stable version, I will put work in to making it easier to compile and run on other machines.  For now, though, the Visual Studio solution and projects point to specific files in specific locations on my machine.  Making this project compile for you will require changes to the projects and solutions.

Finally, this project relies on another, "ExperimentRunner".  Git's submodule capability doesn't seem to fit this situation.  The "ExperimentRunner" folder is expected to exist on the compiling machine.  As mentioned above, the solution and project files make reference to exactly where it exists on my machine.
