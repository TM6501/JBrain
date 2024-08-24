# Rules
CC=gcc
CXX=g++
DEBUGFLAGS=-g -Wall
RELEASEFLAGS=-O2
CPPFLAGS=-std=c++2a -pthread -I /data/CPP_Libraries/yaml-cpp-yaml-cpp-0.7.0/include -I ./CGP_CPP -I ./CGP_CPP_TestApp -I ./CGym -I ./ExperimentRunner -I /data/CPP_Libraries -I /usr/include/python3.10
CPPLIBFLAGS=-L./libraries -L/usr/lib/python3.10/config-3.10-aarch64-linux-gnu -L/usr/lib/aarch64-linux-gnu -L/data/CPP_Libraries/yaml-cpp-yaml-cpp-0.7.0/build -lyaml-cpp -lExperimentRunner -lCGym -lCGP_CPP -lpython3.10 -lcrypt -lpthread -ldl  -lutil -lm -lm
CPPMONOLIBFLAGS=-L/usr/lib/python3.10/config-3.10-aarch64-linux-gnu -L/usr/lib/aarch64-linux-gnu -L/data/CPP_Libraries/yaml-cpp-yaml-cpp-0.7.0/build -lyaml-cpp -lpython3.10 -lcrypt -lpthread -ldl  -lutil -lm -lm
ARFLAGS=rvs

# We make the debug version with the release version every time because
# I don't know enough about makefiles to do it differently:
.cpp.o:
	$(CXX) $(CPPFLAGS) -c $*.cpp -o $*.o
	$(CXX) $(CPPFLAGS) $(DEBUGFLAGS) -c $*.cpp -o $*_d.o

TestApp_d: CGym/AbstractGymEnv.o CGym/DungeonRoomEnv.o CGym/pch.o CGym/SimpleRepeat.o CGP_CPP/AbstractCGPIndividual.o CGP_CPP/AbstractTesterClass.o CGP_CPP/BasicGene.o CGP_CPP/CGPFFANNIndividual.o CGP_CPP/CGPFunctions.o CGP_CPP/CGPGenerator.o CGP_CPP/GeneralCGPSolver.o CGP_CPP/GymTester.o CGP_CPP/GymWorldGenerator.o CGP_CPP/pch.o CGP_CPP/Enums.o ExperimentRunner/Experiment.o ExperimentRunner/MultiExperiment.o ExperimentRunner/pch.o CGP_CPP_TestApp/CGP_CPP_TestApp.o
	$(CXX) $(CPPFLAGS) $(DEBUGFLAGS) -o CGP_CPP_TestApp/testApp_d.out CGym/AbstractGymEnv_d.o CGym/DungeonRoomEnv_d.o CGym/pch_d.o CGym/SimpleRepeat_d.o CGP_CPP/AbstractCGPIndividual_d.o CGP_CPP/AbstractTesterClass_d.o CGP_CPP/BasicGene_d.o CGP_CPP/CGPFFANNIndividual_d.o CGP_CPP/CGPFunctions_d.o CGP_CPP/CGPGenerator_d.o CGP_CPP/GeneralCGPSolver_d.o CGP_CPP/GymTester_d.o CGP_CPP/GymWorldGenerator_d.o CGP_CPP/pch_d.o CGP_CPP/Enums_d.o ExperimentRunner/Experiment_d.o ExperimentRunner/MultiExperiment_d.o ExperimentRunner/pch_d.o CGP_CPP_TestApp/CGP_CPP_TestApp_d.o $(CPPMONOLIBFLAGS)

TestApp_r: CGym/AbstractGymEnv.o CGym/DungeonRoomEnv.o CGym/pch.o CGym/SimpleRepeat.o CGP_CPP/AbstractCGPIndividual.o CGP_CPP/AbstractTesterClass.o CGP_CPP/BasicGene.o CGP_CPP/CGPFFANNIndividual.o CGP_CPP/CGPFunctions.o CGP_CPP/CGPGenerator.o CGP_CPP/GeneralCGPSolver.o CGP_CPP/GymTester.o CGP_CPP/GymWorldGenerator.o CGP_CPP/pch.o CGP_CPP/Enums.o ExperimentRunner/Experiment.o ExperimentRunner/MultiExperiment.o ExperimentRunner/pch.o CGP_CPP_TestApp/CGP_CPP_TestApp.o
	$(CXX) $(CPPFLAGS) $(RELEASEFLAGS) -o CGP_CPP_TestApp/testApp_r.out CGym/AbstractGymEnv.o CGym/DungeonRoomEnv.o CGym/pch.o CGym/SimpleRepeat.o CGP_CPP/AbstractCGPIndividual.o CGP_CPP/AbstractTesterClass.o CGP_CPP/BasicGene.o CGP_CPP/CGPFFANNIndividual.o CGP_CPP/CGPFunctions.o CGP_CPP/CGPGenerator.o CGP_CPP/GeneralCGPSolver.o CGP_CPP/GymTester.o CGP_CPP/GymWorldGenerator.o CGP_CPP/pch.o CGP_CPP/Enums.o ExperimentRunner/Experiment.o ExperimentRunner/MultiExperiment.o ExperimentRunner/pch.o CGP_CPP_TestApp/CGP_CPP_TestApp.o $(CPPMONOLIBFLAGS)

#Clean
clean: 
	rm -f libraries/*.a && rm -f CGP_CPP/*.o && rm -f CGym/*.o && rm -f ExperimentRunner/*.o && rm -f CGP_CPP_TestApp/*.o
