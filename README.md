# centrality

#To setup the eic Library:

#Make your local analysis folder somewhere:

#$ mkdir -p eic/sofware
#$ cd eic/software
#$ git clone https://gitlab.com/eic/eic-smear.git
#$ mkdir build
#$ mkdir install
#$ cd build
#$ export EICDIRECTORY=/Users/carorg/projects/eic/software/install

#$ cmake ../eic-smear/ -DCMAKE_INSTALL_PREFIX=$EICDIRECTORY
#$ make
#$ make install

##Then, the tricky part: Create a file, inside analysis folder, that will setup your "EIC" environment, name it setupEnvironment.C with the following content:

##void setupEnvironment()
##{
##    gInterpreter->AddIncludePath("-I/Users/carorg/projects/eic/software/install/include/");
##    gSystem->Load("/Users/carorg/projects/eic/software/install/lib/libeicsmear.dylib");
##}

##Then, you are ready to work, enter ROOT, and type '.x setupEnvironment.C', that tell ROOT where the headers and the libraries are. You are ready to work, for example, to execute '.x runEICTree.C'. 

##You will always have to set up the environment yourself before executing your jobs.

#To run simulations on BeAGLE i need to used the BNL CLuster.
#To run jobs
#.x setupEnvironment.C every time that exit root
#root -l runEICTree.C++ "Analysis Code to obtain kinematical variables"
#Analysis_ratios.C++ "Analysis code to obtain multiplicity ratios of the differents root files for Target and Deutirium"

