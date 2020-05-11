README:

All files in this repository are required in order to run the code to generate data points that map velocity, vorticity and streamfunction values at different grid locations. 

The documentation for what each file does is written within the source code, do refer to them for help.

Contact: I can be reached at andrewng789@gmail.com for further clarifications.

Input: To vary initial and run time conditions of the flow, these will be done within the ProgramOptions.cpp file.
Output: After successful running of the code, you will be presented with mutiple .txt files with data of velocity, vorticity and streamfunction for all defined timesteps through out the time interval specified OR once convergence is reached, which ever comes first. 

*********
* Linux *
*********

Once you have navigated to the directory containing all source files, including CMakeLists.txt, the following commands are to be inputted: 
__________________________________________________
\.cmake 
make
run FullSolver
__________________________________________________