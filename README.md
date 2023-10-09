# analyticRT
Models for Regge trajectories and scattering amplitudes which satisfy S-matrix constraints.

Compilation of the base library requires only [ROOT](https://root.cern.ch/) (tested with version 6.17 and 6.24) with [*MathMore*](https://root.cern.ch/mathmore-library) and [Boost C++](https://www.boost.org/) (version $\geq$ 1.68) libraries.

##  INSTALLATION
To install clone normally and use:
```bash
cd analyticRT
mkdir build && cd build
cmake ..
cmake --build . --target install
```
This will create the core library `/lib/libANALYTIC.so` with the linkable library as well as ROOT dictionary (.pcm) files. 

Additionally a [scripting executable](./src/cling/analyticRT.cpp) will be installed into `/bin/analyticRT` which short-cuts loading the libraries into an interactive ROOT session and running a .cpp file as an interpreted script.   
This executable requires the environment variable `ANALYTICRT` to be set to the top-level directory in order to find auxilary files. This can be done as such:
```bash
export ANALYTICRT=/path/to/analyticRT # for bash
setenv ANALYTICRT /path/to/analyticRT # for csh
```


##  USAGE

To run a script use the executable described above which mimics a Python-like environment without requiring recompilation when changes are made to model files:
```bash
analyticRT my_script.cpp
```
or add the bin directory to $PATH to call `analyticRT` from any directory. 