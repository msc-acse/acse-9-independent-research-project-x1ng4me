# Parallelization of iSALE-2D Using OpenMP
### **Independent Research Project**

## **Introduction**
iSALE (impact-SALE) is a well-established shock physics code (written in Fortran 95) that has been developed and used for more than two decades to simulate impact cratering. This project is a development of iSALE-2D, which is an extension of source code with refactoring and adding OpenMP.

Since the source code of iSALE-2D is not public, the permission is **required** before installing iSALE and doing further development. The recommended platform for using iSALE is Ubuntu or Mac OS X.

## **Instructions for installing all the necessary iSALE pre-requisites for macOS 10.13 (High Sierra)**
### **Install Xcode from the mac app store**
### **Install the command line tools:**
```
xcode-select --install
```
### **Install XQuartz**
- Download and install XQuartz from here: [XQuartz](http://xquartz.macosforge.org/landing/)  
- You will need to logout and back in again after this.
### **Install Intel Fortran Compiler ifort**
- Download and install ifort from here: [ifort](https://software.intel.com/en-us/fortran-compilers)
### **Install gcc**
- First install `gcc`. As of Nov 2017, this installs gcc 7.2.0, which seems to work fine for iSALE.  
```
sudo port install gcc7
```  
- Tell your mac to use this gcc version and not the Xcode version.  
```
sudo port select --set gcc mp-gcc7
```  
At this point, check you are using the correct version by running `gcc -v` or `ifort -v`. It should tell you that you are using the right version here. If you see that it is the Xcode version, try selecting it again, or restarting the terminal.
### **Install Vtune Amplifier**
Intel VTune Amplifier is a performance analysis tool targeted for users developing serial and multithreaded applications. The tool is delivered as a Performance Profiler with Intel Performance Snapshots and supports local and remote target analysis on the Windows, Linux, and Android platforms. Though you cannot analyze applications running on the macOS systems, you still can install the VTune Amplifier on macOS and analyze remote targets. Therefore, a remote machine with Linux OS system is required.  
- Download and install Vtune Amplifier from here: [Vtune](https://software.intel.com/en-us/vtune)

At this point you should be **ready to install iSALE**.
## **Getting Started**
### **Installing and Compiling of iSALE**
Once the permission is granted, please clone the iSALE root directory from the GitHub.  
The next step of the installation is to configure iSALE, which will attempt to install iSALE and all its associated packages.  
- You might install iSALE by using the following order to meet all requirements:  
```
FC=ifort CC=icc CXX=icpc ./configure --prefix=<your install path> --with-pysaleplot --openmp
```  
- If configuration succeeded, use following command to compile iSALE.  
```  
make  
```  
This command will create all the programs (subprograms and libraries) needed to run iSALE.  
- Note that if you are reconfiguring (e.g., after an update or with different configure options) it is prudent to make clean before recompiling:  
```  
make clean  
make  
```  
- To complete the installation of iSALE it is useful to copy all the files needed to run iSALE to a separate directory. This is done using the command:
```  
make install  
```  
Then the iSALE is fully installed and ready for developments.
### **Installing Development during the IRP**
- Please download the developed source code:  
```  
git clone https://github.com/msc-acse/acse-9-independent-research-project-x1ng4me
```  
- Find the following files in the `/iSALE/src/2D` folder.
```  
advect.F90  
advect_routines.F90
setup_initialize.F90
```
- Replace them with the files in the IRP folder cloned from the GitHub.  
- The use the following commands to reconfiguring iSALE.  
```  
make clean  
make  
make install  
```  
For now, all changes will be installed to under the iSALE folder. The editor for `Fortran 95` is `vim`, you can directly edit or view the code by:  
```  
vim advect.F90 (for example)
```  

## **Running Vtune Amplifier Analysis**
- Please copy the `bash` file to the folders of sample cases. This file will be used to run OpenMP jobs in Vtune Amplifier. You can change how many threads you want to use in this file.
```  
cp runscript.bash <your file path>  
```  
- Vtune Amplifier UI interface has some issues with running multi-threads' jobs. To get the summary of the ruuning status, hot-spots analysis and hyperperformance reports, please run the Vtune Amplifier in the command line with:  
```  
cd <your install directory for the sample cases>  

/opt/intel/vtune_amplifier/bin32/amplxe-cl -collect hotspots -app-working-dir <install directory><sample case name> -- <install directory><sample case name>/runscript.bash  
```  
The Vtune Amplifier will create a file folder named like `r001hs` in your directory for the sample cases.  
Please copy the whole folder from the remote machine to your local computer to check the results.  
```  
scp -r <host-name>:<your install directory for the sample cases>/r000hs .  
```  
Then, in the Vtune Amplifier UI interface, please click `open results` on the top left corner, and locate the folder you copied to your computer, open the `*.amplxe` file.

**Note:** Please make sure copy the whole folder not only the `*.amplxe` file to the local computer.

## **Repository Information**
- Developments - The source code after refactoring and parallelization  
  - Old Version - The previous code files I got at the beginning of the IRP, with comments editing  
  - Mid-term Version - The intermedia version of the code files which finished part of the code refactoring, but before any loops merging and parallelization.  
  - Latest Version - The latest version of the code files for final submission, which includes the parallelization, code refactoring and documentation.  
- Script - The running script which should be used for analysis  
  - runscript.bash - The script for running Vtune Amplifier to create analysis report  
  - vtune_output.txt - The output after each test, which includes the the information of the parallelization, and test status  
- Tests - All tests results during the IRP
  - test.ipynb - The jupyter notebook which includes all tests' output.  

## **Build With**
Fortran 95  
Python 3

## **Author and Course Information**
**Author:** Xianzheng Li  
**CID:** 01611306  
**GitHub Alias:** x1ng4me  
**Course** : ACSE-9 Independent Research Project  
**Supervisor** : Prof. Gareth Collins  
  
## **Acknowledgments**
- Apologies for any inconvenience due to the permission of the source code  
- Thanks to the support from my supervisor Prof. Gareth Collins and Dr. Thomas Davison.  
- Thanks for the advice from Prof. Gerard Gorman.  
