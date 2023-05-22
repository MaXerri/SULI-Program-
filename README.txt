Code from Crystallography projects during my Science Undergraduate Labratory Internship program (Summer 2022).
I initially coded this coversion process in python and then remade the code in C++ while learning the language.  The compiled c++ code ran significantly faster than my original python code. The C++ scipts were the final product of the programalong with my final presentation and report which can be found on my linkden. The C++ code can be found in the folder dc7unsrt_c++_conversions and the original python scripts are in this curent directory.

__________________________________________________________________________________________________________________________________________________________________

dc7unsrt conversions for the ncdist repository 
This additional code was made by Mario A. Xerri for Herbert J Bernstein and Lawrence C. Andrews's ncdist repository
August 2022 
xmario02@yahoo.com 

The research perfromed use the Wigner Seitz representation of the unit to define a 7- dimesnion representation of the cell which includes teh 3 edge lengths, 3 shortest face diagonals and the shortest body diagonals.  Starting from a Niglgi reduced cell in its G6 representation (Andrews and Bernstein, 1988), the Wigner Seitz cell can be constructed.

These scripts convert text file lists of PDB cell parameter entries to their dc7unsrt representation and then recover the Niggli cell from dc7unsrt.  Recovered Nigli-cells are compared to thier original to ensure consistency and each step of the conversion process in shown in an outputted excel spreadsheet. To learn more about the functionality of the ncdist repository look at th README.txt after you download the ncdist package (steps for this are shown below). These scripts use functions from the ncdist repository so it must be downloaded before running converter_main.cpp which has the main function for these scripts.  

_______________________________________________________________________________________

ADDING NEW FILES TO READ

In the main_converter.cpp code, to read a new textfile with cell parameters, change crystal,txt to the name of the file you are going to read in the scipt and make sure itis formatted in the same way as crystal.txt is (this file can be found in within this current directory).  

_________________________________________________________________________________________

RUNNIGN THE SCRIPTS

From the build directory run the following commands:
cmake .
make 
./main_executable

The output csv files will be in the build directory after running the script 

_______________________________________________________________________________________

INSTALLATION

Preliminaries:  You need a development system with cmake, C, C++, Fortran, R, Rcpp, 
RcppParallel and RcppArmadillo installed.

Download the package from https://github.com/yayahjb/ncdist.git

If you downloaded as a zip, unpack the kit.  In any case if you have the source in
ncdist

cd ncdist/build
cmake -DCMAKE_INSTALL_PREFIX=<prefix> ..
make all
make install

will install the kit in <prefix>, specifically, it will 

install ncdist in <prefix>/bin/ncdist
install librcpp_ncdist in <prefix>/lib/librcpp_ncdist.<shared_library_extension>

In most cases you will want <prefix> to be $CCP4
