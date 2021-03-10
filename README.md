# spinContaminationChecker
A program to carry out spin contamination checks, analysis, and related models.

To run a calculation from a Gaussian matrix element file named `matrixFile.mat`, use the command
```
spinContaminationChecker.exe matrixFile.mat
```

An optional _second_ command line argument can be given to set a print flag. If the print flag is set to -1, then DEBUG print mode is turned on. The default print flag value is 0.

An optional _third_ command line arguement can be given to set the number of openMP processes to use (assuming an openMP version has been compiled). By default, this value is set to 1, or serial mode.

As an example, to run a calculation from a Gaussian matrix element file named `matrixFile.mat` at print level 1 with 48 cores one uses the command
```
spinContaminationChecker.exe matrixFile.exe 1 48
```
