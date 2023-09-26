# Templates for the calibration of elemental mercury gas generators
This repository contains two software templates for the calibration of elemental mercury gas generators. 
One template is for the single point calibration, and one for the multi point calibration.
The template for the single point calibration calculates 

* the output ratio between the measurements of the candidate and the measurements of the reference generator
* the corrected candidate concentration
* the complete uncertainty budget associated to the measurement.

The template for the multi point calibration calculates the above listed items for each setpoint plus

* the coefficients of the calibration curve
* the uncertainty associated to the coefficients of the calibration curve.

Details about the equations implemented can be found in Protocol for the *SI-traceable calibration of elemental mercury (Hg<sup>0</sup>) gas generators used in the field*.

## Documentation
The documentation of the templates is included in this repository (see PDF file *documentation_20230921.pdf* for the latest version).
To facilitate future updates of the documentation, the LaTex source file (*SI-Hg documentation.zip*) used to generate the PDF file is also included in this report.

## Requirements
Details about the software requirements necessary to execute the code are listed in the documentation.
In case of modification of the code, please update the listed requirements.

There are no requirements regarding operating system or hardware.

## Developement period
The scripts have been developed in 2022 and in 2023. 
Version 1.0 of the templates has been finalized in Septemeber 2023.

## Software validation
Both the single point and the multi point template have been validated. 
See the validation report for details and results of the validation.

## Automated test execution
There is no automated test execution.

## Software design
Both the single point and the multi point template consist of an Excel file which calls a Python script. 
The data are to be inserted directly in the Excel file. The Python script reads and processes the data. 
The results of the data processing are displayed in new sheets in the Excel file.
For more information, read the documentation.

## Related projects
EMPIR 19RNM03 SI-Hg

## License
This code is provided under the Creative Commons Attribution 4.0 International Public License.
