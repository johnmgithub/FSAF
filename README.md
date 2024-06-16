This is a Java port of parts of the **Fast Subband Adaptive Filtering** process developed by Michael Tsiroulnikov aka Michael Zrull. 
See Michael Zrull (2020) Michael Tsiroulnikov (2022) Fast Subband Adaptive Filtering (FSAF) [MATLAB Central File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/83363-fast-subband-adaptive-filtering-fsaf) and Loudspeakers for AEC: Measurement and Linearization [MATLAB Central File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/117715-loudspeakers-for-aec-measurement-and-linearization">www.mathworks.com/matlabcentral/fileexchange/117715-loudspeakers-for-aec-measurement-and-linearization). 
The MATLAB code is Copyright (C) 2001-2022 MICHAEL ZRULL (TSIROULNIKOV). All rights reserved. 
Direct or indirect commercial derivatives or usage of the code is prohibited without the express written permission of Michael Tsiroulnikov. 
The MATLAB code was released under the GNU General Public License v.3+ [www.gnu.org/licenses](https://www.gnu.org/licenses/).
That license also applies to this Java port.

The Java code is copyright (c) 2024 John Mulcahy, all rights reserved. 

The port does not include any audio interface handling. See the main method in FsafSpkid.java for example usage.

Extended versions of Jama and Jampack are used for Matrix handling. 
The Parallel class of Dave Hale's [Mines Java Toolkit](https://inside.mines.edu/~dhale/jtk/index.html) is used for parallel execution of the filter band processing.
The HebiRobotics [MAT File Lbrary](https://github.com/HebiRobotics/MFL) is used to load MAT files.

There is a complete implementation of measurement using FSAF in [Room EQ Wizard](https://www.roomeqwizard.com), versions V5.40 beta 32 onwards.
