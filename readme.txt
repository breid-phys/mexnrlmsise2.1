A MATLAB MEX wrapper for the NRLMSISE2.1 Empirical Model

A MEX file allows users to execute FORTRAN code from within MATLAB.
This repository adds two MATLAB function, one to interface with the MEX code, 
the other to automatically download the apf107.dat file to provide needed indices
to MSIS. 

nrlmsise2_1.m should automatically compile a MEX file when first run.

See nrlmsise2_1.m for details.

#######################################################################
 MSIS® (NRL-SOF-014-1) SOFTWARE
 NRLMSIS® empirical atmospheric model software. Use is governed by the
 Open Source Academic Research License Agreement contained in the file
 nrlmsis2.1_license.txt, which is part of this software package. BY
 USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
 CONDITIONS OF THE LICENSE.  
#######################################################################

NRLMSIS 2.1 Whole-Atmosphere Empirical Model of Temperature and Neutral Species
  Densities

VERSION HISTORY
  08 MAR 19 Version 1.97 (Beta version)
  26 MAY 20 Version 2.0 (Release version)
  04 APR 22 Version 2.1 (Release version with NO density)

AUTHORS
  John Emmert (john.emmert@nrl.navy.mil)
  Douglas Drob (douglas.drob@nrl.navy.mil)
  McArthur Jones Jr. (mcarthur.jones@nrl.navy.mil)
  
REFERENCE FOR NRLMSIS 2.0
  Emmert, J. T., Drob, D. P., Picone, J. M., Siskind, D. E., Jones, M. Jr., 
  Mlynczak, M. G., et al. (2021). NRLMSIS 2.0: A whole-atmosphere empirical model
  of temperature and neutral species densities. Earth and Space Science, 8, 
  e2020EA001321. https://doi.org/10.1029/2020EA001321

ACKNOWLEDGEMENTS
  This work was supported by the Office of Naval Research and NASA.
