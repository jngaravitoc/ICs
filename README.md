# ICs
This is a hacked version of the initial conditions code GalIC.

The original code can be found at: 

http://www.h-its.org/tap-software-e/galic-code/


This hacked version includes:

1. The correct NFW and Henrquist equivalence conversion following the Appendix of van der Marel (2012) http://adsabs.harvard.edu/abs/2012ApJ...753....8V

2. Now all the quantities (R, M, A etc..) are calculated at the virial radius. 

3. The Hubble constant is set to 70/km/s/Mpc
(set in Gadget h=1.0)

4. The Disk and the Bulge scale length are now input parameters.
