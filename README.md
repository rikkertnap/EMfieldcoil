## Description


Calculates EM field of a coil in quasi-static approximation

### Requirement  

matplotlib, numpy scipy, argparse

### Running

python EMfields.py 


usage: EMfields.py [-h] [--R R] [--N N] [--L L] [--I0 I0] [--f F] [--save]

If no values specified it uses default values, Namely <br>


  --R R     Coil radius (m) default 0.134 m <br>
  --N N     Number of turns default 500 <br>
 --L L      Axial width (m) default 0.02 m <br>
  --I0 I0   Current amplitude (A) default 0.05 A <br>
  --f F     Frequency (Hz) default 200.0 Hz <br>
  --save    Save figures <br>
