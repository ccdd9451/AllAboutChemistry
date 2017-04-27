#!/home/kh621/.installed/anaconda/bin/python

from xilio import write
from subprocess import call, Popen
from pathlib import Path
import sys

BASEDIR = Path('/home/kh621/PeptideSim')

def main():
    singlerun()



def singlerun():
    from xilio import load, dump
    from os import environ
    env = environ
    if not env.get('CUDA_VISIBLE_DEVICES'):
        process_num = input(
            'Please specify the gpu index for calculations:')
        env['CUDA_VISIBLE_DEVICES'] = str(process_num)
        topfile = BASEDIR/'top'
        try:
            top = load(topfile)
        except:
            top = set()

        if process_num in top:
            raise OSError(
                'GPU line {} is busy!'.format(process_num))
        else:
            top.add(process_num)
        dump(topfile, top)
    else:
        process_num = env.get('CUDA_VISIBLE_DEVICES')

    with shelf_with_locker() as shelf:
        workpep = shelf['queued'].pop()
        shelf['running'].append(workpep)

    workdir = BASEDIR/workpep
    pepname = convert_short_to_long(workpep)

    setupfiles(workdir, pepname, env)
    call(['./ambsc'], cwd=workdir, env=env)

    with shelf_with_locker() as shelf:
        shelf['running'].remove(workpep)
        shelf['finished'][workpep] = None
        # Will update here to directly get the protential
        # energy after simulation finishs

    if process_num in load(topfile):
        Popen([__file__], cwd=BASEDIR, env=env)
    sys.exit(0)

from contextlib import contextmanager
@contextmanager
def shelf_with_locker():
    from shelve import open as sopen
    from fcntl import flock, LOCK_EX, LOCK_UN
    shelf_locker = BASEDIR/'slf.lck'
    pep_queue = BASEDIR/'pep_queue'
    try:
        lck = open(shelf_locker, w)
        flock(lck, LOCK_EX)
        shelf = load(pep_queue)
        yield shelf
    except:
        raise
    finally:
        dump(pep_queue, shelf)
        flock(lck, LOCK_UN)
        lck.close()
    return workpep

def setupfiles(directory, pseqs, env):
    import os
    try:
        os.mkdir(directory)
    except OSError as e:
    # Skip mkdir if dir is exist, only print out error
        print(e)
    write(directory/'ambsc', script, executable=True)
    write(directory/'tlsc', tleapfile.format(names=pseqs))
    call(['tleap','-s','-f','tlsc'], cwd=directory, env=env)


Aminonames = (
    'ALA ARG ASN ASP CYS ' +
    'GLU GLN GLY HIS HYP ' +
    'ILE LEU LYS MET PHE ' +
    'PRO GLP SER THR TRP ' +
    'TYR VAL').split()
Aminoshorts = \
    'ARNDCEQGHOILKMFPUSTWYV'
assert len(Aminonames) == len(Aminoshorts)
Aminodict = { s:n for s, n in \
        zip(Aminoshorts, Aminonames)}

def convert_short_to_long(short):
    long = map(lambda x: Aminodict[x],
            list(short))
    long = ' '.join(long)
    return long

tleapfile = """\
source leaprc.ff14SB
foo = sequence {{ {names} }}
set default pbradii mbondi3
saveamberparm foo prmtop inpcrd
quit
"""

script = """\
#!/bin/bash

callcuda(){
pmemd.cuda -O \\
	-i 1in -o 1out \\
	-p prmtop -c inpcrd \\
	-r 1rst -inf 1mdinfo

pmemd.cuda -O \\
	-i 2in -o 2out \\
	-p prmtop -c 1rst \\
	-r 2rst -x 2.mdcrd \\
	-inf 2mdinfo

pmemd.cuda -O \\
	-i 3in -o 3out \\
	-p prmtop -c 2rst \\
	-r 3rst -x 3.mdcrd \\
	-inf 3info
}

expandfiles(){

echo \\
'Minimize
 &cntrl
  imin=1, !Choose a minimization run
  ntx=1, !Read coordinates but not velocities
  maxcyc=2000, !Maximum minimization cycles
  ncyc=1000, !The steepest descent algorithm for the first 0-ncyc cycles
  ntpr=100, !Print to the Amber mdout output file every () cycles
  igb=8, !Generalized Born function is set to 8
  cut=99999.0 !GPU Amber asked to do so, MAKE IT CLEAR SOONER
 /' > 1in

echo \\
'Heating
 &cntrl
  imin=0, !Choose a molecular dynamics (MD) run
  ntx=1, !Read coordinates but not velocities
  irest=0, !Restart the simulation from rst file
  nstlim=10000, !Number of MD steps in run
  dt=0.002, !Time step in picoseconds
  ntf=2, !Setting to not calculate force for SHAKE constrained bonds
  ntc=2, !Enable SHAKE to constrain all bonds involving hydrogen
  tempi=0.0, !Initial thermostat temperature in K
  temp0=300.0, !Final thermostat temperature in K
  ntpr=100, !Every ntpr steps, mdinfo and mdout are printed
  ntwx=100, !Every ntwx steps, the coordinates will be written to the mdcrd file
  ntb=0, !no periodicity is applied and PME is off
  ntp=0, !No pressure scaling
  ntt=3, !defines a Langevin thermostat
  igb=8, !Generalized Born function is set to 8
  gamma_ln=2.0, !Collision frequency in ps−1
  nmropt=1, !NMR restraints and weight changes will be read.
  ig=-1, !Randomize the seed for the pseudo-random number generator
  cut=99999.0 !GPU Amber asked to do so, MAKE IT CLEAR SOONER
 /
&wt type='TEMP0', istep1=0, istep2=9000, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=9001, istep2=10000, value1=300.0, value2=300.0 /
&wt type='END' / ' > 2in

echo \\
'Production
 &cntrl
  imin=0, !Choose a molecular dynamics (MD) run
  ntx=5, !Coordinates and velocities are read from rst file
  irest=1, !Restart the simulation from rst file
  nstlim=5000000, !Number of MD steps in run // 5e6 frames are 10 ps
  dt=0.002, !Time step in picoseconds
  ntf=2, !Setting to not calculate force for SHAKE constrained bonds
  ntc=2, !Enable SHAKE to constrain all bonds involving hydrogen
  temp0=300.0, !Initial thermostat temperature in K
  ntpr=100, !Every ntpr steps, mdinfo and mdout are printed
  ntwx=100, !Every ntwx steps, the coordinates will be written to the mdcrd file
  ntb=0, !no periodicity is applied and PME is off
  ntp=0, !No pressure scaling
  ntt=3, !defines a Langevin thermostat
  igb=8, !Generalized Born function is set to 8
  gamma_ln=2.0, !Collision frequency in ps−1
  ig=-1, !Randomize the seed for the pseudo-random number generator
  cut=99999.0 !GPU Amber asked to do so, MAKE IT CLEAR SOONER
 /' > 3in
}

processdat(){
mkdir Analysis; cd Analysis
process_mdout.perl ../3out
cd ..
}

cleanupfiles(){
mkdir archieved
mv 1* archieved
mv 2* archieved
mv tlsc archieved
}

touch date
echo 'start at: '`date` >>timelog
############ Script starts here #############
#
#
#

expandfiles
callcuda
processdat
cleanupfiles

#
#
#
############ Script running ends ############
echo 'end at: '`date` >>timelog
exit
"""

if __name__ == '__main__':
    main()
