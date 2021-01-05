# ORCHESTRA
Southern Ocean 1/12 NEMO configuration

## Quick Start:

```
git clone git@github.com:jdha/ORCHESTRA.git
./ORCHESTRA/scripts/setup/ORCHESTRA_setup -w $PWD/orch -x $PWD/orch -s $PWD/ORCHESTRA -m archer
cd orch/nemo/CONFIG/CORE2NYF-ORCH0083-LIM3/EXP_ZPS
```
Edit the project code in  `runscript.pbs` then:
```
qsub runscript.pbs
```

### Forcing/Input data:

[ORCHESTRA](http://gws-access.ceda.ac.uk/public/jmmp_collab/ORCHESTRA)

_this is automatically transferred when the setup script is executed_

### Important:

Note that this repository only holds the CORE2 Normal Year Forcing Experiment. There are, at present, two run directories. One for the ZPS vertical coordinates, the other for the SZT hybrid coorinates. At the moment the namelist and compiler options are as close to the ZPS setup as possible. 

This repository is unlike similar ones for NOC-MSM in as much as it contains most of the NEMO code base. We cannot reference the NEMO repository as the revision number used for ORCHESTRA no longer exists. This is due to a rebasing of the SVN trunk onto a parallel development branch.
