#!/bin/bash

TEMPERATURE=$1

if [ $# -eq 0 ]
then
  printf "Usage:\n\tcreate-graphics.sh temperature\n\n"
  exit 1
fi

mkdir -p graphics/ graphics-out

tail -n +4 gamd.log > graphics/headerless-gamd.log
awk 'NR%1==0' graphics/headerless-gamd.log |awk -v TEMP=$TEMPERATURE '{print ($8+$7)/(0.001987*TEMP)" " $2 " " ($8+$7)}' > graphics/weights.dat
cpptraj -p ../data/dip.prmtop  -y output.dcd -i psi-dat-commands.cpptraj
awk '{print $2}' graphics/psi-cpptraj.dat |tail -n +2 >psi.dat

cpptraj -p ../data/dip.prmtop  -y output.dcd -i phi-dat-commands.cpptraj
awk '{print $2}' graphics/phi-cpptraj.dat |tail -n +2 >phi.dat

cpptraj -p ../data/dip.prmtop  -y output.dcd -i phi-psi-commands.cpptraj
awk '{print $2, $3}' graphics/phi-psi-cpptraj.dat |tail -n +2 > phi-psi.dat



PyReweighting-1D.py -input psi.dat -T $TEMPERATURE -cutoff 10 -Xdim -180 180 -disc 6 -Emax 20 -job amdweight_CE -weight graphics/weights.dat | tee -a graphics/reweight-variable-cumulant-expansion-1D.log
mv -v pmf-c1-psi.dat.xvg graphics-out/pmf-psi-reweight-CE1.xvg
mv -v pmf-c2-psi.dat.xvg graphics-out/pmf-psi-reweight-CE2.xvg
mv -v pmf-c3-psi.dat.xvg graphics-out/pmf-psi-reweight-CE3.xvg
mv -v psi.dat graphics/



PyReweighting-1D.py -input phi.dat -T $TEMPERATURE -cutoff 10 -Xdim -180 180 -disc 6 -Emax 20 -job amdweight_CE -weight graphics/weights.dat | tee -a graphics/reweight-variable-cumulant-expansion-1D.log
mv -v pmf-c1-phi.dat.xvg graphics-out/pmf-phi-reweight-CE1.xvg
mv -v pmf-c2-phi.dat.xvg graphics-out/pmf-phi-reweight-CE2.xvg
mv -v pmf-c3-phi.dat.xvg graphics-out/pmf-phi-reweight-CE3.xvg
mv -v phi.dat graphics/


PyReweighting-2D.py -T $TEMPERATURE  -cutoff 10 -input phi-psi.dat -Xdim -180 180 -discX 6 -Ydim -180 180 -discY 6 -Emax 20 -job amdweight_CE -weight graphics/weights.dat | tee -a graphics/reweight-variable-cumulant-expansion-2D.log
mv -v pmf-c1-phi-psi.dat.xvg graphics-out/pmf-2D-phi-psi-reweight-CE1.xvg
mv -v pmf-c2-phi-psi.dat.xvg graphics-out/pmf-2D-phi-psi-reweight-CE2.xvg
mv -v pmf-c3-phi-psi.dat.xvg graphics-out/pmf-2D-phi-psi-reweight-CE3.xvg
mv -v 2D_Free_energy_surface.png graphics-out/pmf-2D-phi-psi-reweight-CE2.png
mv -v phi-psi.dat graphics/

mv -v weights-c1-psi.dat.xvg graphics-out/
mv -v weights-c2-psi.dat.xvg graphics-out/
mv -v weights-c3-psi.dat.xvg graphics-out/
mv -v weights.png graphics-out/weights-psi.png

