#!/bin/sh
#
BINDIR=/home/ktagami/work/ChargedDefect/Sxdefectalign/bin
#
$BINDIR/sxdefectalign --charge 2 --eps 12.88 --vdef ../electrostatic_pot.cube --vref ../../../Pristine/q_0/electrostatic_pot.cube --center 0.125,0.125,0.125 --relative --ecut 270 --qe > Log.1
