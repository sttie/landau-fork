#!/bin/bash

ERAPARH=../era
BUILDCONF=Release64OSX
BUILDPATH=Release64OSX/GNU-MacOSX/libmoons.dylib
SATNAME=nep-sat
SATPATH=nepsat

set -e

echo "Compling .dau to .c"
time racket ${SATNAME}.dau
mv ${SATNAME}.c ./moons/src/
cd moons

echo "Building libmoons"
make CONF=${BUILDCONF}
cd ../../

echo "Coping libmoons.dylib to era/lib64"
cp ./tests/moons/dist/${BUILDPATH} ${ERAPARH}/lib64/

cd ${ERAPARH}/${SATPATH}
echo "Checkout to landau-test branch"
git checkout landau-test

echo "Computing O-C and LSM parameters to compare"
sh ./iteration.sh 

echo "Converting lsm/* to .csv"
cd utils
racket db_to_csv.slon > lsm_db_new.csv

echo "Comparing lsm_db_new.csv with lsm_db_old.csv"
if [[ $(diff lsm_db_old.csv lsm_db_new.csv) == "" ]]
then
  echo "Integration test -- PASSED"
else
  echo "Integration test -- FAILED"
  printf "To open vimdiff press y:"
  read -r
  if [[ $REPLY == "y" ]]
  then
    vimdiff -u NONE lsm_db_old.csv lsm_db_new.csv
  fi
fi
