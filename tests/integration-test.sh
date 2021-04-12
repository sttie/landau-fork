#!/bin/bash

ERAPARH=../era
BUILDCONF=Release64OSX
BUILDPATH=Release64OSX/GNU-MacOSX/libmoons.dylib
SATNAME=nep-sat
SATPATH=nepsat

set -e

echo "Compling .dau to .c"
time racket ${SATNAME}.dau -c -extfl
mv ${SATNAME}.c ./moons/src/
cd moons

echo "Building libmoons"
make CONF=${BUILDCONF}
cd ../../

echo "Coping libmoons.dylib to era/lib64"
cp ./tests/moons/dist/${BUILDPATH} ${ERAPARH}/lib64/

cd ${ERAPARH}/${SATPATH}
ERA_BRANCH=$(git rev-parse --abbrev-ref HEAD)
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
  echo "Integration test C -- PASSED"
else
  echo "Integration test C -- FAILED"
  printf "To open vimdiff press y:"
  read -r
  if [[ $REPLY == "y" ]]
  then
    vimdiff -u NONE lsm_db_old.csv lsm_db_new.csv
  fi
fi

echo "Back to initial era branch"
git checkout ${ERA_BRANCH}
