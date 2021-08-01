echo "integrating spacecraft.dau"
time racket test_spacecraft.rkt > test-spacecraft-out-new.txt
echo "comparing test-spacecraft-out-new.txt with test-spacecraft-out.txt"
if [[ $(diff test-spacecraft-out-new.txt test-spacecraft-out.txt) == "" ]]
then
  echo "Integration test racket -- PASSED"
else
  echo "Integration test racket -- FAILED"
  printf "To open vimdiff press y:"
  read -r
  if [[ $REPLY == "y" ]]
  then
    vimdiff -u NONE test-spacecraft-out-new.txt test-spacecraft-out.txt
  fi
fi
