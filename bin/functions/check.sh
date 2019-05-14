#!/bin/bash

check_dependencies(){
  counter=0
  printf '\n%s\t%20s\n' "DEPENDENCY" "STATUS"
  printf '%s\t%20s\n'   "----------" "------"

  for command in "$@"
  do
    length_command=$(echo $command | wc -m)
    distance_table=$((30 - $length_command))
    distance_expression=$(echo "%${distance_table}s")
    printf '%s' $command
    if ! [ -x "$(which $command 2> /dev/null)" ]; then
      printf $distance_expression
      printf "NOT INSTALLED\n"
      let counter++
    else
      printf $distance_expression
      printf "INSTALLED\n"
    fi
  done

  if [ $counter -gt 0 ]; then
    printf "ERROR: $counter missing dependencies, aborting execution\n" &>> $log_file
    exit 1
  fi
}

check_dockers(){
  counter=0
  printf '\n%s\t%20s\n' "DOCKER" "STATUS"
  printf '%s\t%20s\n' "----------" "------"

  for docker in "$@"
  do
    length_docker=$(echo $docker | wc -m)
    distance_table=$((30 - $length_docker))
    distance_expression=$(echo "%${distance_table}s")
    printf '%s' $docker

    if [[ "$(docker images -q $docker 2> /dev/null)" == "" ]]; then
      printf $distance_expression
      printf "NOT INSTALLED\n"
      let counter++
    else
      printf $distance_expression
      printf "INSTALLED\n"
    fi
  done

  echo $counter
  if [ $counter -gt 0 ]; then
    printf "ERROR: $counter missing dockers, aborting execution\n" &>> $log_file
    exit 1
  fi
}
