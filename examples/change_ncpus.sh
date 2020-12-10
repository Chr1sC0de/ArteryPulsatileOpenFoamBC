#!/bin/bash

# Recursively searches some folder for ncpus and replaces with the desired number of cpus

while getopts ":d:n:" opt; do
  case $opt in
    d) top_dir="$OPTARG"
    ;;
    n) n_cpus="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ -z $top_dir ]
then
    top_dir="."
fi

if [ -z $n_cpus ]
then
    n_cpus="8"
fi

oldstring="ncpus=[0-9]+"
newstring="ncpus=$n_cpus"

stuff="s/$oldstring/$newstring/g"

grep -E -r -l $oldstring $top_dir | xargs sed -i -E $stuff 

oldstring="numberOfSubdomains [0-9]+"
newstring="numberOfSubdomains $n_cpus"

stuff="s/numberOfSubdomains[[:space:]][1-9]+/numberOfSubdomains_$n_cpus/g"

grep -E -r -l "$oldstring" $top_dir | xargs sed -i -E $stuff
grep -E -r -l "numberOfSubdomains_[0-9]+" $top_dir | xargs sed -i -E 's/_/ /g'
