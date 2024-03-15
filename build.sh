#!/bin/bash
scenario="$1"

if [ "$scenario" = "" ];
then
	echo You must enter a scenario name
else
    if [ ! -d "outputs/$scenario" ]; then
   		echo "Creating $scenario directory"
   		mkdir "outputs/$scenario"
	fi
	docker build --build-arg SCENARIO=$scenario -t dsse-$scenario:1.1.0 . -f Dockerfile
    docker volume create --name oedisi_outputs --opt type=none --opt device=${PWD}/outputs/$scenario --opt o=bind
    docker run --rm -it --mount source=oedisi_outputs,target=/home/outputs dsse-$scenario:1.1.0
fi