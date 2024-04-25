#!/bin/bash
if [ ! -d "outputs/custom" ]; then
	echo "Creating outputs/custom directory"
	mkdir -p "outputs/custom"
fi
docker build -t pnnl-dsse-ekf:custom . -f Dockerfile.custom
docker volume create --name oedisi_outputs --opt type=none --opt device=${PWD}/outputs/custom --opt o=bind
docker run --rm -it --mount source=oedisi_outputs,target=/home/outputs pnnl-dsse-ekf:custom