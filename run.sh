#!/bin/bash
scenario="$1"
cd /home/

if [ "$scenario" = "" ];
then
	echo You must enter a scenario name
else
	if [ ! -d "outputs/$scenario" ]; then
   		echo "Creating $scenario directory"
   		mkdir "outputs/$scenario"
	fi
	pkill -9 helics_broker
	pkill -9 python
	oedisi run --runner build_$scenario/system_runner.json
fi