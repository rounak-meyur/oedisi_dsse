#!/bin/bash
if [ ! -d "outputs/custom" ]; then
	echo "Creating outputs directory"
	mkdir -p outputs/custom
fi
pkill -9 helics_broker
pkill -9 python
oedisi build --component-dict scenario/components.json --system scenario/system.json --target-directory build
oedisi run --runner build/system_runner.json