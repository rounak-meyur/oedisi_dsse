# OEDISI DSSE
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE) on the IEEE 123 bus test feeder and several Smart-DS test feeders.

## DOCKER - Standard Scenarios
To build a container with the IEEE 123 bus and the SMART-DS scenarios, use the following command. To simply run the DSSE algorithm for a specific scenario modify the argument `scenario` to one of the following preconfigured settings: 
- ieee123 
- small 
- medium 
- large (requires more than 24GB of memory)
```shell
bash build.sh < scenario >
```
This will build the container `openenergydatainitiative/pnnl-dsse-ekf` with the the scenario name as the tag `< scenario >`. Also, it will create an outputs directory to save the co-simulation outputs. Finally, it will run the container which automatically starts the co-simulation and saves outputs.

**NOTE**: To run the `large` scenario, make sure that the memory limit of Docker container is increased to at least 24GB. This can be changed in Docker Desktop application. Go to Settings -> Resources -> Resource Allocation -> Memory limit and increase the slider to meet the minimum memory requirement of 24GB. 

### Load Container from OpenEDI
```shell
docker pull openenergydatainitiative/pnnl-dsse-ekf:< scenario >
```
### Run Container
First, create an `outputs` folder in the local directory to save the outputs for the DSSE scenario. Then mount this volume while running the Docker container to save the outputs directly in this folder.
```shell
mkdir -p outputs/< scenario >
docker volume create --name oedisi_outputs --opt type=none --opt device=${PWD}/outputs/< scenario > --opt o=bind
docker run --rm -it --mount source=oedisi_outputs,target=/home/outputs openenergydatainitiative/pnnl-dsse-ekf:< scenario >
```

## DOCKER - Custom Scenario
To build a container with a custom scenario, find the folder in the local system which contains the OpenDSS file for the custom scenario. Then, change **line 79** of `Dockerfile.custom` to copy the directory with the OpenDSS files to the container home directory.
```shell
COPY < path to local directory > /home/custom
```
Inside the `scenario/custom/system.json` file, change the configuration parameter `existing_feeder_file` on **line 44** to point to the OpenDSS `Master.dss` file location. Then build the Docker container using the following command. 
```shell
docker build -t pnnl-dsse-ekf:custom . -f Dockerfile.custom
```
This will build the container `pnnl-dsse-ekf` with the tag `custom`. To run the container, run the following command.
```shell
docker run -it --rm --entrypoint bash -w "/home/" pnnl-dsse-ekf:custom
```
Once inside the container, to simply run the DSSE algorithm for a specific scenario modify the argument `scenario` to one of the following preconfigured settings: 
- custom

Outputs are saved inside the outputs directory within a subdirectory named by the scenario name. Replace the argument `scenario` below to execute state estimation for the desired scenario.
```shell
bash run.sh
```
