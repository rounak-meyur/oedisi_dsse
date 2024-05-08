# OEDISI DSSE
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE) on the IEEE 123 bus test feeder and several Smart-DS test feeders.

## DOCKER - Build Standard Scenarios
To build a container with the IEEE 123 bus and the SMART-DS scenarios, use the following command. To simply run the DSSE algorithm for a specific scenario modify the argument `scenario` to one of the following preconfigured settings: 
- ieee123 
- small 
- medium 
- large (requires more than 24GB of memory)
```shell
bash create.sh < scenario >
```
This will build the container `openenergydatainitiative/pnnl-dsse-ekf` with the the scenario name as the tag `< scenario >`. Also, it will create an outputs directory to save the co-simulation outputs. Finally, it will run the container which automatically starts the co-simulation and saves outputs.

## DOCKER - Load and Run Containers from OpenEDI
The containers can be loaded directly from DockerHub for the four pre-defined scenarios as follows:
- [ieee123](https://hub.docker.com/r/openenergydatainitiative/pnnl-dsse-ekf-ieee-123)
- [small](https://hub.docker.com/r/openenergydatainitiative/pnnl-dsse-ekf-sfo-p1u)
- [medium](https://hub.docker.com/r/openenergydatainitiative/pnnl-dsse-ekf-sfo-p6u)
- [large](https://hub.docker.com/r/openenergydatainitiative/pnnl-dsse-ekf-sfo-p9u)

### [IEEE 123 bus feeder](https://hub.docker.com/r/openenergydatainitiative/pnnl-dsse-ekf-ieee-123)
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE) on the IEEE 123 bus test feeder.

#### Load Container
```shell
docker pull openenergydatainitiative/pnnl-dsse-ekf-ieee-123:0.0.0
```
#### Run Container
First, create an `outputs` folder in the local directory to save the outputs for the DSSE scenario. Then mount this volume while running the Docker container to save the outputs directly in this folder.
```shell
mkdir -p "outputs/ieee123"
docker volume create --name oedisi_outputs --opt type=none --opt device=${PWD}/outputs/ieee123 --opt o=bind
docker run --rm -it --mount source=oedisi_outputs,target=/home/outputs openenergydatainitiative/pnnl-dsse-ekf-ieee-123:0.0.0
```

### [Smart-DS Small sized feeder](https://hub.docker.com/r/openenergydatainitiative/pnnl-dsse-ekf-sfo-p1u)
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE) on the small Smart-DS feeder.

#### Load Container
```shell
docker pull openenergydatainitiative/pnnl-dsse-ekf-sfo-p1u:0.0.0
```
#### Run Container
First, create an `outputs` folder in the local directory to save the outputs for the DSSE scenario. Then mount this volume while running the Docker container to save the outputs directly in this folder.
```shell
mkdir -p "outputs/small"
docker volume create --name oedisi_outputs --opt type=none --opt device=${PWD}/outputs/small --opt o=bind
docker run --rm -it --mount source=oedisi_outputs,target=/home/outputs openenergydatainitiative/pnnl-dsse-ekf-sfo-p1u:0.0.0
```

### [Smart-DS Medium sized feeder](https://hub.docker.com/r/openenergydatainitiative/pnnl-dsse-ekf-sfo-p6u)
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE) on the medium Smart-DS feeder.

#### Load Container
```shell
docker pull openenergydatainitiative/pnnl-dsse-ekf-sfo-p6u:0.0.0
```
#### Run Container
First, create an `outputs` folder in the local directory to save the outputs for the DSSE scenario. Then mount this volume while running the Docker container to save the outputs directly in this folder.
```shell
mkdir -p "outputs/medium"
docker volume create --name oedisi_outputs --opt type=none --opt device=${PWD}/outputs/medium --opt o=bind
docker run --rm -it --mount source=oedisi_outputs,target=/home/outputs openenergydatainitiative/pnnl-dsse-ekf-sfo-p6u:0.0.0
```

### [Smart-DS Large sized feeder](https://hub.docker.com/r/openenergydatainitiative/pnnl-dsse-ekf-sfo-p9u)
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE) on the large Smart-DS feeder.

#### Load Container
```shell
docker pull openenergydatainitiative/pnnl-dsse-ekf-sfo-p9u:0.0.0
```
#### Run Container
First, create an outputs folder in the local directory to save the outputs for the DSSE scenario. Then mount this volume while running the Docker container to save the outputs directly in this folder.
```shell
mkdir -p "outputs/large"
docker volume create --name oedisi_outputs --opt type=none --opt device=${PWD}/outputs/large --opt o=bind
docker run --rm -it --mount source=oedisi_outputs,target=/home/outputs openenergydatainitiative/pnnl-dsse-ekf-sfo-p9u:0.0.0
```
**NOTE**: To run the `large` scenario, make sure that the memory limit of Docker container is increased to at least 24GB. This can be changed in Docker Desktop application. Go to Settings -> Resources -> Resource Allocation -> Memory limit and increase the slider to meet the minimum memory requirement of 24GB.

## DOCKER - Custom Scenario
To build a container with a custom scenario, find the folder in the local system which contains the OpenDSS file for the custom scenario. Then, change **line 79** of `Dockerfile.custom` to copy the directory with the OpenDSS files to the container home directory.
```shell
COPY < path to local directory > /home/custom
```
Inside the `scenario/custom/system.json` file, change the configuration parameter `existing_feeder_file` on **line 44** to point to the OpenDSS `Master.dss` file location. Then build and run the container using the `custom.sh`
```shell
bash custom.sh
```
This will build the container `pnnl-dsse-ekf` with the the scenario name as the tag `custom`. Also, it will create an outputs directory to save the co-simulation outputs. Finally, it will run the container which automatically starts the co-simulation and saves outputs.