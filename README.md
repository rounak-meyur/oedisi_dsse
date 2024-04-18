# OEDISI DSSE
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE)

## DOCKER - Standard Scenarios
To build a container with the IEEE 123 bus and the SMART-DS scenarios, use the following command. 
```shell
docker build -t pnnl-dsse-ekf:base . -f Dockerfile
```
This will build the container `pnnl-dsse-ekf` with the tag `base`. To run the container, run the following command.
```shell
docker run -it --rm --entrypoint bash -w "/home/" pnnl-dsse-ekf:base
```
Once inside the container, to simply run the DSSE algorithm for a specific scenario modify the argument `scenario` to one of the following preconfigured settings: 
- ieee123 
- small 
- medium 
- large (requires more than 24GB of memory)
Outputs are saved inside the outputs directory within a subdirectory named by the scenario name. Replace the argument `scenario` below to execute state estimation for the desired scenario.
```shell
bash run.sh < scenario >
```
**NOTE**: To run the `large` scenario, make sure that the memory limit of Docker container is increased to at least 24GB. This can be changed in Docker Desktop application. Go to Settings -> Resources -> Resource Allocation -> Memory limit and increase the slider to meet the minimum memory requirement of 24GB. 

## DOCKER - Custom Scenario
To build a container with a custom scenario, find the folder in the local system which contains the OpenDSS file for the custom scenario. Then, change **line 79** of `Dockerfile.custom` to copy the directory with the OpenDSS files to the container.
```shell
COPY < path to local directory > /home/< name of directory >
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
bash run.sh < scenario >
```
