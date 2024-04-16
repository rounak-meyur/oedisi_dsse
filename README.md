# OEDISI DSSE
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE)

## Docker
To build the container, use the following command. This will build the container `openenergydatainitiative/pnnl-dsse-ekf` with the tag `1.0.0`.
```shell
docker build -t openenergydatainitiative/pnnl-dsse-ekf:1.0.0 . -f Dockerfile
```
To run the container, run the following command.
```shell
docker run -it --rm --entrypoint bash openenergydatainitiative/pnnl-dsse-ekf:1.0.0
```
Once inside the container, to simply run the DSSE algorithm for a specific scenario modify the argument `scenario` to one of the preconfigured settings: ieee123, small, medium or large. Outputs are saved in the outputs directory within a subdirectory named by the scenario name. Replace the `scenario` below to point to the desired scenario.
```shell
bash /home/run.sh < scenario >
```
To run the `large` scenario, make sure that the memory limit of Docker container is increased to at least 24GB. This can be changed in Docker Desktop application. Go to Settings -> Resources -> Resource Allocation -> Memory limit and increase the slider to meet the minimum memory requirement of 24GB. 