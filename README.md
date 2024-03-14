# OEDISI DSSE
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE)

## Docker
To simply run the DOPF algorithm for a specific scenario modify the build arg *SCENARIO* to one of the preconfigured settings: ieee123, small, medium or large. ieee123 is the default scenario. Outputs are saved in the outputs directory.

```shell
docker build --build-arg SCENARIO=ieee123 -t dsse-ieee123:1.1.0 . -f Dockerfile
docker volume create --name oedisi_outputs --opt type=none --opt device=${PWD}/outputs --opt o=bind
docker run --rm -it --mount source=oedisi_outputs,target=/home/outputs dsse-ieee123:1.1.0
```

## Build and Run
The DSSE is run as a C++ executable. The C++ code has been cloned from the [GridAppsD repository](#https://github.com/GRIDAPPSD/gridappsd-state-estimator/tree/OEDISI.1.1). To modify the C++ code, you can edit the source code in the /build/gridappsd-state-estimator/state-estimator/ directory and compile using
```shell
cd /build/gridappsd-state-estimator/
make -C state-estimator
```

Replace the \<scenario\> below to point to the desired scenario folder name

```shell
cd /home/
./run.sh <scenario>
```
