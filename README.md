# OEDISI DSSE
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE)

## Docker
To simply run the DSSE algorithm for a specific scenario modify the build arg *SCENARIO* to one of the preconfigured settings: ieee123, small, medium or large. Outputs are saved in the outputs directory within a subdirectory named by the scenario name.

```shell
bash build.sh ieee123
```

## Build and Run
The DSSE is run as a C++ executable. The C++ code has been cloned from the [GridAppsD repository](#https://github.com/GRIDAPPSD/gridappsd-state-estimator/tree/OEDISI.1.1). To modify the C++ code, you can edit the source code in the /build/gridappsd-state-estimator/state-estimator/ directory and compile using
```shell
cd /build/gridappsd-state-estimator/
make -C state-estimator
```
