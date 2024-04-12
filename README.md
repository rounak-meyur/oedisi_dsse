# OEDISI DSSE
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE)

## Docker
To simply run the DSSE algorithm for a specific scenario modify the build arg *SCENARIO* to one of the preconfigured settings: ieee123, small, medium or large. Outputs are saved in the outputs directory within a subdirectory named by the scenario name.

```shell
bash run.sh ieee123
```

## Build and Run
Change Line 82 of `Dockerfile` to copy the directory with the OpenDSS files to the container.
```shell
COPY < path to local directory > /home/< name of directory >
```
Inside the `system.json` file, change the `existing_feeder_file` to point to the OpenDSS `Master.dss` file location. 

Then build the Docker container using
```shell
docker build -t opendss-test:0.0.1 . -f Dockerfile
```

Once it is built, run the container
```shell
docker run -it --rm --entrypoint bash opendss-test:0.0.1
```
