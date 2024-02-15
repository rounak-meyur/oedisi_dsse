# OEDISI DOPF
Open Energy Data Initiative - Solar Systems Integration Data and Analytics (OEDI-SI) Distribution System State Estimation (DSSE)

## Docker

To simply run the DOPF algorithm for a specific scenario modify the build arg *SCENARIO* to one of the preconfigured settings: small, medium, large, or ieee123. ieee123 is the default scenario. Outputs are saved in the mounted volume to your local directory.

```shell
    docker build -t oedisi-dsse:1.1.0 . -f Dockerfile
    docker run --rm -it --entrypoint bash oedisi-dsse:1.1.0
```
#
 
## Setup

```shell
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Build and Run
Replace the \<scenario\> below to point to the desired scenario folder name

```shell
./run.sh <scenario>
```
