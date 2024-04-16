FROM ubuntu:20.04

# ----------------------------------------------------
# INSTALL HELICS
# ----------------------------------------------------

RUN apt-get update \
    && export DEBIAN_FRONTEND="noninteractive" \
    && export TZ="America/Pacific" \
    && apt install -y \
       libboost-dev \
       libzmq5-dev \
       git \
       cmake-curses-gui \
       build-essential

RUN mkdir /build \
    && cd /build \
    && git clone https://github.com/GMLC-TDC/HELICS \
    && cd HELICS \
    && mkdir build \
    && cd build \
    && cmake -DHELICS_BUILD_CXX_SHARED_LIB=True ../ \
    && make \
    && make install

# ----------------------------------------------------
# INSTALL Activemq c++ extensions
# ----------------------------------------------------
RUN apt install -y m4 \
       wget \
       libaprutil1-dev
       
RUN cd /build \
    && wget http://archive.apache.org/dist/activemq/activemq-cpp/3.9.5/activemq-cpp-library-3.9.5-src.tar.gz \
    && tar -xzf activemq-cpp-library-3.9.5-src.tar.gz \
    && cd activemq-cpp-library-3.9.5 \
    && ./configure \
    && make \
    && make install 

RUN apt install -y liblapack-dev \
       libblas-dev \
       libssl-dev

# ----------------------------------------------------
# INSTALL State Estimator
# ----------------------------------------------------
RUN cd /build \
    && git clone --depth 1 --branch OEDISI.1.1 https://github.com/GRIDAPPSD/gridappsd-state-estimator \
    && cd gridappsd-state-estimator \
    && git clone https://github.com/GRIDAPPSD/SuiteSparse \
    && git clone https://github.com/GRIDAPPSD/json \
    && LD_LIBRARY_PATH=/build/gridappsd-state-estimator/SuiteSparse/lib/ make -C SuiteSparse LAPACK=-llapack BLAS=-lblas \
    && make -C state-estimator \
    && rm -rf .git SuiteSparse/.git json/.git 



# ----------------------------------------------------
# INSTALL Python requirements 
# ----------------------------------------------------
RUN apt install -y python3 \
        python3-pip \
        python-is-python3 \
        sudo \
        vim




WORKDIR /home/
COPY feeder_federate /home/feeder_federate
COPY estimator_federate /home/estimator_federate
COPY measuring_federate /home/measuring_federate
COPY recorder_federate /home/recorder_federate
COPY requirements.txt /home/requirements.txt
COPY scenario /home/scenario
RUN mkdir /home/outputs
COPY run.sh /home/run.sh


RUN pip install -r requirements.txt \
    && rm -rf /root/.cache/pip/wheels


RUN for dir in scenario/*/; \
    do \
    SCENARIO=$(echo $dir | cut -f2 -d /); \
    oedisi build \
    --component-dict scenario/${SCENARIO}/components.json \
    --system scenario/${SCENARIO}/system.json \
    --target-directory build_${SCENARIO}; \
    done

WORKDIR /home/