FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get -y upgrade && apt-get install -y build-essential wget curl unzip cmake python3-pip mono-complete python-is-python3 libqt5network5 libqt5sql5 git

COPY . /root/unbeqonet
WORKDIR /root/unbeqonet

RUN ./compile_and_setup_depencies.sh
ENV PATH="${PATH}://root/unbeqonet/bin"
