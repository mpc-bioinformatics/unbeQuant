#!/bin/bash

# Create Dockerfile, needed for every step
docker build -t unbeqonet:local . -f docker/Dockerfile
