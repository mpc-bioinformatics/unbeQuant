# Alternative LFQ Workflow

This workflow quantifies identified as well as unidentified features! Furthermore, you can run this workflow 
with multiple FDRs and also retreive the identification results.

# Depenencies

Execute `compile_and_setup_dependencies.sh` to set up all needed dependencies. make sure you are in a Python-Environment,
to prevent errors or dependency conflicts.

Also have `java` and `mono` installed. Java is needed nextflow and Mono isused for the ThermoRawFileParser.

## Executing via docker

Create docker:

> docker build -t unbeqonet:local . -f docker/Dockerfile

Then simply execute this workflow by additionally adding `-with-docker` into the cli


## Executing and setting up for nf-cloud

1. Run setup file (generating needed docker images)

2. Copy xxx as a workflow template (which gets you an initial runnable workflow in nf-cloud)

TBD

# How to run

Check out the main workflow. There all required parameters are listed. Each sub workflow contains additional parameters which may be configured.


