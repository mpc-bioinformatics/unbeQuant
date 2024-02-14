# Alternative LFQ Workflow

This workflow quantifies identified as well as unidentified features! Furthermore, you can run this workflow 
with multiple FDRs and also retreive the identification results.

# Depenencies

Execute `compile_and_setup_dependencies.sh` to set up all needed dependencies. make sure you are in a Python-Environment,
to prevent errors or dependency conflicts.

Also have `java` and `mono` installed. Java is needed nextflow and Mono isused for the ThermoRawFileParser.

# How to run

Check out the main workflow. There all required parameters are listed. Each sub workflow contains additional parameters which may be configured.
