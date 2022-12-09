#!/bin/bash

# Build ProtGraphTraverseInt
rm -rf bin/ProtGraphTraverseInt/build
rm -rf bin/protgraphtaverseint

cmake -B bin/ProtGraphTraverseInt/build -S bin/ProtGraphTraverseInt 
# Some machines (like older ubuntus) need to specify the new gcc compiler via:  "-D CMAKE_C_COMPILER=gcc-11 -D CMAKE_CXX_COMPILER=g++-11"
# If you want to limit the variants by X, then please modify the source-code
OLDDIR=$(pwd)
cd bin/ProtGraphTraverseInt/
cmake --build ./build
cp build/protgraphtraverseint ../protgraphtraverseint
cd $OLDDIR



# Build ProtGraphTraverseFloat
rm -rf bin/ProtGraphTraverseFloat/build
rm -rf bin/protgraphtraversefloat

cmake -B bin/ProtGraphTraverseFloat/build -S bin/ProtGraphTraverseFloat 
# Some machines (like older ubuntus) need to specify the new gcc compiler via:  "-D CMAKE_C_COMPILER=gcc-11 -D CMAKE_CXX_COMPILER=g++-11"
# If you want to limit the variants by X, then please modify the source-code

OLDDIR=$(pwd)
cd bin/ProtGraphTraverseFloat/
cmake --build ./build
cp build/protgraphtraversefloat ../protgraphtraversefloat
cd $OLDDIR


# Install ProtGraph
pipenv --rm 
# Install for ProtGraph global_export
pipenv run pip install protgraph==0.3.8
pipenv run pip uninstall -y igraph
pipenv run pip install igraph==0.9.11
pipenv run pip install pip install apsw==3.36.0.post1
pipenv run pip uninstall -y biopython
pipenv run pip install git+https://github.com/biopython/biopython.git@947868c487a12799d51173c5f651a44ecb3fb6fa
pipenv run pip install XlsxWriter==3.0.3
