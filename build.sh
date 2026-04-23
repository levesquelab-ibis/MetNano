#!/bin/bash

# Create directory inside the conda environment if needed
mkdir -p $PREFIX/bin/

# Install main launcher (bash script)
cp MetNano $PREFIX/bin/MetNano
chmod +x $PREFIX/bin/MetNano

# Install Python source code
cp -r src $PREFIX/bin/

# Make utility scripts executable
chmod +x $PREFIX/bin/src/utilities/*.sh
