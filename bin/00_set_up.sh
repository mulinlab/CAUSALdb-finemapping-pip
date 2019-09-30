#!/bin/bash

# conda env create -f finemapcondaenv.yml

# PAINTOR v3.0
wget -c https://github.com/gkichaev/PAINTOR_V3.0/archive/3.0.zip
unzip 3.0.zip
cd PAINTOR_V3.0-3.0
bash install.sh
cd ..
rm 3.0.zip

# CAVIARBF
wget -c https://bitbucket.org/Wenan/caviarbf/get/7e428645be5e.zip
unzip 7e428645be5e.zip
cd Wenan-caviarbf-7e428645be5e
make
cd ..
rm 7e428645be5e.zip

# FINEMAP v1.3.1
wget -c http://www.christianbenner.com/finemap_v1.3.1_x86_64.tgz
tar -zxvf finemap_v1.3.1_x86_64.tgz
rm finemap_v1.3.1_x86_64.tgz