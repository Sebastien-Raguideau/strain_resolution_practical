#!/bin/bash

echo "install repos/tools needed for workshop"

###### seb ll -h need accomodation ######
sed -i 's/ls -alF/ls -lhaF/g' .bashrc 

###### Install miniconda ######
export REPOS=$HOME"/repos"
mkdir -p $REPOS
cd $REPOS
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh

/bin/bash Miniconda3-py38_4.12.0-Linux-x86_64.sh -b -p $REPOS/miniconda3
/home/ubuntu/repos/miniconda3/condabin/conda init
/home/ubuntu/repos/miniconda3/condabin/conda config --set auto_activate_base false

__conda_setup="$('/home/ubuntu/repos/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/ubuntu/repos/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/ubuntu/repos/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/ubuntu/repos/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

###### Install STRONG ######
cd $REPOS

# requirement
sudo apt-get update
sudo apt-get -y install libbz2-dev libreadline-dev cmake g++ zlib1g zlib1g-dev

# clone
git clone --recurse-submodules https://github.com/chrisquince/STRONG.git
cd STRONG
git submodule foreach git pull origin master

# real install
./install_STRONG.sh 
ln -s $REPOS/STRONG/bin/STRONG $REPOS/miniconda3/envs/STRONG/bin/

###### install env for Intro ######
cd $REPOS
wget https://raw.githubusercontent.com/Sebastien-Raguideau/strain_resolution_practical/main/conda_env_Intro.yaml
mamba env create -f conda_env_Intro.yaml
rm conda_env_Intro.yaml


###### Add databases ######
export DATABASE=$HOME"/Databases"
mkdir -p $DATABASE
cd $DATABASE

# rpsblast
wget https://strongtest.s3.climb.ac.uk/rpsblast_cog_db.tar.gz
tar -xvzf rpsblast_cog_db.tar.gz
rm rpsblast_cog_db.tar.gz

# Kraken 
mkdir -p $REPOS/MiniKraken
cd $REPOS/MiniKraken
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20220926.tar.gz
tar -xvzf k2_standard_16gb_20220926.tar.gz

###### Install Bandage ######
wget https://github.com/rrwick/Bandage/releases/download/v0.9.0/Bandage_Ubuntu-x86-64_v0.9.0_AppImage.zip
unzip Bandage_Ubuntu-x86-64_v0.9.0_AppImage.zip
ln -s Bandage_Ubuntu-x86-64_v0.9.0_AppImage /home/ubuntu/repos/miniconda3/bin/Bandage

###### add silly jpg ######
cd 
wget https://raw.githubusercontent.com/Sebastien-Raguideau/strain_resolution_practical/main/Figures/image_you_want_to_copy.jpg
wget https://raw.githubusercontent.com/Sebastien-Raguideau/strain_resolution_practical/main/Figures/image_you_want_to_display.jpg

###### install feh ######
sudo apt -y install feh
sudo apt -y install evince

###### remove script ######
#rm /home/ubuntu/install_vm.sh
