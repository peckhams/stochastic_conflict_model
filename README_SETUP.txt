S.D. Peckham
September 8, 2021

===============================
 Container Setup Instructions
===============================
Choose "Ubuntu" as the base image.

## Let "$" denote "/home/clouseau$".

$ git clone https://github.com/peckhams/stochastic_conflict_model
Cloning into 'stochastic_conflict_model'...
remote: Enumerating objects: 80, done.
remote: Counting objects: 100% (80/80), done.
remote: Compressing objects: 100% (72/72), done.
remote: Total 80 (delta 29), reused 0 (delta 0), pack-reused 0
Unpacking objects: 100% (80/80), 348.65 KiB | 7.42 MiB/s, done.
$

# Install miniconda for Linux
# See:  https://gist.github.com/arose13/fcc1d2d5ad67503ba9842ea64f6bac35
$  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
$  bash ~/miniconda.sh -b -p ~/miniconda 
$  rm ~/miniconda.sh
$  export PATH=~/miniconda/bin:$PATH

$ which conda
/home/clouseau/miniconda/bin/conda

# Create directory for model output files
$ cd
$ mkdir Conflict
$ cd Conflict
$ mkdir Output
$ cd

# Install required Python packages with conda
$ conda update -n base conda
$ conda install numpy
$ conda install pandas
$ conda install matplotlib
$ conda install gdal

# Install the model package using its setup.py
$ cd stochastic_conflict_model
$ pip install -e .

# Run the model with default CFG file
$ python conflict --cfg_file './input_files/conflict.cfg'

# Run the model with another CFG file
$ python conflict --cfg_file './input_files/conflict_Upopcount1.cfg'

# Run the model with another CFG file
$ python conflict --cfg_file './input_files/conflict_Upopcount2.cfg'
  
# Output files are in:  /home/clouseau/Conflict/Output/
# Output files are spatial grid stacks, currently in the simple binary RTS format
    with RTI header file.  (RTS = RiverTools Sequence; RTI = RiverTools information)
# NetCDF format for output files is coming soon.