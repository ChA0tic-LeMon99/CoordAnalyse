# == README for CoordAnalyse.py == #

                                                                               # == Program Description == # 
                                                                                
# CoordAnalyse analyses text file that contains atom names and atomic coordinates. Input file must folow format <atom name> <x-coordinate> <y-coordinate> <z-coordinate> separated by spaces
# A reference table 'reference.dat' is provided. It has the format <atomic number> <chemical symbol> <element name> <atomic mass> <covalent bond radius> 
# Atomic mass and bond radius from reference table are used to calculate weighted centre of gravity and to determine whether atoms are covalently bonded
# Analysis results will be shown in command line and be written to output file(s)                                                                        
                                                                                
                                                                                   # == How To Use == #
                                                                                    
==To start CoordAnalyse from command line== 
# Type: python <space> CoordAnalyse.py <space> <input filename> <space> reference.dat 
# For example: python CoordAnalyse.py testdata.dat reference.dat
# Comment on version of python and what libraries are required

# When prompted for output file, enter an available output filename
# CoordAnalyse will start the analysis. Results will be written to output.dat file

==When asked 'Anything else?'==
# Enter additional commands for more analysis or to exit program

==To calculate angle between two vectors==
# Type 'ang'
# CoordAnalyse will prompt you for reference atom, first atom, and second atom 
# Type the atom name as shown in input file. Atom name is case sensitive
# For example: C5'
# Results from 'ang' command is not written to output file

==To calcuate distance between two atoms==
# Type 'dis'
# CoordAnalyse will prompt you for first and second atom 
# Results from 'dis' command is not written to output file

==To plot atomic coordinates==
# Type 'plot'
# CoordAnalyse will prompt you for plot filename
==Default colour scheme==
# Carbon - Red
# Nitrogen - Green
# Oxygen - Blue
# Hetatom - Magenta
==Default viewing angle==
# Elevation angle: 30
# Azimuthal angle: 70
# The plot will be saved as a .png image 

==To exit program==
# Type 'exit'

                                                                                    # == Output files == #

# A .dat text file will always be created. It records:
# Number of atoms
# Number of atoms for each element
# Maximum and minimum extent of molecule
# Weighted and unweighted centre of gravity
# Which atoms are covalently bonded and the distance between them

# If 'plot' command is used, a .png image file will be created
# It is a 3 dimensional plot that shows atomic positions in x-, y-, and z- direction


