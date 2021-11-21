#!/usr/bin/env python
import sys # Get sys package for system functions
import numpy as np # numpy for arrays
import re # re for regular expression processing
import os.path # os.path to check file existence
import math # math for cosine functions
import matplotlib # matplotlib to plot coordinates
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d # mplot3d to plot 3 dimensional plot
                                                                      
                                                                                    ## Check that command-line arguments are correct                                                                 
if (len(sys.argv) != 3): # Check that there are three command line arguments 
    print ("Invalid input\nEnter 'program name' 'input filename' 'reference filename'") # Print error message when command line arguments less than 3 
    sys.exit (-1) # Exit program
    
if (os.path.isfile(sys.argv[1]) == False): # Checks if input file exists 
    print ("Invalid input filename\nEnter 'program name' 'input filename' 'reference filename'") # Print error message if input file does not exist     
    sys.exit (-1) # Exit program

if (os.path.isfile(sys.argv[2]) == False): # Checks if reference table file exists 
    print ("Invalid reference table filename\nEnter 'program name' 'input filename' 'reference filename'") # Print error message if periodic table file does not exist   
    sys.exit (-1) # Exit program
    
                                                                                 ## Read input files and create output file
                                                                                 
input_name = sys.argv[1] # Define input_name as input file
perio_table = sys.argv[2] # Define _name as referenece table file

def prompt(): # Define a function prompt() that stores output filename
    output_name = input("Enter name for output file\n") # Prompt user for output filename and assigned it to variable output_name
    try:
        f = open("%s.dat" % output_name, "x") # Check if output file already exists        
    except:
        print("Error: Filename already in use")
        prompt() # Use recursion to keep calling prompt() until a valid output filename is entered
    return output_name  # Return the value of valid output filename
    
output_name = prompt()

print ("\nWelcome to CoordAnalyse!\n") # Print welcome message   

data_items = [] # Create list to store every line in input file as a list of words
with open (input_name, "r") as fi: # Read each line from input file
    for line in fi:
        words = line.split() # Splits the current line into a list where each word is a list item
        data_items.append(words) # Add each list of words into the end of data_items list
        
                                                                         ## Global parameters: atom number, coordinate number, atom class

num_atoms = len (data_items) # Assign number of lines in input file to variable num_atoms
coord_num = 3 # Number of atomic coordinates 
atom_list = [] # Create a list to storing objects that belong to class 'atom2_class'
coords = np.zeros(coord_num) # Create an array for atomic coordinates
table_elements = 112 # Number of elements in reference table

class atom2_class: # Create a class to store atom properties: atom name, element, coordinates, atomis mass, and covalent bond radii
    def __init__(self, name, ele2, coords, mass, radii):
        self.name = name
        self.ele2 = ele2
        self.coords = coords
        self.mass = mass
        self.radii = radii                                
        
                                                     ## Fill in the value for property in class atom2_class and check input data for errors
                                                     
error_log = False # Create a variable to record whether an error exists in input file
row = 0  # Create a variable to count row number   
                       
for i in data_items: # Iterate through each list of words in data_items list 
    row = row + 1 # Variable 'row' stores the current row number
    name = i[0] # Assign the first word to variable 'name' 
    ele = re.findall("^[A-Z]{1}[a-z]*", name) # Use regular expression to find alphabets in values stored in variable 'name', then assign alphabets to variable 'ele' 
    ele2 = ''.join(ele) # Converts ele from list item to string
        
    for j in range(coord_num): # Iterate through every atomic coordinate
        try: # Try to assign coordinates in input file to the 'coords' array
            coords[j] = i[j + 1] 
        except: # When an exception is encountered
            print("Invalid coordinate(s) in Row:", row, "\n", name, i[1:j + 2], "\n") # Print invalid coordinate in current row   
            error_log = True # Log error into 'error_log'       
                          
    with open (perio_table, "r") as fi: # Read the lines from reference table file
        count = 0 # Initialise row number
        for line in fi: # Iterate through every line in reference table
            words = line.split()
            if (ele2 == words[1]): # If input file element name matches reference table element
                mass = words[3] # Assign atomic mass value in reference table to 'mass'
                radii = words[4] # Assign bond radii value in reference table to 'radii'
                break # Exit loop when an element match is found
            elif (count == table_elements - 1 and ele2 != words[1]): # If input file element name does not match any reference table element
                print("Invalid element in Row:", row, "\n", name, coords, "\n") # Print invalid element in current row
                error_log = True # Logs error into error_log      
            else:
                count = count + 1 # Count row number      
                                                   
    atm_list = atom2_class(name, ele2, np.copy(coords), mass, radii) # Transfer the value for each 'atom2_class' property to object 'atm_list'
    atom_list.append(atm_list) # Insert object 'atm_list' to the end of 'atom_list' list
    
if(error_log == True): # If error exists in input file, print error message and exit program 
    print("\nPlease check input file for error(s) and try again\n")
    sys.exit(-1)                
        
                                                                            ## Print out number of atoms and write to output file
                                                                            
print ("There are ", len (data_items), " atoms in ", input_name) # Print number of atoms
f = open("%s.dat" % output_name, "a") # Open output file and write number of atoms to the end of file
f.write("There are {} atoms in {}\n\n".format(len (data_items), input_name))

                                                                            ## Count atoms that are the same element      
                                                                         
counted_arr = np.zeros(num_atoms) # Create an array of size number of atoms to record which atom has already been counted
counted_arr[0] = False # Set the first element to false, which means the first atom has not been counted  

print ("\nNumber of each element")
f = open("%s.dat" % output_name, "a")
f.write("\nNumber of element\n")

for i in range (num_atoms): # Iterate through every atom
    if(counted_arr[i] == False): # If the current atom has not been counted before
        count = 1 # Initialise count
        for j in range (i + 1, num_atoms - 1): # Iterate through every atom after atom'i'
            if (atom_list[i].ele2 == atom_list[j].ele2): # If atom'j' element matches atom'i' element
                count = count + 1 # Count the number of atom'j' that are the same element as atom'i'
                counted_arr[j] = True # Set 'j' element to true, which means atom'j' has been counted
            else: # If atom'j' element does not match atom'i' element
                count = count # Count remains unchanged
                if (counted_arr[j] == True): # If atom'j' has been counted before,  atom'j' remains counted
                    counted_arr[j] = True
                else:
                    counted_arr[j] = False # If atom'j' has not been counted before, atom'j' remains uncounted
        print (atom_list[i].ele2, ":", count) # Print and write to output file the number of atoms for each element
        f = open("%s.dat" % output_name, "a")
        f.write("{} : {}\n".format(atom_list[i].ele2, count))      
                                      
                                                                                                                
                                                                          ## Calculate maximun, minimum, and simple centre of gravity                                                                           
maxcoord = np.copy(atom_list[0].coords) # Copy coordinates of first atom to maxcoord and mincoord
mincoord = np.copy(atom_list[0].coords)
    
for i in range(num_atoms): # Iterate through every atom  
    for j in range (coord_num): # Iterate through every atomic coordinate
        if (atom_list[i].coords[j] > maxcoord[j]): # If current coordinate is larger than current maxcoord, copy current coordinate to maxcoord
            maxcoord[j] = np.copy(atom_list[i].coords[j]) 
        if (atom_list[i].coords[j] < mincoord[j]): # If current coordinate is smaller than current mincoord, copy current coordinate to mincoord
            mincoord[j] = np.copy(atom_list[i].coords[j]) 
                    
sumcoord = np.zeros(coord_num) # Array to store sum of each atomic coordinate
for i in range(num_atoms):              
    sumcoord = sumcoord + atom_list[i].coords # Update the sum of coordinates  
cog = sumcoord/num_atoms # Divide sum of coordinates by number of atoms to calculate centre of gravity  
    
print ("\nThe maximum extent of the molecule is ", maxcoord[0], ",", maxcoord[1], ",", maxcoord[2]) # Print and write to output maximum coordinates
f = open("%s.dat" % output_name, "a")
f.write("\nThe maximum extent of the molecule is {}, {}, {}\n".format(maxcoord[0], maxcoord[1], maxcoord[2]))

print ("The minimum extent of the molecule is ", mincoord[0], ",", mincoord[1], ",", mincoord[2], "\n") # Print and write to output minimum coordinates
f = open("%s.dat" % output_name, "a")
f.write("The minimum extent of the molecule is {}, {}, {}\n\n".format(mincoord[0], mincoord[1], mincoord[2]))
 
print ("The simple centre of gravity is ", str(round(cog[0], 3)), ",", str(round(cog[1], 3)), ",", str(round(cog[2], 3))) # Print and write to output cog rounded to 3 decimal places
f = open("%s.dat" % output_name, "a")
f.write("The simple centre of gravity is {}, {}, {}\n".format(str(round(cog[0], 3)), str(round(cog[1], 3)), str(round(cog[2], 3))))

                                                                                      ## Calculate weighted centre of gravity

sum_mass = 0 # Initialise the sum of atomic mass
sum_coord = np.zeros(coord_num) # Create an array to store sum of each atomic coordinate
for i in range (num_atoms): # Iterate through every atom
    sum_mass = sum_mass + float(atom_list[i].mass) # Add current atomic mass to 'sum_mass'
    weighted_coord = np.multiply(atom_list[i].coords, float(atom_list[i].mass)) # Multiply atomic coordinates by atomic mass
    sum_coord = sum_coord + weighted_coord # Add weighted coordinate to sum_coord
    
weight_cog = sum_coord / sum_mass # Divide sum of weighted coordinates by sum of atomic mass  
  
print ("The weighted centre of gravity is ", str(round(weight_cog[0], 3)), ",", str(round(weight_cog[1], 3)), ",", str(round(weight_cog[2], 3)), "\n") 
# Print and write to output weighted cog rounded to 3 decimal places
f = open("%s.dat" % output_name, "a")
f.write("The weighted centre of gravity is {}, {}, {}\n\n".format(str(round(weight_cog[0], 3)), str(round(weight_cog[1], 3)), str(round(weight_cog[2], 3))))
  
                                                                                ## Calculate covalent bond using 1.8 cutoff   
                                                                                
print ("Covalent bond using 1.8 cutoff")   
f = open("%s.dat" % output_name, "a")
f.write("Covalent bond using 1.8 cutoff\n")                                                            
for i in range (num_atoms):
    for j in range (i + 1, num_atoms - 1):
        distance = np.linalg.norm(atom_list[i].coords - atom_list[j].coords) # Calculate distance between two atoms
        if (distance < 1.8):
            print ("Atom ", atom_list[i].name, " is bonded to ", atom_list[j].name, "distance", str(round(distance, 3))) # print atoms bonded to adjacent atom
            f = open("%s.dat" % output_name, "a")
            f.write("Atom {} is bonded to {} distance {}\n".format(atom_list[i].name, atom_list[j].name, str(round(distance, 3))))
            
                                                                               ## Calculate covalent bond using reference table
                                                                               
print ("\nCovalent bond using reference")
f = open("%s.dat" % output_name, "a")
f.write("\nCovalent bond using reference\n")                                                        
for i in range (num_atoms): # Iterate through every atom
    for j in range (i + 1, num_atoms - 1): # Iterate through every atom after atom'i'
        covalent_radii = float(atom_list[i].radii) + float(atom_list[j].radii) # Covalent bond length is the total covalent bond radii of atom'i' and atom'j'
        up_bound = covalent_radii + 0.1  # Set the allowed degree of deviation from ideal bond length 
        low_bound = covalent_radii - 0.1
        distance = np.linalg.norm(atom_list[i].coords - atom_list[j].coords) # Use linear algebra to calculate distance between two atoms
        if (distance < up_bound): # If atomic distance is less than the upper bound of covalent bond length
            print ("Atom ", atom_list[i].name, " is bonded to ", atom_list[j].name, "distance", str(round(distance, 3))) # Print and write to output atom'i' is bonded to atom'j'
            f = open("%s.dat" % output_name, "a")        
            f.write("Atom {} is bonded to {} distance {}\n".format(atom_list[i].name, atom_list[j].name, str(round(distance, 3))))         
                                                                                                                                                                    
                                                                                                          
                                                                             ## Calculate angle between two user-selected bonds                         

def calc_angle2(reference_at, coord_1_at, coord_2_at): # Define a function for calculating angle which accepts three atom names as arguments
    for i in range (num_atoms): # Iterate through every element in 'atom_list' to find object with name that matches atom name from user input 
        if (atom_list[i].name == reference_at):
            reference = atom_list[i].coords  
        elif (atom_list[i].name == coord_1_at):
            coord_1 = atom_list[i].coords
        elif (atom_list[i].name == coord_2_at):
            coord_2 = atom_list[i].coords    
                                                                                   
    vector_1 = coord_1 - reference # Calculate vector 1 and vector 2 by subtracting reference coordinates from atom 1 and atom 2 coordinates
    vector_2 = coord_2 - reference  
    mod_1 = 0 # Initialise magnitudes and dot product
    mod_2 = 0
    dot_product = 0
    
    for i in range (coord_num): # Iterate through every atomic coordinate
        mod_1 = mod_1 + (vector_1[i] ** 2) # Calculate and add coordinate squared
        mod_2 = mod_2 + (vector_2[i] ** 2)
        dot_product = dot_product + (vector_1[i] * vector_2[i]) # Update dot product
        if (i == coord_num - 1): # When the last atomic coordinate is reached
            mod_1 = mod_1 ** 0.5  # Calculate vector magnitude
            mod_2 = mod_2 ** 0.5
            mod_1_2 = mod_1 * mod_2 # Multiply vector 1 and vector 2 magnitudes
            angle = (180 / math.pi) * math.acos(dot_product / mod_1_2) # Calculate inverse cosine using acos(), then converts radian to angle in degree
    print ("\n", reference_at, "--", coord_1_at, "and", reference_at, "--", coord_2_at, "Angle", ":\n", str(round(angle, 2)), "degree\n") # Print the angle between the vectors

                                                           
                                                                              ## Calculate distances between two user-selected atoms
                                                                                   
def calc_distance(coord_1_at, coord_2_at): # Define a function for calculating distance which accepts two atom names as arguments
    for i in range (num_atoms): # Iterate through every object in 'atom_list' to find object with name that matches atom name from user input
        if (atom_list[i].name == coord_1_at): 
            coord_1 = atom_list[i].coords
        elif (atom_list[i].name == coord_2_at):
            coord_2 = atom_list[i].coords                                                                           
    distance = np.linalg.norm(coord_1 - coord_2) # Uses linear algebra to calculate distance between two coordinates
    print ("\n", coord_1_at, "and", coord_2_at, "Distance", ":\n", str(round(distance,3)), "Angstrom\n") # Print atomic distance
    
                                                                                 ## Count C, N, O

def count_ele(element): # Define count_ele() function which takes element name as argument
    count = 0 # Initialise count
    for i in range(num_atoms): # Find number of atoms belonging to the same element as the query element
        if (atom_list[i].ele2 == element):
            count = count + 1
    return count
            

                                                                                     ## Plot atomic coordinates

def plot(): # Define plot() function which plots atomic coordinates
    mol_len = (maxcoord[0] - mincoord[0]) + (maxcoord[1] - mincoord[1])
    mol_height = (maxcoord[2] - mincoord[2])
    plt.rcParams['figure.figsize'] = [(mol_len / 1.8), (mol_height * 1.3)] # Set length and height of plot
    fig = plt.figure("Coord Plot") 
    fig.suptitle(input_name, fontsize=22, fontweight=500) # Set plot title to input filename
    ax = plt.axes(projection='3d') # Plot 3 dimensional axes  
    ax.view_init(30, 70) # Specifies elevation and azimuth angle
    
    ax.set_xlabel('X-axis', fontsize=12, fontweight=500) # Set labels for x-, y-, and z- axis
    ax.set_ylabel('Y-axis', fontsize=12, fontweight=500)
    ax.set_zlabel('Z-axis', fontsize=12, fontweight=500)     
    
    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
  

    C_count = count_ele("C") # Calls count_ele() to calculate number of carbon
    Cx_arr = np.zeros(C_count) # Create array of size number of carbon to store coordinates of all carbon atoms
    Cy_arr = np.zeros(C_count)
    Cz_arr = np.zeros(C_count)
    Carr_count = 0 # Initialise array number
    
    N_count = count_ele("N") # Calls count_ele() to calculate number of nitrogen
    Nx_arr = np.zeros(N_count) # Create array of size number of nitrogen to store coordinates of all nitrogen atoms
    Ny_arr = np.zeros(N_count)
    Nz_arr = np.zeros(N_count)
    Narr_count = 0 # Initialise array number
    
    O_count = count_ele("O") # Calls count_ele() to calculate number of oxygen
    Ox_arr = np.zeros(O_count) # Create array of size number of oxygen to store coordinates of all oxygen atoms
    Oy_arr = np.zeros(O_count)
    Oz_arr = np.zeros(O_count)
    Oarr_count = 0 # Initialise array number
    
    Het_count = num_atoms - C_count - N_count - O_count # Calculate number of hetatom by subtracting total number of C, N, and O from number of atoms
    Hetx_arr = np.zeros(Het_count) # Create array of size number of hetatom to store coordinates of all het atoms
    Hety_arr = np.zeros(Het_count)
    Hetz_arr = np.zeros(Het_count)
    Hetarr_count = 0 # Initialise array number
    
    curr_coords = np.zeros(coord_num) # Create an array to store coordinates of current atom
    for i in range (num_atoms): # Iterate through every atom
        curr_ele = atom_list[i].ele2 # Assign element of current atom to curr_ele
        curr_coords = np.copy(atom_list[i].coords) # Copy current atomic coordinates to curr_coords array
        x = curr_coords[0] 
        y = curr_coords[1]
        z = curr_coords[2]       
         
        ax.text(x, y, z,  '%s' % atom_list[i].name, size=8, fontweight=500, color='k', zorder= 5, position = (x + 0.18, y + 0.18, z + 0.18))  # Label coordinate with atom name
        
        if (curr_ele == "C"): # If current element is carbon, store value of each of its coordinates to the correct index position in the correct carbon coordinate array 
            Cx_arr[Carr_count] = x 
            Cy_arr[Carr_count] = y
            Cz_arr[Carr_count] = z
            Carr_count = Carr_count + 1 # Update the next available array index
                        
        elif (curr_ele == "N"):
            Nx_arr[Narr_count] = x
            Ny_arr[Narr_count] = y
            Nz_arr[Narr_count] = z
            Narr_count = Narr_count + 1            
            
        elif (curr_ele == "O"):
            Ox_arr[Oarr_count] = x
            Oy_arr[Oarr_count] = y
            Oz_arr[Oarr_count] = z
            Oarr_count = Oarr_count + 1 
            
        else: 
            Hetx_arr[Hetarr_count] = x
            Hety_arr[Hetarr_count] = y
            Hetz_arr[Hetarr_count] = z
            Hetarr_count = Hetarr_count + 1 
            
    # Plot the x-, y-, and z- coordinate array of atoms of different element in colour-coded circles, with opacity that decreases with distance to front of plot
    ax.scatter3D(Cx_arr, Cy_arr, Cz_arr, s=150, c="Red", depthshade=True) 
    ax.scatter3D(Nx_arr, Ny_arr, Nz_arr, s=150, c="Green", depthshade=True)
    ax.scatter3D(Ox_arr, Oy_arr, Oz_arr, s=150, c="Blue", depthshade=True)
    ax.scatter3D(Hetx_arr, Hety_arr, Hetz_arr, s=150, c="Magenta", depthshade=True)
            
    for i in range (num_atoms): # Iterate through every atom
        for j in range (i + 1, num_atoms - 1): # Iterate through every atom after atom'i'
            covalent_radii = float(atom_list[i].radii) + float(atom_list[j].radii) # Covalent bond length is the total covalent bond radii of atom'i' and atom'j'
            up_bound = covalent_radii + 0.1  # Set the allowed degree of deviation from ideal bond length 
            low_bound = covalent_radii - 0.1
            distance = np.linalg.norm(atom_list[i].coords - atom_list[j].coords) # Use linear algebra to calculate distance between two atoms
            if (distance < 1.8): # If atomic distance is below 1.8
                x_coord = [atom_list[i].coords[0], atom_list[j].coords[0]]
                y_coord = [atom_list[i].coords[1], atom_list[j].coords[1]] 
                z_coord = [atom_list[i].coords[2], atom_list[j].coords[2]]                    
                ax.plot3D(x_coord, y_coord, z_coord, "gray", linewidth=1.9) # Draw a grey line which represents covalent bond between the two atoms
                
                        
    def picprompt(): # Define a function picprompt() that stores plot filename
        output_name = input("Enter filename for plot output\n") # Prompt user for plot filename and assigned it to variable output_name
        try:
            f = open("%s.png" % output_name, "x") # Check if plot filename already exists        
        except:
            print("Error: Filename already in use")
            prompt() # Use recursion to keep calling picprompt() until a valid filename is entered
        return output_name  # Return the value of valid plot filename    
            
    fig.savefig(picprompt()) # Save plot as a picture with plot filename

                                                                                      ## Custom commands                                                                                

def other(): # Define a function other() that process commands from user
    Other = input("Anything else?\n") # Prompts user for command
    if(Other == "ang"): # Command 'ang' calls calc_angle() function
        reference_at = input("Reference atom:")
        coord_1_at = input("First atom:")
        coord_2_at = input("Second atom:") # Prompts user for atom names
        calc_angle2(reference_at, coord_1_at, coord_2_at)
        other() # Use recursion to call other() until 'exit' is encountered
        
    elif(Other == "dis"): # Command 'dis' calls calc_distance() function
        coord_1_at = input("First atom:")
        coord_2_at = input("Second atom:") # Prompts user for atom names
        calc_distance(coord_1_at, coord_2_at)
        other() # Use recursion to call other() until 'exit' is encountered
        
    elif(Other == "plot"): # Command 'plot' calls plot() function
        plot()
        other()
    
    elif(Other == "exit"): # Command 'exit' quits the program
        sys.exit(0)        
           
    else:
        print("\nSorry, invalid command. Try again") # Print error messgae if command is none of the above
        other() # Use recursion to call other() until a valid command is entered 
                
other()