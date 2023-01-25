#################
#  Distortion Value Visualizer
#    
#   User Guide:
#        1. create your database of info to analyze. All values used in this 
#           program can be found using VESTA. This file needs to be in .csv format.
#           See structural_deviation.csv as an example. Separator labels are fine
#           as long as they are no longer than one cell (the code will skip them).
#        2. run the program. You'll be prompted to enter your file name. Python
#           needs the .csv at the end too, so don't forget that. Make sure you're
#           running python in the correct directory so it can find the file. If
#           you input a file name python can't find, you'll be given an error
#           message and prompted to try again.
#        3. You'll then be prompted with a menu and asked to choose options to
#           display on the x and y axes. After that, you'll be asked to choose
#           values to display using the color and size of the scatterplot
#           markers. These last two values are optional and the code will run
#           fine without them. You may enter both, either, or none.
#        4. The program will then display your graph with appropriate axis
#           labels as well as labels for the color gradient and size legend
#           if appropriate.
#        5. the program will ask if you want to make another graph. Input y if
#           you do, and you will once again be prompted with the menu as before. 
#           If you are done, input n, and the program will end. 
#   
##################

import csv
import matplotlib.pyplot as plt

GRAPH_MENU = '\n\t\t\t\t\tMENU:\n\t\
                a = structure type\n\t\
                b = average bond length (A)\n\t\
                c = polyhedra volume (A^3)\n\t\
                d = distortion index\n\t\
                e = quadratic elongation\n\t\
                f = bond angle variance (deg^2)\n\t\
                g = effective coodination number (ECoN)\n\t\
                h = cation\n\t\
                i = ionic radius (A)\n\t\
                j = cation charge'

def open_file():
    """
    PURPOSE: prompt user for input file name, open, return file pointer. 
        Display error message if file doesn't exist
    takes no input
    Returns: fp -- file pointer
    """
    continue_str = 'yes' #initialize continue str to start loop
    while continue_str != 'no':
        file_name = input("Input file of data to examine: ")
        try:
            fp = open(file_name, 'r')
            continue_str = 'no'
        except FileNotFoundError: 
            #if the file doesn't work, print error message and loop
            print("Error: file does not exist. Please try again.")
            continue_str = 'yes'
    return fp

def read_file(fp):
    '''
    PURPOSE: read file pointer and store distortion and structure data assc
            with compounds to be studied
    fp: file pointer
    returns dist_dict: dictionary with compounds as keys and distortion/
                        structure info in a list the value
    '''
    dist_dict = {}
    
    reader = csv.reader(fp)
    header = next(reader, None) #skip header
    
    for line in reader:
        if len(line[1]) > 0: #if the line isn't blank/header, pull info out
            dist_dict[line[0]] = []    
            cryst_syst = line[1]
            structure = line[3]
            ICSD_num = line[4]
            space_grp = line[5]
            avg_bond_len = line[6]
            poly_vol = line[7]
            dist_index = line[8]
            quad_elong = line[9]
            bond_angle_var = line[10]
            ECoN = line[11]
            cation = line[12]
            cat_rank = line[13]
            radius = line[14]
            charge = line[15]
            
            #add values to list in dict
            dist_dict[line[0]] = [cryst_syst, structure, ICSD_num, \
                    space_grp, avg_bond_len, poly_vol, dist_index, quad_elong,\
                    bond_angle_var, ECoN, cation, cat_rank, radius, charge]
            
        else:
            header = next(reader, None) #skip blank lines or secondary headers

    return dist_dict

def graph_values(dist_dict):
    '''
    PURPOSE: prompt user to tie values to variables (x and y axes, color, and size)
            and pull respective values from the imported data dict into a list
            for each variable
    dist_dict: dictionary of imported data values
    returns:
    x/y/c/s list: list of data for each variable
    x/y/c/s lable: lable of data shown by each variable ('n' if not tied to data)
    '''
    print(GRAPH_MENU)
    #prompt user for three graphed characteristics from these choices
    y = input('Please select an option for the y-axis: ')
    x = input('Please select an option for the x-axis: ')
    print('\nThe next two prompts allow you to show more than two variables at'\
          ' once. Filling these out is optional. Type \'n\' if you do not wish'\
            ' to fill out the field.')
    c = input('Please select an option to be shown by the marker color: ')
    s = input('Please select an option to be shown by the marker size: ')
    
    #translate user choice to index value
    choice_dict = {'a':1, 'b':4, 'c':5, 'd':6, 'e':7, 'f':8, 'g':9, 'h':11,\
                   'i':12, 'j':13}
        
    label_dict = {'a':'Structure Type', 'b':'Average Bond Length (A)', \
                  'c':'Polyhedra Volume (A^3)','d':'Distortion Index', \
                  'e':'Quadratic Elongation', 'f':'Bond Angle Variance (deg^2)',\
                  'g':'Effective Coordination Number (ECoN)', 'h':'Cation',\
                  'i':'Ionic Radius (A)', 'j':'Cation Charge'}   

        
    for key, value in choice_dict.items():
        if y == key:
            y_index = value
            y_label = label_dict[key]
        if x == key:
            x_index = value
            x_label = label_dict[key]
        if c == key:
            c_index = value
            c_label = label_dict[key]
        if s == key:
            s_index = value
            s_label = label_dict[key]
  
    #index ditionary to pull out relevant values
    y_list = []
    x_list = []
    c_list = []
    s_list = []
    
    for key, value in dist_dict.items():
        #turn str to flt and append to correct list
        y_val = float(value[y_index])
        y_list.append(y_val)
        x_val = float(value[x_index])
        x_list.append(x_val)
        if c != 'n':
            c_val = float(value[c_index])
            c_list.append(c_val)
        else:
            c_label = 'no'
        if s != 'n':
            s_val = float(value[s_index]) * 90 
            #scale up size values so points visible
            #the scale value may need to be adjusted depending on what value you 
            #choose to display
            s_list.append(s_val)
        else:
            s_label = 'no'
                      
    return [x_list, y_list, c_list, s_list], x_label, y_label, c_label, s_label
        

########### main code ##########

fp = open_file()

dist_dict = read_file(fp)

#graph user-inputted values

cont_str = 'y' #initialize continue string

while cont_str == 'y': #keep making new graphs while user says continue
    values, x_label, y_label, c_label, s_label = graph_values(dist_dict)
    x_vals = values[0]
    y_vals = values[1]
    c_vals = values[2]
    s_vals = values[3]
    
    #user has selected options for both size and color
    if c_label != 'no' and s_label != 'no':
        sc = plt.scatter(x_vals, y_vals, c = c_vals, s = s_vals, alpha = 0.3, \
                cmap = 'viridis')
    #user has selected an option for color but not size
    elif c_label != 'no' and s_label == 'no':
        sc = plt.scatter(x_vals, y_vals, c = c_vals, alpha = 0.3, \
                cmap = 'viridis')
    #user has selected an option for size but not color
    elif c_label == 'no' and s_label != 'no':
        sc = plt.scatter(x_vals, y_vals, s = s_vals, alpha = 0.3)
    #user has only selected x and y axis options
    elif c_label == 'no' and s_label == 'no':
        sc = plt.scatter(x_vals, y_vals)
    
    if c_label != 'no':
        #add and label colorbar    
        cbar = plt.colorbar()
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(c_label, rotation=270)
    
    if s_label != 'no':
        #add and label size legend
        plt.legend(*sc.legend_elements("sizes", num=6), title = s_label,\
                   fontsize = 'small')
    
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()
    
    
    cont_str = input('Make another graph? y/n: ')