import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from scipy.stats import skew
from scipy.spatial import ConvexHull
import dask
import robocrys as rc

"""
This code is a set of functions that analyze a given pymatgen Structure object in order to extract information about the transition metal-oxygen bond lengths and polyhedra in the structure. The main functions are: 
- get_transition_metal_sites: returns the site numbers of all transition metals in the structure
- _get_metal_bond_lengths: returns all bond lengths between transition metals and oxygen atoms in the structure
- get_average_metal_bond_length: returns the average bond length between transition metals and oxygen atoms in the structure
- get_oxygen_sites: returns the site numbers of all oxygen atoms in the structure
- get_std_polyhedra: returns the standard deviation of the bond lengths for each transition metal in the structure
- get_mean_polyhedra: returns the mean of the bond lengths for each transition metal in the structure
- get_skew_polyhedra: returns the skew of the bond lengths for each transition metal in the structure
- get_distortion_index: returns the distortion index of the bond lengths for each transition metal in the structure, calculated as (length of bond - average length of bond)/average length of bond
- get_poly_oxygens: returns the number of nearest neighbors for each transition metal in the structure

The code also imports the following external libraries:
- numpy: for numerical computations
- pandas: for data manipulation
- matplotlib.pyplot: for plotting data
- pymatgen: for handling and manipulating crystal structures
- scipy.stats: for statistical computations
- scipy.spatial: for calculating geometric properties of sets of points

"""

# write a function that given a pymatgen structure reports all the transition metal-oxygen bond lengths
# a function that given a pymatgen structure reports all the site numbers which contain a transition metal
def get_transition_metal_sites(structure):
    sites = [i.species for i in structure]
    tms = [i for i, site in enumerate(structure) if site.species.elements[0].is_transition_metal \
        or site.species.elements[0].is_metalloid \
            or site.species.elements[0].is_alkali \
                or site.species.elements[0].is_alkaline \
                or site.species.elements[0].is_halogen \
                or site.species.elements[0].is_metal]
    # check if any of the species on each site is a transition metal
    species_per_site = [any([i.is_transition_metal for i in site]) for site in sites]
    return tms

def _get_metal_bond_lengths(structure):
    bond_lengths = structure.distance_matrix
    tms = get_transition_metal_sites(structure)
    oxygens = get_oxygen_sites(structure)
    tms_bond_lengths = [[structure.get_distance(tms[i], ox) for ox in oxygens] for i in range(len(tms))]
    filter_lengths = np.array([[i for i in tms_bond_lengths[i] if i < 2.8] for i in range(len(tms_bond_lengths))])
    return filter_lengths

def get_average_metal_bond_length(structure):
    tms_bond_lengths = _get_metal_bond_lengths(structure)
    tms_bond_lengths =  np.array(tms_bond_lengths).flatten()
    return np.mean(tms_bond_lengths)
    
def get_oxygen_sites(structure):
    sites = [i.species for i in structure]
    # check if any of the species on each site is an oxygen atom
    oxygens = np.array([i for i in range(len(structure)) if structure[i].species.elements[-1].value == 'O'])
    return oxygens

def get_std_polyhedra(structure):
    tms_bond_lengths = _get_metal_bond_lengths(structure)
    return np.array([np.std(i) for i in tms_bond_lengths])

def get_mean_polyhedra(structure):
    tms_bond_lengths = _get_metal_bond_lengths(structure)
    return np.array([np.mean(i) for i in tms_bond_lengths])

def get_skew_polyhedra(structure):
    tms_bond_lengths = _get_metal_bond_lengths(structure)
    return np.array([skew(i) for i in tms_bond_lengths])

def get_distortion_index(structure):
    tms_bond_lengths = _get_metal_bond_lengths(structure)
    avg_bond_length = [np.mean(i) for i in tms_bond_lengths]
    # for each polyhedra calculate the distortion index. (length of bond - average length of bond)/average length of bond
    D = [np.sum([np.abs(i - avg_bond_length[j])/avg_bond_length[j] for i in tms_bond_lengths[j]])/len(tms_bond_lengths[j]) for j in range(len(tms_bond_lengths))]
    return D

def get_poly_oxygens(structure):
    tms = get_transition_metal_sites(structure)
    # use the _get_metal_bond_lengths function to get the number of nearest neighbours for each transition metal
    n_oxy_neighbors = [len(i) for i in _get_metal_bond_lengths(structure)]
    # for each site get the n nearest neighbours for each transition metal
    oxy_neighbor_sites = [structure.get_neighbors(structure[i], n_oxy_neighbors[i]) for i in tms]
    return oxy_neighbor_sites

def get_structure(filename):
    return Structure.from_file(filename)

# determine the polyhedra type by the number of neighbors
def get_polyhedra_types(structure):
    n_neighbors = get_polyhedra_atoms(structure)
    n_neighbors = [len(i)-1 for i in n_neighbors]
    print(n_neighbors)
    struc_types = []
    for n in n_neighbors:
        struc_type = 'octahedral'
        if n == 6:
            struc_type = 'octahedral'
        if n == 5:
            struc_type = 'tetrahedral'
        struc_types.append(struc_type)
    return struc_types

def get_l0(structure):
    tms_bond_lengths = _get_metal_bond_lengths(structure)
    struc_type = np.unique(np.array(get_polyhedra_types(structure)).astype(np.str_))[0]
    print(struc_type)
    sites, vol = prune_poly(structure)
    vol = np.mean(vol)
    if struc_type == 'tetrahedral':
        s = (vol*6*np.sqrt(2))**(1/3)
        h = np.sqrt(2/3)*s
        l0 = h/2
    if struc_type == 'octahedral':
        s = ((vol)*3/np.sqrt(2))**(1/3)
        l0 = np.sqrt(2)*s
    return l0

def quadratic_elongation(structure):
    l0 = np.array(get_l0(structure))
    l = np.array(_get_metal_bond_lengths(structure))
    quad_elongation = np.mean((l/l0)**2)
    return quad_elongation

def get_bond_angles(structure):
    tms = get_transition_metal_sites(structure)
    oxygens = get_oxygen_sites(structure)
    tms_bond_angles = [[structure.get_angle(tms[i], ox, tms[i+1]) for ox in oxygens] for i in range(len(tms)-1)]
    return tms_bond_angles

def bond_angle_variance(structure):
    return 

def get_structure(filename):
    return Structure.from_file(filename)

def get_material_name(structure):
    return structure.composition.reduced_formula

def get_material_structure_list(file_dir):
    materials_list = os.listdir(file_dir)
    materials_list = [i for i in materials_list if i[-4:] == '.cif']
    return [get_structure(file_dir + '/' + i) for i in materials_list]

def get_difference_length(structure):
    difference= [np.max(i) - np.min(i) for i in _get_metal_bond_lengths(Structure.from_file(structure))]
    #difference between max and min of lengths
    return difference

def _get_polyhedra_atoms(structure):
    tms = get_transition_metal_sites(structure)
    polyhedra_list = [structure.get_neighbors_in_shell(structure.cart_coords[i], 1.4, 1.4) for i in tms]
    site_list = [[site.index for site in polyhedra_list[i]] for i in range(len(polyhedra_list))]
    new_site_list = []
    for i, j in enumerate(site_list):
        temp_site_list = j
        temp_site_list.append(tms[i])
        new_site_list.append(temp_site_list)
    return new_site_list

def get_polyhedra_atoms(structure):
    new_site_list = prune_poly(structure)[0]
    return new_site_list

def get_poly_volume(structure, site_list):
    # get the coordinates of the sites
    coords = structure.cart_coords
    # given a set of points use a convex hull to calculate the volume of the polyhedra
    volume1 = [ConvexHull(coords[site_list[i][0:-1], :]).volume for i in range(len(site_list))]
    # import ConvexHull from scipy.spatial
    return volume1

# we are going to make a function that gets the volume of the polyhedra. If the volume of the polyhedra is greater than 10 
# then we need to remove that polyhedra from the list of polyhedra. 
def prune_poly(structure):
    # get the polyhedra atoms
    polyhedra_atoms = _get_polyhedra_atoms(structure)
    # get the volume of the polyhedra
    volume = get_poly_volume(structure, polyhedra_atoms)
    # get the transition metal sites
    tms = get_transition_metal_sites(structure)
    # find the indicies of polyhedra that have a volume greater than 10
    vols2 = []
    tms2 = []
    sites2 = []
    for i, j in enumerate(volume):
        if j < 10 and j > 4:
            vols2.append(j)
            tms2.append(tms[i])
            sites2.append(polyhedra_atoms[i])
    return sites2, vols2, 
    
    
def get_dataframe(file_dir):
    # get the transition metal sites
    structures = get_material_structure_list(file_dir)
    tms = []
    oxygens = []
    bond_lengths = []
    distortion_index = []
    average_bond_length = []
    std_bond_length = []
    mean_bond_length = []
    skew_bond_length = []
    poly_sites = []
    material_name = []
    quad_elongs = []
    grav_caps = []
    files = os.listdir(file_dir)
    for idx, i in enumerate(structures):
        try:
            tms.append(get_transition_metal_sites(i))
            oxygens.append(get_oxygen_sites(i))
            bond_lengths.append(_get_metal_bond_lengths(i))
            distortion_index.append(get_distortion_index(i))
            average_bond_length.append(get_average_metal_bond_length(i))
            std_bond_length.append(get_std_polyhedra(i))
            mean_bond_length.append(get_mean_polyhedra(i))
            skew_bond_length.append(get_skew_polyhedra(i))
            poly_sites.append(get_poly_oxygens(i))
            material_name.append(get_material_name(i))
            #quad_elongs.append(quadratic_elongation(i))
            grav_caps.append(theoretical_gravimetric_capacity(i))
        except:
            print('FAILED TO LOAD', files[i])
        
    df = pd.DataFrame({'material_name': material_name, 'tms': tms, 'oxygens': oxygens, 'bond_lengths': bond_lengths, \
                        'distortion_index': distortion_index, 'average_bond_length': average_bond_length,
                       'std_bond_length': std_bond_length, 'mean_bond_length': mean_bond_length, \
                        'skew_bond_length': skew_bond_length, 'poly_sites': poly_sites, \
                        'theoretical_gravimetric_capacity': grav_caps})
    return df

# write a function to test all the functions we have written so far. 
def unit_test(structure):
    prune_poly(structure)
    get_l0(structure)
    prune_poly(structure), 
    get_poly_volume(structure, 
    prune_poly(structure)[0])
    get_polyhedra_atoms(structure)
    get_transition_metal_sites(structure)
    get_oxygen_sites(structure)
    _get_metal_bond_lengths(structure)
    get_distortion_index(structure)
    get_average_metal_bond_length(structure)
    get_std_polyhedra(structure)
    get_mean_polyhedra(structure)
    get_skew_polyhedra(structure)
    get_poly_oxygens(structure)
    get_material_name(structure)
    return 'unit test passed'

# for each element create a list of the difference between the most common oxidation state and all the other oxidation states of the element. Report 
# the difference

def get_oxidation_state_difference(structure):
    elems = np.array(structure.composition.elements)
    elems = np.array([i for i in elems if i.is_transition_metal == True or \
                        i.is_alkali or i.is_alkaline or i.is_metalloid or i.is_metal or i.is_halogen])
    print(elems)
    d_elems = []
    for i in elems:
        if i.oxidation_states != []:
            d_elem = []
            for j in i.oxidation_states:
                if j >= 0: 
                    d_elem = np.append(d_elem, j - i.common_oxidation_states[0])
            d_elems.append(d_elem)
    return np.abs(np.concatenate(d_elems))

def theoretical_gravimetric_capacity(structure):
    F = 96485
    oxy_diffs = get_oxidation_state_difference(structure)
    molar_mass = structure.composition.reduced_composition.weight
    Q = np.mean([F*i/3.6/molar_mass for i in oxy_diffs])
    return Q

##### SLOPPY CODE ADDITION ####
### summary: new code added with robocrys!
### DATE:/26/2023, we need to reconstruct this to fix the redudancy in the code
### using the new robocrys functions we need update the volume functions as well 
# we need to combine the dataframes into one dataframe but keeping track of the structure they came from
# lets make a staggered dataframe with the structure name as the index
# turn the above into a function that loads a structure and returns a dataframe of the sites

import numpy as np


def clean_pym_struc(structure):
    structure.merge_sites()
    #structure.remove_oxidation_states()
    return structure

def extract_site_df(structure, condensor):
    con_test = condensor.condense_structure(structure)
    site_dict  = con_test['sites']
    # there is a dict of dicts, so we need to convert it to a list of dicts
    sites_dicts = [site_dict[i] for i in site_dict.keys()]
    sites_df = pd.DataFrame(sites_dicts)
    geo2 = [[sites_df['geometry'][i][j] for j in sites_df['geometry'][i].keys()] for i in sites_df['geometry'].keys()]
    nnn2 = [[sites_df['nnn'][i][j] for j in sites_df['nnn'][i].keys()] for i in sites_df['nnn'].keys()]
    sites_df = sites_df.drop(['geometry', 'nnn'], axis=1)
    sites_df = sites_df.assign(geometry=geo2)
    sites_df['polyhedra'] = [i[0] for i in sites_df['geometry']]
    sites_df['polyhedra_likeness'] = [i[1] for i in sites_df['geometry']]
    sites_df = sites_df.drop('geometry', axis=1)
    sites_df['symmetry'] = [con_test['spg_symbol'] for i in range(len(sites_df))]
    sites_df['formula'] = [con_test['formula'] for i in range(len(sites_df))]
    sites_df['crystal_system'] = [con_test['crystal_system'] for i in range(len(sites_df))]
    df = sites_df
    return df

def get_df_robocrys(structure):
    condensor = rc.StructureCondenser()
    con_test = condensor.condense_structure(structure)
    site_dict  = con_test['sites']
    print(con_test)
    # there is a dict of dicts, so we need to convert it to a list of dicts
    df = extract_site_df(structure, condensor)
    return df

def read_df(filename):
    pym_structure = get_structure(filename)
    filename = filename.split('/')[-1]
    #pym_structure = clean_pym_struc(get_structure(filename))
    print('read pymatgen structure for ', filename)
    rc_structure = get_df_robocrys(pym_structure)
    print('read robocrys dataframe for ', filename)
    rc_structure['structure'] = [pym_structure for i in range(len(rc_structure))]
    print('appending data together for ', filename)
    return rc_structure

def fix_robo_df(df):
    # for each row in the dataframe, we need to get the element name and from the element name. The current element name could contain the oxidation state. 
    # we need to remove the oxidation state and add it to the oxidation state column. It is written as 3+ or 3-. in the dataframe make it a float
    # to get the elem name we need to split the string where the first number appears
    split_pos = [[i for i, value in enumerate(list_value) if value.isdigit()][0] for list_value in df['element'].tolist()]

    elems = [df['element'][i][:split_pos[i]] for i in range(len(df['element']))]

    oxidation_states = [df['element'][i][split_pos[i]:] for i in range(len(df['element']))]

    structures = df['structure']
    # if there is a '+' in the name then we can take the numerical value as positive, if there is a '-' then we can take the numerical value as negative. It is possible 
    # that the value is 0, in which case we can just take the value as 0
    oxidation_value = [float(oxidation_states[i][:-1]) if oxidation_states[i][-1] == '+' else -float(oxidation_states[i][:-1]) if oxidation_states[i][-1] == '-' else 0 for i in range(len(oxidation_states))]

    df['center_atom'] = elems

    df['oxidation_state'] = oxidation_value

    #df['mean_nn_distance'] = [np.mean(get_neighbor_bond_lengths(df)) for i in range(len(df['nn']))]

    df['distortion_index'] = get_distortion_index(df)
    df['quadratic_elongation'] = get_quadratic_elongation(df)
    df['mean_nn_distance'] = get_mean_nn_distance(df)
    return df

def get_neighbor_bond_lengths(df):
    # get the octahedral sites
    structure = df['structure'].iloc[0]

    center_atoms = df['center_atom']
    # get all the sites that contain the center atom
    center_sites = [[i for i, value in enumerate(structure.sites) if center_atoms[j] in value.species_string] for j in range(len(center_atoms))]

    # depending on the len of "nn" we will grab the n smallest bond lengths for the oct site row in the distance matrix
    nn_number = [len(df['nn'][i]) for i in range(len(df['nn']))]
    center_rows = [structure.distance_matrix[i,:] for i in center_sites]
    # now we need to get the n smallest bond lengths for each row
    nearest_neighbors = [np.sort(center_rows[i])[:nn_number[i]] for i in range(len(center_rows))]
    #print(nearest_neighbors[0])
    # remove the first value as this is the distance to itself
    nearest_neighbors = [nearest_neighbors[i][0][1:] for i in range(len(nearest_neighbors))]
    #print(nearest_neighbors[0])
    return nearest_neighbors

def get_distortion_index(df):
    nearest_neighbors = get_neighbor_bond_lengths(df)
    bond_lengths = nearest_neighbors
    nn_number = [len(nearest_neighbors[i]) for i in range(len(nearest_neighbors))]
    mean_bond_length = [np.mean(bond_lengths[i]) for i in range(len(bond_lengths))]
    D = [np.sum((bond_lengths[i] - mean_bond_length[i])/mean_bond_length[i])/len(bond_lengths[i]) for i in range(len(bond_lengths))]
    return D

def get_l0(df):
    nearest_neighbors = get_neighbor_bond_lengths(nearest_neighbors)
    bond_lengths = nearest_neighbors
    nn_number = [len(nearest_neighbors[i]) for i in range(len(nearest_neighbors))]
    mean_bond_length = [np.mean(bond_lengths[i]) for i in range(len(bond_lengths))]
    # we will assume that l0 is the mean bond length
    return mean_bond_length

def get_quadratic_elongation(df):
    nearest_neighbors = get_neighbor_bond_lengths(df)
    bond_lengths = nearest_neighbors
    nn_number = [len(nearest_neighbors[i]) for i in range(len(nearest_neighbors))]
    mean_bond_length = [np.mean(bond_lengths[i]) for i in range(len(bond_lengths))]
    Q = [np.sum((bond_lengths[i]/mean_bond_length[i])**2)/len(bond_lengths[i]) for i in range(len(bond_lengths))]
    return Q

def get_mean_nn_distance(df):
    bond_lengths = get_neighbor_bond_lengths(df)
    return [np.mean(bond_lengths[i]) for i in range(len(bond_lengths))]


###### USE THIS ONE TO DIRECTLY READ THE DATAFRAME FROM THE FILE
def read(filename):
    try:
        df = read_df(filename)
        df = fix_robo_df(df)
    except:
        df = pd.DataFrame()
    return df
