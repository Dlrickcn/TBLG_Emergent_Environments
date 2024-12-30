# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 15:22:58 2023
@author: yunus
"""
from ase.io.xyz import read_xyz
from ase.io import read
from pymatgen.core import Lattice, Structure, Molecule
from pymatgen.io.ase import AseAtomsAdaptor
from dscribe.descriptors import SOAP
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler, Normalizer
import numpy as np
import pandas as pd
import os
from ase.visualize import plot
from ase.neighborlist import NeighborList
from ase import Atoms
import itertools

np.set_printoptions(threshold=np.inf)
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

path = 'output'
os.makedirs(path,exist_ok=True)

def read_dos(file):
    energy = np.loadtxt(file,usecols=0)
    dos = np.loadtxt(file,usecols=1)
    return energy, dos

def get_atoms(file_name):
    crystal_vector = np.loadtxt(f"./Atom_Koordinatlari/{file_name}", skiprows=2, usecols=[6,7,8], max_rows=3, dtype=np.float64)
        
    structure = Molecule.from_file(filename=f"./Atom_Koordinatlari/{file_name}")
    atoms = AseAtomsAdaptor.get_atoms(structure)
        
    atoms.cell = crystal_vector
    atoms.set_chemical_symbols(['C']*len(atoms))
    
    pos_new = []
    for k in range(len(atoms)):
        if atoms.positions[k][2]> 31:
            pos_new.append([atoms.positions[k][0],atoms.positions[k][1],31.0])
        else:
            pos_new.append([atoms.positions[k][0], atoms.positions[k][1],30])
    atoms.set_positions(pos_new)
    atoms.pbc = [True,True,False]

    return atoms


def get_dos(file_name,atoms):
    energy_list = []
    pdos_list = []
    for atom in range(len(atoms)):
        energy, dos = read_dos(f'{file_name}/pldos_{atom}.datnorm')
        energy_list.append(energy)
        pdos_list.append(dos)
        
    pdos_list = np.array(pdos_list)
    ldos = pdos_list.sum(axis=0)
    ldos = np.array(ldos)
    pdos_list = np.array(pdos_list)
    energy_list = np.array(energy_list)
    
    return energy_list, pdos_list, ldos

def get_angle(i_pos,neig_pos):
    norm_i = i_pos / np.linalg.norm(i_pos)
    norm_n = neig_pos /  np.linalg.norm(neig_pos)
    angle = np.arccos(np.clip(np.dot(norm_i, norm_n), -1.0, 1.0))    
    return np.degrees(angle)

def get_features(atoms):
    r_cut = 6
    n_max = 10
    l_max = 10
    sigma = 0.5
    species = atoms.get_chemical_symbols()
    
    soap_periodic = SOAP(r_cut=r_cut,n_max=n_max,l_max=l_max,sigma=sigma,rbf='gto',periodic=True,species=species)

    features = soap_periodic.create(atoms)
    rounded_features = np.round(features, decimals=5,)
    
    pca = PCA(n_components=6)
    scaler = StandardScaler()
        
    features_scale = scaler.fit_transform(features)
    pcafeat = np.round(pca.fit_transform(features_scale), decimals=3)
    unique_feature , unique_index , unique_counts = np.unique(rounded_features, axis=0, return_counts=True, return_index= True)
    unique_feature_pca , unique_index_pca , unique_counts_pca = np.unique(np.round(pca.fit_transform(features_scale), decimals=3), axis=0, return_counts=True, return_index= True)


    PC_values = np.arange(pca.n_components_) + 1
    print('PCA_VARIANCE',pca.explained_variance_ratio_)
    print('SUM_PCA_VARIANCE',sum(pca.explained_variance_ratio_))
    plt.plot(PC_values, pca.explained_variance_ratio_, 'o-', linewidth=2, color='blue')
    plt.title('Scree Plot')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained')
    #plt.savefig(f'output/{file_name}_1.png',dpi=50)
    plt.show()
  
    fig, (ax1,ax2) = plt.subplots(2,1)
    fig.set_size_inches(6, 9)
    #plt.figure(dpi=15,figsize=(10,60))
    fig.tight_layout()

    ax1.set_title(f"Unique elements of {file_name} numbers: {len(unique_index)}")
    print(len(atoms))
    plot.plot_atoms(atoms, ax1,show_unit_cell=2,radii=[0.4]*len(atoms),colors=['red']*len(atoms))
    ax2.scatter(unique_feature_pca.T[0],unique_feature_pca.T[1],marker='o',s=200,alpha=0.9,edgecolor='orange',linewidths=5)

    ax1.axis('off')
    ax2.set_xlabel('PCA1')
    ax2.set_ylabel('PCA2')
    #plt.savefig(f'output/{file_name}.png',dpi=50)
    plt.show()
  
    return pcafeat, unique_feature_pca, unique_index, unique_counts

def plot_pca(train, test):
    pca = PCA(n_components=2)
    features = np.concatenate((train,test))
    model = pca.fit(features)
    features_train = model.transform(train).T
    features_test = model.transform(test).T
    plt.plot(features_train[0],features_train[1],'.',c='blue',label='Train')
    plt.plot(features_test[0],features_test[1],'.',c='red',label='Test')
    plt.xlabel('PCA1')
    plt.ylabel('PCA2')
    plt.legend()
    plt.show()    

atoms_list = []
unique_list = []
features_list = []
features_pca_list = []
features_ave_list = []

unique_array_length = []
unique_array_index = []
unique_array_counts = []
unique_array_pca = []

files = os.listdir(f'./Atom_Koordinatlari')[:]
files = [i[:-4] for i in files] #.xyz uzantısını çöpe atmak için.

moire = pd.read_excel("database.xlsx", usecols="A", skiprows=[-1])
lattice = pd.read_excel("database.xlsx", usecols="B", skiprows=[-1])['lattice']
angle = pd.read_excel("database.xlsx", usecols="C", skiprows=[-1])['angle']
unitcell_atom = pd.read_excel("database.xlsx", usecols="D", skiprows=[-1])['unitcell_atom']
supercell = pd.read_excel("database.xlsx", usecols="E", skiprows=[-1])['supercell']

moire_list = []
name_list = []
angle_list = []
unitcell_atom_list = []
supercell_list = []
lattice_list = []
space_groups = []

from ase import spacegroup

for i in range(len(moire['moire'])):
    name = moire['moire'][i].replace(",", "_").replace(" ", "").replace("[","").replace("]","")
    
    if name in files:
    
        name_list.append(name)
        moire_list.append(name)
        file_name = f"{name}.xyz"
        path = f"./Atom_Koordinatlari/{file_name}"
        if os.path.isfile(path) == True:
            atoms = get_atoms(file_name)
            
            space_groups.append(spacegroup.get_spacegroup(atoms).symbol)
            
            features, features_pca,  uniq_index, uniq_counts = get_features(atoms)
            atoms_list.append(atoms)
            unique_array_length.append(len(uniq_index)) # Unique atomların sayısının listesi
            unique_array_index.append(uniq_index)   # Unique atomların indexleri
            unique_array_counts.append(uniq_counts)
            features_list.append(features)
            features_pca_list.append(features_pca)
            angle_list.append(angle[i])
            unitcell_atom_list.append(unitcell_atom[i])
            supercell_list.append(supercell[i])
            lattice_list.append(lattice[i])
    
        else:
            print(path,' yok.')

### Unit Cell Atom Size ile Unique Atomları Çizdirmek için.
fig, ax = plt.subplots(1,1, figsize=(9,10))
ax.plot(unitcell_atom_list,unique_array_length,'o')
ax.set_xlabel('Unit Cell Atom Count', fontsize=12, fontname='Times New Roman')
ax.set_ylabel('Unique Atom Count', fontsize=12, fontname='Times New Roman')
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
plt.savefig('output/unique_vs_unitcell.png')
plt.show()


dict_database = {'moire':moire_list, 'lattice':lattice_list, 'angle':angle_list, 'unitcell_atom':unitcell_atom_list, 'supercell':supercell_list, 'unique_atom':unique_array_length, 'unique_hist':unique_array_counts,'space_group':space_groups}
df = pd.DataFrame.from_dict(dict_database)
df.to_csv('database.csv')
with pd.ExcelWriter("data_new.xlsx") as writer:
    df.to_excel(writer)  
print('Finished calculated PCA.')
