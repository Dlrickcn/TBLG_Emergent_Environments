# -*- coding: utf-8 -*-
"""
Created on 30.12.2024

@author: dilaraickecan, yunusemreokyayli
"""
import numpy as np
import matplotlib.pyplot as plt
import joblib
from ase.io.vasp import read_vasp

data = joblib.load('data_small.pkl')  
freq = data['freq']         

model_ridge = joblib.load('models/model_krr.pkl')
model_pca = joblib.load('models/model_pca.pkl')
model_soap = joblib.load('models/model_soap.pkl')
model_scaler = joblib.load('models/model_scaler.pkl')

atoms = read_vasp('POSCAR_61_31_61_30')
              
pos_new = []
for k in range(len(atoms)):
    if atoms.positions[k][2] > 31:
        pos_new.append([atoms.positions[k][0],atoms.positions[k][1],31.0])
    else:
        pos_new.append([atoms.positions[k][0],atoms.positions[k][1],30])
                    
atoms.set_positions(pos_new)

soap_features = model_soap.create(atoms)
features = model_pca.transform(soap_features)
predicted_ldos = []

for test_feature_atom in features:    
    predicted_ldos.append(model_ridge.predict(test_feature_atom.reshape(1,-1)))

predicted_tdos = np.sum(np.array(predicted_ldos),axis=0).flatten()

plt.title('Predict phDOS')
plt.plot(freq,predicted_tdos)
plt.xlabel('Frequency (THz)')
plt.ylabel('Total DOS (a.u.)')
plt.savefig('predicted_PhDOS.png',dpi=200)
plt.show()
