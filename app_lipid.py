import streamlit as st
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path

#path='/Users/michaelallwright/Documents/data/lipid/model_data/'
path='data/'


def file_selector(folder_path=path):
	filenames = os.listdir(folder_path)
	selected_filename = st.selectbox('Select a file', filenames)
	return os.path.join(folder_path, selected_filename)

def findx(string):

	#based on molecule structure being of the form 'Cer (d18:1/18:1)' returns key variables
    b1=string.find('(')
    d=string.find('/')
    b2=string.find(')')
    c1=string[b1+2:b1+4]
    dc1=string[d-1:d]
    c2=string[d+1:d+3]
    dc2=string[d+4:d+5]
    oh=0
    if re.search('OH',string):
        oh=1
    
    
    
    return c1,dc1,c2,dc2,oh

def data_proc(df):

    df['log_mass']=df['Mass'].apply(lambda x:np.log(x))
    df['mass_squared']=df['Mass'].apply(lambda x:pow(x,2))
    df['mass_sqrt']=df['Mass'].apply(lambda x:pow(x,0.5))
    
    df['s_back_dbl_mass_ratio']=df['spingoid_backbone_dbl_bonds']/df['Mass']
    df['s_back_dbl_sp_carb_ratio']=df['spingoid_backbone_dbl_bonds']/df['spingoid_backbone_carb']
    df['fatty_acyl_dbl_mass_ratio']=df['fatty_acyl_dbl_bonds']/df['Mass']
    df['fatty_acyl_fa_carb_ratio']=df['fatty_acyl_dbl_bonds']/df['fatty_acyl_carb']
    
    one_hot=pd.get_dummies(df["Species"])
    df=pd.concat([df,one_hot],axis=1)

    return df

st.title('ReTime: Lipid Retention Time Predictor')

app_mode = st.sidebar.selectbox("Choose the app mode",
["Show instructions", "Run the app", "Show the source code"])

def read_markdown_file(markdown_file):
	return Path(markdown_file).read_text()

if app_mode=="Show instructions":
	intro_markdown = read_markdown_file("intro.md")
	st.markdown(intro_markdown, unsafe_allow_html=True)

else:	

	st.write('Select File')
	filename = file_selector()

	df=pd.read_csv(filename)

	value_vars=['Cer 3mth mice', 'Cer_diabetes_cohort', 'Ceramide 6mth mice',
	       'Human Hippocamupus Cer', 'Short Cer']
	#df_proc1.fillna(0,inplace=True)
	id_vars=[c for c in df.columns if c not in value_vars]
	df_proc=pd.melt(df,id_vars=id_vars, value_vars=value_vars)
	df_proc.rename(columns={'variable':'Species','value':'Retention Time'},inplace=True)
	df_proc=data_proc(df_proc)	
	df_proc['Species_type']="Ceramide"

	df_proc=df_proc[['Lipid ID', 'spingoid_backbone_carb', 'spingoid_backbone_dbl_bonds',
	       'fatty_acyl_carb', 'fatty_acyl_dbl_bonds', 'Mass',
	       'Precursor Mass [M + H]+', 'Fragment ion 1', 'Fragment ion 2',
	       'Atom mass increase from prior SL with same degrees of unsaturation',
	       'Cer 3mth mice', 'Cer_diabetes_cohort', 'Ceramide 6mth mice',
	       'Human Hippocamupus Cer', 'Short Cer']]

	st.write("Visualise Data")
	st.write(df_proc)

	lipids=list(df_proc['Lipid ID'])
"""
	r=np.zeros(len(lipids))

	for i,l in enumerate(lipids):
		r[i]=st.number_input(l)
"""
