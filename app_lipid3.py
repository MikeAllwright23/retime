###################### STREAMLIT APP FILE TO PREDICT RETENTION TIMES (RETIME) ######################

#This code drives a webapp powered by streamlit and enables users to provide an Excel table with known retention times and masses for 
# a set of Lipids, and based on what is known predict the retention time of other lipids to a certain level of accuracy.
# the program has 6 stages from importing the file through to running the model and exporting a file with predicted retention times
# for the unknown lipids

# CREATED BY
# DATE


import streamlit as st
import numpy as np
import time
import pandas as pd
#import openpyxl
import matplotlib.pyplot as plt
import os
from pathlib import Path
import re
import seaborn as sns

from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score



# where data is stored
path='raw_data/'
predictors=['spingoid_backbone_dbl_bonds','fatty_acyl_dbl_bonds',\
'spingoid_backbone_carb', 'fatty_acyl_carb','log_mass','Mass','mass_squared','mass_sqrt']#'Mass','mass_squared','mass_sqrt',
file_out='new_retenton_time_predictions.csv'


#1 - import raw file produced by users

def file_selector(folder_path=path):
	filenames = os.listdir(folder_path)
	selected_filename = st.selectbox('Select a file', filenames)
	return os.path.join(folder_path, selected_filename)


#2 - calculate all variables
#

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

	mol_series="d18:"+dc1+"/"+"XX:"+dc2
	if re.search('OH',string):
		oh=1
	
	return c1,dc1,c2,dc2,oh,mol_series

def mass_trans(df,verbose=False):

	if verbose:
		st.write('columns available:',dict(zip([c for c in df.columns],[len(c) for c in df.columns])))	

	#creates more variables for modelling
	df['log_mass']=df['Mass'].apply(lambda x:np.log(x))
	df['mass_squared']=df['Mass'].apply(lambda x:pow(x,2))
	df['mass_sqrt']=df['Mass'].apply(lambda x:pow(x,0.5))
	
	return df

#3 - split train and test

def train_test_try(df,preds):
	mask1=(df['Type']=="Train")

	X_train=df.loc[mask1,preds]
	y_train=df.loc[mask1,'Retention Time']

	mask2=(df['Type']=="Predict")

	X_test=df.loc[mask2,preds]
	y_test=df.loc[mask2,'Retention Time']


	#y_test=df['Retention Time'][mask2]
	
	return X_train,y_train,X_test,y_test

scaler = StandardScaler()
def scale_df(X,cols):
	X[cols]=scaler.fit_transform(X[cols])
	return X


def read_markdown_file(markdown_file):
	return Path(markdown_file).read_text()


#4 - model training data

#5 - predict on test data

#6 - output revised csv

st.set_page_config(page_title='ReTimeML',  layout='wide', page_icon=':machine learning:')

st.title('ReTimeML: A Retention Time Predictor for the LC - MS/MS analysis of ceramides and sphingomyelins*')



#OLD and poss useful


app_mode = st.sidebar.selectbox("Choose the app mode",
["Show instructions", "Run the app", "Show the source code","Download Example Data"])


if app_mode=="Show the source code":
	st.write("For source code please contact michael.allwright@sydney.edu.au")

elif app_mode=="Show instructions":
	intro_markdown = read_markdown_file("intro.md")
	st.markdown(intro_markdown, unsafe_allow_html=True)

elif app_mode=="Run the app":

	intro_markdown = read_markdown_file("intro.md")
	st.markdown(intro_markdown, unsafe_allow_html=True)

	uploaded_file = st.file_uploader("Choose a file")
	
	if uploaded_file is not None:

		file_name=re.sub('.xlsx|.csv','',uploaded_file.name)

		#st.write('You selected the following file: '+file_name)

		st.write("Now running the ReTimeML algorithm...")


		try:
			df_all=pd.read_csv(uploaded_file)
			csv=True
		except:
			df_all=pd.read_excel(uploaded_file,sheet_name=None)#,engine="openpyxl")
			csv=False
		
		if csv is True:
			
			df_all=dict({'outfile':df_all})

		#automatically assert Type field if not present:


		#st.write(file_type)
		clf_lass = linear_model.Lasso(alpha=0.05)
		clf_ridge=linear_model.Ridge(alpha=0.1)
		
		for k in df_all.keys():

			#file_out='new_retenton_time_predictions.csv'

			#truncate columns
			
			df=df_all[k]

			
			df.columns=[c.strip() for c in df.columns]

			#mask=(df['Lipid ID'].astype(str).apply(lambda x:len(str(x)))>5)
			#df=df.loc[mask,]

			#st.write(df.columns)

			#st.write(df)

			if "RT" in df.columns:
				df.rename(columns={'RT':'Retention Time'},inplace=True)

			if 'Type' not in df.columns:
				df['Type']='Train'
				mask=pd.isnull(df['Retention Time'])
				df.loc[mask,'Type']='Test'

			df=mass_trans(df)

			dfk=df.copy()
			
			carbs=df['Lipid ID'].astype(str).apply(findx)

			#st.write(df['carbs'])

			for i,c in enumerate(['spingoid_backbone_carb', 'spingoid_backbone_dbl_bonds','fatty_acyl_carb', 'fatty_acyl_dbl_bonds',
							  'OH','mol_series']):
				#st.write(c)
				df[c]=[carb[i] for carb in carbs]
			
			df=df[pd.notnull(df['Lipid ID'])]
			
			
			#st.write(df)
			mask=(df['Type']=="Train")
			df_train=df.loc[mask,]

			preds2=[p for p in predictors if df_train[p].nunique()>1]

			#st.write(preds2)

			

			#st.write(dict(zip(predictors,[df_train[p].unique() for p in predictors])))


			#st.write(df)

			#st.write(df.dtypes)
			
			#st.write(df)
			df=scale_df(df,preds2)

			
			
			X_train,y_train,X_test,y_test=train_test_try(df,preds=preds2)

			# DIFFERENT POTENTIAL MODELS BELOW
			#reg =clf_reg.fit(X_train, y_train)
			ridge=clf_ridge.fit(X_train, y_train)
			#xgb_mod_trained=xgb_mod.fit(X_train, y_train)
			
			lass = clf_lass.fit(X_train, y_train)
			
			
			if df['Lipid ID'].apply(lambda x:x[0:3]=="Cer")[0]:
				#st.write("Cer")
				model_select="ridge"
				df['pred_ret']=ridge.predict(df[preds2])
			
			else:
				df['pred_ret']=lass.predict(df[preds2])
				model_select="lasso"
			

			

			#now show the dumbell chart with predictions followed by the ability to download

			df_out=df[['Lipid ID','mol_series','Type','pred_ret','Retention Time']]

			df_out=pd.merge(dfk[['Lipid ID','Mass']],df_out,on='Lipid ID')

			df_out.columns=['Lipid ID','Mass','mol_series','Type','Predicted Retention Time','Actual Retention Time']

			for c in ['Predicted Retention Time','Actual Retention Time']:

				df_out[c]=df_out[c].apply(lambda x:round(x,2))

			st.subheader("Break down of key molecules")
			#st.write("Here are the retime generated retention times:")
			#st.write(df_out)	

			@st.cache
			def convert_df(df):
				# IMPORTANT: Cache the conversion to prevent computation on every rerun
				return df.to_csv().encode('utf-8')

			csv = convert_df(df_out)
			
			file_out=file_name+' '+k+'.csv'

			#st.write(df_out.head())
			df_out.rename(columns={'Mass':'m/z'},inplace=True)
			df_out.sort_values(by='mol_series',inplace=True)
			fig=sns.lmplot(data=df_out,x='m/z',y='Predicted Retention Time',hue='mol_series',order=2)
			st.pyplot(fig)

			

			st.subheader("Download Output Data")

			st.download_button(
				label="Download data as CSV",
				data=csv,
				file_name=path+file_out,
				mime='text/csv',
			)

			
			st.write(df_out)


		#run models

elif app_mode=="Download Example Data":

	st.header("Example Data Files Download")

	st.write("Download supplementary files to test ReTimeML here:")

	lips4=['Cer (d18:2/16:1)','Cer (d18:2/16:0)','Cer (d18:1/16:1)','Cer (d18:1/16:0)','Cer (d18:0/16:0)','Cer (d18:1/17:0)','Cer (d18:2/18:1)','Cer (d18:2/18:0)','Cer (d18:1/18:1)',
	'Cer (d18:1/18:0)','Cer (d18:0/18:0)','Cer (d18:2/20:1)','Cer (d18:2/20:0)','Cer (d18:1/20:1)','Cer (d18:1/20:0)','Cer (d18:0/20:0)','Cer (d18:2/22:1)',
	'Cer (d18:2/22:0)','Cer (d18:1/22:1)','Cer (d18:1/22:0)','Cer (d18:0/22:0)','Cer (d18:2/24:1)','Cer (d18:2/24:0)','Cer (d18:1/24:1)','Cer (d18:1/24:0)','Cer (d18:0/24:0)']

	mass4=[534.5,536.5,536.5,538.6,540.5,552.7,
	562.5,564.5,564.5,566.6,568.6,590.6,
	592.6,592.6,594.6,596.6,618.6,620.6,
	620.6,622.6,624.6,646.6,648.6,648.6,650.7,652.9]

	retimes4=[np.nan,np.nan,np.nan,15.05,15.27,15.33,np.nan,np.nan,15.05,15.7,15.91,np.nan,np.nan,np.nan,16.28,np.nan,np.nan,np.nan,np.nan,16.94,
		np.nan,np.nan,np.nan,17.08,17.68,17.95]
	

	types4=['Test','Test','Test','Train','Train','Train','Test',
	'Test','Train','Train','Train','Test','Test','Test',
	'Train','Test','Test','Test','Test','Train',
	'Test','Test','Test','Train','Train','Train']

	lips5=['SM (d18:1/12:0)','SM (d18:2/16:1)','SM (d18:2/16:0)','SM (d18:1/16:1)','SM (d18:1/16:0)','SM (d18:0/16:0)','SM (d18:2/18:1)',
	'SM (d18:2/18:0)','SM (d18:1/18:1)','SM (d18:1/18:0)','SM (d18:0/18:0)','SM (d18:2/20:1)','SM (d18:2/20:0)','SM (d18:1/20:1)',
	'SM (d18:1/20:0)','SM (d18:0/20:0)','SM (d18:2/22:1)','SM (d18:2/22:0)','SM (d18:1/22:1)','SM (d18:1/22:0)',
	'SM (d18:0/22:0)','SM (d18:2/24:1)','SM (d18:2/24:0)','SM (d18:1/24:1)','SM (d18:1/24:0)','SM (d18:0/24:0)',]

	mass5=[647.5,699.6,701.6,701.6,703.5,705.6,
	727.6,729.6,729.6,731.6,733.6,755.6,
	757.6,757.6,759.6,761.6,783.7,785.7,
	785.7,787.7,789.7,811.7,813.7,813.7,815.7,817.7,]

	retimes5=[2.73,np.nan,np.nan,3.04,3.9,np.nan,
	np.nan,np.nan,np.nan,4.96,5.58,np.nan,
	np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,
	np.nan,np.nan,np.nan,np.nan,np.nan,9.83,11.1,np.nan,]


	types5=['Train','Test','Test','Train','Train','Test',
	'Test','Test','Test','Train','Train','Test',
	'Test','Test','Test','Test','Test','Test',
	'Test','Test','Test','Test','Test','Train','Train','Test']

	df_supp4=pd.DataFrame({'Lipid ID':lips4,'Mass':mass4,'Retention Time':retimes4,'Type':types4})

	df_supp5=pd.DataFrame({'Lipid ID':lips5,'Mass':mass5,'Retention Time':retimes5,'Type':types5})

	@st.cache
	def convert_df(df):
		# IMPORTANT: Cache the conversion to prevent computation on every rerun
		return df.to_csv().encode('utf-8')

	csv_supp4 = convert_df(df_supp4)
	csv_supp5 = convert_df(df_supp5)

	st.subheader("Download Example Data Files")

	st.download_button(
		label="Download Supplementary File S4",
		data=csv_supp4,
		file_name="Supplementary File S4.csv",
		mime='text/csv',
	)

	st.download_button(
		label="Download Supplementary File S5",
		data=csv_supp5,
		file_name="Supplementary File S5.csv",
		mime='text/csv',
	)



#import data in format lipid ID, Mass and retention time

#when mass is not available, fill in with mass from lipid list which we create

