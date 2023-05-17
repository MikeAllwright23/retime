import streamlit as st
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import re
from sklearn import linear_model
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score



#from sklearn.linear_model import Lasso

#path='/Users/michaelallwright/Documents/data/lipid/model_data/'
path='data/'
predictors=['spingoid_backbone_dbl_bonds','fatty_acyl_dbl_bonds',\
'spingoid_backbone_carb', 'fatty_acyl_carb','mass_squared','mass_sqrt','log_mass','Mass']#'Mass',



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

def mass_trans(df):

	#creates more variables for modelling
	df['log_mass']=df['Mass'].apply(lambda x:np.log(x))
	df['mass_squared']=df['Mass'].apply(lambda x:pow(x,2))
	df['mass_sqrt']=df['Mass'].apply(lambda x:pow(x,0.5))
	
	return df

def train_test_try(df,preds):
    mask1=(df['Train/Test']=="Train")

    X_train=df.loc[mask1,preds]
    y_train=df.loc[mask1,'Retention Time']

    #y_test=df['Retention Time'][mask2]
    
    return X_train,y_train


def dumbell_plot(df,experiment=''):
	# Setup plot size.
	fig, ax = plt.subplots(figsize=(6,9))

	# Create grid 
	# Zorder tells it which layer to put it on. We are setting this to 1 and our data to 2 so the grid is behind the data.
	ax.grid(which="major", axis='both', color='#758D99', alpha=0.6, zorder=1)

	# Remove splines. Can be done one at a time or can slice with a list.
	ax.spines[['top','right','bottom']].set_visible(False)

	# Setup data
	gdp_dumbbell = df.copy()

	# Plot data
	# Plot horizontal lines first
	ax.hlines(y=df['Lipid ID'], xmin=df['lass_ret'], xmax=df['Retention Time'], color='#758D99', zorder=2, linewidth=2, label='_nolegend_', alpha=.8)
	# Plot bubbles next
	ax.scatter(df['lass_ret'], df['Lipid ID'], label='Low Variable', s=60, color='#006BA2', zorder=3)
	ax.scatter(df['Retention Time'], df['Lipid ID'], label='High Variable', s=60, color='#DB444B', zorder=3)

	# Set xlim
	x_min=df['Retention Time'].min()-0.5
	x_max=df['Retention Time'].max()+0.5
	ax.set_xlim(x_min, x_max)

	# Reformat x-axis tick labels
	ax.xaxis.set_tick_params(labeltop=True,      # Put x-axis labels on top
	                         labelbottom=False,  # Set no x-axis labels on bottom
	                         bottom=False,       # Set no ticks on bottom
	                         labelsize=15,       # Set tick label size
	                         pad=-1)             # Lower tick labels a bit

	# Reformat y-axis tick labels
	ax.set_yticklabels(gdp_dumbbell['Lipid ID'],       # Set labels again
	                   ha = 'left')              # Set horizontal alignment to left
	ax.yaxis.set_tick_params(pad=300,            # Pad tick labels so they don't go over y-axis
	                         labelsize=15,       # Set label size
	                         bottom=False)       # Set no ticks on bottom/left

	# Set Legend
	ax.legend(['Predicted Retention Time', 'Actual Retention Time'], loc=(-0.9,1.09), ncol=2, frameon=False, handletextpad=-.1, handleheight=1)

	# Add in line and tag
	ax.plot([-0.59	, .9],                 # Set width of line
	        [1.17, 1.17],                # Set height of line
	        transform=fig.transFigure,   # Set location relative to plot
	        clip_on=False, 
	        color='#E3120B', 
	        linewidth=.6)
	ax.add_patch(plt.Rectangle((-0.59,1.17),               # Set location of rectangle by lower left corder
	                           0.05,                       # Width of rectangle
	                           -0.025,                      # Height of rectangle. Negative so it goes down.
	                           facecolor='#E3120B', 
	                           transform=fig.transFigure, 
	                           clip_on=False, 
	                           linewidth = 0))

	
	# Add in title and subtitle
	ax.text(x=-0.59, y=1.09, s=experiment+': Retention Time predictions (minutes)', transform=fig.transFigure, ha='left', fontsize=16, weight='bold', alpha=.8)
	ax.text(x=-0.59, y=1.04, s="Difference between predicted and actual Retention Time", transform=fig.transFigure, ha='left', fontsize=14, alpha=.8)

	# Set source text
	ax.text(x=-0.59, y=0.04, s="""Source: "Tim Couttas""", transform=fig.transFigure, ha='left', fontsize=9, alpha=.7)


	st.pyplot(fig)

	# Export plot as high resolution PNG
	#plt.savefig('economist_dumbbell.png',    # Set path and filename
	#            dpi = 300,                          # Set dots per inch
	#            bbox_inches="tight",                # Remove extra whitespace around plot
	#            facecolor='white')                  # Set

scaler = StandardScaler()
def scale_df(X,cols):
    X[cols]=scaler.fit_transform(X[cols])
    return X

st.title('ReTime: Lipid Retention Time Predictor')

app_mode = st.sidebar.selectbox("Choose the app mode",
["Show instructions", "Run the app", "Show the source code"])

def read_markdown_file(markdown_file):
	return Path(markdown_file).read_text()

if app_mode=="Show instructions":
	intro_markdown = read_markdown_file("instructions.md")
	st.markdown(intro_markdown, unsafe_allow_html=True)

else:	

	filename = file_selector()
	df=pd.read_csv(filename)

	

	df=mass_trans(df)
	
	df['carbs']=df['Lipid ID'].apply(findx)

	for i,c in enumerate(['spingoid_backbone_carb', 'spingoid_backbone_dbl_bonds','fatty_acyl_carb', 'fatty_acyl_dbl_bonds',
					  'OH']):
		df[c]=df['carbs'].apply(lambda x:x[i])

	species=list(df['Species'].unique())
	spec_sel = st.selectbox('Select Species',species)

	clf_lass = linear_model.Lasso(alpha=0)

	mask=(df["Species"]==spec_sel)
	df2=df.loc[mask,]
	mask=(df2['Train/Test']=="Train")
	df_train=df2.loc[mask,]

	preds2=[p for p in predictors if df2.loc[mask,p].nunique()>1]
    
	df2=scale_df(df2,preds2)
    
	X_train,y_train=train_test_try(df2,preds=preds2)
    #reg =clf_reg.fit(X_train, y_train)
	lass =clf_lass.fit(X_train, y_train)
    #ridge=clf_ridge.fit(X_train, y_train)
    #xgb_mod_trained=xgb_mod.fit(X_train, y_train)
    
	df2['lass_ret']=lass.predict(df2[preds2])

	mask=(df2['Train/Test']=="Train")&pd.notnull(df2['Retention Time'])
	df_test=df2.loc[mask,]

	r2=r2_score(df_test['lass_ret'],df_test['Retention Time'])
	st.subheader("R Squared: "+str(round(r2,4)))


	data_model = st.selectbox('Select Data to display and output',['Validation Only','All'])


	if data_model=='Validation Only':
		dumbell_plot(df_test,experiment=spec_sel)
	else:
		dumbell_plot(df2,experiment=spec_sel)
	plt.show()




	#now show the dumbell chart with predictions followed by the ability to download

	df_out=df2[['Lipid ID','Train/Test','lass_ret','Retention Time']]

	df_out.columns=['Lipid ID','Train/Test','Predicted Retention Time','Actual Retention Time']

	for c in ['Predicted Retention Time','Actual Retention Time']:

		df_out[c]=df_out[c].apply(lambda x:round(x,2))

	st.write(df_out)	

	@st.cache
	def convert_df(df):
	    # IMPORTANT: Cache the conversion to prevent computation on every rerun
	    return df.to_csv().encode('utf-8')

	csv = convert_df(df_out)

	st.download_button(
	    label="Download data as CSV",
	    data=csv,
	    file_name=spec_sel+'_retenton_time_predictions.csv',
	    mime='text/csv',
	)


	#run models


	#import data in format lipid ID, Mass and retention time

	#when mass is not available, fill in with mass from lipid list which we create

