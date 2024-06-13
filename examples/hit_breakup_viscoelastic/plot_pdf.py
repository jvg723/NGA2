import numpy as np
import matplotlib as plt
import matplotlib.ticker as mticker
import seaborn as sns
import pandas as pd

# Import our diameter data
# elvira=pd.DataFrame({'Normalized diameter':np.genfromtxt('elvira/diameter.192.elvira/diameter_1.00000E+01')})
# elvira['Method']='ELVIRA'
# elvira['Weight']=0.16
# lvira=pd.DataFrame({'Normalized diameter':np.genfromtxt('lvira/diameter.192.lvira/diameter_1.00000E+01')})
# lvira['Method']='LVIRA'
# lvira['Weight']=0.135
newt=pd.DataFrame({'Normalized diameter':np.genfromtxt('cases_6112024/01_hb/diameter_2.00000E+01')})
newt['Method']='newt'
newt['Weight']=0.1
beta_05=pd.DataFrame({'Normalized diameter':np.genfromtxt('cases_6112024/02_hb/diameter_2.00000E+01')})
beta_05['Method']='viscr=0.5'
beta_05['Weight']=0.1
beta_005=pd.DataFrame({'Normalized diameter':np.genfromtxt('cases_6112024/03_hb/diameter_2.00000E+01')})
beta_005['Method']='viscr=0.05'
beta_005['Weight']=0.1

data=pd.concat([newt,beta_05,beta_005],ignore_index=True)

# Set plotting style
sns.set_theme()
sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})

# Histogram-style pdf
g=sns.displot(data,bins=30,x='Normalized diameter',hue='Method',weights='Weight',log_scale=True,multiple='dodge',kde=True)
g.set(yscale='log')
g.set(ylabel='Scaled density')
g.set(xbound=[0.026,1])
g.set(ybound=[1e-1,1e3])

# Log-log fitted pdf
# g=sns.kdeplot(data,x='Normalized diameter',hue='Method',weights='Weight',log_scale=[True,True],common_norm=True,bw_adjust=0.8,cut=0)
# g.set(ylabel='Scaled density')
# g.set(xbound=[1e-4,1])
# g.set(ybound=[1e-4,2])

plt.pyplot.show()
