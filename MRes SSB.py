#!/usr/bin/env python
# coding: utf-8

# In[40]:


import numpy as np  
import matplotlib.pyplot as plt  
from scipy.integrate import odeint  
import math
plt.rcParams.update({'font.sans-serif':'Helvetica'})


# In[49]:


# define the function
def PET_degradation(s,t,params):  
    
    km, km0, kdm, kp, kdp, K, n, kcatP, kmP, X = params 
    m_PETase, p_PETase, cPET = s
        
    rate_m_PETase_prod = X * (km + km0)
    rate_p_PETase_prod = X * (kp * m_PETase)
    
    rate_m_PETase_loss = kdm * m_PETase
    rate_p_PETase_loss = kdp * p_PETase
    
    dm_PETase = rate_m_PETase_prod - rate_m_PETase_loss
    dp_PETase = rate_p_PETase_prod - rate_p_PETase_loss
    
    dcPET = - (kcatP * p_PETase * c_PET)/(kmP + c_PET)
    
    ds = [dm_PETase, dp_PETase, dcPET]
    
    return ds  


# In[50]:


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set(xlabel = 'Time (mins)')
ax.set(ylabel = 'PETase expression (copies per cell)')

km = 30
km0 = 0.03
kdm = 0.3466
kp = 6.931
kdp = 0.06931
K = 40
n = 2
kcatP = 5.9
kmP = 4.6

#intitial condtions
m_PETase0 = 5
p_PETase0 = 0
c_PET = 5

# set time observations
t_max = 1000
t_obs = np.linspace(0,t_max,t_max+1) 

X_vals = [1, 2]

for X in X_vals:
    
    params = [ km, km0, kdm, kp, kdp, K, n, kcatP, kmP, X ]
    s0 = [ m_PETase0, p_PETase0, c_PET0 ]

    # run simulation
    s_obs = odeint(PET_degradation,s0,t_obs,args=(params,))  

    m_PETase_obs = s_obs[:,0]
    p_PETase_obs = s_obs[:,1]
    cPET_obs = s_obs[:,2]
    
    ax.plot(t_obs, p_PETase_obs, label = f"Burden growth factor = {X}X")
ax.legend()
plt.show()
fig.savefig('PETase expression over time', dpi = 300, format = 'png')


# In[51]:


# create figure and plot results

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set(xlabel = 'Time (mins)')
ax.set(ylabel = 'PET concentration (mM)')

km = 30
km0 = 0.03
kdm = 0.3466
kp = 6.931
kdp = 0.06931
K = 40
n = 2
kcatP = 5.9
kmP = 4.6

#intitial condtions
m_PETase0 = 5
p_PETase0 = 0
c_PET = 5

# set time observations
t_max = 1000
t_obs = np.linspace(0,t_max,t_max+1) 

X_vals = [1, 2]

for X in X_vals:
    
    params = [ km, km0, kdm, kp, kdp, K, n, kcatP, kmP, X ]
    s0 = [ m_PETase0, p_PETase0, c_PET0 ]

    # run simulation
    s_obs = odeint(PET_degradation,s0,t_obs,args=(params,))  

    m_PETase_obs = s_obs[:,0]
    p_PETase_obs = s_obs[:,1]
    cPET_obs = s_obs[:,2]
    
    ax.plot(t_obs, cPET_obs, label = f"Burden growth factor = {X}X")
ax.legend()
plt.show()
fig.savefig('PET degradation over time', dpi = 300, format = 'png')


# In[ ]:




