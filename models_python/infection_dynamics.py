import os
import sys
sys.path.append('..')
from lib.mppaper import *
from lib.simulation_objects import *
import lib.mpsetup as mpsetup
import shutil
from scipy.integrate import solve_ivp
import matplotlib


dir_io=sys.argv[1] # directory with input files

dir_out_tot='{inp}/data'.format(inp=dir_io) # directory with output plots
dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots

if not os.path.exists(dir_out_plots_tot):
    os.makedirs(dir_out_plots_tot) 
    
      

if not os.path.exists(dir_out_tot):
    os.makedirs(dir_out_tot) 
    
      
# LOAD PARAMETERS AND INITIALIZE SIMULATION
    
param_file='{inp}/params.txt'.format(inp=dir_io)

pars = {}

                    # ~ echo "# 1 <p>  2 <mu>  3 <K_D> 4 <I> 5 <phi> 6 <omega>  7 <model>    "  >> $input


params = np.genfromtxt(param_file, dtype="f8, |S10, |S10, |S15, |S10, f8, |U20, |U20", names=['p','mu','K_D', 'I', 'phi', 'omega', 'model', 'IC_mode'], encoding='UTF-8')

model  =str(params['model']) # cocktail, singlephage
IC_mode  =str(params['IC_mode']) # rightphage, wrongphage, randomized_right
print(model)
print('r1r2' in model)


pars['model']  =model
pars['IC_mode']  =IC_mode


pars['p']  =float(params['p'])
pars['mu'] =float(params['mu'])
pars['K_D']=float(params['K_D'])
pars['I']  =float(params['I'])
pars['omega']  =float(params['omega'])
pars['phi']  =float(params['phi'])

print(params)
print((params.shape))
print(pars)
print(isinstance(pars, dict))

print((pars['phi']))
    
#bacterial growth parameters
pars['r']      = 1. # growth rate
pars['K_bact'] =  10.**10. # bacts carrying capacity
pars['d']      = pars['r']/pars['K_bact'] # density dependent death rate per individual
print((pars['d']))
    
#phage growth parameters
pars['lambda'] = 100. #burst size
#bacteria-immune interaction parameters
pars['kappa'] = 8.2 * 10.**-8.   # immune killing rate     

S_eq=pars['omega']/(pars['lambda']* pars['phi'])
print(S_eq)

P_eq = (pars['r'] * (1 - S_eq/pars['K_bact']) - pars['kappa']*pars['I'] /(1 + S_eq/pars['K_D']))/pars['phi']
print(P_eq)

file_out='{inp}/allparameters.json'.format(inp=dir_io)

with open(file_out,'w') as f_handle:
    json.dump(pars, f_handle)
    

dt_int = 0.01 # hours, maximum integration time step

tfin=200 # end time in hours
tf=tfin # end time in hours

S    = []
R    = []
P1   = []
P2   = []

if IC_mode == 'rightphage':  
    
      
    P1_0 = 4.*pars['K_bact'] # High IC should be preferred, otherwise always cleared abobe I_bif
    P1_0 = min([P_eq*10.,P1_0])
    P1_0 = max([P1_0,100.])
    P2_0 = 0.
    S_0  = 0.1*pars['K_bact']
    S_0  = min([S_0, S_eq*10.])
    S_0  = max([S_0,100.])
    R_0  = 10.
    
elif IC_mode == 'rightphage_low':  
    
    
    P1_0 = 10**7. #
    P2_0 = 0.
    S_0  = 4*10**8.
    R_0  = 10.
    
else:
    print("IC not implemented")
    sys.exit()

if  'cocktail' in model:
    P1_0 = P1_0/2.
    P2_0 = P1_0
   
    
   
if model == 'nophage':
    S_0  = 1.*pars['K_bact']
    R_0  = 0.
    P1_0 = 0.
    P2_0 = 0.
   

data_fin = np.array([S_0, R_0, P1_0, P2_0 ])

if  'r1r2' in model:
    R2_0 = R_0
   
    data_fin = np.array([S_0, R_0,R2_0, P1_0, P2_0 ])
 
file_out='{inp}/initial_conditions.dat'.format(inp=dir_out_tot)

with open(file_out,'a') as f_handle:
    np.savetxt(f_handle, data_fin.T, fmt='%40.15f',  header=" S \t R \t P1 \t P2")
    f_handle.write("\n")

#integrate equations

derivatives_pass = derivatives(pars)

bact_ext.terminal = True
# ~ bact_ext.terminal = False
bact_ext.direction = -1
bact_ext_r1r2.terminal = True
bact_ext_r1r2.direction = -1


if  'r1r2' in model:
    sol = solve_ivp(derivatives_pass, [0, tf], [S_0, R_0, R2_0, P1_0, P2_0],  method='RK45', dense_output=True, events=bact_ext_r1r2, max_step=dt_int, rtol= 1e-8, atol=1e-8)

else:
    sol = solve_ivp(derivatives_pass, [0, tf], [S_0, R_0, P1_0, P2_0],  method='RK45', dense_output=True, events=bact_ext, max_step=dt_int, rtol= 1e-8, atol=1e-8)


print((sol.t))

tf=int(sol.t[-1])

print(tf)

#evaluate solution

t_eval= np.linspace(0, tf, tf+1, endpoint=True)
z = sol.sol(t_eval)

print((t_eval.shape))
print((z.shape))
print(t_eval)

if  'r1r2' in model:
   
    S, R, R2, P1, P2 = z

else:
    S, R, P1, P2 = z
    R2 = np.zeros_like(R)

print((S.shape))

# store output 

data_fin = np.array([t_eval,S, R, R2, P1, P2 ])
 
file_out='{inp}/infection_time_series.dat'.format(inp=dir_out_tot)

with open(file_out,'a') as f_handle:
    np.savetxt(f_handle, data_fin.T, fmt='%40.15f',  header="time \t S \t R1 \t R2 \t P1 \t P2")
    f_handle.write("\n")

timelast=sol.t[-1]
B_tot=S+R+R2

print((timelast, tfin)) 
print((S.size, B_tot.size , S.size/2))

if timelast<tfin - 24.:
    bact_tot  = 0
    R_tot     = 0
    phage_tot = 0
    
    bact_mean  = 0
    R_mean     = 0
    phage_mean = 0
    bact_min  = 0
    phage_min = 0
    
    bact_fin       =0
    phage_fin      =0
    frac_resist_fin=0
    
else:
    
    timeidx= np.argmax(t_eval>tfin - 24.)
    bact_tot =  np.sum(B_tot[timeidx:])
    R_tot =  np.sum(R[timeidx:] + R2[timeidx:])
    phage_tot = np.sum(P1[timeidx:]+P2[timeidx:])
    
    bact_mean =  np.mean(B_tot[timeidx:])
    R_mean =  np.mean(R[timeidx:] + R2[timeidx:])
    phage_mean = np.mean(P1[timeidx:]+P2[timeidx:])
    
    bact_min =  np.amin(B_tot[timeidx:])
    phage_min =  np.amin(P1[timeidx:]+P2[timeidx:])
   
    bact_fin       =  B_tot[-1]
    phage_fin      =  P1[-1]+P2[-1]
    
    frac_resist_fin=float(R[-1] + R2[-1])/bact_fin


              
data_out = np.array([bact_mean, R_mean, phage_mean, timelast, bact_fin, phage_fin, frac_resist_fin, bact_min, phage_min, bact_tot, R_tot, phage_tot])

print(bact_min)


file_out='{inp}/global_features.txt'.format(inp=dir_out_tot)

with open(file_out,'a') as f_handle:
    np.savetxt(f_handle, data_out, fmt='%15.15f', newline=" ")
    f_handle.write("\n")


    
    
# plot time series
      
thisfigsize = figsize
thisfigsize[1] *= 0.75



fig = plt.figure(figsize=(thisfigsize[0]*3./3, thisfigsize[0]/2.))
grid = gridspec.GridSpec(1, 1, left=0.15, right=0.7, top=0.91, bottom=0.22, wspace=0.4, hspace=0.35)
labeled_axes = []
    
ax = plt.Subplot(fig, grid[0, 0])
fig.add_subplot(ax)
labeled_axes.append(ax)

ax.plot(t_eval, S , linestyle='-', color='r', label='S')
ax.plot(t_eval, R , linestyle='-', color='k', label='R')
ax.plot(t_eval, R2 , linestyle='-', color='grey', label='R2')
ax.plot(t_eval, P1, linestyle='-', color='g', label='P1')
ax.plot(t_eval, P2, linestyle='-', color='b', label='P2')

ax.set_xlabel('time (h)')
ax.set_ylabel('density')
ax.set_yscale('log')
ax.set_ylim(ymin=10**-5)
ax.xaxis.labelpad = axis_labelpad
ax.yaxis.labelpad = axis_labelpad
ax.legend(frameon=False, ncol=1, columnspacing=0.5, handletextpad=0.2,
    loc='upper right', bbox_to_anchor=(1.5, 1.18))

mpsetup.despine(ax) 
    
    
#### finish figure ####
labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
xycoords=('axes fraction'), fontweight = 'bold')
out_file='{out}/pop_densities_vs_time.pdf'.format(out=dir_out_plots_tot)
fig.savefig(out_file)






