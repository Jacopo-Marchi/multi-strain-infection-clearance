import os
import sys
sys.path.append('..')
from lib.mppaper import *
from lib.simulation_objects import *
import lib.mpsetup as mpsetup
import shutil
from scipy.integrate import solve_ivp
import matplotlib
import matplotlib.gridspec as gridspec



dir_io=sys.argv[1] # directory with input files

dir_out_tot='{inp}/data'.format(inp=dir_io) # directory with output plots
dir_out_plots_tot='{inp}/plots'.format(inp=dir_io) # directory with output plots
#dir_out_plots = dir_out_plots.replace(".", "d")

if not os.path.exists(dir_out_plots_tot):
    os.makedirs(dir_out_plots_tot) 
    
      

if not os.path.exists(dir_out_tot):
    os.makedirs(dir_out_tot) 
    
      
# LOAD PARAMETERS AND INITIALIZE SIMULATION
    
param_file='{inp}/params.txt'.format(inp=dir_io)

pars = {}

                    # ~ echo "# 1 <p>  2 <mu>  3 <K_D> 4 <I> 5 <phi> 6 <omega>  7 <model>    "  >> $input


params = np.genfromtxt(param_file, dtype="f8, |S10, |S10, |S15, |S10, f8, |U30, |U20", names=['p','mu','K_D', 'I', 'phi', 'omega', 'model', 'IC_mode'], encoding='UTF-8')

model  =str(params['model']) # cocktail, singlephage
IC_mode  =str(params['IC_mode']) # rightphage, wrongphage, randomized_right
print(model)
print('r1r2' in model)
print('a' in 'abaco')

print(IC_mode)
print(np.__version__)


pars['model']  =model
pars['IC_mode']  =IC_mode


sig=1
pars['p']  =float(params['p'])
if pars['p'] ==-1:
    pars['model']  =model= 'cocktail_nores_r1r2'
elif pars['p'] <0:
    pars['model']  =model= 'cocktail_diffphis_r1r2'
    sig=-pars['p']
    pars['p']  =1.




pars['mu'] =float(params['mu'])
pars['K_D']=float(params['K_D'])
pars['I']  =float(params['I'])
pars['omega']  =float(params['omega'])
pars['phi']  =float(params['phi'])
pars['phi2']  =sig * pars['phi']
# ~ pars['phi2']  = pars['phi']/5.

print(params)
print((params.shape))
print((pars['phi']))
print((pars['I'])) 
    
    

pars['K_N']   = 1e7
pars['alpha'] = 0.97

#bacterial growth parameters
pars['r']      = 0.75 # growth rate
pars['K_bact'] =  10.**10. # bacts carrying capacity
pars['d']      = pars['r']/pars['K_bact'] # density dependent death rate per individual
print((pars['d']))
    
#phage growth parameters
pars['lambda'] = 100. #burst size
#bacteria-immune interaction parameters
pars['kappa'] = 8.2 * 10.**-8.   # immune killing rate     
    
print((pars['I'] * pars['kappa'])) 

    
pars['eta']  =2.
pars['L']  =10
L=pars['L']
pars['P_c']  = 1.5*10.**7. # phage saturation density 


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
    # ~ S_0  = 0.04*pars['K_bact']
    S_0  = 4*10**8.
    R_0  = 10.
    
else:
    print("IC not implemented")
    sys.exit()

if  'cocktail' in model:
    P1_0 = P1_0/2.
    P2_0 = P1_0

if  model == 'single_p0d5_r1r2':
    P2_0 = P1_0
    P1_0 = 0.



if  'nores' in model:
    R_0=0



if  'r1r2' in model:
    R2_0 = R_0
   
if model == 'nophage':
    S_0  = 1.*pars['K_bact']
    R_0  = 0.
    P1_0 = 0.
    P2_0 = 0.
   


I_0 = min([2.7*10.**6., pars['I']/5.])    

Ei_1_0= np.zeros(L)
Ei_2_0= np.zeros(L)
    
in_conds=np.zeros(2*L+5)
in_conds[:5] = S_0, R_0, P1_0, P2_0, I_0
in_conds[5:5+L] = Ei_1_0
in_conds[5+L:5+2*(L)] = Ei_2_0

data_fin = np.array([S_0, R_0, P1_0, P2_0 ])

if  'r1r2' in model:
 
    in_conds=np.zeros(2*L+6)
    in_conds[:6] = S_0, R_0, R2_0,  P1_0, P2_0, I_0
    in_conds[6:6+L] = Ei_1_0
    in_conds[6+L:6+2*(L)] = Ei_2_0
    data_fin = np.array([S_0, R_0, R2_0, P1_0, P2_0 ])
    
 

if  'infsink' in model:
 
    in_conds=np.zeros(5*L+6)
    in_conds[:6] = S_0, R_0, R2_0,  P1_0, P2_0, I_0
    in_conds[6:6+L] = Ei_1_0
    in_conds[6+L:6+2*(L)] = Ei_2_0
    in_conds[6+2*L:6+3*(L)] = Ei_2_0
    in_conds[6+3*L:6+4*(L)] = Ei_2_0
    in_conds[6+4*L:6+5*(L)] = Ei_2_0
    data_fin = np.array([S_0, R_0, R2_0, P1_0, P2_0 ])
    
    
 
file_out='{inp}/initial_conditions.dat'.format(inp=dir_out_tot)

with open(file_out,'a') as f_handle:
    np.savetxt(f_handle, data_fin.T, fmt='%40.15f',  header=" S \t R \t P1 \t P2")
    f_handle.write("\n")

#integrate equations

derivatives_pass = derivatives_real(pars)
bact_ext_real_pass = bact_ext_gen(pars)

bact_ext_real_pass.terminal = True
bact_ext_real_pass.direction = -1
sol = solve_ivp(derivatives_pass, [0, tf], in_conds,  method='RK45', dense_output=True, events=bact_ext_real_pass, max_step=dt_int, rtol= 1e-8, atol=1e-8)

print((sol.t))

tf=int(sol.t[-1])

print(tf)

#evaluate solution

t_eval= np.linspace(0, tf, tf+1, endpoint=True)
z = sol.sol(t_eval)

print((t_eval.shape))
print((z.shape))
print(t_eval)


E_res_i_1=  np.zeros_like(L)
E_res_i_2=  np.zeros_like(L)
E_res2_i_2= np.zeros_like(L)



if  'infsink' in model:
    S, R, R2, P1, P2,Is  = z[:6]
    Ei_1= z[6:6+L]
    Ei_2= z[6+L:6+2*(L)]
    E_res_i_1= z[6+2*(L):6+3*(L)]
    E_res_i_2= z[6+3*(L):6+4*(L)]
    E_res2_i_2= z[6+4*(L):6+5*(L)]

elif  'r1r2' in model:
    S, R, R2, P1, P2, Is = z[:6]

    
    Ei_1= z[6:6+L]
    Ei_2= z[6+L:6+2*(L)]

else:
    S, R, P1, P2, Is = z[:5]
    
    
    Ei_1= z[5:5+L]
    Ei_2= z[5+L:5+2*(L)]
    
    R2 = np.zeros_like(R)
    
    

B_tot=S+R + R2 +  np.sum(Ei_1 + Ei_2 + E_res_i_1 + E_res_i_2+ E_res2_i_2, axis=0 )

print(("sum infected"))
print((S))
print((np.sum(Ei_1 + Ei_2, axis=0 )))
print((np.sum(E_res_i_1 + E_res_i_2+ E_res2_i_2, axis=0 )))


print((B_tot.shape))
print((Ei_1.shape))
print((E_res_i_1.shape))
print((S.shape))
print((S.shape))

# store output 

data_fin = np.array([t_eval,S, R, R2, P1, P2, B_tot])
print((data_fin.shape))

 
file_out='{inp}/infection_time_series.dat'.format(inp=dir_out_tot)

with open(file_out,'a') as f_handle:
    np.savetxt(f_handle, data_fin.T, fmt='%40.15f',  header="time \t S \t R1\t R2 \t P1 \t P2")
    f_handle.write("\n")


timelast=sol.t[-1]

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



print(bact_tot)
print(phage_tot)

              
data_out = np.array([bact_mean, R_mean, phage_mean, timelast, bact_fin, phage_fin, frac_resist_fin, bact_min, phage_min, bact_tot, R_tot, phage_tot])

print(bact_min)


file_out='{inp}/global_features.txt'.format(inp=dir_out_tot)

with open(file_out,'a') as f_handle:
    np.savetxt(f_handle, data_out, fmt='%15.15f', newline=" ")
    f_handle.write("\n")


    
    
# plot time series
      
thisfigsize = figsize
thisfigsize[1] *= 0.75



# ~ fig = plt.figure(figsize=(thisfigsize[0]*3./3, thisfigsize[0]/2.))
fig = plt.figure(figsize=(8,4))
grid = gridspec.GridSpec(1, 1, left=0.12, right=0.9, top=0.91, bottom=0.12, wspace=0.4, hspace=0.35)
labeled_axes = []
    
ax = plt.Subplot(fig, grid[0, 0])
fig.add_subplot(ax)
labeled_axes.append(ax)

ax.plot(t_eval, S , linestyle='-', color='r', label='Susceptible bacteria')
ax.plot(t_eval, R , linestyle='-', color='k', label='Resistant bacteria 1')
ax.plot(t_eval, R2 , linestyle='-', color='grey', label='Resistant bacteria 2')
ax.plot(t_eval, P1, linestyle='-', color='g', label='Phage 1')
ax.plot(t_eval, P2, linestyle='-', color='b', label='Phage 2')
ax.plot(t_eval, Is, linestyle='-', color='m', label='Immune system')

       
ax.set_xlabel('Time (h)')
ax.set_ylabel(r"Density (g$^{-1}$)")
ax.set_xlim(0,200)
ax.set_ylim(10,5e13)
ax.set_yscale('log')
ax.legend(loc='upper right',frameon=False)
ax.xaxis.labelpad = axis_labelpad
ax.yaxis.labelpad = axis_labelpad
# ~ ax.legend(frameon=False, ncol=1, columnspacing=0.5, handletextpad=0.2,
    # ~ loc='upper right', bbox_to_anchor=(1.5, 1.18))

mpsetup.despine(ax) 
    
    
#### finish figure ####
labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
xycoords=('axes fraction'), fontweight = 'bold')
#    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
#mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
out_file='{out}/pop_densities_vs_time.pdf'.format(out=dir_out_plots_tot)
#    print out_file
fig.savefig(out_file)
#    plt.show()
#    print out_file






