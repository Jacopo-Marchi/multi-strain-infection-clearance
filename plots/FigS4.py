
import os

import sys
sys.path.append('..')
from lib.simulation_objects import *





dir_out_plots_tot='../../multistrain_paper_plots' # directory with output plots

if not os.path.exists(dir_out_plots_tot):
    os.makedirs(dir_out_plots_tot) 
    
      
# LOAD PARAMETERS AND INITIALIZE SIMULATION
    
pars = {}


model='cocktail_diffphis_r1r2' # 
IC_mode='rightphage_low' # 
print(model)
print('r1r2' in model)
print('a' in 'abaco')

print(IC_mode)
print(np.__version__)


pars['model']  =model
pars['IC_mode']  =IC_mode


sig=0.2
pars['p']  =1.
pars['mu'] =4.3e-8
pars['K_D']=4.1e7
pars['I']  =1e4
pars['omega']  =0.07
pars['phi']  =5e-7
pars['phi2']  =sig * pars['phi']
# ~ pars['phi2']  = pars['phi']/5.

    

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



dt_int = 0.01 # hours, maximum integration time step

tfin=tf=200 # end time in hours
tf=tfin # end time in hours

S    = []
R    = []
P1   = []
P2   = []

    
P1_0 = 10**7. #
P2_0 = 0.
# ~ S_0  = 0.04*pars['K_bact']
S_0  = 4*10**8.
R_0  = 10.
    

if  'cocktail' in model:
    P1_0 = P1_0/2.
    P2_0 = P1_0



if  'r1r2' in model:
    R2_0 = R_0
   


I_0 = min([2.7*10.**6., pars['I']/5.])    

Ei_1_0= np.zeros(L)
Ei_2_0= np.zeros(L)
    
 
in_conds=np.zeros(2*L+6)
in_conds[:6] = S_0, R_0, R2_0,  P1_0, P2_0, I_0
in_conds[6:6+L] = Ei_1_0
in_conds[6+L:6+2*(L)] = Ei_2_0




 
#integrate equations

derivatives_pass = derivatives_real(pars)
bact_ext_real_pass = bact_ext_gen(pars)

bact_ext_real_pass.terminal = True
bact_ext_real_pass.direction = -1
sol = solve_ivp(derivatives_pass, [0, tf], in_conds,  method='RK45', dense_output=True, events=bact_ext_real_pass, max_step=dt_int, rtol= 1e-8, atol=1e-8)

print((sol.t))

print(tf)

#evaluate solution
# ~ tf=int(sol.t[-1])


t_eval= np.linspace(0, int(sol.t[-1]), int(sol.t[-1])+1, endpoint=True)

   
# ~ if sol.t[-1] < tf:
        
    # ~ z = sol.sol(sol.t[-1])
    
    # ~ print((z.shape))
    

    # ~ S, R, R2, P1, P2, Is = z[:6]

    # ~ in_conds=np.zeros(2*L+6)
    # ~ in_conds[:6] = 0, 0, 0,  P1, P2, Is
    # ~ in_conds[6:6+L] = Ei_1_0
    # ~ in_conds[6+L:6+2*(L)] = Ei_2_0
    
   
    # ~ sol2 = solve_ivp(derivatives_pass,  [sol.t[-1], tf], in_conds,  method='RK45', dense_output=True, max_step=dt_int, rtol= 1e-8, atol=1e-8)

    
        
    # ~ print(tf)
    # ~ print((sol.t[-1]))
    
    # ~ #evaluate solution
    
    # ~ idev=np.searchsorted(t_eval, sol.t[-1])
    
    # ~ t_eval = np.insert(t_eval, idev, sol.t[-1])
    
    # ~ z = np.hstack((sol.sol(t_eval[t_eval<=sol.t[-1]]),sol2.sol(t_eval[t_eval>sol.t[-1]])))
    # ~ print((z.shape))


# ~ else:

z = sol.sol(t_eval)

print((t_eval.shape))
print((z.shape))
print(t_eval)


E_res_i_1=  np.zeros_like(L)
E_res_i_2=  np.zeros_like(L)
E_res2_i_2= np.zeros_like(L)


if  'r1r2' in model:
    S, R, R2, P1, P2, Is = z[:6]

    
    Ei_1= z[6:6+L]
    Ei_2= z[6+L:6+2*(L)]

else:
    S, R, P1, P2, Is = z[:5]
    
    
    Ei_1= z[5:5+L]
    Ei_2= z[5+L:5+2*(L)]
    
    R2 = np.zeros_like(R)
    
    

B_tot=S+R + R2 +  np.sum(Ei_1 + Ei_2 + E_res_i_1 + E_res_i_2+ E_res2_i_2, axis=0 )

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
# ~ ax.set_xlim(0,100)
ax.set_ylim(0.4,5e14)
ax.set_yscale('log')
# ~ ax.legend(loc='upper right',frameon=False)
ax.xaxis.labelpad = axis_labelpad
ax.yaxis.labelpad = axis_labelpad
ax.legend(frameon=False, ncol=2, columnspacing=0.5, handletextpad=0.2,
    loc='upper right', bbox_to_anchor=(1.05, 1.1))

mpsetup.despine(ax) 
    
    
#### finish figure ####
labeldict = dict(labelstyle=r'{\sf \textbf{%s}}', fontsize='medium',
xycoords=('axes fraction'), fontweight = 'bold')
#    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(-0.2, 0.95), **labeldict)
#mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(-0.3, 0.95), **labeldict)
out_file='{out}/FigS4.pdf'.format(out=dir_out_plots_tot)
#    print out_file
fig.savefig(out_file)
#    plt.show()
#    print out_file



