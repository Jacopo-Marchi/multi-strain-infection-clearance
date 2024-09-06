import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

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
import matplotlib.transforms as transforms



# ~ rcparams['font.size'] =  14
# ~ matplotlib.rcParams.update(rcparams)


dir_out_plots_tot='../../multistrain_paper_plots' # directory with output plots
#dir_out_plots = dir_out_plots.replace(".", "d")

if not os.path.exists(dir_out_plots_tot):
    os.makedirs(dir_out_plots_tot) 
    
      
pars = {}


file_params='../../models_exploration_Iphi_extthr_FIN/singlephage/short_grid_scan_fctphagepars_rightphage/summary_analysis_p_0d/allparameters.json'

with open(file_params,'r') as f_handle:
    pars = json.loads(f_handle.read())


model    =pars['model'] 
IC_mode  =pars['IC_mode']


pars['d']      = pars['r']/pars['K_bact'] # density dependent death rate per individual
      

pars['phi']  =1e-11

Is =[15e6, 29e6, 29e6]

titles=[r'$I<I_0$', r'$I>I_0$', r'$I>I_0$, no phage']


fig, axes = plt.subplots(1,len(Is),figsize=(12,4))
labeled_axes = []
for i in range(len(Is)):

    ax = axes[i]
    labeled_axes.append(ax)

    pars['I']  =Is[i]
    
    S_eq=pars['omega']/(pars['lambda']* pars['phi'])
    print(S_eq)
    
    P_eq = (pars['r'] * (1 - S_eq/pars['K_bact']) - pars['kappa']*pars['I'] /(1 + S_eq/pars['K_D']))/pars['phi']
    print(P_eq)
    
    B_UI=(pars['K_bact'] - pars['K_D'])/2. - np.sqrt(((pars['K_bact'] + pars['K_D'])**2 )/4. - pars['K_bact'] * pars['K_D']*pars['kappa']*pars['I']/pars['r'])

    dt_int = 0.01 # hours, maximum integration time step
    
    tfin=200 # end time in hours
    tf=tfin # end time in hours
    
    S    = []
    R    = []
    P1   = []
    P2   = []
    if IC_mode == 'rightphage':  # Higher IC for visualization purposes
        
        P1_0 = 4*pars['K_bact'] 
        P1_0 = min([P_eq*10.,P1_0])
        P1_0 = max([P1_0,100.])
        P2_0 = 0.
        S_0  = 0.1*pars['K_bact']
        S_0  = min([S_0, S_eq*10.])
        S_0  = max([S_0,100.])
        R_0  = 1.
        
    elif IC_mode == 'rightphage_low':  
   
        P1_0 = 10**7. #
        P2_0 = 0.
        S_0  = 4*10**8.
        R_0  = 10.
    
    else:
        print("IC not implemented")
        sys.exit()
    
    if model == 'cocktail':
        P1_0 = P1_0/2.
        P2_0 = P1_0
       
    if model == 'nophage':
        S_0  = 1.*pars['K_bact']
        R_0  = 0.
        P1_0 = 0.
        P2_0 = 0.
    if i == 2:
        P1_0 = 0.
        P2_0 = 0.
       
    # Run the model
         
    derivatives_pass = derivatives(pars)
    
    bact_ext.terminal = True
    # ~ bact_ext.terminal = False
    bact_ext.direction = -1
    sol = solve_ivp(derivatives_pass, [0, tf], [S_0, R_0, P1_0, P2_0],  method='RK45', dense_output=True, events=bact_ext, max_step=dt_int, rtol= 1e-8, atol=1e-8)
    t_eval= np.linspace(0, tf, tf+1, endpoint=True)
        
    if sol.t[-1] < tf:
            
        z = sol.sol(sol.t[-1])
        
        print((z.shape))
        
        S, R, P1, P2 = z 
        
        sol2 = solve_ivp(derivatives_pass, [sol.t[-1], tf], [0,0, P1, P2],  method='RK45', dense_output=True, max_step=dt_int, rtol= 1e-8, atol=1e-8)
        
        
        print(tf)
        print((sol.t[-1]))
        
        #evaluate solution
        
        idev=np.searchsorted(t_eval, sol.t[-1])
        
        t_eval = np.insert(t_eval, idev, sol.t[-1])
        
        z = np.hstack((sol.sol(t_eval[t_eval<=sol.t[-1]]),sol2.sol(t_eval[t_eval>sol.t[-1]])))
        print((z.shape))
        
    
    
    else:
            
        
        #evaluate solution
        
        z = sol.sol(t_eval)
    
    print((t_eval.shape))
    print((z.shape))
    
    S, R, P1, P2 = z


    
    #Plot the results
    labels=['','','']
       
    if i==0:
        
        labels=['Susceptible bacteria','Resistant bacteria','Phage']
    
    ax.plot(t_eval, S , linestyle='-', color='r', label=labels[0])
    ax.plot(t_eval, R , linestyle='-', color='k', label=labels[1])
    # ~ if i<2:
    ax.plot(t_eval, P1, linestyle='-', color='g', label=labels[2])
    # ~ ax.plot(t_eval, P2, linestyle='-', color='b', label='P2')
    
    ax.set_xlabel('Time (h)')
    ax.set_ylabel(r"Density (g$^{-1}$)", labelpad=20)
    ax.set_xlim(0,100)
    ax.set_ylim(10,5e13)
    ax.set_yscale('log')
    ax.legend(frameon=False, columnspacing=0.5, handletextpad=0.2,
      loc='upper right', bbox_to_anchor=(1.05, 1.05))
    ax.axhline(y=B_UI, color='k', linestyle='--', linewidth=2)
    trans = transforms.blended_transform_factory(
        ax.get_yticklabels()[0].get_transform(), ax.transData)
    ax.text(-0.015,B_UI + 0.2*B_UI, r'$B^U_I$', color="k", transform=trans, 
            ha="center", va="top")

    

    ax.set_title('{title}'.format(title=titles[i]))

    
    # ~ plt.setp(ax.spines.values(),linewidth=1.5)
    ax.tick_params(labelsize=14,direction='in',width=1.5)
    # ~ ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    
    if i<2:
        # ~ ax.set_ylabel("")
        # ~ frame1 = plt.gca()
        # ~ frame1.axes.get_yaxis().set_visible(False)
        # ~ ax.get_yaxis().set_visible(False)
        ax.set_yticklabels([])
        
    if i>0:
        ax.set_ylabel("")
 
out_file='{out}/Fig2.pdf'.format(out=dir_out_plots_tot)

fig.tight_layout()
# ~ fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

fig.savefig(out_file,bbox_inches='tight')

