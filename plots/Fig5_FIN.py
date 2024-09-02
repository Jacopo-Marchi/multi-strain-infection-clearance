

import os

import sys
sys.path.append('..')
from lib.simulation_objects import *




def plot_summary_diagram(phi, I, array_z , zlabel):
    
    # ~ rcparams['font.size'] =  14
    # ~ matplotlib.rcParams.update(rcparams)

    array_x=phi
    xlabel=r'$\phi (g / (h \ \rm{PFU}))$'

    array_y=I
    ylabel=r'$I (\rm{cells} / g)$'


    print((array_z.size, array_x.size))
    
    

    print((array_z.size, array_x.size))
    
    fig = plt.figure(figsize=(5,4))
    grid = gridspec.GridSpec(1, 1, left=0.15, right=0.93, top=0.95, bottom=0.15,
             wspace=0.4, hspace=0.35)
    labeled_axes = []
    ax = plt.Subplot(fig, grid[0, 0])
    fig.add_subplot(ax)
    labeled_axes.append(ax)
    
    phase_diagr_all(array_x, array_y, array_z,ax)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
        
    ax.plot(phi_finer,I_LB,color='b',linestyle='--', linewidth=3)
    ax.plot(phi_finer,I_LB2,color='g',linestyle='--', linewidth=2)
    ax.plot(phi_finer,I_LB3,color='gray',linestyle='--', linewidth=2)

    
    I_roach    =[2.4e7, 2.7e6]
    phi_roach  =[5.4e-8, 5.4e-8]
        
    ax.plot(phi_roach, I_roach, marker='D', markersize=4, color='b',linestyle='')

    trans = transforms.blended_transform_factory(ax.transData, ax.get_xticklabels()[0].get_transform())
    ax.annotate(r'$I_s$', xy=(np.unique(phi)[1], I_LB[1]), xycoords='data',xytext=(0, 5), textcoords='offset pixels', color='b')    #
    
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    
    out_file='{out}/Fig5.pdf'.format(out=dir_out_plots_tot)

    
    
    fig.savefig(out_file)#,bbox_inches='tight'


dir_io='../../models_exploration_Iphi_realpars_FIN/singlephage/short_grid_scan_fctphagepars_rightphage_low/summary_analysis_p_0d/' # directory with input files

dir_out_plots_tot='../../multistrain_paper_plots'.format(inp=dir_io) # directory with output plots

if not os.path.exists(dir_out_plots_tot):
    os.makedirs(dir_out_plots_tot) 
    

     
file_in='{inp}/obs_fct_phagepars.txt'.format(inp=dir_io) 
file_in_IC='{inp}/IC_fct_phagepars.txt'.format(inp=dir_io) 

data=[]
with open(file_in) as f:
    lines=f.readlines()
    for line in lines:
        if not line.startswith("#"):
            line = bytes(line, 'utf-8').decode('utf-8', 'ignore')
            myarray = np.fromstring(line, dtype=float, sep=' ')
            data.append(myarray)
            #print(myarray)
    print((len(data)))
    
data = np.asarray(data)


        
data_IC=[]
with open(file_in_IC) as f:
    lines=f.readlines()
    for line in lines:
        if not line.startswith("#"):
            line = bytes(line, 'utf-8').decode('utf-8', 'ignore')
            myarray = np.fromstring(line, dtype=float, sep=' ')
            data_IC.append(myarray)
            #print(myarray)
    print((len(data_IC)))
data_IC = np.asarray(data_IC)
print((data_IC.shape))



I                    =    data[:,0 ]       
mu                   =    data[0,1 ]       
p                    =    data[0,2 ]       
K_D                  =    data[0,3 ]
phi                  =    data[:,4 ]
omega                =    data[0,5 ]
bact_tot             =    data[:,6 ]
R_tot                =    data[:,7 ]
phage_tot            =    data[:,8 ]
timelast             =    data[:,9 ]
bact_fin             =    data[:,10]
phage_fin            =    data[:,11]
frac_resist_fin      =    data[:,12]
bact_min             =    data[:,13]
phage_min            =    data[:,14]

print(I)
print(bact_tot)




S_0            =    data_IC[0,6]
P0            =    data_IC[0,7]

# ~ S_0  = 4*10**8.



# ~ # other simulation parameters 
# ~ r      = 0.75 # growth rate
# ~ K_bact =  10.**10. # bacts carrying capacity
# ~ d     = r/K_bact # density dependent death rate per individual
# ~ lambd = 100. #burst size
# ~ kappa = 8.2 * 10.**-8.   # immune killing rate     
# ~ P_c= 1.5*10.**7.

# other simulation parameters, these are the same for all the grid
pars = {}

file_params='{inp}/allparameters.json'.format(inp=dir_io)

with open(file_params,'r') as f_handle:
    pars = json.loads(f_handle.read())
    
r      = pars['r'] # growth rate
K_bact =  pars['K_bact'] # bacts carrying capacity
lambd = pars['lambda'] #burst size
kappa = pars['kappa']   # immune killing rate    
P_c = pars['P_c']   # immune killing rate    


d     = r/K_bact # density dependent death rate per individual 
print(d)

I_bif= r/kappa

phi_finer=np.geomspace(1e-9, 1e-6, num=1000)
I_suigreatb0= (r - phi_finer*P_c)/kappa * (1 + S_0/np.unique(K_D))
I_LB = np.maximum(I_suigreatb0, I_bif)

I_suigreatb02= (r - phi_finer*P_c)/kappa * (1 + S_0*3/np.unique(K_D))
I_LB2 = np.maximum(I_suigreatb02, I_bif)

I_suigreatb03= (r - phi_finer*P_c)/kappa * (1 + (S_0/3)/np.unique(K_D))
I_LB3 = np.maximum(I_suigreatb03, I_bif)






plot_summary_diagram(phi, I,    bact_tot , 'number of bacteria')
