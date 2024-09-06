
import os

import sys
sys.path.append('..')
from lib.simulation_objects import *




def Trascendental_phi_iterative_Pdyn(phi_it, r, K_bact, omega, lambd   ):
    
    #~ return 0
    count =0
    
    if np.isfinite(phi_it):
        
        # ~ S0_it= K_bact/10
        S0_it= 0.1*omega/(lambd*phi_it)
        if phi_it> (100*omega/(lambd*K_bact)):
            # ~ S0_it = 10*omega/(lambd*phi_it)
            S0_it= 0.1*omega/(lambd*phi_it)

        P0_it=  2*K_bact
        
        count+=1
        
        phi_next= ((np.log(S0_it/10))*omega/r + 2  + (np.log(phi_it*P0_it/r)))*r/(2*P0_it)
        print((np.log(phi_it*P0_it/r)))
        
        if 5*(1- omega/(lambd*K_bact*phi_it))*r/phi_it < P0_it:
            
            print ("next piece fct")
            
            P0_it=  5*(1- omega/(lambd*K_bact*phi_it))*r/phi_it
            
            phi_next=10*omega/(lambd*K_bact* (8 - (np.log(S0_it/10)*omega/r) - np.log(phi_it*P0_it/r)))
        
        
        
        print (count, phi_it, phi_next)
        
        if np.abs( phi_it - phi_next) < 1e-13:
            print("found root")
            
            print (count, phi_it, phi_next)

            return phi_next
            
            print (count, phi_it, phi_next)

        
        else:
            return Trascendental_phi_iterative_Pdyn(phi_next, r, K_bact, omega, lambd)
            
        return
            




def Trascendental_phi_Newton_Pdyn(phi_it, r, K_bact, omega, lambd   ):
    
    
    for i in range(1000):
            
        #~ return 0
        count =0
        
            
        if np.isfinite(phi_it):
            
            # ~ S0_it= K_bact/10
            S0_it= 0.1*omega/(lambd*phi_it)
            P0_it=  2*K_bact
            
            count+=1
        
            # ~ f_prime = 1/phi_it - (2*P0_it)/r
            f_prime = 1/phi_it - (2*P0_it)/r - omega/(r*phi_it)

            if phi_it> (100*omega/(lambd*K_bact)):
                # ~ S0_it = 10*omega/(lambd*phi_it)
                S0_it= 0.1*omega/(lambd*phi_it)
                
                f_prime = 1/phi_it - (2*P0_it)/r - omega/(r*phi_it)
            
            f = (np.log(S0_it/10))*omega/r + 2  + (np.log(phi_it*P0_it/r)) - phi_it *(2*P0_it)/r

            phi_next = phi_it - f/f_prime
           
            if  5*(1- omega/(lambd*K_bact*phi_it))*r/phi_it < P0_it: # this piece sends phi below lower bound
                
                print ("next piece fct")
                
                P0_it=  5*(1- omega/(lambd*K_bact*phi_it))*r/phi_it
                
                f_prime = - omega/(r*phi_it)   + omega/(lambd*K_bact*phi_it*phi_it)*(1/(1-omega/(lambd*K_bact*phi_it)) - 10.)

               
                f_second =  omega/(r*phi_it*phi_it)   - 2*omega/(lambd*K_bact*phi_it*phi_it*phi_it)*( 1/(1-(omega/(lambd*K_bact*phi_it))) - 10.) - ( 1/((1-(omega/(lambd*K_bact*phi_it)))**2) )*(omega/(lambd*K_bact*phi_it*phi_it))**2
                
                # ~ phi_next=10*omega/(lambd*K_bact* (9 - ((np.log(S0_it/10)/np.log(phi_it*P0_it/r))*omega/r)))
              
            
                f = (np.log(S0_it/10))*omega/r + 2  + (np.log(phi_it*P0_it/r)) - phi_it *(2*P0_it)/r
                
                phi_next = phi_it - 2*f*f_prime/(2*(f_prime**2) -  f*f_second )

                print((f_second))
                print((2*f*f_prime/(2*(f_prime**2) -  f*f_second )))
            

              
                
            
            
            print((np.log(phi_it*P0_it/r)))
            print((P0_it))
            print((phi_it*P0_it/r))
            print((f))
            print((f_prime))
            

            
            print (count, phi_it, phi_next)

            
                       
            if phi_next*P0_it/r<1:
                print("Solution below lower bound")
                sys.exit()
            
    
            if np.abs( phi_it - phi_next) < 1e-13:
                print("found root")
                
                print (count, phi_it, phi_next)
    
                return phi_next
                
                print (count, phi_it, phi_next)
    
        
    
            phi_it = phi_next

    # If the solution doesn't converge within the maximum number of iterations
    raise ValueError("Solution did not converge.")

    
        
         
        
        
            




    
def plot_summary_diagram(phi, I, array_z , zlabel, pars):

    array_x=phi
    xlabel=r'$\phi (g / (h \ \rm{PFU}))$'

    array_y=I
    ylabel=r'$I (\rm{cells} / g)$'
   

    print((array_z.size, array_x.size))
    
    
    
    
    
    fig = plt.figure(figsize=(12,8), constrained_layout=False)
    grid = fig.add_gridspec(nrows=3,ncols=6, left=0.07, right=0.95, top=0.96, bottom=0.08,
             wspace=0.3, hspace=0.4)
    ax = fig.add_subplot(grid[0:-1, 1:-1])
    
    
    
    labeled_axes = []
    labeled_axes.append(ax)
    
    phase_diagr_all(array_x, array_y, array_z,ax)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    
    ax.plot(np.unique(phi),np.maximum(I_JTB, I_bif),color='b',linestyle='--', linewidth=3)

    trans = transforms.blended_transform_factory(ax.transData, ax.get_xticklabels()[0].get_transform())
    plt.vlines(x=phi_trasc_Pdyn, ymin= 0, ymax=I_bif,color='g',linestyle='--', linewidth=3)
    ax.annotate(r'$\phi_c$', xy=(phi_trasc_Pdyn, 1e4), xycoords='data',xytext=(8, -1), textcoords='offset pixels', color='g')    #

    ax.annotate(r'$I_s$', xy=(np.unique(phi)[1], I_JTB[1]), xycoords='data',xytext=(0, 20), textcoords='offset pixels', color='b')    #
    
    I_fig3 =[20e6, 30e6, 1e6]
    phi_fig3  =[7e-12, 7e-12, 4e-9]
    # ~ markers  =[r"$\rm{B}$", r"$\rm{C}$", r"$\rm{D}$"]
    markers  =["B", "C","D"]
    for i, m in enumerate(markers):
        # ~ ax.plot(phi_fig3[i], I_fig3[i], marker=m, markersize=7, color='k',linestyle='')
        ax.annotate(m, xy=(phi_fig3[i], I_fig3[i]), xycoords='data',xytext=(0, 0), textcoords='offset pixels', color='k',fontsize=12,weight='bold', ha='center', va='center')    #,fontsize=10

    ax.set_xscale("log")
    ax.set_yscale("log")
    
    
    # PLOT EXAMPLE DYNAMICS
    # PLOT EXAMPLE DYNAMICS
    
    model    =pars['model'] 
    IC_mode  =pars['IC_mode'] 
    
    
    for i in range(len(I_fig3)):
        
        # ~ ax = fig.add_subplot(3,3,7+i)
        ax = fig.add_subplot(grid[2, 2*i:2*i+2])


        # ~ ax = axes[i]
        labeled_axes.append(ax)
    
        pars['I']   =I_fig3[i]
        pars['phi'] =phi_fig3[i]
        
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
           
        
        
        
        if  'r1r2' in model:
            R2_0 = R_0
           
         
           
        # Run the model
        derivatives_pass = derivatives(pars)
        
        bact_ext.terminal = True
        bact_ext.direction = -1
        bact_ext_r1r2.terminal = True
        bact_ext_r1r2.direction = -1
        
        
        if  'r1r2' in model:
            sol = solve_ivp(derivatives_pass, [0, tf], [S_0, R_0, R2_0, P1_0, P2_0],  method='RK45', dense_output=True, events=bact_ext_r1r2, max_step=dt_int, rtol= 1e-8, atol=1e-8)
        
        else:
            sol = solve_ivp(derivatives_pass, [0, tf], [S_0, R_0, P1_0, P2_0],  method='RK45', dense_output=True, events=bact_ext, max_step=dt_int, rtol= 1e-8, atol=1e-8)

        t_eval= np.linspace(0, tf, tf+1, endpoint=True)
            
        if sol.t[-1] < tf:
                
            z = sol.sol(sol.t[-1])
            
            print((z.shape))
            
    
            
            if  'r1r2' in model:
               
                S, R, R2, P1, P2 = z
                sol2 = solve_ivp(derivatives_pass,  [sol.t[-1], tf], [0, 0, 0, P1, P2],  method='RK45', dense_output=True, max_step=dt_int, rtol= 1e-8, atol=1e-8)

            
            else:
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
        
        
        if  'r1r2' in model:
           
            S, R, R2, P1, P2 = z
        
        else:
            S, R, P1, P2 = z
            R2 = np.zeros_like(R)    
                
        #Plot the results
           
        ax.plot(t_eval, S , linestyle='-', color='r', label='Susceptible bacteria')
        ax.plot(t_eval, R , linestyle='-', color='k', label='Resistant bacteria 1')
        # ~ ax.plot(t_eval, R , linestyle='-', color='grey', label='Resistant bacteria 2')
        ax.plot(t_eval, P1, linestyle='-', color='g', label='Phage 1')
        # ~ ax.plot(t_eval, P2, linestyle='-', color='b', label='Phage 2')
        
        ax.set_xlabel('Time (h)')
        ax.set_ylabel(r"Density (g$^{-1}$)", labelpad=20)
        ax.set_xlim(0,100)
        ax.set_ylim(10,5e13)
        ax.set_yscale('log')
        if i==2:
            ax.legend(loc='upper right',frameon=False)
        
        if B_UI>0:
            ax.axhline(y=B_UI, color='k', linestyle='--', linewidth=2)
            trans = transforms.blended_transform_factory(
                ax.get_yticklabels()[0].get_transform(), ax.transData)
            ax.text(-0.02,B_UI-B_UI*0.2, r'$B^U_I$', color="k", transform=trans, 
                ha="center", va="center")
    
        
    
        ax.tick_params(direction='in',width=1.5)
        
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
            
        
    labeldict = dict(labelstyle=r'%s', fontsize=20,
                 xycoords=('axes fraction'), fontweight = 'bold')
    mpsetup.label_axes([labeled_axes[0]], labels='A', xy=(0.05,  0.89), **labeldict)
    mpsetup.label_axes([labeled_axes[1]], labels='B', xy=(0.05,  0.89), **labeldict)
    mpsetup.label_axes([labeled_axes[2]], labels='C', xy=(0.05,  0.89), **labeldict)
    mpsetup.label_axes([labeled_axes[3]], labels='D', xy=(0.05,  0.89), **labeldict)
    
    out_file='{out}/Fig4.pdf'.format(out=dir_out_plots_tot)
    # ~ fig.tight_layout()
    
    
    fig.savefig(out_file, dpi = 300)#,bbox_inches='tight'


#### subfigure A ####
## import data


model='cocktail_r1r2'    

dir_io='../../models_exploration_Iphi_extthr_FIN/{model}/short_grid_scan_fctphagepars_rightphage/summary_analysis_p_1d/'.format(model=model)  # directory with input files

dir_out_plots_tot='../../multistrain_paper_plots'.format(inp=dir_io, model=model) # directory with output plots

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
print((data.shape))


# ~ data_IC=[]
# ~ with open(file_in_IC) as f:
    # ~ lines=f.readlines()
    # ~ for line in lines:
        # ~ if not line.startswith("#"):
            # ~ line = bytes(line, 'utf-8').decode('utf-8', 'ignore')
            # ~ myarray = np.fromstring(line, dtype=float, sep=' ')
            # ~ data_IC.append(myarray)
            # ~ #print(myarray)
    # ~ print((len(data_IC)))
# ~ data_IC = np.asarray(data_IC)
# ~ print((data_IC.shape))



        
# ~ print data

# ~ echo "# 1 I " " 2 mu " " 3 p" " 4 K_D "  " 5 phi "  " 6 omega " " 7 bact_tot "  " 8 R_tot "  " 9 phage_tot "  " 10 timelast "  " 11 bact_fin "   " 12 phage_fin "  " 13 frac_resist_fin "  " 14 bact_min "  " 15 phage_min "   > $file_out


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

# other simulation parameters, these are the same for all the grid
pars = {}

file_params='{inp}/allparameters.json'.format(inp=dir_io)

with open(file_params,'r') as f_handle:
    pars = json.loads(f_handle.read())
    
r      = pars['r'] # growth rate
K_bact =  pars['K_bact'] # bacts carrying capacity
lambd = pars['lambda'] #burst size
kappa = pars['kappa']   # immune killing rate    


d     = r/K_bact # density dependent death rate per individual 
print(d)

I_JTB= (1 + omega/(lambd*np.unique(phi)*K_D))**2. * d*K_D/kappa
I_bif= r/kappa
print((I_JTB, I_bif))
S_eq=omega/(lambd*np.unique(phi))
print(S_eq)
I_cocktail=(1 + 2*omega/(lambd*np.unique(phi)*K_D))**2. * (d+ r*np.unique(mu)/S_eq)*K_D/kappa
print(I_cocktail)

B_M = K_D*(np.sqrt(kappa*I*K_bact/(r*K_D))-1) 

phitrans=np.unique(omega)/(lambd*B_M)



# ~ S0            =    data_IC[:,6]
# ~ P0            =    data_IC[:,7]



# solve trascendental



phi_trasc_Pdyn=Trascendental_phi_iterative_Pdyn(r/(2*K_bact)*10, r, K_bact, omega, lambd   )
phi_trasc_Pdyn_Newt=Trascendental_phi_Newton_Pdyn(r/(2*K_bact), r, K_bact, omega, lambd   )
    
    
print(phi_trasc_Pdyn)
print(phi_trasc_Pdyn_Newt)

    
# ~ sys.exit()
    


plot_summary_diagram(phi, I,    bact_tot , 'number of bacteria', pars)
