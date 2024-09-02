import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy.integrate import solve_ivp
import json
from scipy.interpolate import interpn
import matplotlib.cm as cm
import shutil
import matplotlib.colors as colors
# ~ from matplotlib.mlab import griddata
from scipy.interpolate import griddata
from lib.mppaper import *
import lib.mpsetup as mpsetup


import matplotlib.transforms as transforms
import cmasher as cmr


# ~ def phage_saturation(P):
    
    

def gridint(x, y, z, resX=1000, resY=1000):
    "Convert 3 column data to matplotlib grid"
    xi = np.geomspace(min(x), max(x), resX)
    yi = np.geomspace(min(y), max(y-0.01), resY)
    
    print((min(xi), min(yi)))
    print((max(xi), max(yi), max(y)))
    Z = interpn((x, y), z, (xi[None,:], yi[:,None]), method='nearest')
    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z
    
    
def phase_diagr_all(array_x, array_y, array_z,ax):    
    ext_thr =10
    
    if  np.any(array_z>ext_thr):
    
        mincol=1e4
        if mincol <=ext_thr:
            mincol=ext_thr 
        maxcol=np.amax(array_z)
        
        norm=matplotlib.colors.LogNorm(vmin=mincol, vmax=maxcol)
        
        xx = np.reshape(array_x, (np.unique(array_x).size,np.unique(array_y).size))
        yy = np.reshape(array_y, (np.unique(array_x).size,np.unique(array_y).size))
        
        # ~ print(xx[0,:])
        # ~ print(xx[:,0])
        # ~ print(yy[0,:])
        # ~ sys.exit()
        
        sortrows=np.argsort(xx[:,0])
        sortcols=np.argsort(yy[0,:])
        xx=xx[sortrows,:]
        xx=xx[:,sortcols]
        yy=yy[sortrows,:]
        yy=yy[:,sortcols]
        # ~ print(xx[0,:])
        # ~ print(xx[:,0])
        # ~ print(yy[0,:])
        # ~ sys.exit()
  
        
        arrplot =  np.reshape(array_z, (np.unique(array_x).size,np.unique(array_y).size))
        arrplot=arrplot[sortrows,:]
        arrplot=arrplot[:,sortcols]
        
        X, Y, Z = gridint(np.unique(array_x), np.unique(array_y), arrplot)
        
        cmap = cmr.get_sub_cmap('Reds', 0.3, 0.7)
        heatmap = ax.pcolormesh(xx, yy, arrplot, cmap=cmap, norm=norm,linewidth=0,rasterized=True,shading='auto') # faster than pcolor , antialiased=False shading='gouraud'
        # ~ heatmap = ax.pcolormesh(X, Y, Z, cmap=cmap, norm=norm,linewidth=0,rasterized=True) # faster than pcolor , antialiased=False shading='gouraud'

        
        
   
        cbar = plt.colorbar(heatmap) 
        cbar.set_label(r"Bacteria density $(\rm{CFU} / g)$", rotation=270, labelpad=+14)    

    ax.set_ylim(ymin=np.amin(array_y)*0.8, ymax=np.amax(array_y)*1.5)
    ax.set_xlim(xmin=np.amin(array_x)*0.8, xmax=np.amax(array_x)*1.5)
    
    
    ax.xaxis.labelpad = axis_labelpad
    ax.yaxis.labelpad = axis_labelpad
    
    # ~ return ax
    
    
    

def derivatives(pars_pass):
    pars=pars_pass


    p      =  pars['p']  
    mu     =  pars['mu'] 
    K_D    =  pars['K_D']
    I      =  pars['I']  
    phi    = pars['phi']
    r      = pars['r']       # growth rate
    # ~ K_bact = pars['K_bact']  # bacts carrying capacity
    d      = pars['d']       # density dependent death rate per individual
    omega  = pars['omega']     # decay rate
    lambd  = pars['lambda']    # burst size
    kappa  = pars['kappa']     # immune killing rate     
        
            
    def derivatives_fin(t, z):
        S, R, P1, P2 = z
        P2=0
        
        
            
            
        dSdt= r * S - d*S*(S+R) + mu*r*(R-S) - (kappa*I*S/(1 + (S+R)/K_D)) - phi*P1*S 
        dRdt= r * R - d*R*(S+R) + mu*r*(S-R) - (kappa*I*R/(1 + (S+R)/K_D)) - p*phi*P1*R
        dP1dt= lambd*phi*P1*S+ lambd*p*phi*P1*R - omega*P1
        dP2dt= 0
    
        
        return [dSdt, dRdt,dP1dt, dP2dt]
        
    def derivatives_fin_cocktail(t, z):
        S, R, P1, P2 = z
        
            
    
        dSdt= r * S - d*S*(S+R) + mu*r*(R-S) - (kappa*I*S/(1 + (S+R)/K_D)) - phi*P1*S 
        dRdt= r * R - d*R*(S+R) + mu*r*(S-R) - (kappa*I*R/(1 + (S+R)/K_D)) - (1-p)*phi*P2*R - p*phi*P1*R
        dP1dt= lambd*phi*P1*S + lambd*p*phi*P1*R  - omega*P1
        dP2dt= (1-p)*lambd*phi*P2*R  - omega*P2
    
    
        
        return [dSdt, dRdt,dP1dt, dP2dt]
        
        
        
    
    def derivatives_fin_cocktail_r1r2(t, z):
    
        # ~ print z.shape
    
        S, R1, R2, P1, P2 = z
        
        
        B_tot=S+R1 + R2 # total number of bacteria
        
        
        dSdt= r * S - d*S*(B_tot) + mu*r*(R1 +R2-S) - (kappa*I*S/(1 + (B_tot)/K_D)) - phi*P1*S - phi*P2*S 
        dR1dt= r * R1 - d*R1*(B_tot) + mu*r*(S-R1) - (kappa*I*R1/(1 + (B_tot)/K_D))  - p*phi*P2*R1
        dR2dt= r * R2 - d*R2*(B_tot) + mu*r*(S-R2) - (kappa*I*R2/(1 + (B_tot)/K_D)) - (1-p)*phi*P2*R2 - phi*P1*R2
        dP1dt= lambd*phi*P1*S + lambd*phi*P1*R2  - omega*P1
        dP2dt= lambd*phi*P2*S +  (1-p)*lambd*phi*P2*R2 + p*lambd*phi*P2*R1  - omega*P2
        
        
        
        
        return [dSdt, dR1dt,dR2dt,dP1dt, dP2dt]
        
    model=pars['model']
    print(model)
    if model == 'singlephage' or model == 'nophage':
        return derivatives_fin
    elif model == 'cocktail':
        return derivatives_fin_cocktail
    elif model == 'cocktail_r1r2':
        return derivatives_fin_cocktail_r1r2
    else:
        print("model not implemented")
        sys.exit()
    
    
    
def bact_ext(t, z): return z[0] + z[1] - 10. # threshold at 10 /mL. If the lung volume is 100 ul that corresponds to 1 cell. Perhaps it's 400 uL
def bact_ext_r1r2(t, z): return z[0] + z[1]+ z[2] - 10. # threshold at 10 /mL. If the lung volume is 100 ul that corresponds to 1 cell. Perhaps it's 400 uL




def derivatives_real(pars_pass): # more realistic model with infected states and phage-saturation
    pars=pars_pass

    p      =  pars['p']  
    mu     =  pars['mu'] 
    # ~ print mu
    # ~ print pars['mu'] 
    K_D    =  pars['K_D']
    I      =  pars['I']  
    phi    = pars['phi']
    r      = pars['r']       # growth rate
    # ~ K_bact = pars['K_bact']  # bacts carrying capacity
    d      = pars['d']       # density dependent death rate per individual
    omega  = pars['omega']     # decay rate
    lambd  = pars['lambda']    # burst size
    kappa  = pars['kappa']     # immune killing rate     
    L    = pars['L']     # infection stages     
    P_c  = pars['P_c']
    eta= pars['eta'] # infection advancement rate     

    alpha      =  pars['alpha']  
    K_N      =  pars['K_N']  


        
    
    def bact_dyn(same, neighbors, tot_active, phage_effect, immune_killing):
        return r * same - d*same*(tot_active) + mu*r*(neighbors - same)  - phage_effect*same - (immune_killing*same/(1 + (tot_active)/K_D))
        
    def infection_stages_dyn(same, susceptibles, adv_rate, phage_effect):
    
        dEi_dt= np.zeros_like(same)
        dEi_dt[0]= phage_effect*(susceptibles ) - adv_rate*same[0]
        dEi_dt[1:]= adv_rate*same[:-1] - adv_rate*same[1:]
        
    
        return dEi_dt
    
    def phage_dyn(same, phage_effect, burst_rate, infected):
        return burst_rate - phage_effect*(infected)- omega*same
    
    def derivatives_fin_phagesat(t, z):
    
        # ~ print z.shape
    
        S, R, P1, P2, Is  = z[:5]
        Ei_1= z[5:5+L]
        Ei_2= np.zeros_like(z[5+L:5+2*(L)])
        
        P2=0
        
        B_tot=S+R +  np.sum(Ei_1 + Ei_2  ) # total number of bacteria
        
        FP1 = P1/(1. + (P1 + P2)/P_c)
        FP2 = P2/(1. + (P1 + P2)/P_c)
        
        # ~ if t< t_inj:
            # ~ FP1 = 0
            # ~ FP2 = 0
            
        dSdt= bact_dyn(S, R , B_tot, phi*FP1 , kappa*I)
        dRdt= bact_dyn(R, S , B_tot, p*phi*FP1, kappa*I)
        dP1dt= phage_dyn(P1, phi*FP1, lambd*eta*L*(Ei_1[-1]), S + p*R )
        dP2dt= 0
        
        dIdt= alpha * Is *(1- Is/I) * B_tot /(B_tot + K_N)   
        # ~ if t< t_inj:
            # ~ dP1dt = 0
            # ~ dP2dt = 0
        
        dEi_1dt= infection_stages_dyn(Ei_1, S + p*R  , eta*L, phi*FP1)
        dEi_2dt= np.zeros_like(z[4+L:4+2*(L)])

        ders=np.zeros_like(z)
        ders[:5] = dSdt, dRdt, dP1dt, dP2dt, dIdt
        ders[5:5+L] = dEi_1dt
        ders[5+L:5+2*(L)] = dEi_2dt
       
        
        return ders
        
    
    def derivatives_fin_phagesat_cocktail(t, z):
    
        # ~ print z.shape
    
        S, R, P1, P2,Is  = z[:5]
        Ei_1= z[5:5+L]
        Ei_2= z[5+L:5+2*(L)]
        
        B_tot=S+R +  np.sum(Ei_1 + Ei_2  ) # total number of bacteria
        
        FP1 = P1/(1. + (P1 + P2)/P_c)
        FP2 = P2/(1. + (P1 + P2)/P_c)
        
        # ~ if t< t_inj:
            # ~ FP1 = 0
            # ~ FP2 = 0
            
        dSdt= bact_dyn(S, R , B_tot, phi*FP1 , kappa*I)
        dRdt= bact_dyn(R, S , B_tot, p*phi*FP1 + (1-p)*phi*FP2, kappa*I)
        dP1dt= phage_dyn(P1, phi*FP1, lambd*eta*L*(Ei_1[-1]), S + p*R )
        dP2dt= phage_dyn(P2, phi*FP2, lambd*eta*L*(Ei_2[-1]), (1-p)*R )
        dIdt= alpha * Is *(1- Is/I) * B_tot /(B_tot + K_N)
        # ~ if t< t_inj:
            # ~ dP1dt = 0
            # ~ dP2dt = 0
        
        dEi_1dt= infection_stages_dyn(Ei_1, S + p*R  , eta*L, phi*FP1)
        dEi_2dt= infection_stages_dyn(Ei_2, (1-p)*R  , eta*L, phi*FP2)

        ders=np.zeros_like(z)
        ders[:5] = dSdt, dRdt, dP1dt, dP2dt, dIdt
        ders[5:5+L] = dEi_1dt
        ders[5+L:5+2*(L)] = dEi_2dt
       
        
        return ders
        
    
    
    def derivatives_fin_phagesat_cocktail_r1r2(t, z):
    
        # ~ print z.shape
    
        S, R1, R2, P1, P2,Is  = z[:6]
        Ei_1= z[6:6+L]
        Ei_2= z[6+L:6+2*(L)]
        
        # ~ print(Ei_1.shape)
        # ~ print(Ei_2.shape)
        
        B_tot=S+R1 + R2 +  np.sum(Ei_1 + Ei_2  ) # total number of bacteria
        
        FP1 = P1/(1. + (P1 + P2)/P_c)
        FP2 = P2/(1. + (P1 + P2)/P_c)
        
        # ~ if t< t_inj:
            # ~ FP1 = 0
            # ~ FP2 = 0
            
        dSdt= bact_dyn(S, R1+R2 , B_tot, phi*FP1 + phi*FP2 , kappa*I)
        dR1dt= bact_dyn(R1, S  , B_tot, p*phi*FP2, kappa*I)
        dR2dt= bact_dyn(R2, S  , B_tot,  phi*FP1 +(1-p)*phi*FP2 , kappa*I)
        dP1dt= phage_dyn(P1, phi*FP1, lambd*eta*L*(Ei_1[-1]), S + R2 )
        dP2dt= phage_dyn(P2, phi*FP2, lambd*eta*L*(Ei_2[-1]), S +p*R1 + (1-p)*R2 )
        dIdt= alpha * Is *(1- Is/I) * B_tot /(B_tot + K_N)
        # ~ if t< t_inj:
            # ~ dP1dt = 0
            # ~ dP2dt = 0
        
        dEi_1dt= infection_stages_dyn(Ei_1, S + R2  , eta*L, phi*FP1)
        dEi_2dt= infection_stages_dyn(Ei_2, S + p*R1 + (1-p)*R2  , eta*L, phi*FP2)

        ders=np.zeros_like(z)
        ders[:6] = dSdt, dR1dt, dR2dt, dP1dt, dP2dt, dIdt
        ders[6:6+L] = dEi_1dt
        ders[6+L:6+2*(L)] = dEi_2dt
       
        
        return ders
    
    
    def derivatives_fin_phagesat_cocktail_nores(t, z):
    
        # ~ print z.shape
    
        S, R1, R2, P1, P2,Is  = z[:6]
        R1=R2=0
        Ei_1= z[6:6+L]
        Ei_2= z[6+L:6+2*(L)]
        
        # ~ print(Ei_1.shape)
        # ~ print(Ei_2.shape)
        
        B_tot=S+R1 + R2 +  np.sum(Ei_1 + Ei_2  ) # total number of bacteria
        
        FP1 = P1/(1. + (P1 + P2)/P_c)
        FP2 = P2/(1. + (P1 + P2)/P_c)
        
        # ~ if t< t_inj:
            # ~ FP1 = 0
            # ~ FP2 = 0
            
        dSdt= bact_dyn(S, R1+R2 , B_tot, phi*FP1 + phi*FP2 , kappa*I)
        dR1dt= 0
        dR2dt= 0
        dP1dt= phage_dyn(P1, phi*FP1, lambd*eta*L*(Ei_1[-1]), S + R2 )
        dP2dt= phage_dyn(P2, phi*FP2, lambd*eta*L*(Ei_2[-1]), S +p*R1 + (1-p)*R2 )
        dIdt= alpha * Is *(1- Is/I) * B_tot /(B_tot + K_N)
        # ~ if t< t_inj:
            # ~ dP1dt = 0
            # ~ dP2dt = 0
        
        dEi_1dt= infection_stages_dyn(Ei_1, S + R2  , eta*L, phi*FP1)
        dEi_2dt= infection_stages_dyn(Ei_2, S + p*R1 + (1-p)*R2  , eta*L, phi*FP2)

        ders=np.zeros_like(z)
        ders[:6] = dSdt, dR1dt, dR2dt, dP1dt, dP2dt, dIdt
        ders[6:6+L] = dEi_1dt
        ders[6+L:6+2*(L)] = dEi_2dt
       
        
        return ders
    
    def derivatives_fin_phagesat_cocktail_r1r2_infsink(t, z):
    
        # ~ print z.shape
    
        S, R1, R2, P1, P2,Is  = z[:6]
        Ei_1= z[6:6+L]
        Ei_2= z[6+L:6+2*(L)]
        E_res_i_1= z[6+2*(L):6+3*(L)]
        E_res_i_2= z[6+3*(L):6+4*(L)]
        E_res2_i_2= z[6+4*(L):6+5*(L)]
        
        
        # ~ print(Ei_1.shape)
        # ~ print(Ei_2.shape)
        
        B_tot=S+R1 + R2 +  np.sum(Ei_1 + Ei_2 + E_res_i_1 + E_res_i_2+ E_res2_i_2 ) # total number of bacteria
        
        FP1 = P1/(1. + (P1 + P2)/P_c)
        FP2 = P2/(1. + (P1 + P2)/P_c)
        
        # ~ if t< t_inj:
            # ~ FP1 = 0
            # ~ FP2 = 0
            
        dSdt= bact_dyn(S, R1+R2 , B_tot, phi*FP1 + phi*FP2 , kappa*I)
        dR1dt= bact_dyn(R1, S  , B_tot, p*phi*FP2, kappa*I)
        dR2dt= bact_dyn(R2, S  , B_tot,  phi*FP1 +(1-p)*phi*FP2 , kappa*I)
        dP1dt= phage_dyn(P1, phi*FP1, lambd*eta*L*(Ei_1[-1] + E_res_i_1[-1]), S + R2 + np.sum(Ei_1 + Ei_2 + E_res_i_1 + E_res2_i_2) )
        dP2dt= phage_dyn(P2, phi*FP2, lambd*eta*L*(Ei_2[-1] + E_res_i_2[-1] + E_res2_i_2[-1]), S +p*R1 + (1-p)*R2 + np.sum(Ei_1 + Ei_2 + (1-p)*E_res_i_1 + (1-p)*E_res2_i_2 + p*E_res_i_2) )
        
        
        
        
        dIdt= alpha * Is *(1- Is/I) * B_tot /(B_tot + K_N)
        # ~ if t< t_inj:
            # ~ dP1dt = 0
            # ~ dP2dt = 0
        
        dEi_1dt= infection_stages_dyn(Ei_1, S   , eta*L, phi*FP1)
        dEi_2dt= infection_stages_dyn(Ei_2, S   , eta*L, phi*FP2)
    
        dE_res_i_1dt= infection_stages_dyn(E_res_i_1, R2 , eta*L, phi*FP1) 
        dE_res_i_2dt= infection_stages_dyn(E_res_i_2, p*R1 , eta*L, phi*FP2) 
        dE_res2_i_2dt= infection_stages_dyn(E_res2_i_2, (1-p)*R2 , eta*L, phi*FP2) 
       
        

        ders=np.zeros_like(z)
        ders[:6] = dSdt, dR1dt, dR2dt, dP1dt, dP2dt, dIdt
        ders[6:6+L] = dEi_1dt
        ders[6+L:6+2*(L)] = dEi_2dt
        ders[6+2*(L):6+3*(L)] = dE_res_i_1dt
        ders[6+3*(L):6+4*(L)] = dE_res_i_2dt
        ders[6+4*(L):6+5*(L)] = dE_res2_i_2dt
        
        return ders
        
    
    def derivatives_fin_phagesat_cocktail_diffphis(t, z):
        phi2 = pars['phi2'] 

        
        # ~ print z.shape
    
        S, R1, R2, P1, P2,Is  = z[:6]
        Ei_1= z[6:6+L]
        Ei_2= z[6+L:6+2*(L)]
        
        # ~ print(Ei_1.shape)
        # ~ print(Ei_2.shape)
        
        B_tot=S+R1 + R2 +  np.sum(Ei_1 + Ei_2  ) # total number of bacteria
        
        FP1 = P1/(1. + (P1 + P2)/P_c)
        FP2 = P2/(1. + (P1 + P2)/P_c)
        
        # ~ if t< t_inj:
            # ~ FP1 = 0
            # ~ FP2 = 0
            
        dSdt= bact_dyn(S, R1+R2 , B_tot, phi*FP1 + phi2*FP2 , kappa*I)
        dR1dt= bact_dyn(R1, S  , B_tot, p*phi2*FP2, kappa*I)
        dR2dt= bact_dyn(R2, S  , B_tot,  phi*FP1 +(1-p)*phi2*FP2 , kappa*I)
        dP1dt= phage_dyn(P1, phi*FP1, lambd*eta*L*(Ei_1[-1]), S + R2 )
        dP2dt= phage_dyn(P2, phi2*FP2, lambd*eta*L*(Ei_2[-1]), S +p*R1 + (1-p)*R2 )
        dIdt= alpha * Is *(1- Is/I) * B_tot /(B_tot + K_N)
        # ~ if t< t_inj:
            # ~ dP1dt = 0
            # ~ dP2dt = 0
        
        dEi_1dt= infection_stages_dyn(Ei_1, S + R2  , eta*L, phi*FP1)
        dEi_2dt= infection_stages_dyn(Ei_2, S + p*R1 + (1-p)*R2  , eta*L, phi2*FP2)

        ders=np.zeros_like(z)
        ders[:6] = dSdt, dR1dt, dR2dt, dP1dt, dP2dt, dIdt
        ders[6:6+L] = dEi_1dt
        ders[6+L:6+2*(L)] = dEi_2dt
       
        
        return ders
        
    model=pars['model']
    print(model)
    if model == 'singlephage' or model == 'nophage':
        return derivatives_fin_phagesat
    elif model == 'cocktail':
        return derivatives_fin_phagesat_cocktail
    elif model == 'cocktail_r1r2' or model == 'single_p0d5_r1r2':
        return derivatives_fin_phagesat_cocktail_r1r2
    elif model == 'cocktail_r1r2_infsink':
        return derivatives_fin_phagesat_cocktail_r1r2_infsink
    elif model == 'cocktail_diffphis_r1r2':
        return derivatives_fin_phagesat_cocktail_diffphis
    elif model == 'cocktail_nores_r1r2':
        return derivatives_fin_phagesat_cocktail_nores
    else:
        print("model not implemented")
        sys.exit()
    
    
    
 
def bact_ext_gen(pars_pass): # more realistic model with infected states and phage-saturation
    pars=pars_pass

    p      =  pars['p']  
    mu     =  pars['mu'] 
    # ~ print mu
    # ~ print pars['mu'] 
    K_D    =  pars['K_D']
    I      =  pars['I']  
    phi    = pars['phi']
    r      = pars['r']       # growth rate
    # ~ K_bact = pars['K_bact']  # bacts carrying capacity
    d      = pars['d']       # density dependent death rate per individual
    omega  = pars['omega']     # decay rate
    lambd  = pars['lambda']    # burst size
    kappa  = pars['kappa']     # immune killing rate     
    L    = pars['L']     # infection stages     
    P_c  = pars['P_c']
    eta= pars['eta'] # infection advancement rate     

    alpha      =  pars['alpha']  
    K_N      =  pars['K_N']  
    
    def bact_ext_real(t, z): 
        
            
        S, R, P1, P2, Is = z[:5]
        
        
        Ei_1= z[5:5+L]
        Ei_2= z[5+L:5+2*(L)]
        
        B_tot=S+R +  np.sum(Ei_1 + Ei_2 , axis=0) # total number of bacteria
        
        
        
        return B_tot - 10. # threshold at 10 /mL. If the lung volume is 100 ul that corresponds to 1 cell. Perhaps it's 400 uL
 
        
    def bact_ext_r1r2_real(t, z): 
        
            
        S, R, R2, P1, P2, Is = z[:6]
        
        
        Ei_1= z[6:6+L]
        Ei_2= z[6+L:6+2*(L)]
        
        B_tot=S+R + R2 +   np.sum(Ei_1 + Ei_2 , axis=0) # total number of bacteria
        
        
        
        return B_tot - 10. # threshold at 10 /mL. If the lung volume is 100 ul that corresponds to 1 cell. Perhaps it's 400 uL


    
    def bact_ext_r1r2_infsink(t, z): 
        
            
        S, R1, R2, P1, P2,Is  = z[:6]
        Ei_1= z[6:6+L]
        Ei_2= z[6+L:6+2*(L)]
        E_res_i_1= z[6+2*(L):6+3*(L)]
        E_res_i_2= z[6+3*(L):6+4*(L)]
        E_res2_i_2= z[6+4*(L):6+5*(L)]
        
        
        # ~ print(Ei_1.shape)
        # ~ print(Ei_2.shape)
        
        B_tot=S+R1 + R2 +  np.sum(Ei_1 + Ei_2 + E_res_i_1 + E_res_i_2+ E_res2_i_2, axis=0 ) # total number of bacteria
        
        
        return B_tot - 10. # threshold at 10 /mL. If the lung volume is 100 ul that corresponds to 1 cell. Perhaps it's 400 uL
    



       
    
    model=pars['model']
    print(model)
    if model == 'singlephage' or model == 'nophage':
        return bact_ext_real
    elif model == 'cocktail':
        return bact_ext_real
    elif model == 'cocktail_r1r2' or model == 'single_p0d5_r1r2' or model == 'cocktail_nores_r1r2':
        return bact_ext_r1r2_real
    elif model == 'cocktail_r1r2_infsink':
        return bact_ext_r1r2_infsink
    elif model == 'cocktail_diffphis_r1r2':
        return bact_ext_r1r2_real
    else:
        print("model not implemented")
        sys.exit()


    
