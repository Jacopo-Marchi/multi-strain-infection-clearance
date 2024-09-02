#!/bin/bash

#GENERATES FIG 5 OF THE MANUSCRIPT

if (( $# != 0 ))
then

    #~ echo "Usage $0, <input_file>  "
    echo "Usage $0,   "
	
else

    
    #~ model='cocktail_r1r2' # 
    #~ model='cocktail_diffphis_r1r2' # 
    #~ model='single_p0d5_r1r2' #
    model='singlephage' # 
    
    IC_mode='rightphage_low' # 

    data_dir_tot=~/Documents/phage_therapy/Jacopo_multistrain_sinergy/models_exploration_Iphi_realpars_FIN/${model}/


 
    if [ ! -d $data_dir_tot ]; then
        mkdir $data_dir_tot
    fi        
    data_dir_tot=${data_dir_tot}/short_grid_scan_fctphagepars_${IC_mode}
    data_dir=${data_dir_tot}/



 
    if [ ! -d $data_dir_tot ]; then
        mkdir $data_dir_tot
    fi        
    py_script_dir='../'plots/ # scripts directory
    script_dir='./' # scripts directory
    program_dir='../'models_python/ # scripts directory
    
    
    curr_dir=`pwd`
    count_tot=0

  

    for p in  0. # 
    do

       
        pstr_tmp="p_${p}" #
                    
        pstr=`echo $pstr_tmp | sed 's/\./d/g' | sed 's/\s\s*/_/g'` # substitutes dots with d, and spaces and tabs with _
        
        dir_out=${data_dir}/summary_analysis_${pstr}/
        file_out=${dir_out}/obs_fct_phagepars.txt
        file_out_IC=${dir_out}/IC_fct_phagepars.txt
        
        rm -fr $dir_out
        mkdir $dir_out
        echo $file_out
        
    
        echo "# 1 I " " 2 mu " " 3 p" " 4 K_D "  " 5 phi "  " 6 omega " " 7 bact_tot "  " 8 R_tot "  " 9 phage_tot "  " 10 timelast "  " 11 bact_fin "   " 12 phage_fin "  " 13 frac_resist_fin "  " 14 bact_min "  " 15 phage_min "   > $file_out
    
    

      
    
        

    #~ for phi in   '2e-7'  '3e-7'  '5e-7' 
    for phi in   '4e-9'  '7e-9'   '1e-8'    '2e-8'    '3e-8'   '4e-8'   '5e-8'  '7e-8' '1e-7'  '4e-7'  '7e-7' '1e-6' '2e-7'  '3e-7'  '5e-7' 
    #~ for phi in   '4e-9'  
    do




    for omega in   '0.07'  # for realparams
    do
    
        for mu in   '4.3e-08'  # '1e-07'    '1e-05'     '1e-03'    '1e-01'  
        do
        
        
            for K_D in  '4.1e07'  # for realparams
            do
                 
                for I in     '1e04'  '1e05'   '1e06'   '2e06'  '5e06' '8e06'   '10e06'    '15e06' '20e06' '30e06'  '40e06'   '60e06'    '80e06'    '100e06'  '2e08' 
                #~ for I in    '1e04'

                do
            
        
        
        
                    param_dir_tmp="p_${p}_mu_${mu}_K_D_${K_D}_I_${I}_phi_${phi}_omega_${omega}" #
                    
                    param_dir_d=`echo $param_dir_tmp | sed 's/\./d/g' | sed 's/\s\s*/_/g'` # substitutes dots with d, and spaces and tabs with _
                    
                    param_dir=`echo $param_dir_d | sed 's/d0\+_/d_/g' | sed 's/d\(0*[123456789][123456789]*\)0\+_/d\1_/g'`
                    
                    
                    data_dir_fin=${data_dir_tot}/${param_dir}/
                    
                    # RUN MODEL, SKIP IF ALREADY GENERATED SYNTH DATA

                    cd $program_dir
                    
                    rm -r $data_dir_fin

                    if [ ! -d $data_dir_fin ]; then
                        mkdir $data_dir_fin
                    fi
                    
                    
                    input=${data_dir_fin}/'params.txt'
                    
                    echo "$data_dir_fin"
                
        
                    rm $input
                    
                    echo "# 1 <p>  2 <mu>  3 <K_D> 4 <I> 5 <phi> 6 <omega>  7 <model>  8 <IC_mode>    "  >> $input
                    echo "$p "  " $mu " " $K_D " " $I " " $phi " " $omega " " $model "  " $IC_mode "   >> $input
                    
                    python infection_dynamics_realmodel.py "$data_dir_fin"
                    
                    
                    
                    # COLLECT OBSERVABLES

                    file_feat=${data_dir_fin}/data/global_features.txt
                    file_IC=${data_dir_fin}/data/initial_conditions.dat
                    
                    
                    
                #~ data_out = np.array([num_bact_tot, avg_bact_tot, var_bact_tot, num_phage_tot, avg_phage_tot, var_phage_tot, avg_IS, var_IS, timelast, timefirst2, timefirst5, timefirst8, num_bact1, num_bact2, num_bact5, num_bact8, avg_num_bact_strain, avg_num_phage_strain, avg_entr_bact, avg_entr_phage])
                
                    line_features=`head -1 $file_feat`
                    
                
                    bact_tot=`echo $line_features | awk '{ print $1 }'`        
                    R_tot=`echo $line_features | awk '{ print $2 }'`   
                    phage_tot=`echo $line_features | awk '{ print $3 }'`   
                    timelast=`echo $line_features | awk '{ print $4 }'`   
                    bact_fin=`echo $line_features | awk '{ print $5 }'`   
                    phage_fin=`echo $line_features | awk '{ print $6 }'`   
                    frac_resist_fin=`echo $line_features | awk '{ print $7 }'`   
                    bact_min=`echo $line_features | awk '{ print $8 }'`   
                    phage_min=`echo $line_features | awk '{ print $9 }'`   
                    
                    S0=`sed '2q;d' $file_IC`  
                    P0=`sed '5q;d' $file_IC`  
                 
                 
                    echo "$I" "$mu" "$p" "$K_D" "$phi" "$omega" "$bact_tot"  "$R_tot"  "$phage_tot"  "$timelast"  "$bact_fin"   "$phage_fin"  "$frac_resist_fin"  "$bact_min"  "$phage_min"  >> $file_out 
                   
                    echo "$I" "$mu" "$p" "$K_D" "$phi" "$omega" "$S0"  "$P0"   >> $file_out_IC 
                    
                    
                    
                    cd $curr_dir
                    echo " "
                    echo "p = "$p" mu = "$mu" K_D = "$K_D" I = "$I" model = "$model", run "$count_tot
                    ((count_tot++))
	    done
	    
	done
		
    done
    done
    done
    
    
        file_params="${data_dir_fin}/allparameters.json"
        file_params_new="${dir_out}/allparameters.json"
        
        cp $file_params $file_params_new

        cd $py_script_dir
        
        # Final phase diagram
        python Fig5_FIN.py 
        
        
        cd $curr_dir    
    done

fi
