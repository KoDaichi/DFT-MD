# [begin]
# name = ## specify the name by regexp
# aliases = 
# entry_type =
# has_table = 
# val_type\d* =
# choice\d* =
# description\d* =
# val_range\d* = 
# default_value =
# columns =
# max_nrow =
# importance\d* =
# restriction\d* =
# invalid_choice\d* =
# check_existence\d* =
# function\d* =
# required_attribute =    ## for tabular data
# recommended_attribute = ## for tabular data
# all_columns_required =  ## for tabular data
# [end]

# the 'control' block

[begin]
name = control$
entry_type = block
has_table = false
description = control the overall condition of the calculation from this block.
[end]

[begin]
name = control\.positron$
entry_type = primitive
val_type   = string
choice     = bulk,defect,off
description = set this variable in order to perform positron analysis
default_value = off
[end]

[begin]
name = control\.nfstopcheck$
entry_type = primitive
val_type   = int
description = specify the minimum number of SCF iterations which must be performed before the calculation can be terminated gracefully
default_value = 1
[end]

[begin]
name = control\.condition$
entry_type = primitive
val_type   = string __or__ int
choice     = preparation,initial,continuation,automatic,fixed_charge,fixed_charge_continuation,-3,-2,-1,0,1,2,3
description = specify the 'condition' of your calculation.
default_value = automatic
check_existence=F_CHGT,F_ZAJ,F_CNTN,F_CNTN_BIN __if__ control.condition __eq__ continuation __or__ 1 __or__ fixed_charge_continuation __or__ 3
check_existence1=F_CHGT __if__ control.condition __eq__ fixed_charge __or__ 2
check_existence2=F_OCCMAT __if__ control.condition __eq__ fixed_charge __or__ 2 && accuracy.hubbard.sw_hubbard __eq__ True
check_existence3=F_CNTN_BIN_PAW __if__ control.condition __eq__ fixed_charge __or__ 2 && accuracy.paw __eq__ True
check_existence4=F_HYBRIDINFO __if__ control.condition __eq__ fixed_charge __or__ 2 && accuracy.hybrid_functional.sw_hybrid_functional __eq__ True
check_existence5=F_SCF_ZAJ __if__ control.condition __eq__ fixed_charge __or__ 2 && accuracy.hybrid_functional.sw_hybrid_functional __eq__ True
[end]

[begin]
name = control\.cpumax$
entry_type = primitive
val_type = float
unit_type = longtime
description = enter the maximum elapsed time of the current simulation in units of 'long time'
default_value = 86400
function = check_cpumax
[end]

[begin]
name = control\.precision_wffile$
entry_type = primitive
val_type = string
choice=double_precision,single_precision
[end]

[begin]
name = control\.cachesize$
entry_type = primitive
val_type = int
description = enter the optimum cache size (may or may not increase the performance of your calculation)
default_value=256
[end]

[begin]
name = control\.max_iteration$
entry_type = primitive
val_type = int
description = enter the maximum number of SCF iterations. 
default_value = 10000
[end]

[begin]
name = control\.max_scf_iteration$
entry_type = primitive
val_type = int
description = enter the maximum number of SCF iterations per md step. 
[end]

[begin]
name = control\.max_mdstep$
entry_type = primitive
val_type = int
description = enter the maximum number of md steps.
[end]

[begin]
name = control\.driver$
entry_type = primitive
val_type = string
choice=general,neb,constraints,meta_dynamics,scdft,sc_dft
description = enter the 'driver' for the ionic evolution, ie, neb, meta_dynamics and others
default_value=general
[end]

[begin]
name = control\.multiple_replica_mode$
entry_type = primitive
val_type = bool
default_value=off
[end]

[begin]
name = control\.multiple_replica_max_iteration$
entry_type = primitive
val_type = int
val_range=1,
description = enter the maximum number of updates in an NEB simulation
default_value=100
[end]

[begin]
name=control\.sw_ekzaj$
entry_type=primitive
val_type=bool
description = when set to 'on', PHASE will output the WFs in ekcal and epsmain format
default_value=off
invalid_choice =True __if__ control.condition __eq__ fixed_charge __or__ fixed_charge_continuation __or__ 2 __or__ 3 __or__ -3 || accuracy.ksampling.method __eq__ __undefined__ __or__ mesh __or__ monk __or__ file __or__ directin
[end]

[begin]
name=control\.sw_dipole_correction$
entry_type=primitive
val_type=bool
default_value=off
[end]

[begin]
name=control\.sw_screening_correction$
entry_type=primitive
val_type=bool
default_value=off
[end]

[begin]
name=control\.sw_vdw_correction$
entry_type=primitive
val_type=bool
default_value=off
check_existence = F_DFTD3PAR __if__ accuracy.vdw_method __eq__ dft-d3 __or__ dftd3
[end]

[begin]
name=control\.sw_fef$
entry_type=primitive
val_type=bool
default_value=off
[end]

[begin]
name=control\.number_of_blocksize$
aliases=control\.blocksize$,control\.nblocksize$,control\.nb$,control\.nblocksize_dgemm$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_mgs$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_betar$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_betar_dot_wfs_nlmta$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_betar_dot_wfs_npe$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_vnonlocal$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_submat$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_submat_latter$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_force$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_rspace_betar$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.nblocksize_rspace_v$
entry_type=primitive
val_type=int
range=1,
[end]

[begin]
name=control\.sw_cif_output$
entry_type=primitive
val_type=bool
default_value=off
[end]

[begin]
name=control\.fixed_charge_option$
entry_type=block
has_table=false
[end]

[begin]
name=control\.fixed_charge_option\.kparallel$
entry_type=primitive
val_type=string
choice=one_by_one,all_at_once
[end]

[begin]
name=control\.sw_serial_fft$
entry_type=primitive
val_type=bool
[end]

[begin]
name=control\.sw_communicator_for_chg$
entry_type=primitive
val_type=bool
[end]

[begin]
name=control\.checkpoint_file$
entry_type=block
has_table=false
[end]

[begin]
name=control\.checkpoint_file\.iteration$
entry_type=primitive
val_type=int
[end]

[begin]
name=control\.checkpoint_file\.iteration_ionic$
entry_type=primitive
val_type=int
[end]

[begin]
name=control\.checkpoint_file\.iteration_unitcell$
entry_type=primitive
val_type=int
[end]

[begin]
name=control\.checkpoint_file\.iteration_neb$
entry_type=primitive
val_type=int
[end]

[begin]
name=control\.checkpoint_file\.iteration_reac$
entry_type=primitive
val_type=int
[end]

[begin]
name=control\.checkpoint_file\.cputime$
entry_type=primitive
val_type=float
unit_type=longtime
[end]

#the 'accuracy' block

[begin]
name = accuracy$
entry_type = block
has_table = false
description = configure the accuracy of the current simulation
[end]

[begin]
name=accuracy\.cutoff_wf$
aliases=accuracy\.cke_wavefunctions$,accuracy\.cke_wf$,accuracy\.cutoff_energy_for_wavefunctions$
entry_type = primitive
val_type = float
unit_type = energy
description = specify the cutoff energy for the WFs.
importance = required
[end]

[begin]
name=accuracy\.cutoff_cd$
aliases=accuracy\.cke_chargedensity$,accuracy\.cke_cd$,accuracy\.cutoff_energy_for_chargedensity$
entry_type = primitive
val_type = float
unit_type = energy
description = specify the cutoff energy for the augmentation charge. 4xcutoff_wf for norm-conserving pp, larger value for the uspp.
importance = required
function=check_cutoff_cd __if__ accuracy.cutoff_wf __eq__ __defined__
[end]

[begin]
name=accuracy\.num_bands$
entry_type=primitive
val_type=int
description=enter the number of bands. you must specify a number which is greater than the number of valence electrons.
importance=recommended
val_range=1,
function = check_num_bands
[end]

[begin]
name=accuracy\.ksampling$
entry_type=block
has_table=false
description = specify the k-point sampling under this block.
[end]

[begin]
name=accuracy\.ksampling\.method$
entry_type=primitive
val_type=string
choice=monk,mesh,gamma,file,directin
description = enter the method for the k-point sampling.
default_value=monk
restriction = mesh __if__ accuracy.smearing.method __eq__ tetrahedral __or__ t __or__ improved_tetrahedron __or__ tetrahedron
check_existence = F_KPOINT __if__ accuracy.ksampling.method __eq__ file
[end]

[begin]
name=accuracy\.ksampling\.mesh$
entry_type=block
has_table=false
function = check_k_mesh
importance=recommended __if__ accuracy.ksampling.method __eq__ __undefined__ __or__ monk && accuracy.ksampling.mp_index __eq__ __undefined__
importance1=recommended __if__ accuracy.ksampling.method __eq__ mesh
[end]

[begin]
name=accuracy\.ksampling\.mesh\.nx$
entry_type=primitive
val_type=int
description= the number of k-mesh with respect to the first recirpocal vector
val_range=1,
default_value=4
[end]

[begin]
name=accuracy\.ksampling\.mesh\.ny$
entry_type=primitive
val_type=int
description= the number of k-mesh with respect to the second recirpocal vector
val_range=1,
default_value=4
[end]

[begin]
name=accuracy\.ksampling\.mesh\.nz$
entry_type=primitive
val_type=int
val_range=1,
description= the number of k-mesh with respect to the third recirpocal vector
default_value=4
[end]

[begin]
name=accuracy\.ksampling\.mp_index$
entry_type=block
has_table=false
importance=recommended __if__ accuracy.ksampling.method __eq__ __undefined__ __or__ monk && accuracy.ksampling.mesh __eq__ __undefined__
[end]

[begin]
name=accuracy\.ksampling\.mp_index\.n1$
entry_type=primitive
val_type=int
description= the number of k-mesh with respect to the first recirpocal vector
val_range=1,
default_value=4
[end]

[begin]
name=accuracy\.ksampling\.mp_index\.n2$
entry_type=primitive
val_type=int
description= the number of k-mesh with respect to the second recirpocal vector
val_range=1,
default_value=4
[end]

[begin]
name=accuracy\.ksampling\.mp_index\.n3$
entry_type=primitive
val_type=int
description= the number of k-mesh with respect to the third recirpocal vector
val_range=1,
default_value=4
[end]

[begin]
name=accuracy\.ksampling\.kshift$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.ksampling\.kshift\.k1$
entry_type=primitive
val_type=float
val_range=0,0.5
description=enter the amount of 'shift' from the Gamma point with respect the 1st reciprocal vector
[end]

[begin]
name=accuracy\.ksampling\.kshift\.k2$
entry_type=primitive
val_type=float
description=
val_range=0,0.5
description=enter the amount of 'shift' from the Gamma point with respect the 2nd reciprocal vector 
[end]

[begin]
name=accuracy\.ksampling\.kshift\.k3$
entry_type=primitive
val_type=float
description=
val_range=0,0.5
description=enter the amount of 'shift' from the Gamma point with respect the 3rd reciprocal vector 
[end]

[begin]
name=accuracy\.ksampling\.kpoints$
entry_type=block
has_table=true
[end]

[begin]
name=accuracy\.ksampling\.kpoints\.table$
entry_type=table
columns=kx,ky,kz,denom,weight
val_type=float,float,float,float,int
importance=required __if__ accuracy.ksampling.method __eq__ directin
[end]

[begin]
name=accuracy\.ksampling\.base_reduction_for_gamma$
entry_type=primitive
val_type=bool
description=set this variable to 'on' in order to reduce calculations at the Gamma point
#restriction = False __if__ accuracy.hybrid_functional.sw_hybrid_functional __eq__ True
default_value=on
[end]

[begin]
name=accuracy\.ksampling\.base_symmetrization_for_gamma$
entry_type=primitive
val_type=bool
description=set this variable to 'on' in order to reduce calculations at the Gamma point
#restriction = False __if__ accuracy.hybrid_functional.sw_hybrid_functional __eq__ True
default_value=on
[end]

[begin]
name=accuracy\.smearing$
entry_type=block
has_table=false
[end]

name=accuracy\.smearing.\tetrahedron$
entry_type=block
has_table=false
[end]

name=accuracy\.smearing.\tetrahedron\.sw_correction$
entry_type=primitive
entry_type=bool
[end]

name=accuracy\.smearing.\tetrahedron\.dimension$
entry_type=primitive
entry_type=int
[end]

[begin]
name=accuracy\.smearing\.method$
entry_type=primitive
val_type=string
choice=parabolic,tetrahedral,tetrahedron,improved_tetrahedron,cold,t,fermi_dirac
description=enter the method for smearing. mind that you must specify either tetrahedral or improved_tetrahedron in order to output the DOS by the tetrahedron method
default_value=parabolic
importance = required __if__ postprocessing.dos.method __eq__ t __or__ tetrahedron __or__ tetrahedral
invalid_choice = __undefined__ __or__ parabolic __or__ cold __if__ postprocessing.dos.method __eq__ tetrahedron __or__ t __or__ tetrahedral && postprocessing.dos.sw_dos __eq__ True
[end]

[begin]
name=accuracy\.electronic_temp$
entry_type=primitive
val_type=float
unit_type=temperature
[end]

[begin]
name=accuracy\.smearing\.width$
entry_type=primitive
val_type=float
unit_type=energy
description=enter the smearing width. 
default_value=0.001 hartree
[end]

[begin]
name=accuracy\.xctype$
entry_type=primitive
val_type=string
choice=ldapw91,ggapbe,ggapbex,ggapw91,revpbe,rpbe,pbesol,vdwdf,vdwdf2,vdwdf-c09x,vdwdf2-c09x,vdwdf-optpbe,vdwdf-optb86b,vdwdf2-b86r,vdwdf-cx
function=check_xctype __if__ accuracy.xctype __eq__ __defined__
[end]

[begin]
name=accuracy\.vdwdf$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.vdwdf\.mode$
entry_type=primitive
val_type=string
choice=oneshot,scf
[end]

[begin]
name=accuracy\.vdwdf\.ndel$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.vdwdf\.nphid$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.vdwdf\.nr12$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.vdwdf\.maxk$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdwdf\.r12max$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=accuracy\.vdwdf\.dq$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdwdf\.lambda$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdwdf\.q0cut$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdwdf\.q0min$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdwdf\.ds$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdwdf\.na_gl$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.vdwdf\.a1$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdwdf\.a2$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdwdf\.eval_kernel_by_interpolation$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.vdwdf\.save_memory_mode$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.vdwdf\.sw_use_wugygi_method$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.vdwdf\.vdwdf_version$
entry_type=primitive
val_type=int
range=1,2
[end]

[begin]
name=accuracy\.vdwdf\.exchange_pot_type$
entry_type=primitive
val_type=string
choice=revpbe,pw86r,c09x,optpbe,optb86b,b86r,lvpw86r
[end]

[begin]
name=accuracy\.hybrid_functional$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_hybrid_functional$
entry_type=primitive
val_type=bool
default_value=off
importance  = required __if__ accuracy.hybrid_functional.sw_exchange_only __eq__ True
#restriction = __undefined__ __or__ False __if__ accuracy.ksampling.method __eq__ __undefined__ __or__ mesh __or__ monk __or__ file __or__ directin
[end]

[begin]
name=accuracy\.hybrid_functional\.functional_type$
entry_type=primitive
val_type=string
choice=hse06,pbe0,hf
[end]

[begin]
name=accuracy\.hybrid_functional\.alpha$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.hybrid_functional\.omega$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_exchange_only$
entry_type=primitive
val_type=bool
default_value=off
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_screened_exchange$
entry_type=primitive
val_type=bool
default_value=off
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_singular_correction$
entry_type=primitive
val_type=bool
default_value=off
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_memory_reduction$
entry_type=primitive
val_type=bool
default_value=off
[end]

[begin]
name=accuracy\.hybrid_functional\.reduction_factor$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.hybrid_functional\.reduction_factor\.f1$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.hybrid_functional\.reduction_factor\.f2$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.hybrid_functional\.reduction_factor\.f3$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.hybrid_functional\.omega$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.hybrid_functional\.omega_pbe$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_rspace$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_rspace_dgm$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_eval_vexx$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_retard_eigval_evaluation$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_precalculate$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_change_axis$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.hybrid_functional\.charge_mesh$
aliases=accuracy\.hybrid_functional\.accuracy$
entry_type=primitive
val_type=string
choice=exact,fine,moderate,coarse
[end]

[begin]
name=accuracy\.hybrid_functional\.cutoff_wf_for_exx$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=accuracy\.hybrid_functional\.sw_output_hybrid_info$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.scf_convergence$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.scf_convergence\.delta_total_energy$
entry_type=primitive
val_type=float
unit_type=energy
#importance = recommended __if__ control.condition __eq__ __undefined__ __or__ initial __or__ continuation __or__ automatic __or__ -1 __or__ 0 __or__ 1
description=enter the convergence criteria for the SCF loop
default_value=1e-9 hartree
function=check_scf_convergence_criteria __if__ control.condition __eq__ __undefined__ __or__ initial __or__ continuation __or__ automatic __or__ -1 __or__ 0 __or__ 1
[end]

[begin]
name=accuracy\.scf_convergence\.succession$
entry_type=primitive
val_type=int
description = specify the number of consecutive times the convergence criterion must be met in order for the calculation to be deemed converged.
[end]

[begin]
name=accuracy\.force_convergence$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.force_convergence\.max_force$
aliases=accuracy\.force_convergence\.delta_force$
entry_type=primitive
val_type=float
unit_type=force
description=enter the convergence criteria for the structural optimization loop
default_value=1e-3 hartree/bohr
[end]

[begin]
name=accuracy\.force_convergence\.delta_force$
entry_type=primitive
val_type=float
unit_type=force
description=enter the convergence criteria for the structural optimization loop
default_value=1e-3 hartree/bohr
[end]

[begin]
name=accuracy\.force_convergence\.output_molecular_forcmx$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.initial_wavefunctions$
entry_type=primitive
val_type=string
choice=random_numbers,matrix_diagon,file,atomic_orbital,atomic_orbitals
description=enter the method for initial WF generation
default_value=random_numbers
#importance   = recommended __if__ accuracy.hubbard.sw_hubbard __eq__ True
#importance1  = required __if__ accuracy.hybrid_functional.sw_hybrid_functional __eq__ True
#restriction = file __if__ accuracy.hybrid_functional.sw_hybrid_functional __eq__ True
check_existence = F_ZAJ __if__ accuracy.initial_wavefunctions __eq__ file
[end]

[begin]
name=accuracy\.initial_charge_density$
entry_type=primitive
val_type=string
choice=gauss,atomic_charge_density,file
description=enter the method for initial charge density generation
default_value=gauss
#importance  = required __if__ accuracy.hybrid_functional.sw_hybrid_functional __eq__ True
#restriction = file __if__ accuracy.hybrid_functional.sw_hybrid_functional __eq__ True
check_existence = F_CHGT __if__ accuracy.initial_charge_density __eq__ file
[end]

[begin]
name=accuracy\.matrix_diagon$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.matrix_diagon\.cutoff_wf$
aliases=accuracy\.matrix_diagon\.cke_initial_matdiagon$,accuracy\.matrix_diagon\.cke_wavefunctions$,accuracy\.matrix_diagon\.cutoff_energy_for_wavefunctions$,accuracy\.matrix_diagon\.cke_wf$
entry_type=primitive
val_type=float
unit_type=energy
description=enter the cutoff-energy used exclusively for initial WF generation.
importance=recommended __if__ accuracy.initial_wavefunctions __eq__ matrix_diagon
[end]

[begin]
name=accuracy\.hubbard$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.hubbard\.sw_hubbard$
entry_type=primitive
val_type=bool
description = set this variable to 'on' in order to use the DFT+U method
default_value=off
[end]

[begin]
name=accuracy\.hubbard\.projectors$
entry_type=block
has_table=true
[end]

[begin]
name=accuracy\.hubbard\.projectors\.table$
entry_type=table
columns=no,ueff
val_type=int,float
importance = required __if__ accuracy.hubbard.sw_hubbard __eq__ True
function = projector_exists __if__ accuracy.hubbard.sw_hubbard __eq__ True
[end]

[begin]
name=accuracy\.projector_list$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.projector_list\.projectors$
entry_type=block
has_table=true
[end]

[begin]
name=accuracy\.projector_list\.projectors\.table$
entry_type=table
columns=no,group,radius,l,t
val_type=int,int,float,int,int
importance = required __if__ accuracy.hubbard.sw_hubbard __eq__ True || postprocessing.pdos.sw_pdos __eq__ True
choice4=0,1,2,3
choice5=1,2
function = check_proj_group
all_columns_required=false
#required_attribute = no __if__ accuracy.hubbard.sw_hubbard __eq__ True
required_attribute=group 
required_attribute1=radius
required_attribute2=l
[end]

[begin]
name=accuracy\.cutoff_pwf$
entry_type=primitive
val_type=float
unit_type=energy
importance = required __if__ control.positron __eq__ __defined__
[end]

[begin]
name=accuracy\.positron_convergence$
entry_type=block
has_table=false
importance=recommended __if__ control.positron __eq__ __defined__
[end]

[begin]
name=accuracy\.positron_convergence\.num_extra_bands$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.positron_convergence\.delta_eigenvalue$
entry_type=primitive
val_type=float
unit_type=energy
importance=recommended __if__ control.positron __eq__ __defined__
[end]

[begin]
name=accuracy\.positron_convergence\.succession$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.positron_convergence\.num_max_iteration$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.positron_convergence\.dtim$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.positron_convergence\.epsilon_ele$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.ek_convergence$
entry_type=block
has_table=false
description=specify the convergence criterion for fixed-charge calculations performed by the ekcal/epsmain program. it is strongly recommended to configure this block in such cases.
importance = recommended __if__ control.condition __eq__ fixed_charge __or__ fixed_charge_continuation __or__ 2 __or__ 3 __or__ -3
[end]

[begin]
name=accuracy\.ek_convergence\.num_max_iteration$
entry_type=primitive
val_type=int
description = specify the max. number of iterations for solving the WFs at fixed charge
[end]

[begin]
name=accuracy\.ek_convergence\.delta_eigenvalue$
entry_type=primitive
val_type=float
unit_type=energy
description = the convergence criteria for solving the WFs at fixed charge
#importance=recommended __if__ control.condition __eq__ fixed_charge __or__ fixed_charge_continuation __or__ 2 __or__ 3 __or__ -3
[end]

[begin]
name=accuracy\.ek_convergence\.succession$
entry_type=primitive
val_type=int
description = specify the number of consecutive times the convergence criteria has to be met in order for the calculation to be deemed converged.
[end]

[begin]
name=accuracy\.ek_convergence\.sw_eval_eig_diff$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.ek_convergence\.num_extra_bands$
entry_type=primitive
val_type=int
default_value=2
[end]

[begin]
name=accuracy\.max_force_trans$
entry_type=primitive
val_type=float
unit_type=force
[end]

[begin]
name=accuracy\.max_torque$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=accuracy\.precalculation$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.precalculation\.nel_ylm$
entry_type=primitive
val_type=int
default_value=9
[end]

#dipole
[begin]
name=accuracy\.dipole_correction$
entry_type=block
has_table=false
importance=required __if__ control.sw_dipole_correction __eq__ True
[end]

[begin]
name=accuracy\.dipole_correction\.direction$
entry_type=primitive
val_type=int
default_value=0
val_range=0,3
[end]

[begin]
name=accuracy\.dipole_correction\.division$
entry_type=primitive
val_type=int
default_value=100
val_range=1,
[end]

[begin]
name=accuracy\.dipole_correction\.vacuum$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.dipole_correction\.vacuum\.rx$
entry_type=primitive
val_type=float
default_value=0.0
val_range=-1.0,1.0
[end]

[begin]
name=accuracy\.dipole_correction\.vacuum\.ry$
entry_type=primitive
val_type=float
default_value=0.0
val_range=-1.0,1.0
[end]

[begin]
name=accuracy\.dipole_correction\.vacuum\.rz$
entry_type=primitive
val_type=float
default_value=0.0
val_range=-1.0,1.0
[end]

[begin]
name=accuracy\.dipole_correction\.electric_field$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.dipole_correction\.electric_field\.ex$
entry_type=primitive
val_type=float
default_value=0.0
val_range=-1.0,1.0
[end]

[begin]
name=accuracy\.dipole_correction\.electric_field\.ey$
entry_type=primitive
val_type=float
default_value=0.0
val_range=-1.0,1.0
[end]

[begin]
name=accuracy\.dipole_correction\.electric_field\.ez$
entry_type=primitive
val_type=float
default_value=0.0
val_range=-1.0,1.0
[end]

#screening correction
[begin]
name=accuracy\.screening_correction$
entry_type=block
has_table=false
importance=required __if__ control.sw_screening_correction __eq__ True
[end]
[begin]
name=accuracy\.screening_correction\.alpha$
entry_type=primitive
val_type=float
default_value=1
[end]

#fermi surface
[begin]
name=accuracy\.fermi_surface$
entry_type=block
has_table=false
[end]
[begin]
name=accuracy\.fermi_surface\.sw_write_bxsf_file$
entry_type=primitive
val_type=bool
default_value=false
[end]

#FEF
[begin]
name=accuracy\.fef$
entry_type=block
has_table=false
importance=required __if__ control.sw_fef __eq__ True
[end]

[begin]
name=accuracy\.fef\.ex$
entry_type=primitive
val_type=float
unit_type=force
[end]

[begin]
name=accuracy\.fef\.ey$
entry_type=primitive
val_type=float
unit_type=force
[end]

[begin]
name=accuracy\.fef\.ez$
entry_type=primitive
val_type=float
unit_type=force
[end]

#vdw
[begin]
name=accuracy\.vdw_method$
entry_type=primitive
val_type=string
choice=williams,grimme,dft-d2,dft-d3,dftd3
importance=recommended __if__ control.sw_vdw_correction __eq__ True
[end]

[begin]
name=accuracy\.vdw_radius$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdw_scaling_factor$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdw_scaling_factor_r$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.vdw_damping_factor_r$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.paw$
aliases=accuracy\.paw_switch
entry_type=primitive
val_type=bool
restriction = True __if__ accuracy.spinorbit.mode __eq__ pawpot
restriction1 = False __if__ accuracy.hybrid_functional.sw_hybrid_functional __eq__ True
[end]

[begin]
name=accuracy\.spinorbit$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.spinorbit\.mode$
entry_type=primitive
val_type=string
choice=pawpot,projector,neglected,builtin,zeff,read_from_pp
[end]

[begin]
name=accuracy\.sw_manual_occupation$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.occupation$
entry_type=block
has_table=true
[end]

[begin]
name=accuracy\.occupation\.table$
entry_type=table
all_columns_required=false
columns=band_index,occ,occ_up,occ_down
val_type=int,float,float,float
[end]

# the structure block
[begin]
name = structure$
entry_type = block
has_table = false
description=specify the model of interest under this block.
[end]

[begin]
name=structure\.method$
entry_type=primitive
val_type=string
choice=directin,file
[end]

[begin]
name=structure\.file$
entry_type=block
has_table=false
[end]

[begin]
name=structure\.file\.filetype$
entry_type=primitive
val_type=string
choice=phase0_input,phase0_output
[end]

[begin]
name=structure\.file\.frame$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure\.unit_cell$
entry_type=block
has_table=false
description=
importance=required
[end]

[begin]
name=structure\.unit_cell_type$
entry_type=primitive
val_type=string
choice=bravais,primitive
[end]

[begin]
name=structure\.unit_cell\.a_vector$
entry_type = vector
num_elements = 3
val_type = float,float,float
description=the a-vector
[end]

[begin]
name=structure\.unit_cell\.b_vector$
entry_type = vector
num_elements = 3
val_type = float,float,float
description=the b-vector
[end]

[begin]
name=structure\.unit_cell\.c_vector$
entry_type = vector
num_elements = 3
val_type = float,float,float
description=the c-vector
[end]

[begin]
name=structure\.unit_cell\.a$
entry_type=primitive
val_type=float
unit_type=length
restriction = __undefined__  __if__ structure.unit_cell_type __eq__ primitive
[end]

[begin]
name=structure\.unit_cell\.b$
entry_type=primitive
val_type=float
unit_type=length
restriction = __undefined__  __if__ structure.unit_cell_type __eq__ primitive
[end]

[begin]
name=structure\.unit_cell\.c$
entry_type=primitive
val_type=float
unit_type=length
restriction = __undefined__  __if__ structure.unit_cell_type __eq__ primitive
[end]

[begin]
name=structure\.unit_cell\.alpha$
entry_type=primitive
val_type=float
unit_type=angle
restriction = __undefined__  __if__ structure.unit_cell_type __eq__ primitive
[end]

[begin]
name=structure\.unit_cell\.beta$
entry_type=primitive
val_type=float
unit_type=angle
restriction = __undefined__  __if__ structure.unit_cell_type __eq__ primitive
[end]

[begin]
name=structure\.unit_cell\.gamma$
entry_type=primitive
val_type=float
unit_type=angle
restriction = __undefined__  __if__ structure.unit_cell_type __eq__ primitive
[end]

[begin]
name=structure\.symmetry$
entry_type = block
has_table = false
description=
[end]

[begin]
name=structure\.symmetry\.method$
entry_type=primitive
val_type=string
choice=manual,automatic
function=check_inversion __if__ structure.symmetry.method __eq__ automatic
[end]

[begin]
name=structure\.symmetry\.crystal_structure$
aliases=structure\.symmetry\.crystal$
entry_type=primitive
val_type=string
choice=diamond,hexagonal,fcc,bcc,simple_cubic,facecentered_cubic,trigonal,hcp,body_centered
[end]

[begin]
name=structure\.symmetry\.tspace$
entry_type=block
has_table=false
[end]

[begin]
name=structure\.symmetry\.tspace\.lattice_system$
aliases=structure\.symmetry\.tspace\.system$
entry_type=primitive
val_type=string
choice = primitive,simple,p,s,1,facecentered,fcc,f,2,bodycentered,bcc,b,3,bottomcentered,basecentered,onefacecentered,bot,ba,o,4,rhombohedral,trigonal,r,t,-1,hexagonal,h,0
#importance = recommended __if__ structure.symmetry.method __eq__ automatic
[end]

[begin]
name=structure\.symmetry\.tspace\.num_generators$
entry_type=primitive
val_type=int
val_range=1,3
[end]

[begin]
name=structure\.symmetry\.tspace\.generators$
entry_type=block
has_table=true
[end]

[begin]
name=structure\.symmetry\.tspace\.generators\.table$
entry_type=table
columns=rotation,tx,ty,tz
val_type=string,fract,fract,fract
max_nrow=3
[end]

[begin]
name=structure\.symmetry\.tspace\.af_generator$
entry_type=block
has_table=true
[end]

[begin]
name=structure\.symmetry\.tspace\.af_generator\.table$
entry_type=table
columns=rotation,tx,ty,tz
val_type=string,fract,fract,fract
importance = required __if__ structure.magnetic_state __eq__ af __or__ antiferro
[end]

[begin]
name=structure\.symmetry\.sw_inversion$
entry_type = primitive
val_type = bool
description=set this variable to 'on' if the system of interest has inversion symmetry.
function = check_inversion __if__ structure.symmetry.sw_inversion __eq__ True
[end]

[begin]
name=structure\.magnetic_state$
entry_type = primitive
val_type = string
choice = para,ferro,antiferro,af,noncollinear
description = 
[end]

[begin]
name=structure\.ferromagnetic_state$
entry_type = block
has_table=false
description = 
[end]

[begin]
name=structure\.ferromagnetic_state\.sw_fix_total_spin$
entry_type=primitive
val_type=bool
description=
restriction = __undefined__ __or__ False __if__ structure.magnetic_state __eq__ __undefined__ __or__ para __or__ af __or__ antiferro
[end]

[begin]
name=structure\.ferromagnetic_state\.total_spin$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=structure\.ferromagnetic_state\.spin_fix_period$
entry_type=primitive
val_type=string __or__ int
description=
[end]

[begin]
name=structure\.atom_list$
entry_type = block
has_table=false
imortance=required
[end]

[begin]
name=structure\.atom_list\.coordinate_system$
entry_type=primitive
val_type=string
choice=cartesian,internal,pucv
description=specify cartesian if you are specifing the atomic coordinates in cartesian coordinates, internal if you are specifing the atomic coordinates in fractional coordinates.
[end]

[begin]
name=structure\.atom_list\.atoms$
entry_type=block
has_table=true
description=
[end]

[begin]
importance=required
name=structure\.atom_list\.atoms\.table$
entry_type=table
columns=no,element,rx,ry,rz,mobile,num_layer,weight,aldos,thermo_group,proj_group,molecule,berry,born,region_group,mobilex,mobiley,mobilez,vx,vy,vz,vdw,number
val_type=int,string,float,float,float,bool,int,int,bool,int,int,int,bool,bool,int,bool,bool,bool,float,float,float,string,int
choice8=1,2
description=enter the atomic cooordinates (and associated attributes) from this table.
description1=the ID of the atom
description2=element name
description3=the x-coordinate 
description4=the y-coordinate 
description5=the z-coordinate 
description6=set this value to 'on' if the corresponding atom should be 'mobile' during ionic relaxation or molecular dynamics
description7=
description8=
description9=
description10=
description11=
description12=
all_columns_required=false
required_attribute = element
required_attribute1 = rx
required_attribute2 = ry
required_attribute3 = rz
required_attribute4 = num_layer __if__ postprocessing.ldos.sw_layerdos __eq__ True && postprocessing.ldos.layerdos.slicing_way __eq__ by_atomic_positions
required_attribute6 = proj_group __if__ postprocessing.pdos.sw_pdos __eq__ True || accuracy.hubbard.sw_hubbard __eq__ True
required_attribute7 = thermo_group __if__ structure_evolution.method __eq__ temperature_control
required_attribute8 = aldos __if__ postprocessing.ldos.sw_aldos __eq__ True && postprocessing.ldos.aldos.naldos_from __eq__ __undefined__
required_attribute9 = aldos __if__ postprocessing.ldos.sw_aldos __eq__ True && postprocessing.ldos.aldos.naldos_to __eq__ __undefined__
required_attribute = molecule __if__ control.driver __eq__ rigid_body
function=valid_proj_group_assigned __if__ postprocessing.pdos.sw_pdos __eq__ True || accuracy.hubbard.sw_hubbard __eq__ True
function1 = check_num_layer __if__ postprocessing.ldos.sw_layerdos __eq__ True && postprocessing.ldos.layerdos.slicing_way __eq__ by_atomic_positions
function2 = check_fract __if__ structure.atom_list.coordinate_system __eq__ internal __or__ __undefined__
#function3 = check_bond_length
[end]

[begin]
name=structure\.element_list$
entry_type=block
has_table=true
[end]

[begin]
name=structure\.element_list\.table$
entry_type=table
columns=id,no,element,atomicnumber,mass,zeta,deviation,dev,standard_deviation,qex,theta,phi,mx,my,mz,moment
val_type=int,int,string,float,float,float,float,float,float,float,float,float,float,float,float,float
importance=required
function = all_elements_are_defined
function1 = check_mass __if__ structure_evolution.method __eq__ velocity_verlet __or__ temperature_control __or__ velocity_scaling || phonon.sw_phonon __eq__ True
function2 = check_ppfiles
max_nrow=16
all_columns_required=false
required_attribute=element
required_attribute1=atomicnumber
required_attribute2=mass __if__ structure_evolution.method __eq__ velocity_verlet __or__ temperature_control __or__ velocity_scaling || phonon.sw_phonon __eq__ True
recommended_attribute=zeta __if__ structure.magnetic_state __eq__ ferro __or__ af __or__ antiferro && control.condition __eq__ __undefined__ __or__ initial __or__ continuation __or__ automatic __or__ -1 __or__ 0 __or__ 1
[end]

[begin]
name=structure\.strain$
entry_type=block
has_table=false
[end]

[begin]
name=structure\.strain\.sw_strained_cell$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure\.strain\.e11$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.strain\.e22$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.strain\.e33$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.strain\.e21$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.strain\.e12$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.strain\.e31$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.strain\.e13$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.strain\.e32$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.strain\.e23$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.region\d+$
entry_type=block
has_table=false
[end]

[begin]
name=structure\.region\d+\.type$
entry_type=primitive
val_type=string
choice=cylinder,box
[end]

[begin]
name=structure\.region\d+\.region_group$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure\.region\d+\.radius$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.cylx$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.cyly$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.xmax$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.xmin$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.ymax$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.ymin$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.zmax$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.zmin$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.orientation$
entry_type=primitive
val_type=int
range=1,3
[end]

[begin]
name=structure\.region\d+\.sigma$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.region\d+\.epsilon$
entry_type=primitive
val_type=energy
unit_type=length
[end]

[begin]
name=structure\.region\d+\.tally$
entry_type=primitive
val_type=bool
[end]

# wavefunction solover and charge density mixing

[begin]
name = wavefunction_solver$
aliases = wf_solver$
entry_type = block
has_table=false
[end]

[begin]
name = wavefunction_solver\.solvers$
aliases = wf_solver\.solvers$
entry_type = block
has_table=true
[end]

[begin]
name = wavefunction_solver\.solvers\.table$
aliases = wf_solver\.solvers\.table$
entry_type = table
columns =id,no,sol,till_n,dts,dte,itr,var,prec,cmix,submat
val_type=int,int,string,int,float,float,int,string,bool,int,bool
choice3=matrixdiagon,lmmsd,lm+msd,rmm,rmm2p,rmm3,davidson,msd,sd,cg,pkosugi,pdavidson,mddavidson,mdkosugi
choice8=tanh,linear
##importance = recommended
Function=check_consistency_of_wavefunction_solver
all_columns_required=false
recommended_attribute=sol
#recommended_attribute1=cmix
#recommended_attribute2=submat
recommended_attribute1=till_n
[end]

[begin]
name = wavefunction_solver\.line_minimization$
aliases=wavefunction_solver\.lineminimization$,wf_solver\.line_minimization,wf_solver\.lineminimization$
entry_type = block
has_table = false
[end]

[begin]
name = wavefunction_solver\.line_minimization\.dt_lower_critical$
aliases = wavefunction_solver\.lineminimization\.dt_lower_critical$,wf_solver\.line_minimization\.dt_lower_critical$,wf_solver\.lineminimization\.dt_lower_critical$
entry_type = primitive
val_type = float
unit_type = none
[end]

[begin]
name = wavefunction_solver\.line_minimization\.dt_lower_factor$
aliases = wavefunction_solver\.lineminimization\.dt_lower_factor$,wf_solver\.line_minimization\.dt_lower_factor$,wf_solver\.lineminimization\.dt_lower_factor$
entry_type = primitive
val_type = float
unit_type = none
[end]

[begin]
name = wavefunction_solver\.line_minimization\.dt_upper_factor$
aliases = wavefunction_solver\.lineminimization\.dt_upper_factor$,wf_solver\.line_minimization\.dt_upper_factor$,wf_solver\.lineminimization\.dt_upper_factor$
entry_type = primitive
val_type = float
val_range=1,
unit_type = none
[end]

[begin]
name = wavefunction_solver\.line_minimization\.dt_upper_critical$
aliases = wavefunction_solver\.lineminimization\.dt_upper_critical$,wf_solver\.line_minimization\.dt_upper_critical$,wf_solver\.lineminimization\.dt_upper_critical$
entry_type = primitive
val_type = float
unit_type = none
[end]

[begin]
name = wavefunction_solver\.line_minimization\.delta_lmdenom$
aliases = wavefunction_solver\.lineminimization\.delta_lmdenom$,wf_solver\.line_minimization\.delta_lmdenom$,wf_solver\.lineminimization\.delta_lmdenom$
entry_type = primitive
val_type = float
unit_type = none
[end]

[begin]
name = wavefunction_solver\.rmm$
aliases = wf_solver\.rmm$
entry_type = block
has_table = false
[end]

[begin]
name = wavefunction_solver\.rmm\.edelta_change_to_rmm$
aliases = wf_solver\.rmm\.edelta_change_to_rmm$
entry_type = primitive
val_type = float
unit_type = energy
[end]

[begin]
name = wavefunction_solver\.rmm\.imgsrmm$
aliases = wf_solver\.rmm\.imgsrmm$
entry_type = primitive
val_type = int
[end]

[begin]
name = wavefunction_solver\.rmm\.save_memory_mode$
aliases = wf_solver\.rmm\.save_memory_mode$
entry_type = primitive
val_type = bool
[end]

[begin]
name=wavefunction_solver\.rmm\.rr_critical_value$
aliases=wf_solver\.rmm\.rr_critical_value$
entry_type=primitive
val_type=float
unit_type=none
[end]

[begin]
name=wavefunction_solver\.subspace_rotation$
aliases=wavefunction_solver\.submat$,wf_solver\.subspace_rotation$,wf_solver\.submat$
entry_type=block
has_table=false
[end]

[begin]
name=wavefunction_solver\.subspace_rotation\.subspace_matrix_size$
aliases=wavefunction_solver\.submat\.subspace_matrix_size$,wf_solver\.submat\.subspace_matrix_size$,wf_solver\.subspace_rotation\.subspace_matrix_size$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.subspace_rotation\.damping_factor$
aliases=wavefunction_solver\.submat\.damping_factor$,wf_solver\.subspace_rotation\.damping_factor$,wf_solver\.submat\.damping_factor$
entry_type=primitive
val_type=float
unit_type=none
val_range=0,1
[end]

[begin]
name=wavefunction_solver\.subspace_rotation\.period$
aliases=wavefunction_solver\.submat\.period$,wf_solver\.subspace_rotation\.period$,wf_solver\.submat\.period$
entry_type=primitive
val_type=int
default_value=1
[end]

[begin]
name=wavefunction_solver\.subspace_rotation\.critical_ratio$
aliases=wavefunction_solver\.submat\.critical_ratio$,wf_solver\.subspace_rotation\.critical_ratio$,wf_solver\.submat\.critical_ratio$
entry_type=primitive
val_type=float
unit_type=none
[end]

[begin]
name=wavefunction_solver\.subspace_rotation\.before_renewal$
aliases=wavefunction_solver\.submat\.before_renewal$,wf_solver\.subspace_rotation\.before_renewal$,wf_solver\.submat\.before_renewal$
entry_type=primitive
val_type=bool
[end]

[begin]
name=wavefunction_solver\.subspace_rotation\.scalapack$
aliases=wavefunction_solver\.submat\.scalapack$,wf_solver\.submat\.scalapack$,wf_solver\.subspace_rotation\.scalapack$
entry_type=block
has_table=false
[end]

[begin]
name=wavefunction_solver\.subspace_rotation\.scalapack\.sw_scalapack$
aliases=wavefunction_solver\.submat\.scalapack\.sw_scalapack$,wf_solver\.submat\.scalapack\.sw_scalapack$,wf_solver\.subspace_rotation\.scalapack\.sw_scalapack$
entry_type=primitive
val_type=bool
[end]

[begin]
name=wavefunction_solver\.davidson$
aliases=wf_solver\.davidson$
entry_type=block
has_table=false
[end]

[begin]
name=wavefunction_solver\.davidson\.max_iter_david$
aliases=wf_solver\.davidson\.max_iter_david$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.davidson\.ndavid$
aliases=wf_solver\.davidson\.ndavid$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.davidson\.max_subspace_size$
aliases=wf_solver\.davidson\.max_subspace_size$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.mddavidson$
aliases=wf_solver\.mddavidson$
entry_type=block
has_table=false
[end]

[begin]
name=wavefunction_solver\.mddavidson\.max_iter_mddavid$
aliases=wf_solver\.mddavidson\.max_iter_mddavid$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.mddavidson\.npartition_mddavid$
aliases=wf_solver\.mddavidson\.npartition_mddavid$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.mddavidson\.delta_eig_occup$
aliases=wf_solver\.mddavidson\.delta_eig_occup$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.mddavidson\.delta_eig_empty$
aliases=wf_solver\.mddavidson\.delta_eig_empty$
entry_type=primitive
val_type=bool
[end]

[begin]
name=wavefunction_solver\.mddavidson\.eps_mddavid$
aliases=wf_solver\.mddavidson\.eps_mddavid$
entry_type=primitive
val_type=bool
[end]

[begin]
name=wavefunction_solver\.mddavidson\.eps_residual_mddavid$
aliases=wf_solver\.mddavidson\.eps_residual_mddavid$
entry_type=primitive
val_type=bool
[end]

[begin]
name=wavefunction_solver\.mdkosugi$
aliases=wf_solver\.mdkosugi$
entry_type=block
has_table=false
[end]

[begin]
name=wavefunction_solver\.mdkosugi\.max_iter_mdkosugi$
aliases=wf_solver\.mdkosugi\.max_iter_mdkosugi$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.mdkosugi\.npartition_mdkosugi$
aliases=wf_solver\.mdkosugi\.npartition_mdkosugi$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.mdkosugi\.delta_eig_occup_mdkosugi$
aliases=wf_solver\.mdkosugi\.delta_eig_occup_mdkosugi$
entry_type=primitive
val_type=int
[end]

[begin]
name=wavefunction_solver\.mdkosugi\.delta_eig_empty_mdkosugi$
aliases=wf_solver\.mdkosugi\.delta_eig_empty_mdkosugi$
entry_type=primitive
val_type=bool
[end]

[begin]
name=wavefunction_solver\.mdkosugi\.eps_mdkosugi$
aliases=wf_solver\.mdkosugi\.eps_mdkosugi$
entry_type=primitive
val_type=bool
[end]

[begin]
name=wavefunction_solver\.mdkosugi\.eps_residual_mdkosugi$
aliases=wf_solver\.mdkosugi\.eps_residual_mdkosugi$
entry_type=primitive
val_type=bool
[end]

[begin]
name = charge_mixing$
entry_type = block
has_table = false
[end]

[begin]
name = charge_mixing\.method$
entry_type = primitive
val_type=string
choice=simple,broyden2,pulay
[end]

[begin]
name = charge_mixing\.rmx$
entry_type = primitive
val_type=float
[end]

[begin]
name = charge_mixing\.istr$
entry_type = primitive
val_type=int
val_range=1,
[end]

[begin]
name = charge_mixing\.nbxmix$
entry_type = primitive
val_type=int
val_range=1,
[end]

[begin]
name = charge_mixing\.sw_mix_charge_hardpart$
entry_type = primitive
val_type=bool
[end]

[begin]
name = charge_mixing\.sw_mix_imaginary_hardpart$
entry_type = primitive
val_type=bool
[end]

[begin]
name = charge_mixing\.sw_mix_bothspins_sametime$
entry_type = primitive
val_type=bool
[end]

[begin]
name = charge_mixing\.sw_mix_occ_matrix$
entry_type = primitive
val_type=bool
[end]

[begin]
name = charge_mixing\.num_mixing_methods$
entry_type = primitive
val_type=int
[end]

[begin]
name = charge_mixing\.sw_recomposing$
entry_type = primitive
val_type=bool
[end]

[begin]
name = charge_mixing\.spin_density_mixfactor$
entry_type=primitive
val_type=float
[end]

[begin]
name = charge_mixing\.mixing_methods$
entry_type = block
has_table = true
[end]

[begin]
name = charge_mixing\.mixing_methods\.table$
entry_type = table
columns =id,no,method,rmxs,rmxe,itr,var,prec,istr,nbmix,nbxmix,update
val_type=int,int,string,float,float,int,string,bool,int,int,int,string
choice3=simple,linear,broyden2,broyden,pulay,dfp
choice7=tanh,linear
choice12=renew,anew
#importance = recommended __if__ control.condition __eq__ __undefined__ __or__ initial __or__ continuation __or__ automatic
val_range4=0,1
val_range5=0,1
function = check_charge_mixing
all_columns_required=false
recommended_attribute=method
recommended_attribute1=rmxs
#recommended_attribute2=rmxe
#recommended_attribute3=prec
[end]

[begin]
name=charge_mixing\.charge_preconditioning$
entry_type=block
has_table=false
[end]

[begin]
name=charge_mixing\.charge_preconditioning\.amix$
entry_type=primitive
val_type=float
unit_type=none
val_range=0,1
[end]

[begin]
name=charge_mixing\.charge_preconditioning\.bmix$
entry_type=primitive
val_type=float
unit_type=none
[end]

[begin]
name=charge_mixing\.charge_preconditioning\.metric_ratio$
entry_type=primitive
val_type=float
unit_type=none
[end]

[begin]
name=charge_mixing\.charge_preconditioning\.sw_precon_diff$
entry_type=primitive
val_type=bool
[end]

[begin]
name=charge_mixing\.charge_preconditioning\.sw_metric_diff$
entry_type=primitive
val_type=bool
[end]

[begin]
name=charge_mixing\.spin_density$
entry_type = block
has_table = false
[end]

[begin]
name=charge_mixing\.spin_density\.sw_apply_metric$
entry_type = primitive
val_type=bool
[end]

[begin]
name=charge_mixing\.spin_density\.sw_apply_precon$
entry_type = primitive
val_type=bool
[end]

[begin]
name=charge_mixing\.spin_density\.sw_force_simple_mixing$
entry_type = primitive
val_type=bool
[end]

[begin]
name=charge_mixing\.spin_density\.sw_mix_bothspins_sametime$
entry_type = primitive
val_type=bool
[end]

# structure_evolution block
[begin]
name=structure_evolution$
entry_type=block
has_table=false
[end]

[begin]
name=structure_evolution\.method$
entry_type=primitive
val_type=string
choice=quench,gdiis,cg,velocity_verlet,temperature_control,damp,velocity_scaling,bfgs,cg2,temperature_pressure_control,pressure_temperature_control,pressure_control
invalid_choice = damp __or__ velocity_scaling __if__ control.driver __eq__ __undefined__ __or__ neb __or__ general
function = check_mobile __if__ control.condition __eq__ __undefined__ __or__ initial __or__ continuation __or__ automatic __or__ -1 __or__ 0 __or__ 1
[end]

[begin]
name=structure_evolution\.dt$
entry_type=primitive
val_type=float
unit_type=time
importance = recommended __if__ structure_evolution.method __eq__ quench __or__ velocity_verlet __or__ temperature_control __or__ damp __or__ velocity_scaling
function = check_dt __if__ structure_evolution.method __eq__ velocity_verlet __or__ temperature_control __or__ velocity_scaling
[end]

[begin]
name=structure_evolution\.stress$
entry_type=block
has_table=false
[end]

[begin]
name=structure_evolution\.stress\.sw_stress$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.stress\.iconstpw$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.stress\.sw_smear_ke$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.stress\.e0$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=structure_evolution\.stress\.sigma$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=structure_evolution\.stress\.a$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=structure_evolution\.stress\.sw_stress_correction$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.stress\.delta_ecut$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=structure_evolution\.n_md_step$
entry_type=primitive
val_type=int
importance = recommended __if__ control.driver __eq__ constraint __and__ structure_evolution.method __eq__ temperature_control __or__ velocity_verlet __or__ velocity_scaling
[end]

[begin]
name=structure_evolution\.gdiis$
entry_type=block
has_table=false
[end]

[begin]
name=structure_evolution\.gdiis\.gdiis_box_size$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.gdiis\.gdiis_hownew$
aliases=structure_evolution\.gdiis\.gdiis_update$
entry_type=primitive
val_type=string
choice=renew,anew
[end]

[begin]
name=structure_evolution\.gdiis\.c_forc2gdiis$
entry_type=primitive
val_type=float
unit_type=force
#importance = recommended __if__ structure_evolution.method __eq__ gdiis __or__ bfgs
[end]

[begin]
name=structure_evolution\.gdiis\.c_iteration2gdiis$
entry_type=primitive
val_type=int
#importance = recommended __if__ structure_evolution.method __eq__ gdiis
[end]

[begin]
name=structure_evolution\.gdiis\.initial_method$
entry_type=primitive
val_type=string
choice=quench,cg,sd,cg2
#importance = recommended __if__ structure_evolution.method __eq__ gdiis __or__ bfgs
[end]

[begin]
name=structure_evolution\.gdiis\.sw_correct_eigenvalue$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.temperature_control$
entry_type=block
has_table=false
[end]

[begin]
name=structure_evolution\.temperature_control\.method$
entry_type=primitive
val_type=string
choice=nose,nose_hoover,velocity_scaling
[end]

[begin]
name=structure_evolution\.temperature_control\.thermostat$
entry_type=block
has_table=true
[end]

[begin]
name=structure_evolution\.temperature_control\.thermostat\.table$
entry_type=table
columns = id,no,qmass,temp,tempi,tempf,till_n,tdamp
val_type = int,int,float,float,float,float,int,float
importance = required __if__ structure_evolution.method __eq__ temperature_control __or__ velocity_scaling __or__ temperature_pressure_control __or__ pressure_temperature_control
#function = check_qmass __if__ structure_evolution.method __eq__ temperature_control
all_columns_required=false
required_attribute = temp __if__ structure_evolution.temperature_control.sw_temperature_profile __eq__ False
required_attribute1 = tempi __if__ structure_evolution.temperature_control.sw_temperature_profile __eq__ True
required_attribute2 = tempf __if__ structure_evolution.temperature_control.sw_temperature_profile __eq__ True
required_attribute3 = till_n __if__ structure_evolution.temperature_control.sw_temperature_profile __eq__ True
#required_attribute1 = qmass __if__ structure_evolution.method __eq__ temperature_control
[end]

[begin]
name=structure_evolution\.temperature_control\.num_thermostat$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.temperature_control\.num_chain$
entry_type=primitive
val_type=int
val_range=1,
[end]

[begin]
name=structure_evolution\.temperature_control\.set_initial_velocity$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.temperature_control\.sw_temperature_profile$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.temperature_control\.tdamp$
entry_type=primitive
val_type=time
[end]

[begin]
name=structure_evolution\.pressure_control$
entry_type=block
has_table=false
importance=required __if__ structure_evolution.method __eq__ pressure_control __or__ temperature_pressure_control __or __ pressure_temperature_control
[end]

[begin]
name=structure_evolution\.pressure_control\.mass_baro$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.pressure_control\.pressure$
entry_type=primitive
val_type=float
unit_type=pressure
[end]

[begin]
name=structure_evolution\.pressure_control\.m11$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.pressure_control\.m22$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.pressure_control\.m33$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.pressure_control\.m12$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.pressure_control\.m13$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.pressure_control\.m23$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.damp$
entry_type=block
has_table=false
[end]

[begin]
name=structure_evolution\.damp\.resample_damping_parameters$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.damp\.automatic_dt$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.damp\.resample_period$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.damp\.div$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.damp\.dt_max$
entry_type=primitive
val_type=float
unit_type=time
[end]

[begin]
name=structure_evolution\.damp\.upper_limit_for_dt$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.molecule$
entry_type=block
has_table=false
[end]

[begin]
name=structure_evolution\.molecule\.dt_rotation\d*$
entry_type=primitive
val_type=float
unit_type=time
[end]

[begin]
name=structure_evolution\.molecule\.dt_translation\d*$
entry_type=primitive
val_type=float
unit_type=time
[end]

[begin]
name=structure_evolution\.molecule\.take_pbc_into_account\d*$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice$
entry_type=block
has_table=false
[end]

[begin]
name=structure_evolution\.lattice\.sw_optimize_lattice$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.sw_optimize_coordinates_once$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.sw_rebuild_pws$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.nhistory$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure_evolution\.lattice\.delta$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.max_stress$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.stress_convergence$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.sw_uniform$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.external_stress$
entry_type=block
has_table=false
[end]

[begin]
name=structure_evolution\.lattice\.external_stress\.s11$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.external_stress\.s22$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.external_stress\.s33$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.external_stress\.s12$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.external_stress\.s21$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.external_stress\.s13$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.external_stress\.s31$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.external_stress\.s23$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.external_stress\.s32$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure_evolution\.lattice\.fix_angle_alpha$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.fix_angle_beta$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.fix_angle_gamma$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.fix_length_a$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.fix_length_b$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.fix_length_c$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.sw_neglect_stress_offdiag$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.sw_read_nfchgt_prev_cell$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.sw_read_nfzaj_prev_cell$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.sw_interpolate_charge$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.lattice\.sw_interpolate_wfs$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.predictor$
entry_type=block
has_table=false
[end]

[begin]
name=structure_evolution\.predictor\.sw_charge_predictor$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.predictor\.sw_wf_predictor$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.predictor\.sw_extrapolate_charge$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure_evolution\.predictor\.rms_threshold$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=structure_evolution\.keep_symmetry_strict$
entry_type=primitive
val_type=bool
[end]

#postprocessing
[begin]
name=postprocessing$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.dos$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.dos\.sw_dos$
entry_type=primitive
val_type=bool
restriction = True __if__ postprocessing.ldos.sw_aldos __eq__ True || postprocessing.ldos.sw_layerdos __eq__ True || postprocessing.pdos.sw_pdos __eq__ True
[end]

[begin]
name=postprocessing\.dos\.sw_dos_gaussdistrib$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.dos\.dos_subroutine$
entry_type=primitive
val_type=int
choice=3,4,5
[end]

[begin]
name=postprocessing\.dos\.method$
entry_type=primitive
val_type=string
choice=gaussian,tetrahedron,tetrahedral,g,t
function=check_tetrahedral_dos __if__ postprocessing.dos.method __eq__ tetrahedral __or__ t
restriction = gaussian __or__ g __if__ postprocessing.dos.sw_dos_gaussdistrib __eq__ True
default_value = gaussian
[end]

[begin]
name=postprocessing\.dos\.deltae_dos$
aliases=postprocessing\.dos\.deltae$,postprocessing\.dos\.deltae_dos_gaussd$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=postprocessing\.dos\.variance$
aliases=postprocessing\.dos\.variance_dos_gaussd$,postprocessing\.dos\.variance_gaussd$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=postprocessing\.dos\.nwd_dos_window_width$
entry_type=primitive
val_type=int
[end]

[begin]
name=postprocessing\.dos\.energy_unit$
entry_type=primitive
val_type=int
choice=1,2
default_value=1
[end]

[begin]
name=postprocessing\.ldos$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.ldos\.sw_aldos$
entry_type=primitive
val_type=bool
invalid_choice = True __if__ control.condition __eq__ fixed_charge __or__ fixed_charge_continuation __or__ 2 __or__ 3 __or__ -3
[end]

[begin]
name=postprocessing\.ldos\.aldos$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.ldos\.aldos\.crtdst$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=postprocessing\.ldos\.aldos\.naldos_from$
entry_type=primitive
val_type=int
val_range=1,
[end]

[begin]
name=postprocessing\.ldos\.aldos\.naldos_to$
entry_type=primitive
val_type=int
val_range=1,
[end]

[begin]
name=postprocessing\.ldos\.sw_layerdos$
entry_type=primitive
val_type=bool
invalid_choice = True __if__ control.condition __eq__ fixed_charge __or__ fixed_charge_continuation __or__ 2 __or__ 3 __or__ -3
[end]

[begin]
name=postprocessing\.ldos\.layerdos$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.ldos\.layerdos\.slicing_way$
entry_type=primitive
val_type=string
choice=by_atomic_positions,regular_intervals
[end]

[begin]
name=postprocessing\.ldos\.layerdos\.deltaz$
entry_type=primitive
val_type=float
unit_type=length
importance = recommended __if__ postprocessing.ldos.layerdos.slicing_way __eq__ regular_intervals
[end]

[begin]
name=postprocessing\.ldos\.layerdos\.normal_axis$
entry_type=primitive
val_type=int
val_range=1,3
[end]

[begin]
name=postprocessing\.ldos\.layerdos\.crtdst$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=postprocessing\.ldos\.layerdos\.sw_ignore_inversion$
entry_type=primitive
val_type=bool
invalid_choice =True __if__ postprocessing.ldos.sw_layerdos __eq__ False __or__ __undefined__ || structure.symmetry.sw_inversion __eq__ False __or__ __undefined__
[end]

[begin]
name=postprocessing\.pdos$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.pdos\.sw_pdos$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.pdos\.sw_orb_popu$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.charge$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.charge\.sw_charge_rspace$
entry_type=primitive
val_type=bool
description=set this variable to 'on' in order to output the valecne charge-density in real space
default_value=off
[end]

[begin]
name=postprocessing\.charge\.filetype$
entry_type=primitive
val_type=string
choice=cube,density_only
description=specify the 'filetype' for the charge density output. it is strongly recommended to set this value to 'cube'.
default_value=density_only
importance=recommended __if__ postprocessing.charge.sw_charge_rspace __eq__ True
[end]

[begin]
name=postprocessing\.charge\.title$
entry_type=primitive
val_type=string
[end]

[begin]
name=postprocessing\.charge\.partial_charge$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.charge\.partial_charge\.sw_partial_charge$
entry_type=primitive
val_type=bool
description = set this variable to 'on' in order to output the partial charge. when this variable is set to 'on', it is recommended to configure the four related variables.
default_value=off
#restriction = False __if__ postprocessing.charge.sw_charge_rspace __eq__ __undefined__ __or__ False
invalid_choice = True __if__ postprocessing.charge.sw_charge_rspace __eq__ __undefined__ __or__ False
[end]

[begin]
name=postprocessing\.charge\.partial_charge\.erange_min$
entry_type=primitive
val_type=float
unit_type=energy
importance=recommended __if__ postprocessing.charge.partial_charge.sw_partial_charge __eq__ True
[end]

[begin]
name=postprocessing\.charge\.partial_charge\.erange_max$
entry_type=primitive
val_type=float
unit_type=energy
importance=recommended __if__ postprocessing.charge.partial_charge.sw_partial_charge __eq__ True
[end]

[begin]
name=postprocessing\.charge\.partial_charge\.erange_delta$
entry_type=primitive
val_type=float
unit_type=energy
importance=recommended __if__ postprocessing.charge.partial_charge.sw_partial_charge __eq__ True
[end]

[begin]
name=postprocessing\.charge\.partial_charge\.partial_charge_filetype$
aliases=postprocessing\.charge\.partial_charge\.filetype$
entry_type=primitive
val_type=string
choice=individual,separate
importance=recommended __if__ postprocessing.charge.partial_charge.sw_partial_charge __eq__ True
[end]

[begin]
name=postprocessing\.stm$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.stm\.sw_stm$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.stm\.sw_deficit_charge$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.wf$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.wf\.sw_wf_rspace$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.wf\.filetype$
entry_type=primitive
val_type=string
choice=cube,density_only
importance=recommended __if__ postprocessing.wf.sw_wf_rspace __eq__ True
[end]

[begin]
name=postprocessing\.wf\.eigenvalue$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.wf\.eigenvalue\.eigmin$
entry_type=primitive
val_type=float
unit_type=energy
importance=recommended __if__ postprocessing.wf.sw_wf_rspace __eq__ True
[end]

[begin]
name=postprocessing\.wf\.eigenvalue\.eigmax$
entry_type=primitive
val_type=float
unit_type=energy
importance=recommended __if__ postprocessing.wf.sw_wf_rspace __eq__ True
[end]

[begin]
name=postprocessing\.wannier$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.wannier\.sw_wannier$
entry_type=primitive
val_type=bool
#restriction=__undefined__ __or__ False __if__ accuracy.ksampling.method __eq__ __undefined__ __or__ mesh __or__ monk __or__ file __or__ directin
invalid_choice= True __if__ accuracy.ksampling.method __eq__ __undefined__ __or__ mesh __or__ monk __or__ file __or__ directin
[end]

[begin]
name=postprocessing\.wannier\.eps_grad$
entry_type=primitive
val_type=float
[end]

[begin]
name=postprocessing\.wannier\.dt$
entry_type=primitive
val_type=float
[end]

[begin]
name=postprocessing\.wannier\.max_iteration$
entry_type=primitive
val_type=int
[end]

[begin]
name=postprocessing\.wannier\.filetype$
entry_type=primitive
val_type=string
choice=cube,density_only
description=specify the 'filetype' for the wanier function output. it is strongly recommended to set this value to 'cube'.
importance=recommended __if__ postprocessing.wannier.sw_wannier __eq__ True
[end]

[begin]
name=postprocessing\.wannier\.sw_random_wannier$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.wannier\.sw_continue$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.wannier\.sw_wannier90$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.wannier\.seedname$
entry_type=primitive
val_type=string
[end]

[begin]
name=postprocessing\.wannier\.nb_wan90$
entry_type=primitive
val_type=int
[end]

[begin]
name=postprocessing\.wannier\.sw_use_hartpart_wan90$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.wannier\.sw_write_unk_file$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.wannier\.spin_component_wan90$
entry_type=primitive
val_type=int
[end]

[begin]
name=postprocessing\.polarization$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.polarization\.sw_bp_property$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.polarization\.property$
entry_type=primitive
val_type=string
choice=effective_charge,polarization,piezoelectric_const
[end]

[begin]
name=postprocessing\.frequency$
entry_type=primitive
val_type=int
[end]

[begin]
name=postprocessing\.sw_band_symmetry_analysis$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.workfunc$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.workfunc\.sw_workfunc$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.workfunc\.sw_add_xc_to_vloc$
entry_type=primitive
val_type=bool
[end]

#print level

[begin]
name=printoutlevel$
aliases=printlevel$
entry_type=block
has_table=false
[end]

[begin]
name=printoutlevel\.[ipri]*base$
aliases=printlevel\.[ipri]*base$
entry_type=primitive
val_type=int
choice=0,1,2,3
#importance=recommended
[end]

[begin]
name=printoutlevel\.[ipri]*parallel_debug$
aliases=printlevel\.[ipri]*parallel_debug$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*paradeb$
aliases=printlevel\.[ipri]*paradeb$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*velocity$
aliases=printlevel\.[ipri]*velocity$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*chargedensity$
aliases=printlevel\.[ipri]*chargedensity$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*wf$
aliases=printlevel\.[ipri]*wf$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*negativecharge$
aliases=printlevel\.[ipri]*negativecharge$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*force$
aliases=printlevel\.[ipri]*force$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*ekzaj$
aliases=printlevel\.[ipri]*ekzaj$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*chargemixing$
aliases=printlevel\.[ipri]*chargemixing$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*occup$
aliases=printlevel\.[ipri]*occup$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*md$
aliases=printlevel\.[ipri]*md$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*vloc$
aliases=printlevel\.[ipri]*vloc$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*fftmap$
aliases=printlevel\.[ipri]*fftmap$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*ipribetar$
aliases=printlevel\.[ipri]*ipribetar$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*ipriberry$
aliases=printlevel\.[ipri]*ipriberry$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*ipripao$
aliases=printlevel\.[ipri]*ipripao$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*dos$
aliases=printlevel\.[ipri]*dos$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*chargedensity$
aliases=printlevel\.[ipri]*chargedensity$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*parallel_debug$
aliases=printlevel\.[ipri]*parallel_debug$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*pulay$
aliases=printlevel\.[ipri]*pulay$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*pp$
aliases=printlevel\.[ipri]*pp$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*timing$
aliases=printlevel\.[ipri]*timing$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*solver$
aliases=printlevel\.[ipri]*solver$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*evdff$
aliases=printlevel\.[ipri]*evdff$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*rmm$
aliases=printlevel\.[ipri]*rmm$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*snl$
aliases=printlevel\.[ipri]*snl$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*gdiis$
aliases=printlevel\.[ipri]*gdiis$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*eigenvalue$
aliases=printlevel\.[ipri]*eigenvalue$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*spg$
aliases=printlevel\.[ipri]*spg$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*kp$
aliases=printlevel\.[ipri]*kp$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*matdiagon$
aliases=printlevel\.[ipri]*matdiagon$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*vlhxcq$
aliases=printlevel\.[ipri]*vlhxcq$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*phig$
aliases=printlevel\.[ipri]*phig$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*phonon$
aliases=printlevel\.[ipri]*phonon$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*totalcharge$
aliases=printlevel\.[ipri]*totalcharge$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*submat$
aliases=printlevel\.[ipri]*submat$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*strcfctr$
aliases=printlevel\.[ipri]*strcfctr$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*parallel$
aliases=printlevel\.[ipri]*parallel$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*inputfile$
aliases=printlevel\.[ipri]*inputfile$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*davidson$
aliases=printlevel\.[ipri]*davidson$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*jobstatus$
aliases=printlevel\.[ipri]*jobstatus$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*rb$
aliases=printlevel\.iprirb$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.[ipri]*predictor$
aliases=printlevel\.ipripredictor$
entry_type=primitive
val_type=int
choice=0,1,2,3
[end]

[begin]
name=printoutlevel\.jobstatus_option$
aliases=printlevel\.jobstatus_option$
entry_type=block
has_table=false
[end]

[begin]
name=printoutlevel\.jobstatus_option\.jobstatus_format$
aliases=printlevel\.jobstatus_option\.jobstatus_format$
entry_type=primitive
val_type=string
choice=tag,tag_line,table
[end]

[begin]
name=printoutlevel\.jobstatus_option\.jobstatus_series$
aliases=printlevel\.jobstatus_option\.jobstatus_series$
entry_type=primitive
val_type=bool
[end]

[begin]
name=printoutlevel\.timing_option$
aliases=printlevel\.timing_option$
entry_type=block
has_table=false
[end]

[begin]
name=printoutlevel\.timing_option\.num_subroutines$
aliases=printlevel\.timing_option\.num_subroutines$
entry_type=primitive
val_type=int
[end]

[begin]
name=printoutlevel\.timing_option\.cputime_diff$
aliases=printlevel\.timing_option\.cputime_diff$
entry_type=primitive
val_type=float
[end]

[begin]
name=printoutlevel\.timing_option\.sw_timing_2ndlevel$
aliases=printlevel\.timing_option\.sw_timing_2ndlevel$
entry_type=primitive
val_type=bool
[end]

[begin]
name=printoutlevel\.timing_option\.sw_flatten$
aliases=printlevel\.timing_option\.sw_flatten$
entry_type=primitive
val_type=bool
[end]

[begin]
name=printoutlevel\.timing_option\.sw_firstlevel_only$
aliases=printlevel\.timing_option\.sw_firstlevel_only$
entry_type=primitive
val_type=bool
[end]

[begin]
name=printoutlevel\.timing_option\.sw_details$
aliases=printlevel\.timing_option\.sw_details$
entry_type=primitive
val_type=bool
[end]

[begin]
name=printoutlevel\.timing_option\.measure_count_limit$
aliases=printlevel\.timing_option\.measure_count_limit$
entry_type=primitive
val_type=int
[end]

[begin]
name=printoutlevel\.negativecharge_option$
aliases=printlevel\.negativecharge_option$
entry_type=block
has_table=false
[end]

[begin]
name=printoutlevel\.negativecharge_option\.max_warnings$
aliases=printlevel\.negativecharge_option\.max_warnings$
entry_type=primitive
val_type=int
[end]

[begin]
name=printoutlevel\.n_fermi_vicinity$
aliases=printlevel\.n_fermi_vicinity$
entry_type=primitive
val_type=int
[end]

[begin]
name=printoutlevel\.magmom$
aliases=printlevel\.magmom$
entry_type=primitive
val_type=int
[end]

#phonon related 

[begin]
name=phonon$
entry_type=block
has_table=false
[end]

[begin]
name=phonon\.sw_phonon$
entry_type=primitive
val_type=bool
invalid_choice = True __if__ control.condition __eq__ fixed_charge __or__ fixed_charge_continuation __or__ 3 __or__ -3 __or__ 1
check_existence = F_FORCE __if__ phonon.sw_phonon __eq__ True && phonon.sw_calc_force __eq__ False
[end]

[begin]
name=phonon\.sw_calc_force$
entry_type=primitive
val_type=bool
[end]

[begin]
name=phonon\.sw_vibrational_modes$
entry_type=primitive
val_type=bool
[end]

[begin]
name=phonon\.displacement$
entry_type = primitive
val_type = float
unit_type = length
val_range=0,
[end]

[begin]
name=phonon\.force_calc$
entry_type=block
has_table=false
[end]

[begin]
name=phonon\.force_calc\.start$
entry_type=primitive
val_type=int
[end]

[begin]
name=phonon\.force_calc\.end$
entry_type=primitive
val_type=int
[end]

[begin]
name=phonon\.norder$
entry_type=primitive
val_type=int
[end]

[begin]
name=phonon\.sw_polynomial_fit$
entry_type=primitive
val_type=bool
[end]

[begin]
name=phonon\.point_group$
entry_type=primitive
val_type=string
[end]

[begin]
name=phonon\.method$
aliases=phonon\.calc_type$
entry_type=primitive
val_type=string
choice=zone_center,dos,band
default_value=zone_center
[end]

[begin]
name=phonon\.lattice$
entry_type=block
has_table=false
[end]

[begin]
name=phonon\.lattice\.l1$
entry_type=primitive
val_type=int
range=1,
importance=required __if__ phonon.method __eq__ dos __or__ band
[end]

[begin]
name=phonon\.lattice\.l2$
entry_type=primitive
val_type=int
range=1,
importance=required __if__ phonon.method __eq__ dos __or__ band
[end]

[begin]
name=phonon\.lattice\.l3$
entry_type=primitive
val_type=int
range=1,
importance=required __if__ phonon.method __eq__ dos __or__ band
[end]

[begin]
name=phonon\.special_points$
entry_type=block
has_table=true
[end]

[begin]
name=phonon\.special_points\.table$
entry_type=table
columns=no,k1,k2,k3
val_type=int,float,float,float
#importance = required __if__ phonon.method __eq__ band
[end]

[begin]
name=phonon\.symmetry_lines$
entry_type=block
has_table=true
[end]

[begin]
name=phonon\.symmetry_lines\.table$
entry_type=table
columns=k1,k2,num_division
val_type=int,int,int
#importance = required __if__ phonon.method __eq__ band
[end]

[begin]
name=phonon\.dos$
entry_type=block
has_table=false
[end]

[begin]
name=phonon\.dos\.mesh$
entry_type=block
has_table=false
[end]

[begin]
name=phonon\.dos\.mesh\.nx$
entry_type=primitive
val_type=int
importance = recommended __if__ phonon.method __eq__ dos
[end]

[begin]
name=phonon\.dos\.mesh\.ny$
entry_type=primitive
val_type=int
importance = recommended __if__ phonon.method __eq__ dos
[end]

[begin]
name=phonon\.dos\.mesh\.nz$
entry_type=primitive
val_type=int
importance = recommended __if__ phonon.method __eq__ dos
[end]

[begin]
name=phonon\.dos\.deltae$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=phonon\.sw_lo_to_splitting$
entry_type=primitive
val_type=bool
[end]

[begin]
name=phonon\.electronic_dielectric_constant$
entry_type=block
has_table=false
[end]

[begin]
name=phonon\.electronic_dielectric_constant\.exx$
entry_type=primitive
val_type=float
[end]

[begin]
name=phonon\.electronic_dielectric_constant\.eyy$
entry_type=primitive
val_type=float
[end]

[begin]
name=phonon\.electronic_dielectric_constant\.ezz$
entry_type=primitive
val_type=float
[end]

[begin]
name=phonon\.electronic_dielectric_constant\.exy$
entry_type=primitive
val_type=float
[end]

[begin]
name=phonon\.electronic_dielectric_constant\.eyz$
entry_type=primitive
val_type=float
[end]

[begin]
name=phonon\.electronic_dielectric_constant\.ezx$
entry_type=primitive
val_type=float
[end]

[begin]
name=phonon\.k_vector$
entry_type=block
has_table=false
[end]

[begin]
name=phonon\.k_vector\.kx$
entry_type=primitive
val_type=float
[end]

[begin]
name=phonon\.k_vector\.ky$
entry_type=primitive
val_type=float
[end]

[begin]
name=phonon\.k_vector\.kz$
entry_type=primitive
val_type=float
[end]

[begin]
name=phonon\.sw_lattice_dielectric_tensor$
entry_type=primitive
val_type=bool
[end]

[begin]
name=phonon\.sw_dielectric_function$
entry_type=primitive
val_type=bool
[end]

[begin]
name=phonon\.energy_range$
entry_type=block
has_table=false
[end]

[begin]
name=phonon\.energy_range\.min_energy$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=phonon\.energy_range\.max_energy$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=phonon\.energy_range\.division_number$
entry_type=primitive
val_type=int
[end]

#neb related

[begin]
name=multiple_replica$
entry_type=block
has_table=false
[end]

[begin]
name=multiple_replica\.method$
entry_type=primitive
val_type=string
choice=nudged_elastic_band_method
[end]

[begin]
name=multiple_replica\.accuracy$
entry_type=block
has_table=false
[end]

[begin]
name=multiple_replica\.accuracy\.dt\d*$
entry_type=primitive
val_type=float
[end]

[begin]
name=multiple_replica\.accuracy\.neb_time_integral$
entry_type=primitive
val_type=string
choice=quench,steepest_descent,12,2
[end]

[begin]
name=multiple_replica\.accuracy\.penalty_function$
entry_type=primitive
val_type=bool
[end]

[begin]
name=multiple_replica\.accuracy\.neb_convergence_condition$
entry_type=primitive
val_type=int __or__ string
choice=1,2,3,4,5,energy_e,phase_force,neb_force,force_at_transition_state,phase_force_normal
[end]

[begin]
name=multiple_replica\.accuracy\.neb_convergence_threshold$
entry_type=primitive
val_type=float
[end]

[begin]
name=multiple_replica\.constraint$
entry_type=block
has_table=false
[end]

[begin]
name=multiple_replica\.constraint\.ci_neb$
entry_type=primitive
val_type=bool
[end]

[begin]
name=multiple_replica\.constraint\.sp_k_init$
entry_type=primitive
val_type=float
[end]

[begin]
name=multiple_replica\.constraint\.sp_k_min$
entry_type=primitive
val_type=float
[end]

[begin]
name=multiple_replica\.constraint\.sp_k_max$
entry_type=primitive
val_type=float
[end]

[begin]
name=multiple_replica\.constraint\.sp_k_variable$
entry_type=primitive
val_type=bool
[end]

[begin]
name=multiple_replica\.structure$
entry_type=block
has_table=false
[end]

[begin]
name=multiple_replica\.structure\.number_of_replicas$
entry_type=primitive
val_type=int
val_range=1,
[end]

[begin]
name=multiple_replica\.structure\.replicas$
entry_type=block
has_table=true
[end]

[begin]
name=multiple_replica\.structure\.replicas\.table$
entry_type=table
columns=replica_number,howtogive_coordinates,end0,end1
val_type=int,string,int,int
choice2=proportional,file,from_endpoints,directin
choice3=0,-1
choice4=0,-1
importance=required __if__ control.driver __eq__ neb || control.multiple_replica_mode __eq__ True
all_columns_required=false
[end]

[begin]
name=multiple_replica\.structure\.endpoint_images$
entry_type=primitive
val_type=string
choice=directin,file
[end]

[begin]
name=multiple_replica\.structure\.howtogive_coordinates$
entry_type=primitive
val_type=string
choice=from_endpoint_images
[end]

[begin]
name=multiple_replica\.structure\.atom_list_end0$
entry_type=block
has_table=false
[end]

[begin]
name=multiple_replica\.structure\.atom_list_end0\.coordinate_system$
entry_type=primitive
val_type=string
choice=cartesian,internal,pucv
[end]

[begin]
name=multiple_replica\.structure\.atom_list_end0\.atoms$
entry_type=block
has_table=true
[end]

[begin]
name=multiple_replica\.structure\.atom_list_end0\.atoms\.table$
entry_type=table
columns=element,rx,ry,rz
val_type=string,float,float,float
importance=required __if__ control.driver __eq__ neb || control.multiple_replica_mode __eq__ True
all_columns_required=true
[end]

[begin]
name=multiple_replica\.structure\.atom_list_end1$
entry_type=block
has_table=false
[end]

[begin]
name=multiple_replica\.structure\.atom_list_end1\.coordinate_system$
entry_type=primitive
val_type=string
choice=cartesian,internal,pucv
[end]

[begin]
name=multiple_replica\.structure\.atom_list_end1\.atoms$
entry_type=block
has_table=true
[end]

[begin]
name=multiple_replica\.structure\.atom_list_end1\.atoms.table$
entry_type=table
columns=element,rx,ry,rz
val_type=string,float,float,float
importance=required __if__ control.driver __eq__ neb || control.multiple_replica_mode __eq__ True
all_columns_required=true
[end]

[begin]
name=multiple_replica\.structure\.atom_list_image\d+$
entry_type=block
has_table=false
[end]

[begin]
name=multiple_replica\.structure\.atom_list_image\d+\.coordinate_system$
entry_type=primitive
val_type=string
choice=cartesian,internal,pucv
[end]

[begin]
name=multiple_replica\.structure\.atom_list_image\d+\.atoms$
entry_type=block
has_table=true
[end]

[begin]
name=multiple_replica\.structure\.atom_list_image\d+\.atoms\.table$
entry_type=table
columns=element,rx,ry,rz
val_type=string,float,float,float
all_columns_required=true
[end]

[begin]
name=multiple_replica\.take_pbc_into_account$
entry_type=primitive
val_type=bool
[end]

[begin]
name=multiple_replica\.local_tangent$
entry_type=primitive
val_type=int
val_range=0,1
[end]

[begin]
name=multiple_replica\.end0_energy$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=multiple_replica\.end1_energy$
entry_type=primitive
val_type=float
unit_type=energy
[end]

#optimization with constraints and the blue-moon ensemble approach
[begin]
name=structure\.constrainable\d+$
entry_type=block
has_table=false
description=configure constraints under this block.
[end]

[begin]
name=structure\.constrainable\d+\.type$
entry_type=primitive
val_type = string
choice=bond_length,bond_angle,dihedral_angle,pos_vec,distance_from_pos,plane,center_of_mass,bond_length_diff,coordination_number,user_defined,bond_angle_diff
description=the type of constraint, ie, bond_length, bond_angle, dihedral_angle and etc.
[end]

[begin]
name=structure\.constrainable\d+\.atom\d+$
entry_type=primitive
val_type = int
description=the id of the atom which participates in the target constraint
[end]

[begin]
name=structure\.constrainable\d+\.mobile$
entry_type=primitive
val_type = bool
[end]

[begin]
name=structure\.constrainable\d+\.monitor$
entry_type=primitive
val_type = bool
[end]

[begin]
name=structure\.constrainable\d+\.reaction_coordinate$
entry_type=block
has_table=false
[end]

[begin]
name=structure\.constrainable\d+\.reaction_coordinate\.sw_reaction_coordinate$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure\.constrainable\d+\.reaction_coordinate\.init_value$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.constrainable\d+\.reaction_coordinate\.increment$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.constrainable\d+\.reaction_coordinate\.final_value$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.constrainable\d+\.plane$
entry_type=block
has_table=false
[end]

[begin]
name=structure\.constrainable\d+\.plane\.normx$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.constrainable\d+\.plane\.normy$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.constrainable\d+\.plane\.normz$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.constrainable\d+\.coordination_number$
entry_type=block
has_table=false
[end]

[begin]
name=structure\.constrainable\d+\.coordination_number\.kappa$
entry_type=primitive
val_type=float
[end]

[begin]
name=structure\.constrainable\d+\.coordination_number\.kappa_inv$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.constrainable\d+\.coordination_number\.rc$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=structure\.constrainable\d+\.distance_from_pos\.posx$
entry_type=primitive
val_type=int
unit_type=length
[end]

[begin]
name=structure\.constrainable\d+\.distance_from_pos\.posy$
entry_type=primitive
val_type=int
unit_type=length
[end]

[begin]
name=structure\.constrainable\d+\.distance_from_pos\.posz$
entry_type=primitive
val_type=int
unit_type=length
[end]

[begin]
name=structure\.reac_coord_generation$
entry_type=primitive
val_type=string
choice=via_file,via_input
[end]

[begin]
name=thermodynamic_integration$
entry_type=block
has_table=false
[end]

[begin]
name=thermodynamic_integration\.nsteps$
entry_type=primitive
val_type=int
val_range=1,
[end]

[begin]
name=thermodynamic_integration\.nequib$
entry_type=primitive
val_type=int
val_range=1,
[end]

[begin]
name=thermodynamic_integration\.istart_reac_coords$
entry_type=primitive
val_type=int
val_range=1,
[end]

[begin]
name=thermodynamic_integration\.nreac_coords$
entry_type=primitive
val_type=int
val_range=1,
[end]

[begin]
name=thermodynamic_integration\.nsample$
entry_type=primitive
val_type=int
val_range=1,
[end]

[begin]
name=thermodynamic_integration\.smooth$
entry_type=primitive
val_type=bool
[end]

[begin]
name=thermodynamic_integration\.basedir$
entry_type=primitive
val_type=string
[end]

#meta-dynamics

[begin]
name=meta_dynamics$
entry_type=block
has_table=false
description=enter meta-dynamics related parameters under this block
[end]

[begin]
name=meta_dynamics\.meta_dynamics_type$
entry_type=primitive
val_type=string
choice=bias_only,bias_and_fictitious,bias_generation
description = specify the 'type' of the meta-dynamics simulation to be performed. by specifing bias_only, the bias potential will be directly added to the system, while by specifing bias_and_fictitious, the meta-dynamics will be peformed through the coordinates and velocity of a fictitous particle. when bias_generation is specified, meta-dynamics will not be performed, only the generation of the bias potential will be done.
[end]

[begin]
name=meta_dynamics\.max_bias_update$
entry_type=primitive
val_type=int
description = enter the maximum number of bias potential update the program should perform. if a negative value is specified, the program will not terminate according to this parameter.
default_value=-1
[end]

[begin]
name=meta_dynamics\.extensive_output$
entry_type=primitive
val_type=bool
[end]

[begin]
name=meta_dynamics\.output_per_rank$
entry_type=primitive
val_type=bool
[end]

[begin]
name=meta_dynamics\.output_cvar_every_step$
entry_type=primitive
val_type=bool
[end]

[begin]
name=meta_dynamics\.continuation_strategy$
entry_type=block
has_table=false
[end]

[begin]
name=meta_dynamics\.continuation_strategy\.randomize_velocity$
entry_type=primitive
val_type=bool
[end]

[begin]
name=meta_dynamics\.continuation_strategy\.scale_velocity$
entry_type=primitive
val_type=bool
[end]

[begin]
name=meta_dynamics\.continuation_strategy\.configuration_from_input$
entry_type=primitive
val_type=bool
[end]

[begin]
name=meta_dynamics\.continuation_strategy\.velocity_scaling_factor$
entry_type=primitive
val_type=float
[end]

[begin]
name=meta_dynamics\.collective_variable$
entry_type=block
has_table=false
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\.mass$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\.k$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\.delta_s$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\.mass_thermo$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\.target_ke$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\.control_velocity$
entry_type=primitive
val_type=bool
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\d+$
entry_type=block
has_table=false
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.type$
entry_type=primitive
val_type = string
choice=bond_length,bond_angle,dihedral_angle,pos_vec,distance_from_pos,plane,center_of_mass,bond_length_diff,coordination_number,user_defined,bond_angle_diff
description=the type of the collective variable, ie, bond_length, bond_angle, dihedral_angle and etc.
importance = required __if__ control.driver __eq__ meta_dynamics
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.atom\d+$
entry_type=primitive
val_type = int
description=the id of the atom which participates for the target collective variable
importance = required __if__ control.driver __eq__ meta_dynamics
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.plane$
entry_type=block
has_table=false
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.plane\.normx$
entry_type=primitive
val_type=float
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.plane\.normy$
entry_type=primitive
val_type=float
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.plane\.normz$
entry_type=primitive
val_type=float
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.coordination_number$
entry_type=block
has_table=false
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.coordination_number\.kappa$
entry_type=primitive
val_type=float
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.coordination_number\.kappa_inv$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.coordination_number\.rc$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.distance_from_pos\.posx$
entry_type=primitive
val_type=int
unit_type=length
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.distance_from_pos\.posy$
entry_type=primitive
val_type=int
unit_type=length
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.distance_from_pos\.posz$
entry_type=primitive
val_type=int
unit_type=length
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.dihedral$
entry_type=block
has_table=false
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.dihedral.abs$
entry_type=primitive
val_type=bool
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.mass$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.k$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.mass_thermo$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.target_ke$
entry_type=primitive
val_type=float
description=
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.control_velocity$
entry_type=primitive
val_type=bool
description=
[end]
[begin]
name=meta_dynamics\.collective_variable\d+\.delta_s$
entry_type=primitive
val_type=float
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.smin$
entry_type=primitive
val_type=float
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.smax$
entry_type=primitive
val_type=float
[end]

[begin]
name=meta_dynamics\.collective_variable\d+\.ds$
entry_type=primitive
val_type=float
[end]

[begin]
name=meta_dynamics\.bias_potential$
entry_type=block
has_table=false
[end]

[begin]
name=meta_dynamics\.bias_potential\.height$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=meta_dynamics\.bias_potential\.update_frequency$
entry_type=primitive
val_type=int
[end]

[begin]
name=meta_dynamics\.bias_potential\.output_frequency$
entry_type=primitive
val_type=int
[end]

#epsmain-related

[begin]
name=control\.use_additional_projector
entry_type=primitive
val_type=bool
restriction = True __if__ epsilon.transition_moment.type __eq__ ks
[end]

[begin]
name=epsilon$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.sw_epsilon$
entry_type=primitive
val_type=bool
restriction = False __if__ control.condition __eq__ initial __or__ continuation __or__ 0 __or__ 1
[end]

[begin]
name=epsilon\.crystal_type$
entry_type=primitive
val_type=string
choice=single,poly
[end]

[begin]
name=epsilon\.fermi_energy$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.fermi_energy\.read_efermi$
entry_type=primitive
val_type=bool
[end]

[begin]
name=epsilon\.fermi_energy\.efermi$
entry_type=primitive
val_type=float
importance=required __if__ epsilon.fermi_energy.read_efermi __eq__ True
[end]

[begin]
name=epsilon\.photon$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.photon\.polar$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.photon\.polar\.ux$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.photon\.polar\.uy$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.photon\.polar\.uz$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.photon\.pointing$
aliases=epsilon\.photon\.poynting$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.photon\.pointing\.px$
aliases=epsilon\.photon\.poynting\.px$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.photon\.pointing\.py$
aliases=epsilon\.photon\.poynting\.py$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.photon\.pointing\.pz$
aliases=epsilon\.photon\.poynting\.pz$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.photon\.energy$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.photon\.energy\.low$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=epsilon\.photon\.energy\.high$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=epsilon\.photon\.energy\.step$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=epsilon\.transition_moment$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.transition_moment\.type$
entry_type=primitive
val_type=string
choice=l,rn,ks,mks
[end]


[begin]
name=epsilon\.transition_moment\.delq$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.transition_moment\.symmetry$
entry_type=primitive
val_type=bool
[end]

[begin]
name=epsilon\.transition_moment\.band_i$
entry_type=primitive
val_type=int
[end]

[begin]
name=epsilon\.transition_moment\.band_f$
entry_type=primitive
val_type=int
[end]

[begin]
name=epsilon\.mass$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.mass\.sw_mass$
entry_type=primitive
val_type=bool
[end]

[begin]
name=epsilon\.mass\.shift$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.mass\.point$
entry_type=primitive
val_type=string
choice=band_edge,input
[end]

[begin]
name=epsilon\.mass\.ik$
entry_type=primitive
val_type=int
[end]

[begin]
name=epsilon\.mass\.ib$
entry_type=primitive
val_type=int
[end]

[begin]
name=epsilon\.mass\.direction$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.mass\.direction\.nx$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.mass\.direction\.ny$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.mass\.direction\.nz$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.bz_integration$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.bz_integration\.method$
entry_type=primitive
val_type=string
choice=parabolic,gaussian,tetrahedron,p,g,t
function=check_tetrahedral_dos __if__ epsilon.bz_integration.method __eq__ tetrahedron __or__ t
[end]

[begin]
name=epsilon\.bz_integration\.spin$
entry_type=primitive
val_type=string
choice=both,major,minor
[end]

[begin]
name=epsilon\.band_gap_correction$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.band_gap_correction\.scissor_operator$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=epsilon\.drude_term$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.drude_term\.drude$
entry_type=primitive
val_type=bool
[end]

[begin]
name=epsilon\.drude_term\.effective_mass$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.drude_term\.damping_factor$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=epsilon\.drude_term\.plasma_frequency$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=epsilon\.drude_term\.conductivity$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.ipriepsilon$
entry_type=primitive
val_type=int
val_range=0,1,2,3
[end]

[begin]
name=epsilon\.nonlinear_optics$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.nonlinear_optics\.process$
entry_type=primitive
val_type=string
choice=off,shg,thg
[end]

[begin]
name=epsilon\.nonlinear_optics\.band$
entry_type=primitive
val_type=string
choice=all,inter,intra
[end]

[begin]
name=epsilon\.nonlinear_optics\.excitation$
entry_type=primitive
val_type=string
choice=all,omega,2omega,3omega
[end]

[begin]
name=epsilon\.nonlinear_optics\.term$
entry_type=primitive
val_type=string
choice=all,omega,2omega,3omega
[end]

[begin]
name=epsilon\.nonlinear_optics\.dres_cut_off$
entry_type=primitive
val_type=float
[end]

[begin]
name=epsilon\.nonlinear_optics\.smearing_fact$
entry_type=primitive
val_type=string
choice=resonance,off_resonance
[end]

[begin]
name=epsilon\.nonlinear_optics\.double_resonance$
aliases=epsilon\.nonlinear_optics\.dble_resonance$
entry_type=block
has_table=false
[end]

[begin]
name=epsilon\.nonlinear_optics\.double_resonance\.method$
aliases=epsilon\.nonlinear_optics\.dble_resonance\.method$
entry_type=primitive
val_type=string
choice=omit,damping
[end]

[begin]
name=epsilon\.nonlinear_optics\.double_resonance\.cut_off$
aliases=epsilon\.nonlinear_optics\.dble_resonance\.cut_off$
entry_type=primitive
val_type=float
unit_type=energy
[end]

#uvsor-berry-phonon

[begin]
name=berry_phase$
entry_type=block
has_table=false
[end]

[begin]
name=berry_phase\.sw_berry_phase$
entry_type=primitive
val_type=bool
[end]

[begin]
name=berry_phase\.g_index$
entry_type=primitive
val_type=int
val_range=1,3
[end]

[begin]
name=berry_phase\.mesh$
entry_type=block
has_table=false
[end]

[begin]
name=berry_phase\.mesh\.n1$
entry_type=primitive
val_type=int
[end]

[begin]
name=berry_phase\.mesh\.n2$
entry_type=primitive
val_type=int
[end]

[begin]
name=berry_phase\.mesh\.j$
entry_type=primitive
val_type=int
[end]

[begin]
name=structure\.atom_list\.displacement$
entry_type=block
has_table=false
[end]

[begin]
name=structure\.atom_list\.displacement\.sw_displace_atom$
entry_type=primitive
val_type=bool
[end]

[begin]
name=structure\.atom_list\.displacement\.displaced_atom$
entry_type=primitive
val_type=int
importance=required __if__ atom_list.displacement.sw_displace_atom __eq__ True
[end]

[begin]
name=structure\.atom_list\.displacement\.ux$
entry_type=primitive
val_type=float
unit_type=length
importance=required __if__ atom_list.displacement.sw_displace_atom __eq__ True
[end]

[begin]
name=structure\.atom_list\.displacement\.uy$
entry_type=primitive
val_type=float
unit_type=length
importance=required __if__ atom_list.displacement.sw_displace_atom __eq__ True
[end]

[begin]
name=structure\.atom_list\.displacement\.uz$
entry_type=primitive
val_type=float
unit_type=length
importance=required __if__ atom_list.displacement.sw_displace_atom __eq__ True
[end]

#LR-TDDFT
[begin]
name=spectrum$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.type$
entry_type=primitive
val_type=string
choice=optics,pacs,eels,ixss
[end]

[begin]
name=spectrum\.longwavelimit$
entry_type=primitive
val_type=bool
[end]

[begin]
name=spectrum\.momentum_transfer$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.momentum_transfer\.deltaq$
entry_type=primitive
val_type=float
[end]

[begin]
name=spectrum\.momentum_transfer\.nx$
entry_type=primitive
val_type=float
[end]

[begin]
name=spectrum\.momentum_transfer\.ny$
entry_type=primitive
val_type=float
[end]

[begin]
name=spectrum\.momentum_transfer\.nz$
entry_type=primitive
val_type=float
[end]

[begin]
name=spectrum\.tddft$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.tddft\.sw_tddft$
entry_type=primitive
val_type=bool
[end]

[begin]
name=spectrum\.tddft\.solver$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.tddft\.solver\.equation$
entry_type=primitive
val_type=string
choice=dyson,bs
[end]

[begin]
name=spectrum\.tddft\.solver\.sw_noda$
entry_type=primitive
val_type=bool
[end]

[begin]
name=spectrum\.tddft\.xc_kernel$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.tddft\.xc_kernel\.kernel_type$
entry_type=primitive
val_type=string
choice=rpa,alda-g,alda-r,lrc
[end]

[begin]
name=spectrum\.tddft\.xc_kernel\.lrc_alpha$
entry_type=primitive
val_type=float
[end]

[begin]
name=spectrum\.tddft\.coulomb_kernel$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.tddft\.coulomb_kernel\.sw_nlf$
entry_type=primitive
val_type=bool
[end]

[begin]
name=spectrum\.tddft\.expansion$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.tddft\.expansion\.numgvec$
entry_type=primitive
val_type=int
[end]

[begin]
name=spectrum\.energy$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.energy\.low$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=spectrum\.energy\.high$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=spectrum\.energy\.step$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=spectrum\.bz_integration$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.bz_integration\.method$
entry_type=primitive
val_type=string
choice=lorentzian,l,gaussiang,tetrahedron,t
[end]

[begin]
name=spectrum\.bz_integration\.width$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=spectrum\.band_gap_correction$
entry_type=block
has_table=false
[end]

[begin]
name=spectrum\.band_gap_correction\.scissor_operator$
entry_type=primitive
val_type=float
unit_type=energy
[end]

#ESM
[begin]
name=accuracy\.esm$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.esm\.sw_esm$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.esm\.z1$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=accuracy\.esm\.w$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=accuracy\.esm\.fix_ef$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=accuracy\.esm\.add_elec$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.esm\.z_wall$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=accuracy\.esm\.bar_height$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=accuracy\.esm\.bar_width$
entry_type=primitive
val_type=float
unit_type=length
[end]

[begin]
name=accuracy\.esm\.electric_field$
entry_type=primitive
val_type=float
unit_type=force
[end]

[begin]
name=accuracy\.esm\.nosmooth$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.esm\.external_potential$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.esm\.gps$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.esm\.gpe$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.esm\.bc$
entry_type=primitive
val_type=string
choice=bare,pe1,pe2
[end]

#FCP
[begin]
name=accuracy\.fcp$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.fcp\.sw_fcp$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.fcp\.mu$
entry_type=primitive
val_type=float
unit_type=energy
[end]

[begin]
name=accuracy\.fcp\.mass$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.fcp\.temperature$
entry_type=primitive
val_type=float
unit_type=temperature
[end]

[begin]
name=accuracy\.fcp\.relax_step$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.fcp\.relax_crit$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.fcp\.tot_charge_first$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.fcp\.tot_charge_last$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.fcp\.qmass$
entry_type=primitive
val_type=float
[end]

#paw one center integral
[begin]
name=paw_one_center_integral$
entry_type=block
has_table=false
[end]

[begin]
name=paw_one_center_integral\.element_list$
entry_type=block
has_table=true
[end]

[begin]
name=paw_one_center_integral\.element_list\.table$
entry_type=table
columns=element,surface_integral_method
val_type=string,string
choice2=gl,sphericalharmonicsexpansion,sphericalharmonics,sphex,gausslegendre
[end]

#metagga
[begin]
name=accuracy\.metagga$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.metagga\.ekin_density_type$
entry_type=primitive
val_type=int
[end]

[begin]
name=accuracy\.metagga\.sw_rspace_ekin_density$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.metagga\.sw_calc_ekindens_hardpart$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.metagga\.sw_add_ekin_hard_on_gspace$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.metagga\.val_c_tb09$
entry_type=primitive
val_type=float
[end]

#scdft
[begin]
name = control\.max_scfdft_iteration$
entry_type = primitive
val_type = int
description = enter the maximum number of SC-DFT iterations per md step. 
[end]

[begin]
name=accuracy\.scdft$
aliases=accuracy\.sc_dft$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.scdft\.delta_epsilon$
aliases=accuracy\.sc_dft$\.delta_epsilon$
entry_type=primitive
val_type=float
[end]

#RTTDDFT
[begin]
name=postprocessing\.rttddft$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.rttddft\.sw_rttddft$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.rttddft\.time_step_max$
entry_type=primitive
val_type=int
[end]

[begin]
name=postprocessing\.rttddft\.time_step_delta$
entry_type=primitive
val_type=float
unit_type=time
[end]

[begin]
name=postprocessing\.rttddft\.propagator_method$
entry_type=primitive
val_type=int
range=1,2
[end]

[begin]
name=postprocessing\.rttddft\.propagator_order$
entry_type=primitive
val_type=int
[end]

[begin]
name=postprocessing\.rttddft\.ext_ie_elec$
entry_type=primitive
val_type=int
[end]

[begin]
name=postprocessing\.rttddft\.ext_ie_hole$
entry_type=primitive
val_type=int
[end]

[begin]
name=postprocessing\.rttddft\.ext_pulse_epsilon$
entry_type=primitive
val_type=float
[end]

[begin]
name=postprocessing\.rttddft\.ext_pulse_kx$
entry_type=primitive
val_type=float
[end]

[begin]
name=postprocessing\.rttddft\.ext_pulse_ky$
entry_type=primitive
val_type=float
[end]

[begin]
name=postprocessing\.rttddft\.ext_pulse_kz$
entry_type=primitive
val_type=float
[end]

#nonlocal in rs
[begin]
name=accuracy\.nonlocal_potential$
entry_type=block
has_table=false
[end]

[begin]
name=accuracy\.nonlocal_potential\.sw_rspace$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.nonlocal_potential\.sw_rspace_v$
entry_type=primitive
val_type=bool
[end]

[begin]
name=accuracy\.nonlocal_potential\.r0_factor$
entry_type=primitive
val_type=float
range=1,
[end]

[begin]
name=accuracy\.nonlocal_potential\.gamma_factor$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.nonlocal_potential\.dq$
entry_type=primitive
val_type=float
[end]

[begin]
name=accuracy\.nonlocal_potential\.projector_optimization$
entry_type=primitive
val_type=string
choice=mask_function,prefitting
[end]

#BoltzTraP
[begin]
name=postprocessing\.boltztrap$
entry_type=block
has_table=false
[end]

[begin]
name=postprocessing\.boltztrap\.sw_boltztrap$
entry_type=primitive
val_type=bool
[end]

[begin]
name=postprocessing\.boltztrap\.prefix$
entry_type=primitive
val_type=string
[end]

[begin]
name=postprocessing\.boltztrap\.header$
entry_type=primitive
val_type=string
[end]

[begin]
name=postprocessing\.boltztrap\.version$
entry_type=primitive
val_type=string
choice=1,2
[end]

