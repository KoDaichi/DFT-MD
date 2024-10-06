module m_timer
  implicit none
  integer, parameter :: dp = kind(0d0)
  integer, parameter :: n_tim = 1750
  integer :: iunit1 = 91
  integer :: iunit2 = 92
#ifdef USE_MEM
  integer :: unitmw = 93
  integer :: unitmr = 94
#endif
!--
!--
  integer, parameter :: n_ite = 10
!--
#ifdef USE_WTIME
  real(dp) :: t_now
  real(dp) :: timers(n_tim,0:n_ite), tsta(n_tim), ttot(3)
  real(dp) :: timers_all(n_tim)
#ifdef USE_MEM
  real(dp) :: t_sta
#endif
#else
  integer(8) :: t_now
  integer(8) :: timers(n_tim,0:n_ite), tsta(n_tim), ttot(3)
  integer(8) :: timers_all(n_tim)
#ifdef USE_MEM
  integer(8) :: t_sta
#endif
#endif
  integer  :: ncount(n_tim,0:n_ite), ncnt_e(n_tim,0:n_ite)
  integer  :: ncount_all(n_tim)
  character(50) :: name(n_tim)
  logical  :: tflag = .false.
  character(len=3) :: ite_num
  character(len=5) :: rnk_num
  character(len=80) :: fname
  character(len=80) :: dname
  character(len=5)  :: crank
  integer :: myrank  
  integer :: orank = 0

end module m_timer

!-----------------------------------------------------------------------
subroutine timer_init
  use m_timer
  implicit none
  include 'mpif.h'
  integer :: ierr
  real(4) :: rdum, rarray(2)
!  real(4), external :: etime
  real(4) :: etime

  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

#ifdef USE_MEM
  if (myrank .eq. 0) then
    open(unitmw,file='mem.txt',status='replace')
    write(unitmw,*) "Kukan Times[s] Virtual_Memory[B] Resident_Set_Size[Page]"
#ifdef USE_WTIME
    t_sta=MPI_WTIME()
#else
    call getclkreg(t_sta)
#endif
  endif
#endif

#ifdef USE_KTRACE
      integer ktrace_target(n_tim)
      common /c_ktrace_target/ ktrace_target

      ! init
      ktrace_target(1:n_tim) = 0
      ! set ktarce target
      ktrace_target(1:n_tim) = 1
!$omp parallel
      call ktrace_enter
!$omp end parallel

#endif

  ! real, user, sys of total
  rdum = etime(rarray)
#ifdef USE_WTIME
  ttot(1) = MPI_WTIME()  ! real
#else
  call getclkreg(ttot(1))! real
#endif
  ttot(2) = rarray(1)    ! user
  ttot(3) = rarray(2)    ! sys

  timers = 0.0d0
  tsta   = 0.0d0
  ncount = 0
  ncnt_e = 0

  name = '-'
  ! set timer name
! name(000) = '---------1---------2---------3---------4---------5---------6----'
  name( 1) =  'Kukan_1                                 '
  name( 2) =  'Kukan_2                                 '
  name( 3) =  'Kukan_3                                 '
  name( 4) =  'Kukan_4                                 '
  name( 5) =  'Kukan_5                                 '
  name( 6) =  'Kukan_6                                 '
  name( 7) =  'Kukan_7                                 '
!fj mod 20110819
! name( 8) =  'Kukan_8_latter_half_&_kukan_9           '
  name( 8) =  'Kukan_8_latter_half                     '
  name( 9) =  'kukan_9                                 '
!fj mod 20110819
  name(10) =  'Kukan_10                                '
  name(11) =  'Kukan_11                                '
  name(12) =  '  Solver_SD_or_MSD_(_or_CG_)            '
  name(13) =  '  Solver_lmSD_or_lmMSD_(_or_CG_)        '
  name(14) =  '  Solver_Submat                         '
  name(15) =  '  Solver_MATRIXDIAGON                   '
  name(16) =  'SCF_Loop                                '
  name(17) =  'Initialization                          '
  name(18) =  '  Solver_RMM                            '
  name(19) =  'Forces                                  '
  name(20) =  'Move_Ions                                '
!fj mod 20110906
  name(1498) = 'Parallel_Region_Overhead                 '
  name(1499) = 'Parallel_Region_Overhead                 '
  name(1500) = 'kukan_9_not_scalapack                    '
!fj mod 20110906
!----- FJ sokutei
  name(21) = 'Initialization'
  name(22) = 'InputData_Analysis'
  name(23) = 'Preparation'
  name(24) = 'Preparation_for_mpi'
  name(25) = 'PseudoPotential_Construction'
  name(26) = 'Ewald_and_Structure_Factor'
  name(27) = 'Initial_Electronic_Structure'
  name(28) = 'INI_1_Gvec'
  name(29) = 'INI_2_PP'
  name(30) = 'INI_3_WF'
  name(31) = 'INI_4_CS'
  name(32) = 'm_ES_modified_gram_schmidt' ! kukan 4,5
  name(33) = 'm_ES_energy_eigen_values'  ! kukan 6
  name(34) = 'm_XC_cal_potential'  !kukan 7
  name(35) = 'm_ESlhxc_potential'  !kukan 11
  name(36) = 'm_ESiVQ_integrate_VlhxcQlm'  !kukan 11
!-----
  name(37) = 'PARA_2Dto3D'  !interface
  name(38) = 'PARA_3Dto2D'  !interface
  name(39) = 'WriteDownData_onto_Files'
!-----
  name(40) = 'm_ESIW_by_randomnumbers' 
  name(41) = 'm_pwBS_generate_G_vectors'
  name(42) = 'count_Gvectors_in_spheres'
  name(43) = 'get_Gvectors_in_a_gmaxp_sphere'
  name(44) = 'm_pwBS_generate_G_vectors_3D:allgather'
!---
! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- MATDIAGON / SUBROUTINE
  name(1201) = 'm_ESmat_solve_Hx_eq_eSx_3D'
  name(1202) = 'wd_hsmat'
  name(1203) = 'solve_Hx_eq_ex_LAPACK'
  name(1204) = 'calphase'
  name(1205) = 'nlmtrx'
  name(1206) = 'ovmtrx'
  name(1207) = 'cholde'
  name(1208) = 'amtrxkinetic'
  name(1209) = 'prjvlhxc_l2vlhxc_t'
  name(1210) = 'amtrxlocalpt'
  name(1211) = 'lhl'
  name(1212) = 'lx'
  name(1213) = 'phase_mult'
  name(1214) = 'cp_zaj_mat_to_zaj_l'
  name(1470) = 'make_index_band_for_scalapack'
!---- MATDIAGON add scalapack
  name(1471) = 'solve_Hx_eq_ex_ScaLAPACK'
  name(1472) = 'scalapack_setup_3D'
  name(1473) = 'make_index_band_for_scalapack_md'
  name(1474) = 'trans_scalapack                         '
  name(1475) = 'trans_scalapack_r                       '
  name(1476) = 'm_ESmat_solve_Hx_eq_eSx_3D:1            '
  name(1477) = 'm_ESmat_solve_Hx_eq_eSx_3D:2            '
  name(1478) = 'm_ESmat_solve_Hx_eq_eSx_3D:3            '
  name(1479) = 'solve_Hx_eq_ex_ScaLAPACK:PDSYGVX        '
  name(1480) = 'solve_Hx_eq_ex_ScaLAPACK:PDSYEVX        '
  name(1481) = 'solve_Hx_eq_ex_ScaLAPACK:DO_1           '
  name(1482) = 'solve_Hx_eq_ex_ScaLAPACK:PZHEGVX        '
  name(1483) = 'solve_Hx_eq_ex_ScaLAPACK:PZHEEVX        '
  name(1484) = 'solve_Hx_eq_ex_ScaLAPACK:DO_2           '
  name(1485) = 'trans_scalapack:DO_1                    '
  name(1486) = 'trans_scalapack:DO_2                    '
  name(1487) = 'trans_scalapack:DO_3                    '
  name(1488) = 'trans_scalapack:DO_4                    '
  name(1489) = 'trans_scalapack:isend/irecv_1           '
  name(1490) = 'trans_scalapack_r:DO_1                  '
  name(1491) = 'trans_scalapack_r:DO_2                  '
  name(1492) = 'trans_scalapack_r:DO_3                  '
  name(1493) = 'trans_scalapack_r:DO_4                  '
  name(1494) = 'trans_scalapack_r:isend/irecv_1         '
! -- Initialization / SUBROUTINE
  name(1215) = 'm_ES_betar_dot_WFs_3D'
  name(1216) = 'cpzaj_l_to_psi_ri'
  name(1217) = 'set_msize'
  name(1218) = 'pre_lmta_k_blk'
  name(1219) = 'post_lmta_k_blk'
  name(1220) = 'm_pwBS_calc_length_of_G_3D'
  name(1221) = 'm_pwBS_for_each_WF'
  name(1222) = 'm_PP_vanderbilt_type_3D'
  name(1223) = 'm_PP_local_part_3D'
  name(1224) = 'm_NLP_betar_dot_PWs_3D'
  name(1225) = 'm_CD_cp_chgq_to_chgqo_3D'
  name(1226) = 'm_CD_initial_CD_by_Gauss_func_3D'
  name(1227) = 'm_CS_gnrt_symmetry_operations'
  name(1228) = 'm_IS_symm_check_of_pos'
  name(1229) = 'm_Parallel_init_mpi_elec_3D'
  name(1230) = 'make_index_band_3D'
  name(1231) = 'make_index_band_for_Gdiv_3D'
  name(1232) = 'm_Parallel_mpi_fft_box'
  name(1233) = 'm_Parallel_mpi_fft_box_cd'
  name(1234) = 'm_Parallel_init_mpi_kngp_B_3D'
  name(1235) = 'm_Parallel_init_mpi_nbmx'
  name(1236) = 'im_Parallel_init_mpi_gga'
  name(1237) = 'm_Parallel_init_mpi_snl_3D'
  name(1238) = 'm_Parallel_init_mpi_atm_3D'
  name(1239) = 'm_Parallel_init_mpi_atm2_3D'
  name(1240) = 'm_Parallel_init_mpi_atm_B_3D'
  name(1241) = 'm_Parallel_init_mpi_mix_3D'
  name(1242) = 'm_Parallel_init_mpi_iba_3D'
  name(1243) = 'm_Parallel_wf_onto_fft_3D'
  name(1244) = 'm_Parallel_fft_onto_wf_3D'
  name(1245) = 'm_Parallel_fft_onto_chgq_3D'
  name(1246) = 'm_Parallel_chgq_onto_fftcd_3D'
  name(1247) = 'm_Parallel_fftcd_onto_chgq_3D'
  name(1248) = 'm_IS_structure_factor_3D'
  name(1249) = 'structure_factor1'
  name(1250) = 'structure_factor2'
  name(1251) = 'm_IS_ewald_3D'
  name(1252) = 'decide_rxyz_size'
  name(1253) = 'substitute_rxyz'
  name(1254) = 'decide_alf'
  name(1255) = 'decide_newldg'
  name(1256) = 'cpspac'
  name(1257) = 'set_ewald_parameters'
  name(1258) = 'ewald_Rspace_summation'
  name(1259) = 'ewald_Gspace_summation'
  name(1260) = 'get_zsum'
  name(1261) = 'add_exp_G2_zsum'
  name(1262) = 'ewald_force_Gspace_summation'

! -- MATDIAGON / COMMUNICATION
  name(1263) = 'solve_Hx_eq_ex_LAPACK:allreduce1'
  name(1264) = 'solve_Hx_eq_ex_LAPACK:allreduce2'
  name(1265) = 'nlmtrx:allreduce'
  name(1266) = 'ovmtrx:allreduce_1'
  name(1267) = 'ovmtrx:allreduce_2'
  name(1268) = 'cholde:allreduce_1'
  name(1269) = 'cholde:allreduce_2'
  name(1270) = 'cholde:allreduce_3'
  name(1271) = 'cholde:allreduce_4'
  name(1272) = 'cholde:allreduce_5'
  name(1273) = 'cholde:allreduce_6'
  name(1274) = 'cholde:allreduce_7'
  name(1275) = 'cholde:allreduce_8'
  name(1276) = 'prjvlhxc_l2vlhxc_t:allreduce'
  name(1277) = 'amtrxlocalpt:allreduce'
  name(1278) = 'lhl:allreduce_1'
  name(1279) = 'lhl:allreduce_2'
  name(1280) = 'lhl:allreduce_3'
  name(1281) = 'lhl:allreduce_4'
  name(1282) = 'lhl:allreduce_5'
  name(1283) = 'lhl:allreduce_6'
  name(1284) = 'lx:allreduce_1'
  name(1285) = 'lx:allreduce_2'
! -- MATDIAGON / DO
  name(1286) = 'solve_Hx_eq_ex_LAPACK:DO_1'
  name(1287) = 'solve_Hx_eq_ex_LAPACK:DO_2'
  name(1288) = 'solve_Hx_eq_ex_LAPACK:DO_3'
  name(1289) = 'solve_Hx_eq_ex_LAPACK:DO_4'
  name(1290) = 'solve_Hx_eq_ex_LAPACK:DO_5'
  name(1291) = 'calphase:DO_1'
  name(1292) = 'nlmtrx:DO_1'
  name(1293) = 'nlmtrx:DO_2'
  name(1294) = 'nlmtrx:DO_3'
  name(1295) = 'nlmtrx:DO_4'
  name(1296) = 'nlmtrx:DO_5'
  name(1297) = 'nlmtrx:DO_6'
  name(1298) = 'ovmtrx:DO_1'
  name(1299) = 'ovmtrx:DO_2'
  name(1300) = 'ovmtrx:DO_3'
  name(1301) = 'ovmtrx:DO_4'
  name(1302) = 'ovmtrx:DO_5'
  name(1303) = 'ovmtrx:DO_6'
  name(1304) = 'ovmtrx:DO_7'
  name(1305) = 'cholde:DO_1'
  name(1306) = 'cholde:DO_2'
  name(1307) = 'cholde:DO_3'
  name(1308) = 'cholde:DO_4'
  name(1309) = 'cholde:DO_5'
  name(1310) = 'cholde:DO_6'
  name(1311) = 'cholde:DO_7'
  name(1312) = 'cholde:DO_8'
  name(1313) = 'cholde:DO_9'
  name(1314) = 'cholde:DO_10'
  name(1315) = 'cholde:DO_11'
  name(1316) = 'cholde:DO_12'
  name(1317) = 'cholde:DO_13'
  name(1318) = 'cholde:DO_14'
  name(1319) = 'amtrxkinetic:DO_1'
  name(1320) = 'prjvlhxc_l2vlhxc_t:DO_1'
  name(1321) = 'amtrxlocalpt:DO_1'
  name(1322) = 'amtrxlocalpt:DO_2'
  name(1323) = 'amtrxlocalpt:DO_3'
  name(1324) = 'amtrxlocalpt:DO_4'
  name(1325) = 'lhl:DO_1'
  name(1326) = 'lhl:DO_2'
  name(1327) = 'lhl:DO_3'
  name(1328) = 'lhl:DO_3'
  name(1329) = 'lhl:DO_5'
  name(1330) = 'lhl:DO_6'
  name(1331) = 'lhl:DO_7'
  name(1332) = 'lx:DO_1'
  name(1333) = 'lx:DO_2'
  name(1334) = 'lx:DO_3'
  name(1335) = 'lx:DO_4'
  name(1336) = 'lx:DO_5'
  name(1337) = 'phase_mult:DO_1'
  name(1338) = 'phase_mult:DO_2'
  name(1339) = 'phase_mult:DO_3'
  name(1340) = 'cp_zaj_mat_to_zaj_l:DO_1'
  name(1341) = 'cp_zaj_mat_to_zaj_l:DO_2'
  name(1342) = 'cp_zaj_mat_to_zaj_l:DO_3'
  name(1343) = 'cp_zaj_mat_to_zaj_l:DO_4'
  name(1344) = 'cp_zaj_mat_to_zaj_l:DO_5'
! -- Initialization / DO    # TIMER_INIDO
  name(1345) = 'm_ESIW_by_randomnumbers_3D:DO_1'
  name(1346) = 'm_ESIW_by_randomnumbers_3D:DO_2'
  name(1347) = 'get_Gvectors_in_a_gmaxp_sphere:DO_1'
  name(1348) = 'get_Gvectors_in_a_gmaxp_sphere:DO_2'
  name(1349) = 'get_Gvectors_in_a_gmaxp_sphere:DO_3'
  name(1350) = 'get_Gvectors_in_a_gmaxp_sphere:DO_4'
  name(1351) = 'count_Gvectors_in_spheres:DO_1'
! -- Initialization / COMMUNICATION    # TIMER_INICOMM
  name(1352) = 'm_Parallel_wf_onto_fft_3D:isend/irecv_1'
  name(1353) = 'm_Parallel_wf_onto_fft_3D:allgater_1'
  name(1354) = 'm_Parallel_wf_onto_fft_3D:isend/irecv_2'
  name(1355) = 'm_Parallel_wf_onto_fft_3D:alltoall_1'
  name(1356) = 'm_Parallel_fft_onto_wf_3D:isend/irecv_1'
  name(1357) = 'm_Parallel_fft_onto_wf_3D:alltoall_1'
  name(1358) = 'm_Parallel_fft_onto_chgq_3D:isend/irecv_1'
  name(1359) = 'm_Parallel_fft_onto_chgq_3D:alltoall_1'
  name(1360) = 'm_Parallel_chgq_onto_fftcd_3D:isend/irecv_1'
  name(1701) = 'm_Parallel_chgq_onto_fftcd_3D:allgather_1'
  name(1702) = 'm_Parallel_chgq_onto_fftcd_3D:isend/irecv_1'
  name(1703) = 'm_Parallel_chgq_onto_fftcd_3D:alltoall_1'
  name(1704) = 'm_Parallel_fftcd_onto_chgq_3D:isend/irecv_1'
  name(1705) = 'm_Parallel_fftcd_onto_chgq_3D:alltoall_1'

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- Finalization / SUBROUTINE
  name(1361) = 'm_PP_rd_PP_parameters_3D'
  name(1362) = 'rd_kgp_array'
  name(1363) = 'rd_kgp_array_p'
  name(1364) = 'bcast_nfcntn_bin'
  name(1365) = 'm_PP_wd_PP_parameters_3D'
  name(1366) = 'wd_kgp_array'
  name(1367) = 'wd_kgp_array_p'
  name(1368) = 'm_NLP_rd_snl_3D'
  name(1369) = 'm_NLP_wd_snl_3D'
  name(1370) = 'm_ESIO_rd_EigenValues_etc'
  name(1371) = 'm_ESIO_wd_EigenValues_etc_3D'
  name(1372) = 'm_ESIO_rd_WFs'
  name(1373) = 'm_ESIO_wd_WFs_3D'
  name(1374) = 'm_CD_rd_chgq_3D'
  name(1375) = 'm_CD_wd_chgq_3D'
  name(1376) = 'm_ESIO_rd_Efermi'
  name(1377) = 'm_ESIO_wd_efermi'
  name(1378) = 'm_ESIO_wd_EigenValues'
  name(1379) = 'm_ESIO_wd_WFs_standardout_3D'
! -- Finalization / COMMUNICATION & DO     #TIMER_IODO TIMER_IOCOMM
  name(1380) = 'rd_kgp_array:READ_1'
  name(1381) = 'rd_kgp_array:bcast_1'
  name(1382) = 'rd_kgp_array:READ_2'
  name(1383) = 'rd_kgp_array:bcast_2'
  name(1384) = 'rd_kgp_array:DO_1'
  name(1385) = 'rd_kgp_array_p:READ_1'
  name(1386) = 'rd_kgp_array_p:DO_1'
  name(1387) = 'bcast_nfcntn_bin:bcast'
  name(1388) = 'wd_kgp_array:allreduce_1'
  name(1389) = 'wd_kgp_array:WRITE_1'
  name(1390) = 'wd_kgp_array:DO_1'
  name(1391) = 'wd_kgp_array:allreduce_2'
  name(1392) = 'wd_kgp_array:WRITE_2'
  name(1393) = 'wd_kgp_array_p:DO_1'
  name(1394) = 'wd_kgp_array_p:WRITE_1'
  name(1395) = 'm_NLP_rd_snl_3D:DO_1'
  name(1396) = 'm_NLP_rd_snl_3D:READ_1'
  name(1397) = 'm_NLP_rd_snl_3D:bcast_1'
  name(1398) = 'm_NLP_rd_snl_3D:DO_2'
  name(1399) = 'm_NLP_wd_snl_3D:WRITE_1'
  name(1400) = 'm_NLP_wd_snl_3D:DO_1'
  name(1401) = 'm_NLP_wd_snl_3D:allreduce_1'
  name(1402) = 'm_NLP_wd_snl_3D:DO_2'
  name(1403) = 'm_NLP_wd_snl_3D:allreduce_2'
  name(1404) = 'm_NLP_wd_snl_3D:WRITE_2'
  name(1405) = 'm_ESIO_rd_EigenValues_etc:READ_1'
  name(1406) = 'm_ESIO_rd_EigenValues_etc:DO_1'
  name(1407) = 'm_ESIO_rd_EigenValues_etc:READ_2'
  name(1408) = 'm_ESIO_rd_EigenValues_etc:DO_2'
  name(1409) = 'm_ESIO_rd_EigenValues_etc:READ_3'
  name(1410) = 'm_ESIO_rd_EigenValues_etc:bcast_1'
  name(1411) = 'm_ESIO_rd_EigenValues_etc:DO_3'
  name(1412) = 'm_ESIO_rd_EigenValues_etc:READ_4'
  name(1413) = 'm_ESIO_rd_EigenValues_etc:bcast_2'
  name(1414) = 'm_ESIO_rd_EigenValues_etc:DO_4'
  name(1415) = 'm_ESIO_wd_EigenValues_etc:DO_1'
  name(1416) = 'm_ESIO_wd_EigenValues_etc:WRITE_1'
  name(1417) = 'm_ESIO_wd_EigenValues_etc:DO_2'
  name(1418) = 'm_ESIO_wd_EigenValues_etc:WRITE_2'
  name(1419) = 'm_ESIO_wd_EigenValues_etc:DO_3'
  name(1420) = 'm_ESIO_wd_EigenValues_etc:WRITE_3'
  name(1421) = 'm_ESIO_wd_EigenValues_etc:DO_4'
  name(1422) = 'm_ESIO_wd_EigenValues_etc:WRITE_4'
  name(1423) = 'm_ESIO_wd_EigenValues_etc:DO_5'
  name(1424) = 'm_ESIO_wd_EigenValues_etc:allreduce_1'
  name(1425) = 'm_ESIO_wd_EigenValues_etc:WRITE_5'
  name(1426) = 'm_ESIO_wd_EigenValues_etc:DO_6'
  name(1427) = 'm_ESIO_wd_EigenValues_etc:allreduce_2'
  name(1428) = 'm_ESIO_wd_EigenValues_etc:WRITE_6'
  name(1429) = 'm_ESIO_wd_EigenValues_etc:DO_7'
  name(1430) = 'm_ESIO_wd_EigenValues_etc:allreduce_3'
  name(1431) = 'm_ESIO_wd_EigenValues_etc:WRITE_7'
  name(1432) = 'm_ESIO_wd_EigenValues_etc:DO_8'
  name(1433) = 'm_ESIO_wd_EigenValues_etc:allreduce_4'
  name(1434) = 'm_ESIO_wd_EigenValues_etc:WRITE_8'
  name(1435) = 'm_ESIO_rd_WFs:DO_1'
  name(1436) = 'm_ESIO_rd_WFs:READ_1'
  name(1437) = 'm_ESIO_rd_WFs:DO_2'
  name(1438) = 'm_ESIO_rd_WFs:READ_2'
  name(1439) = 'm_ESIO_rd_WFs:bcast'
  name(1440) = 'm_ESIO_wd_WFs:DO_1'
  name(1441) = 'm_ESIO_wd_WFs:WRITE_1'
  name(1442) = 'm_ESIO_wd_WFs:DO_2'
  name(1443) = 'm_ESIO_wd_WFs:allreduce'
  name(1444) = 'm_ESIO_wd_WFs:WRITE_2'
  name(1445) = 'm_CD_rd_chgq_3D:READ_1'
  name(1446) = 'm_CD_rd_chgq_3D:DO_1'
  name(1447) = 'm_CD_rd_chgq_3D:READ_2'
  name(1448) = 'm_CD_rd_chgq_3D:READ_3'
  name(1449) = 'm_CD_rd_chgq_3D:bcast'
  name(1450) = 'm_CD_rd_chgq_3D:DO_2'
  name(1451) = 'm_CD_rd_chgq_3D:READ_4'
  name(1452) = 'm_CD_wd_chgq_3D:DO_1'
  name(1453) = 'm_CD_wd_chgq_3D:WRITE_1'
  name(1454) = 'm_CD_wd_chgq_3D:WRITE_2'
  name(1455) = 'm_CD_wd_chgq_3D:allreduce'
  name(1456) = 'm_CD_wd_chgq_3D:WRITE_3'
  name(1457) = 'm_CD_wd_chgq_3D:WRITE_4'
  name(1458) = 'm_ESIO_rd_Efermi:READ_1'
  name(1459) = 'm_ESIO_rd_Efermi:bcast_1'
  name(1460) = 'm_ESIO_wd_Efermi:WRITE_1'
  name(1461) = 'm_ESIO_wd_EigenValues_put_kpart:DO'
  name(1462) = 'm_ESIO_wd_EigenValues_put_kpart:allreduce'
  name(1463) = 'm_ESIO_wd_EigenValues:WRITE_1'
!!$  name(1464) = 'm_ESIO_wd_EigenValues_3D:DO_2'
!!$  name(1465) = 'm_ESIO_wd_EigenValues_3D:allreduce_2'
!!$  name(1466) = 'm_ESIO_wd_EigenValues_3D:WRITE_2'
  name(1467) = 'm_ESIO_wd_WFs_standardout_3D:DO_1'
  name(1468) = 'm_ESIO_wd_WFs_standardout_3D:allreduce_1'
  name(1469) = 'm_ESIO_wd_WFs_standardout_3D:WRITE_1'

  name(1490) = 'FORCE                                   '

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_1 / SUBROUTINE / ID= 101 -> 110
  name(101) =  'm_ES_Vlocal_in_Rspace_3D                '
  name(102) =  'cnt_vlhxc_to_fft_box_3D                 '
  name(103) =  'mp_vlhxc_l_onto_afft_3D                 '
  name(104) =  'm_FFT_WF_3D                             '
  name(105) =  'm_FFT_Direct_3D                         '
  name(106) =  'm_FFT_Inverse_3D                        '
  name(107) =  'm_FFT_Vlocal_W_3D                       '
  name(108) =  'm_FFT_W_Vlocal_W_3D                     '
  name(109) =  'm_FFT_CD_Direct_3D                      '
  name(110) =  'm_FFT_CD_Inverse_3D                     '

! -- KUKAN_1 / COMMUNICATION & DGEMM & DO / ID= 111 -> 189
  name(111) =  'cnt_vlhxc_to_fft_box_3D:isend/irecv_1   '
  name(112) =  'cnt_vlhxc_to_fft_box_3D:DO_1            '
  name(113) =  'cnt_vlhxc_to_fft_box_3D:isend/irecv_2   '
  name(114) =  'cnt_vlhxc_to_fft_box_3D:DO_2            '
  name(115) =  'mp_vlhxc_l_onto_afft_3D:DO_1            '
  name(116) =  'mp_vlhxc_l_onto_afft_3D:isend/irecv_1   '
  name(117) =  'mp_vlhxc_l_onto_afft_3D:DO_2            '
  name(118) =  'm_FFT_Direct_3D:Initial                 '
  name(119) =  'm_FFT_Direct_3D:Initial_FFTW            '
  name(120) =  'm_FFT_Direct_3D:Allocate                '
  name(121) =  'm_FFT_Direct_3D:DO_4                    '
  name(122) =  'm_FFT_Direct_3D:DO_5                    '
  name(123) =  'm_FFT_Direct_3D:alltoall                '
  name(124) =  'm_FFT_Direct_3D:DO_6                    '
  name(125) =  'm_FFT_Direct_3D:DO_7                    '
  name(126) =  'm_FFT_Direct_3D:DO_8                    '
  name(127) =  'm_FFT_Direct_3D:alltoall                '
  name(128) =  'm_FFT_Direct_3D:DO_9                    '
  name(129) =  'm_FFT_Direct_3D:DO_10                   '
  name(130) =  'm_FFT_Direct_3D:DO_11                   '
  name(131) =  'm_FFT_Direct_3D:DO_12                   '
  name(132) =  'm_FFT_Direct_3D:allgatherv              '
  name(133) =  'm_FFT_Direct_3D:DO_13                   '
  name(134) =  'm_FFT_Inverse_3D:Initial                '
  name(135) =  'm_FFT_Inverse_3D:Initial_FFTW           '
  name(136) =  'm_FFT_Inverse_3D:DO_3                   '
  name(137) =  'm_FFT_Inverse_3D:alltoall               '
  name(138) =  'm_FFT_Inverse_3D:DO_4                   '
  name(139) =  'm_FFT_Inverse_3D:DO_5                   '
  name(140) =  'm_FFT_Inverse_3D:DO_6                   '
  name(141) =  'm_FFT_Inverse_3D:DO_7                   '
  name(142) =  'm_FFT_Inverse_3D:alltoall               '
  name(143) =  'm_FFT_Inverse_3D:DO_7                   '
  name(144) =  'm_FFT_Inverse_3D:DO_8                   '
  name(145) =  'm_FFT_Inverse_3D:DO_9                   '
  name(146) =  'm_FFT_Inverse_3D:DO_10                  '
  name(147) =  'm_FFT_Inverse_3D:DO_11                  '
  name(148) =  'm_FFT_Inverse_3D:DO_12                  '
  name(149) =  'm_FFT_Inverse_3D:allgatherv             '
  name(150) =  'm_FFT_Inverse_3D:DO_13                  '
  name(151) =  'm_FFT_Vlocal_W_3D:DO_1                  '
  name(152) =  'm_FFT_W_Vlocal_W_3D:DO_1                '
  name(153) =  'm_FFT_W_Vlocal_W_3D:DO_2                '
  name(154) =  'm_FFT_W_Vlocal_W_3D:DO_3                '
  name(155) =  'm_FFT_W_Vlocal_W_3D:DO_4                '
  name(156) =  'm_FFT_W_Vlocal_W_3D:DO_5                '
  name(157) =  'm_FFT_CD_Direct_3D:Initial              '
  name(158) =  'm_FFT_CD_Direct_3D:Allocate             '
  name(159) =  'm_FFT_CD_Direct_3D:DO_3                 '
  name(160) =  'm_FFT_CD_Direct_3D:DO_4                 '
  name(161) =  'm_FFT_CD_Direct_3D:DO_5                 '
  name(162) =  'm_FFT_CD_Direct_3D:alltoall             '
  name(163) =  'm_FFT_CD_Direct_3D:DO_6                 '
  name(164) =  'm_FFT_CD_Direct_3D:DO_7                 '
  name(165) =  'm_FFT_CD_Direct_3D:DO_8                 '
  name(166) =  'm_FFT_CD_Direct_3D:alltoall             '
  name(167) =  'm_FFT_CD_Direct_3D:DO_9                 '
  name(168) =  'm_FFT_CD_Direct_3D:DO_10                '
  name(169) =  'm_FFT_CD_Direct_3D:DO_11                '
  name(170) =  'm_FFT_CD_Direct_3D:DO_12                '
  name(171) =  'm_FFT_CD_Direct_3D:allgatherv           '
  name(172) =  'm_FFT_CD_Direct_3D:DO_13                '
  name(173) =  'm_FFT_CD_Inverse_3D:Initial             '
  name(174) =  'm_FFT_CD_Inverse_3D:Allocate            '
  name(175) =  'm_FFT_CD_Inverse_3D:DO_3                '
  name(176) =  'm_FFT_CD_Inverse_3D:alltoall            '
  name(177) =  'm_FFT_CD_Inverse_3D:DO_4                '
  name(178) =  'm_FFT_CD_Inverse_3D:DO_5                '
  name(179) =  'm_FFT_CD_Inverse_3D:DO_6                '
  name(180) =  'm_FFT_CD_Inverse_3D:DO_7                '
  name(181) =  'm_FFT_CD_Inverse_3D:alltoall            '
  name(182) =  'm_FFT_CD_Inverse_3D:DO_8                '
  name(183) =  'm_FFT_CD_Inverse_3D:DO_9                '
  name(184) =  'm_FFT_CD_Inverse_3D:DO_10               '
  name(185) =  'm_FFT_CD_Inverse_3D:DO_11               '
  name(186) =  'm_FFT_CD_Inverse_3D:DO_12               '
  name(187) =  'm_FFT_CD_Inverse_3D:DO_13               '
  name(188) =  'm_FFT_CD_Inverse_3D:allgatherv          '
  name(189) =  'm_FFT_CD_Inverse_3D:DO_14               '
  name(266) =  'cnt_vlhxc_to_fft_box_3D:allgather_1     '
  name(267) =  'cnt_vlhxc_to_fft_box_3D:alltoall_1      '
  name(268) =  'mp_vlhxc_l_onto_afft_3D:alltoallv_1     '
! -- KUKAN_1 / FFT / ID= 190         
  name(190)=   'kukan_1:FFT_INVERSE                     '
! -- KUKAN_1 / DO / ID= 251 -> 189
  name(251) =  'm_FFT_Direct_3D:DO_14                   '
  name(252) =  'm_FFT_Direct_3D:DO_15                   '
  name(253) =  'm_FFT_Direct_3D:DO_16                   '
  name(254) =  'm_FFT_Inverse_3D:DO_14                  '
  name(255) =  'm_FFT_Inverse_3D:DO_15                  '
  name(256) =  'm_FFT_Inverse_3D:DO_16                  '
  name(257) =  'm_FFT_CD_Direct_3D:DO_15                '
  name(258) =  'm_FFT_CD_Direct_3D:DO_16                '
  name(259) =  'm_FFT_CD_Direct_3D:DO_17                '
  name(260) =  'm_FFT_CD_Inverse_3D:DO_15               '
  name(261) =  'm_FFT_CD_Inverse_3D:DO_16               '
  name(262) =  'm_FFT_CD_Inverse_3D:DO_17               '
  name(263) =  'm_FFT_Inverse_3D:FFTW                   '
  name(264) =  'm_FFT_Inverse_3D:Allocate               '

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_2 / SUBROUTINE / ID= 201 -> 208
  name(201) =  'm_ES_Vnonlocal_W_3D                     '
  name(202) =  'pre_m_ES_Vnonlocal_W                    '
  name(203) =  'm_ES_gather_f_3d_to_2d                  '
  name(204) =  'calc_phase_blk                          '
  name(205) =  'Vnonlocal_W_part_with_blk_3D            '
  name(206) =  'add_vnlph_l_with_eko_blk_3D             '
  name(207) =  'Vnonlocal_W_part_without_blk_3D         '
  name(208) =  'add_vnlph_l_without_eko_blk_3D          '
! -- KUKAN_2 / COMMUNICATION & DGEMM & DO / ID= 209 -> 231
  name(209) =  'pre_m_ES_Vnonlocal_W:DO_1               '
  name(210) =  'm_ES_gather_f_3d_to_2d:DO_1             '
  name(211) =  'm_ES_gather_f_3d_to_2d:DO_2             '
  name(212) =  'm_ES_gather_f_3d_to_2d:allgatherv_1     '
  name(213) =  'calc_phase_blk:DO_1                     '
  name(214) =  'Vnonlocal_W_part_with_blk_3D:DO_1       '
  name(215) =  'Vnonlocal_W_part_with_blk_3D:DO_2       '
  name(216) =  'Vnonlocal_W_part_with_blk_3D:DO_3       '
  name(217) =  'add_vnlph_l_with_eko_blk_3D:DGEMM_1     '
  name(218) =  'add_vnlph_l_with_eko_blk_3D:DGEMM_2     '
  name(219) =  'add_vnlph_l_with_eko_blk_3D:DGEMM_3     '
  name(220) =  'add_vnlph_l_with_eko_blk_3D:DGEMM_4     '
  name(221) =  'add_vnlph_l_with_eko_blk_3D:DGEMM_5     '
  name(222) =  'add_vnlph_l_with_eko_blk_3D:DGEMM_6     '
  name(223) =  'add_vnlph_l_with_eko_blk_3D:DO_1        '
  name(224) =  'add_vnlph_l_with_eko_blk_3D:DO_2        '
  name(225) =  'add_vnlph_l_with_eko_blk_3D:DO_3        '
  name(226) =  'Vnonlocal_W_part_without_blk_3D:DO_1    '
  name(227) =  'Vnonlocal_W_part_without_blk_3D:DO_2    '
  name(228) =  'Vnonlocal_W_part_without_blk_3D:DO_3    '
  name(229) =  'add_vnlph_l_without_eko_blk_3D:DGEMM_1  '
  name(230) =  'add_vnlph_l_without_eko_blk_3D:DGEMM_2  '
  name(231) =  'add_vnlph_l_without_eko_blk_3D:DGEMM_3  '

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_3 / SUBROUTINE / ID= 301 -> 309
  name(301) =  'evolve_WFs_in_MSD_direction_3D          '
  name(302) =  'm_ES_decide_precon_factor_3D            '
  name(303) =  'kinetic_energy_3D                       '
  name(304) =  'm_ES_WF_in_Rspace_3D                    '
  name(305) =  'Vnonlocal_Diagonal_part_3D              '
  name(306) =  'Vnonlocal_D_norm_conserve_case          '
  name(307) =  'Vnonlocal_D_vanderbilt_case             '
  name(308) =  'map_fft_to_WF_3D                        '
  name(309) =  'modified_steepest_descent_3D            '
! -- KUKAN_3 / COMMUNICATION & DGEMM & DO / ID= 310 -> 326
  name(310) =  'm_ES_decide_precon_factor_3D:DO_1       '
  name(311) =  'kinetic_energy_3D:DO_1                  '
  name(312) =  'kinetic_energy_3D:allreduce_1           '
  name(313) =  'm_ES_WF_in_Rspace_3D:DO_1               '
  name(314) =  'm_ES_WF_in_Rspace_3D:isend/irecv_1      '
  name(315) =  'm_ES_WF_in_Rspace_3D:DO_2               '
  name(316) =  'Vnonlocal_Diagonal_part_3D:DO_1         '
  name(317) =  'Vnonlocal_D_norm_conserve_case:DO_1     '
  name(318) =  'Vnonlocal_D_norm_conserve_case:DO_2     '
  name(319) =  'Vnonlocal_D_vanderbilt_case:DO_1        '
  name(320) =  'Vnonlocal_D_vanderbilt_case:allreduce_1 '
  name(321) =  'Vnonlocal_D_vanderbilt_case:DO_2        '
  name(322) =  'Vnonlocal_D_vanderbilt_case:DO_3        '
  name(323) =  'map_fft_to_WF_3D:isend/irecv_1          '
  name(324) =  'map_fft_to_WF_3D:DO_1                   '
  name(325) =  'map_fft_to_WF_3D:DO_2                   '
  name(326) =  'modified_steepest_descent_3D:DO_1       '
  name(329) =  'm_ES_WF_in_Rspace_3D:alltoallv_1        '
  name(330) =  'map_fft_to_WF_3D:alltoallv_1            '
! -- KUKAN_8 / FFT / ID= 327 -> 328
  name(327)=  'kukan_3:FFT_INVERSE                     '
  name(328)=  'kukan_3:FFT_DIRECT                      '

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_4 / SUBROUTINE / ID= 401 -> 410
  name(401) = 'm_ES_betar_dot_WFs_4_each_k_3D'
  name(402) = 'cpzaj_l_to_psi_ri'
  name(403) = 'set_msize'
  name(404) = 'pre_lmta_k_blk'
  name(405) = 'post_lmta_k_blk'
  name(406) = 'm_ES_betar_dot_WFs_4_lmta_k_blk_3D'
  name(407) = 'G_dot_R_mult_snl_blk'
  name(408) = 'betar_dot_WFs_core_blk'
  name(409) = 'betar_dot_WFs_core2_blk'
  name(410) = 'G_dot_R_map_blk_3D'
! -- KUKAN_4 / DGEMM / ID= 411 -> 418 
  name(411) = 'betar_dot_WFs_core_blk:DGEMM_1'
  name(412) = 'betar_dot_WFs_core_blk:DGEMM_2'
  name(413) = 'betar_dot_WFs_core_blk:DGEMM_3'
  name(414) = 'betar_dot_WFs_core_blk:DGEMM_4'
  name(415) = 'betar_dot_WFs_core2_blk:DGEMM_1'
  name(416) = 'betar_dot_WFs_core2_blk:DGEMM_2'
  name(417) = 'betar_dot_WFs_core2_blk:DGEMM_3'
  name(418) = 'betar_dot_WFs_core2_blk:DGEMM_4'
! -- KUKAN_4 / COMMUNICATION / ID= 419 -> 422
  name(419) = 'm_ES_betar_dot_WFs_4_each_k_3D:reduce_scatter'
  name(420) = 'm_ES_betar_dot_WFs_4_each_k_3D:DO_p'
  name(421) = 'm_ES_betar_dot_WFs_4_each_k_3D:allgatherv'
  name(422) = 'm_ES_betar_dot_WFs_4_each_k_3D:DO_unp'
! -- KUKAN_4 / DO / ID= 423 -> 426
  name(423) = 'cpzaj_l_to_psi_ri:DO_1'
  name(424) = 'post_lmta_k_blk:DO_1'
  name(425) = 'G_dot_R_mult_snl_blk:DO_1'
  name(426) = 'G_dot_R_mult_snl_blk:DO_2'

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_5 / SUBROUTINE / ID= 501 -> 514
  name(501) = 'm_ES_MGS_4_each_k'
  name(502) = 'm_ES_MGS_4_each_k_3D'
  name(503) = 'mgs_4_each_k_G_3D'
  name(504) = 'WSW_t'
  name(505) = 'cp_bpr2bprtw'
  name(506) = 'normalize_bp_and_psi_t'
  name(507) = 'W1SW2_t_r'
  name(508) = 'modify_bp_and_psi_t_r'
  name(509) = 'W1SW2_t_r_block'
  name(510) = 'modify_bp_and_psi_t_r_block'
  name(511) = 'm_ES_F_transpose_r_3D'
  name(512) = 'm_ES_F_transpose_back_r_3D'
  name(513) = 'm_ES_W_transpose_r_3D'
  name(514) = 'm_ES_W_transpose_back_r_3D'
! -- KUKAN_5 / DGEMM / ID= 515 -> 526
  name(515) = 'W1SW2_t_r_block:DGEMM_1'
  name(516) = 'W1SW2_t_r_block:DGEMM_2'
  name(517) = 'W1SW2_t_r_block:DGEMM_3'
  name(518) = 'W1SW2_t_r_block:DGEMM_4'
  name(519) = 'W1SW2_t_r_block:DGEMM_5'
  name(520) = 'W1SW2_t_r_block:DGEMM_6'
  name(521) = 'modify_bp_and_psi_t_r_block:DGEMM_1'
  name(522) = 'modify_bp_and_psi_t_r_block:DGEMM_2'
  name(523) = 'modify_bp_and_psi_t_r_block:DGEMM_3'
  name(524) = 'modify_bp_and_psi_t_r_block:DGEMM_4'
  name(525) = 'modify_bp_and_psi_t_r_block:DGEMM_5'
  name(526) = 'modify_bp_and_psi_t_r_block:DGEMM_6'
! -- KUKAN_5 / COMMUNICATION / ID= 527 -> 539,582 -> 590
  name(527) = 'mgs_4_each_k_G_3D:DO_p'
  name(528) = 'mgs_4_each_k_G_3D:DO_p'
  name(529) = 'mgs_4_each_k_G_3D:bcast_1'
  name(530) = 'mgs_4_each_k_G_3D:DO_unp'
  name(531) = 'mgs_4_each_k_G_3D:bcast_2'
  name(532) = 'mgs_4_each_k_G_3D:DO_unp_1'
  name(533) = 'mgs_4_each_k_G_3D:DO_unp_2'
  name(534) = 'mgs_4_each_k_G_3D:DO_p'
  name(535) = 'mgs_4_each_k_G_3D:allreduce'
  name(536) = 'WSW_t:allreduce'
  name(537) = 'W1SW2_t_r:DO_p'
  name(538) = 'W1SW2_t_r:allreduce'
  name(539) = 'W1SW2_t_r:DO_unp'
  name(582) = 'm_ES_F_transpose_r_3D:isend/irecv_1'
  name(583) = 'm_ES_F_transpose_r_3D:alltoall_1'
! -- KUKAN_5 / DO / ID= 540 -> 581
  name(540) = 'mgs_4_each_k_G_3D:DO_1'
  name(541) = 'mgs_4_each_k_G_3D:DO_2'
  name(542) = 'mgs_4_each_k_G_3D:DO_3'
  name(543) = 'mgs_4_each_k_G_3D:DO_4'
  name(544) = 'mgs_4_each_k_G_3D:DO_5'
  name(545) = 'mgs_4_each_k_G_3D:DO_6'
  name(546) = 'mgs_4_each_k_G_3D:DO_7'
  name(547) = 'mgs_4_each_k_G_3D:DO_8'
  name(548) = 'WSW_t:DO_1'
  name(549) = 'WSW_t:DO_2'
  name(550) = 'WSW_t:DO_3'
  name(551) = 'WSW_t:DO_4'
  name(552) = 'cp_bpr2bprtw:DO_1'
  name(553) = 'cp_bpr2bprtw:DO_2'
  name(554) = 'normalize_bp_and_psi_t:DO_1'
  name(555) = 'normalize_bp_and_psi_t:DO_2'
  name(556) = 'W1SW2_t_r:DO_1'
  name(557) = 'W1SW2_t_r:DO_2'
  name(558) = 'W1SW2_t_r:DO_3'
  name(559) = 'W1SW2_t_r:DO_4'
  name(560) = 'W1SW2_t_r:DO_5'
  name(561) = 'W1SW2_t_r:DO_6'
  name(562) = 'W1SW2_t_r:DO_7'
  name(563) = 'W1SW2_t_r:DO_8'
  name(564) = 'modify_bp_and_psi_t_r:DO_1'
  name(565) = 'modify_bp_and_psi_t_r:DO_2'
  name(566) = 'modify_bp_and_psi_t_r:DO_3'
  name(567) = 'modify_bp_and_psi_t_r:DO_4'
  name(568) = 'modify_bp_and_psi_t_r:DO_5'
  name(569) = 'modify_bp_and_psi_t_r:DO_6'
  name(570) = 'W1SW2_t_r_block:DO_1'
  name(571) = 'W1SW2_t_r_block:DO_2'
  name(572) = 'W1SW2_t_r_block:DO_3'
  name(573) = 'm_ES_F_transpose_r_3D:DO_1'
  name(574) = 'm_ES_F_transpose_r_3D:DO_2'
  name(575) = 'm_ES_F_transpose_back_r_3D:DO_1'
  name(576) = 'm_ES_F_transpose_back_r_3D:DO_2'
  name(577) = 'm_ES_W_transpose_r_3D:DO_1'
  name(578) = 'm_ES_W_transpose_r_3D:DO_2'
  name(579) = 'm_ES_W_transpose_r_3D:DO_3'
  name(580) = 'm_ES_W_transpose_back_r_3D:DO_1'
  name(581) = 'm_ES_W_transpose_back_r_3D:DO_2'
  name(585) = 'gram_schmidt_recursive'
  name(586) = 'gram_schmidt_angle'
  name(587) = 'gram_schmidt_block'
  name(588) = 'gram_schmidt_block:DO_1'
  name(589) = 'gram_schmidt_block:DO_2'
  name(590) = 'gram_schmidt_block:DO_3'
  name(591) = 'gram_schmidt_block:allreduce'
  name(592) = 'gram_schmidt_block:DO_4'
! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_6 / SUBROUTINE / ID= 601 -> 607
  name(601) =  'm_ES_eigen_values_for_each_k_3D         '
  name(602) =  'W_T_W                                   '
  name(603) =  'W_Vnonlocal_W                           '
  name(604) =  'm_ES_sort_eigen_values_3D               '
  name(605) =  'heap_sorting                            '
  name(606) =  'cp_eigen_values_for_af                  '
  name(607) =  'expand_neordr_and_nrvf_ordr             '
! -- KUKAN_6 / COMMUNICATION & DGEMM & DO / ID= 608 -> 635
  name(608) =  'm_ES_eigen_values_for_each_k_3D:alreduce'
  name(609) =  'W_T_W:DO_1                              '
  name(610) =  'W_Vnonlocal_W:DO_1                      '
  name(611) =  'W_Vnonlocal_W:DO_2                      '
  name(612) =  'm_ES_sort_eigen_values_3D:DO_1          '
  name(613) =  'm_ES_sort_eigen_values_3D:DO_2          '
  name(614) =  'm_ES_sort_eigen_values_3D:allgather_1   '
  name(615) =  'm_ES_sort_eigen_values_3D:DO_3          '
  name(616) =  'm_ES_sort_eigen_values_3D:DO_4          '
  name(617) =  'm_ES_sort_eigen_values_3D:DO_5          '
  name(618) =  'm_ES_sort_eigen_values_3D:DO_6          '
  name(619) =  'm_ES_sort_eigen_values_3D:DO_7          '
  name(620) =  'm_ES_sort_eigen_values_3D:DO_8          '
  name(621) =  'm_ES_sort_eigen_values_3D:DO_9          '
  name(622) =  'm_ES_sort_eigen_values_3D:DO_10         '
  name(623) =  'm_ES_sort_eigen_values_3D:DO_11         '
  name(624) =  'm_ES_sort_eigen_values_3D:DO_12         '
  name(625) =  'm_ES_sort_eigen_values_3D:DO_13         '
  name(626) =  'm_ES_sort_eigen_values_3D:DO_14         '
  name(627) =  'm_ES_sort_eigen_values_3D:DO_15         '
  name(628) =  'm_ES_sort_eigen_values_3D:DO_16         '
  name(629) =  'm_ES_sort_eigen_values_3D:DO_17         '
  name(630) =  'm_ES_sort_eigen_values_3D:DO_18         '
  name(631) =  'm_ES_sort_eigen_values_3D:DO_19         '
  name(632) =  'heap_sorting:DO_1                       '
  name(633) =  'heap_sorting:DO_2                       '
  name(634) =  'cp_eigen_values_for_af:DO_1             '
  name(635) =  'expand_neordr_and_nrvf_ordr:DO_1        '
! -- KUKAN_6 / FFT / ID= 608 -> 635
  name(636) =  'kukan_6:FFT_INVERSE                     '
  name(637) =  'kukan_6:FFT_DIRECT                      '

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_7 / SUBROUTINE / ID= 701 -> 799
  name(701) =  'ChargeDensity_Construction_3D           '
  name(702) =  'FermiEnergyLevel_3D                     '
  name(703) =  'm_ESoc_fermi_parabolic_3D               '
  name(704) =  'get_occup_l_and_tot                     '
  name(705) =  'width2                                  '
  name(706) =  'get_entropic_term                       '
  name(707) =  'get_entropy                             '
  name(708) =  'get_total_spin0                         '
  name(709) =  'm_ESoc_fermi_tetrahedron_3D             '
  name(710) =  'fermi1                                  '
  name(711) =  'nsdos2                                  '
  name(712) =  'fermi2                                  '
  name(713) =  'nstt3i                                  '
  name(714) =  'nstt2i                                  '
  name(715) =  'm_ESoc_fermi_ColdSmearing_3D            '
  name(716) =  'm_CD_softpart_3D                        '
  name(717) =  'add_occupied_densities                  '
  name(718) =  'map_fft_to_chgq_3D                      '
  name(719) =  'substitute_CD_for_chgq                  '
  name(720) =  'm_CD_hardpart_3D                        '
  name(721) =  'summation_of_ff_3D                      '
  name(722) =  'add_hardpart_to_chgq_l_3D               '
  name(723) =  'm_PP_set_index_arrays1                  '
  name(724) =  'm_PP_set_index_arrays2                  '
  name(725) =  'substitute_qitgred                      '
  name(726) =  'substitute_ylmred                       '
  name(727) =  'calc_phase_div                          '
  name(728) =  'sum_hsr_dot_gauntc0                     '
  name(729) =  'add_hardpart_to_chgq_l_div0             '
  name(730) =  'add_hardpart_to_chgq_l_div1             '
  name(731) =  'charge_average_3D                       '
  name(732) =  'cp_chgq_by_ngpt                         '
  name(733) =  'm_CD_conversion_check_3D                '
  name(734) =  'm_CD_hardpart_hsr_3D                    '
  name(735) =  'Renewal_of_OccMat_3D                    '
  name(736) =  'm_CD_den_mat_3D                         '
  name(737) =  'm_OP_occ_mat_ylm                        '
  name(738) =  'm_OP_occ_mat_ao_kt                      '
  name(739) =  'm_TE_total_energy                       '
  name(740) =  'get_band_energy_3D                      '
  name(741) =  'get_xc_and_HE_of_old_CD_3D              '
  name(742) =  'get_local_potential_energy_3D           '
  name(743) =  'get_hubbard_energy_3D                   '
  name(744) =  'get_nonlocal_potential_energy_3D        '
  name(745) =  'get_hartree_energy_3D                   '
  name(746) =  'get_dipole_energy_3D                    '
  name(747) =  'get_kinetic_energy                      '
  name(748) =  'get_entropic_term                       '
  name(749) =  'sumup_all_energies                      '
  name(750) =  'm_XC_cal_potential_3D                   '
  name(751) =  'xc_allocate_3D                          '
  name(752) =  'set_f2or1_3D                            '
  name(753) =  'map_charge_onto_a_fftcd_box             '
  name(754) =  'check_of_negative_CD_3D                 '
  name(755) =  'cp_afft_to_chgrhr_3D                    '
  name(756) =  'check_lmn_even_3D                       '
  name(757) =  'ggaxcp_diff_3D                          '
  name(758) =  'ggaxcp0_3D                              '
  name(759) =  'abs_grad_rho_up_down_total_3D           '
  name(760) =  'g_xyz_chden_l_3D                        '
  name(761) =  'boundary_zero_into_afft_3D              '
  name(762) =  'cp_afft_to_cgrad_rho_3D                 '
  name(763) =  'add_sq_afft_to_grad_rho_3D              '
  name(764) =  'ex_ggapw91_3D                           '
  name(765) =  'cr_ggapw91_3D                           '
  name(766) =  'ex_ggapbe_3D                            '
  name(767) =  'cr_ggapbe_3D                            '
  name(768) =  'xclda_3D                                '
  name(769) =  'ggabek_3D                               '
  name(770) =  'ggaprd_3D                               '
  name(771) =  'ggaxcp_3D                               '
  name(772) =  'dFxc_over_ddgradrho_3D                  '
  name(773) =  'dFxc_dgradrho_dot_gradrho2_3D           '
  name(774) =  'G_xyz_afft_3D                           '
  name(775) =  'add_negative_afft_to_grad_rho_3D        '
  name(776) =  'finally_gga_xc_pot_3D                   '
  name(777) =  'xcpotf_3D                               '
  name(778) =  'xcpotf_wigner_3D                        '
  name(779) =  'xcpotf_pzold_3D                         '
  name(780) =  'xcpotf_xalfa_3D                         '
  name(781) =  'xcpotf_pz_3D                            '
  name(782) =  'xcpotf_vwn_3D                          '
  name(783) =  'xcpotf_mjw_bh_gl_3D                     '
  name(784) =  'cpafft_3D                               '
  name(785) =  'map_fftcd_onto_charge                   '
  name(786) =  'rhos_diff                               '
  name(787) =  'rhopc_diff_3D                           '
  name(788) =  'rhoh_diff_3D                            '
  name(789) =  'even_case_3D                            '
  name(790) =  'odd_case_3D                             '
  name(791) =  'real_case_3D                            '
  name(792) =  'complex_case_3D                         '
  name(793) =  'dgrhodh1_3D                             '
  name(794) =  'stress_correlation_part_3D              '
  name(795) =  'sum_s_gga12                             '
  name(796) =  'dgrhodh2_3D                             '
  name(797) =  'stress_exchange_part_3D                 '
  name(798) =  'map_drhodh1_3D                          '
  name(799) =  'map_drhodh2_3D                          '
! -- KUKAN_7 / COMMUNICATION & DGEMM & DO / ID= 800 -> 900, 1001 -> 1049
  name(800) =  'm_ESoc_fermi_parabolic_3D:DO_1          '
  name(801) =  'm_ESoc_fermi_parabolic_3D:allreduce_1   '
  name(802) =  'm_ESoc_fermi_parabolic_3D:allgather_1   '
  name(803) =  'm_ESoc_fermi_parabolic_3D:DO_2          '
  name(804) =  'get_occup_l_and_tot:DO_1                '
  name(805) =  'get_entropic_term:DO_1                  '
  name(806) =  'm_ESoc_fermi_tetrahedron_3D:DO_1        '
  name(807) =  'm_ESoc_fermi_tetrahedron_3D:allreduce_1 '
  name(808) =  'm_ESoc_fermi_tetrahedron_3D:DO_2        '
  name(809) =  'm_ESoc_fermi_tetrahedron_3D:DO_3        '
  name(810) =  'm_ESoc_fermi_tetrahedron_3D:DO_4        '
  name(811) =  'fermi1:DO_1                             '
  name(812) =  'fermi1:DO_2                             '
  name(813) =  'fermi1:DO_3                             '
  name(814) =  'nsdos2:DO_1                             '
  name(815) =  'fermi2:DO_1                             '
  name(816) =  'nstt3i:DO_1                             '
  name(817) =  'nstt3i:DO_2                             '
  name(818) =  'nstt3i:DO_3                             '
  name(819) =  'nstt2i:DO_1                             '
  name(820) =  'nstt2i:DO_20                            '
  name(821) =  'm_ESoc_fermi_ColdSmearing_3D:DO_1       '
  name(822) =  'm_ESoc_fermi_ColdSmearing_3D:allreduce_1'
  name(823) =  'm_ESoc_fermi_ColdSmearing_3D:DO_2       '
  name(824) =  'm_ESoc_fermi_ColdSmearing_3D:DO_3       '
  name(825) =  'm_CD_softpart_3D:allreduce_1            '
  name(826) =  'add_occupied_densities:DO_1             '
  name(827) =  'map_fft_to_chgq_3D:isend/irecv_1        '
  name(828) =  'map_fft_to_chgq_3D:DO_1                 '
  name(829) =  'map_fft_to_chgq_3D:DO_2                 '
  name(830) =  'substitute_CD_for_chgq:DO_1             '
  name(831) =  'summation_of_ff_3D:DO_1                 '
  name(832) =  'summation_of_ff_3D:allreduce_1          '
  name(833) =  'summation_of_ff_3D:allreduce_2          '
  name(834) =  'add_hardpart_to_chgq_l_3D:DO_1          '
  name(835) =  'add_hardpart_to_chgq_l_3D:allgather_1   '
  name(836) =  'add_hardpart_to_chgq_l_3D:DO_2          '
  name(837) =  'add_hardpart_to_chgq_l_3D:allreduce_1   '
  name(838) =  'm_PP_set_index_arrays1:DO_1             '
  name(839) =  'm_PP_set_index_arrays2:DO_1             '
  name(840) =  'substitute_qitgred:DO_1                 '
  name(841) =  'substitute_ylmred:DO_1                  '
  name(842) =  'substitute_ylmred:DO_2                  '
  name(843) =  'calc_phase_div:DO_1                     '
  name(844) =  'sum_hsr_dot_gauntc0:DO_1                '
  name(845) =  'sum_hsr_dot_gauntc0:DO_2                '
  name(846) =  'add_hardpart_to_chgq_l_div0:DO_1        '
  name(847) =  'add_hardpart_to_chgq_l_div1:DO_1        '
  name(848) =  'charge_average_3D:DO_1                  '
  name(849) =  'charge_average_3D:allgather_1           '
  name(850) =  'charge_average_3D:DO_2                  '
  name(851) =  'charge_average_3D:DO_3                  '
  name(852) =  'charge_average_3D:DO_4                  '
  name(853) =  'cp_chgq_by_ngpt:isend/irecv_1           '
  name(854) =  'cp_chgq_by_ngpt:DO_1                    '
  name(855) =  'cp_chgq_by_ngpt:DO_2                    '
  name(856) =  'cp_chgq_by_ngpt:DO_3                    '
  name(857) =  'm_CD_conversion_check_3D:DO_1           '
  name(858) =  'm_CD_conversion_check_3D:allreduce_1    '
  name(859) =  'm_CD_conversion_check_3D:DO_2           '
  name(860) =  'm_CD_conversion_check_3D:bcast_1        '
  name(861) =  'm_OP_occ_mat_ylm:DO_1                   '
  name(862) =  'm_OP_occ_mat_ao_kt:DO_1                 '
  name(863) =  'get_band_energy_3D:DO_1                 '
  name(864) =  'get_band_energy_3D:allreduce_1          '
  name(865) =  'get_xc_and_HE_of_old_CD_3D:DO_1         '
  name(866) =  'get_xc_and_HE_of_old_CD_3D:DO_2         '
  name(867) =  'get_xc_and_HE_of_old_CD_3D:allreduce_1  '
  name(868) =  'get_local_potential_energy_3D:DO_1      '
  name(869) =  'get_local_potential_energy_3D:allreduce '
  name(870) =  'get_nonlocal_potential_energy_3D:DO_1   '
  name(871) =  'get_nonlocal_potential_energy_3D:allredu'
  name(872) =  'get_hartree_energy_3D:DO_1              '
  name(873) =  'get_hartree_energy_3D:allreduce_1       '
  name(874) =  'get_dipole_energy_3D:DO_1               '
  name(875) =  'get_dipole_energy_3D:allreduce_1        '
  name(876) =  'xc_allocate_3D:DO_1                     '
  name(877) =  'set_f2or1_3D:DO_1                       '
  name(878) =  'set_f2or1_3D:DO_2                       '
  name(879) =  'map_charge_onto_a_fftcd_box:DO_1        '
  name(880) =  'map_charge_onto_a_fftcd_box:DO_2        '
  name(881) =  'map_charge_onto_a_fftcd_box:DO_3        '
  name(882) =  'map_charge_onto_a_fftcd_box:isend/irecv1'
  name(883) =  'map_charge_onto_a_fftcd_box:DO_4        '
  name(884) =  'check_of_negative_CD_3D:DO_1            '
  name(885) =  'check_of_negative_CD_3D:DO_2            '
  name(886) =  'check_of_negative_CD_3D:DO_3            '
  name(887) =  'check_of_negative_CD_3D:allreduce_1     '
  name(888) =  'cp_afft_to_chgrhr_3D:DO_1               '
  name(889) =  'check_lmn_even_3D:DO_1                  '
  name(890) =  'ggaxcp_diff_3D:allreduce_1              '
  name(891) =  'ggaxcp_diff_3D:allreduce_2              '
  name(892) =  'abs_grad_rho_up_down_total_3D:DO_1      '
  name(893) =  'abs_grad_rho_up_down_total_3D:DO_2      '
  name(894) =  'g_xyz_chden_l_3D:DO_1                   '
  name(895) =  'boundary_zero_into_afft_3D:DO_1         '
  name(896) =  'cp_afft_to_cgrad_rho_3D:DO_1            '
  name(897) =  'add_sq_afft_to_grad_rho_3D:DO_1         '
  name(898) =  'ex_ggapw91_3D:DO_1                      '
  name(899) =  'cr_ggapw91_3D:DO_1                      '
  name(900) =  'ex_ggapbe_3D:DO_1                       '
  name(1001)=  'cr_ggapbe_3D:DO_1                       '
  name(1002)=  'xclda_3D:DO_1                           '
  name(1003)=  'ggabek_3D:DO_1                          '
  name(1004)=  'ggaprd_3D:DO_1                          '
  name(1005)=  'ggaxcp_3D:allreduce_1                   '
  name(1006)=  'ggaxcp_3D:allreduce_2                   '
  name(1007)=  'dFxc_dgradrho_dot_gradrho2_3D:DO_1      '
  name(1008)=  'G_xyz_afft_3D:DO_1                      '
  name(1009)=  'add_negative_afft_to_grad_rho_3D:DO_1   '
  name(1010)=  'finally_gga_xc_pot_3D:DO_1              '
  name(1011)=  'xcpotf_3D:allreduce_1                   '
  name(1012)=  'xcpotf_3D:allreduce_2                   '
  name(1013)=  'xcpotf_wigner_3D:DO_1                   '
  name(1014)=  'xcpotf_wigner_3D:DO_2                   '
  name(1015)=  'xcpotf_pzold_3D:DO_1                    '
  name(1016)=  'xcpotf_xalfa_3D:DO_1                    '
  name(1017)=  'xcpotf_pz_3D:DO_1                       '
  name(1018)=  'xcpotf_vwn_3D:DO_1                     '
  name(1019)=  'xcpotf_mjw_bh_gl_3D:DO_1                '
  name(1020)=  'cpafft_3D:DO_1                          '
  name(1021)=  'map_fftcd_onto_charge:isend/irecv_1     '
  name(1022)=  'map_fftcd_onto_charge:DO_1              '
  name(1023)=  'map_fftcd_onto_charge:DO_2              '
  name(1024)=  'map_fftcd_onto_charge:DO_3              '
  name(1025)=  'rhos_diff:DO_1                          '
  name(1026)=  'rhopc_diff_3D:DO_1                      '
  name(1027)=  'even_case_3D:DO_1                       '
  name(1028)=  'odd_case_3D:DO_1                        '
  name(1029)=  'real_case_3D:DO_1                       '
  name(1030)=  'real_case_3D:DO_2                       '
  name(1031)=  'complex_case_3D:DO_1                    '
  name(1032)=  'dgrhodh1_3D:DO_1                        '
  name(1033)=  'stress_correlation_part_3D:DO_1         '
  name(1034)=  'sum_s_gga12:allreduce_1                 '
  name(1035)=  'sum_s_gga12:allreduce_2                 '
  name(1036)=  'sum_s_gga12:allreduce_3                 '
  name(1037)=  'sum_s_gga12:allreduce_4                 '
  name(1038)=  'sum_s_gga12:send/recv_1                 '
  name(1039)=  'sum_s_gga12:send/recv_2                 '
  name(1040)=  'dgrhodh2_3D:DO_1                        '
  name(1041)=  'dgrhodh2_3D:DO_2                        '
  name(1042)=  'stress_exchange_part_3D:DO_1            '
  name(1043)=  'stress_exchange_part_3D:DO_2            '
  name(1044)=  'map_drhodh1_3D:DO_1                     '
  name(1045)=  'map_drhodh1_3D:allreduce_1              '
  name(1046)=  'map_drhodh1_3D:allreduce_2              '
  name(1047)=  'map_drhodh2_3D:DO_1                     '
  name(1048)=  'map_drhodh2_3D:allreduce_1              '
  name(1049)=  'map_drhodh2_3D:allreduce_2              '
  name(693)=  'map_fft_to_chgq_3D:alltoallv_1           '
  name(694)=  'cp_chgq_by_ngpt:alltoall_1(List)         '
  name(695)=  'map_fftcd_onto_charge:alltoallv_1        '
  name(696)=  'map_charge_onto_a_fftcd_box:alltoallv_1  '
! -- KUKAN_7 / FFT / ID= 1050 -> 1051
  name(691)=  'kukan_7:FFT_INVERSE                     '
  name(692)=  'kukan_7:FFT_DIRECT                      '

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_8&9 / SUBROUTINE / ID= 901 -> 913
  name(901) =  'evolve_WFs_in_subspace_3D               '
  name(902) =  'set_col_partition_3D                    '
  name(903) =  'scalapack_setup_3D                      '
  name(904) =  'set_hmat_3D                             '
  name(905) =  'gather_zmat_all                         '
  name(906) =  'pdsyev_driver_3D                        '
  name(907) =  'pzheev_driver_3D                        '
  name(908) =  'dsyev_driver                            '
  name(909) =  'zheev_driver                            '
  name(910) =  'hobsvw_driver                           '
  name(911) =  'chobsd_driver                           '
  name(912) =  'subspace_rotation_real_3D               '
  name(913) =  'subspace_rotation_imag_3D               '
!fj add 20110819
  name(963) =  'set_hmat_3D_scl                         '
  name(964) =  'trans_scalapack                         '
  name(965) =  'trans_scalapack_r                       '
!fj add 20110819
!fj add 20110901
  name(988) =  'eig_send                                '
!fj add 20110901
! -- KUKAN_8&9 / COMMUNICATION & DGEMM & DO / ID= 914 -> 959
  name(914) =  'evolve_WFs_in_subspace_3D:allreduce_1   '
  name(915) =  'evolve_WFs_in_subspace_3D:isend/irecv_1 '
  name(916) =  'evolve_WFs_in_subspace_3D:DO_1          '
  name(917) =  'evolve_WFs_in_subspace_3D:DO_2          '
  name(918) =  'evolve_WFs_in_subspace_3D:DO_3          '
  name(919) =  'evolve_WFs_in_subspace_3D:DO_4          '
  name(920) =  'evolve_WFs_in_subspace_3D:DO_5          '
  name(921) =  'evolve_WFs_in_subspace_3D:DGEMM_1       '
  name(922) =  'evolve_WFs_in_subspace_3D:DGEMM_2       '
  name(923) =  'evolve_WFs_in_subspace_3D:DO_6          '
  name(924) =  'evolve_WFs_in_subspace_3D:DGEMM_3       '
  name(925) =  'evolve_WFs_in_subspace_3D:DO_7          '
  name(926) =  'evolve_WFs_in_subspace_3D:DO_8          '
  name(927) =  'evolve_WFs_in_subspace_3D:DO_9          '
  name(928) =  'evolve_WFs_in_subspace_3D:DO_10         '
  name(929) =  'evolve_WFs_in_subspace_3D:DO_11         '
  name(930) =  'evolve_WFs_in_subspace_3D:DO_12         '
  name(931) =  'evolve_WFs_in_subspace_3D:DO_13         '
  name(932) =  'evolve_WFs_in_subspace_3D:DO_14         '
  name(933) =  'evolve_WFs_in_subspace_3D:allreduce_2   '
  name(934) =  'evolve_WFs_in_subspace_3D:DO_15         '
  name(935) =  'evolve_WFs_in_subspace_3D:DO_16         '
  name(936) =  'evolve_WFs_in_subspace_3D:DO_17         '
  name(937) =  'evolve_WFs_in_subspace_3D:DO_18         '
  name(938) =  'evolve_WFs_in_subspace_3D:DO_19         '
  name(939) =  'evolve_WFs_in_subspace_3D:DO_20         '
  name(940) =  'set_col_partition_3D:DO_1               '
  name(941) =  'set_col_partition_3D:DO_2               '
  name(942) =  'make_usermap                            '
  name(943) =  'set_hmat_3D:DO_1                        '
  name(944) =  'set_hmat_3D:reduce_scatter_1            '
  name(945) =  'set_hmat_3D:DO_2                        '
  name(946) =  'set_hmat_3D:reduce_scatter_2            '
  name(947) =  'set_hmat_3D:DO_3                        '
  name(948) =  'gather_zmat_all:DO_1                    '
  name(949) =  'gather_zmat_all:allgather               '
  name(950) =  'gather_zmat_all:DO_2                    '
  name(951) =  'gather_zmat_all:DO_3                    '
  name(952) =  'gather_zmat_all:DO_4                    '
  name(953) =  'subspace_rotation_real_3D:DGEMM_1       '
  name(954) =  'subspace_rotation_real_3D:DO_1          '
  name(955) =  'subspace_rotation_imag_3D:DGEMM_1       '
  name(956) =  'subspace_rotation_imag_3D:DO_2          '
  name(957) =  'subspace_rotation_imag_3D:DO_3          '
  name(958) =  'subspace_rotation_imag_3D:DGEMM_2       '
  name(959) =  'subspace_rotation_imag_3D:DO_4          '
  name(960) =  'evolve_WFs_in_subspace_3D:DO_21         '
! -- KUKAN_8 / FFT / ID= 961 -> 962
  name(961)=  'kukan_8:FFT_INVERSE                     '
  name(962)=  'kukan_8:FFT_DIRECT                      '

!fj add 20110819
  name(966) =  'set_hmat_3D_scl:DO_1                    '
  name(967) =  'set_hmat_3D_scl:reduce_scatter_1        '
  name(968) =  'set_hmat_3D_scl:DO_2                    '
  name(969) =  'set_hmat_3D_scl:reduce_scatter_2        '
  name(970) =  'set_hmat_3D_scl:DO_3                    '
  name(971) =  'trans_scalapack:DO_1                    '
  name(972) =  'trans_scalapack:DO_2                    '
  name(973) =  'trans_scalapack:DO_3                    '
  name(974) =  'trans_scalapack:DO_4                    '
  name(975) =  'trans_scalapack:isend/irecv_1           '
  name(976) =  'trans_scalapack_r:DO_1                  '
  name(977) =  'trans_scalapack_r:DO_2                  '
  name(978) =  'trans_scalapack_r:DO_3                  '
  name(979) =  'trans_scalapack_r:DO_4                  '
  name(980) =  'trans_scalapack_r:isend/irecv_1         '
  name(981) =  'evolve_WFs_in_subspace_3D:DO_21         '
  name(982) =  'evolve_WFs_in_subspace_3D:allgather     '
  name(983) =  'evolve_WFs_in_subspace_3D:DO_22         '
  name(984) =  'evolve_WFs_in_subspace_3D:DO_23         '
  name(985) =  'evolve_WFs_in_subspace_3D:DO_24         '
  name(986) =  'evolve_WFs_in_subspace_3D:DO_25         '
  name(987) =  'evolve_WFs_in_subspace_3D:allgather     '
  name(992) =  'trans_scalapack:alltoall_1(List)        '
  name(993) =  'trans_scalapack_r:alltoall_1(List)      '
  name(994) =  'evolve_WFs_in_subspace_3D:alltoallv_1   '
!fj add 20110820
!fj add 20110901
  name(989) =  'eigsend:send                            '
  name(990) =  'eigsend:recv                            '
  name(991) =  'dsyev_driver  dsyev || dsyevd           '
!fj add 20110901

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_11 / SUBROUTINE / ID= 1101 -> 1145
  name(1101) = 'ChargeDensity_Mixing_3D                 '
  name(1102) = 'm_CD_prepare_precon_3D                  '
  name(1103) = 'm_CD_simple_mixing_3D                   '
  name(1104) = 'alloc_chgqstore_recompose_chgq_3D       '
  name(1105) = 'precon_4_charge_mix_3D                  '
  name(1106) = 'compose_chgq_dealloc_chgqstore_3D       '
  name(1107) = 'm_CD_mix_broyden1_3D                    '
  name(1108) = 'simple_mix1_3D                          '
  name(1109) = 'mix_broyden_alloc3_3D                   '
  name(1110) = 'dF_F_d0_u_v_and_dd_3D                   '
  name(1111) = 'set_ncrspd_mxiter_etc                   '
  name(1112) = 'rotate_cmix_arrays                      '
  name(1113) = 'renew_u_br_3D                           '
  name(1114) = 'mult1s_3D                               '
  name(1115) = 'mult1s5_3D                              '
  name(1116) = 'subtr_j_th_term_3D                      '
  name(1117) = 'renew_v_3D                              '
  name(1118) = 'renew_d_br_3D                           '
  name(1119) = 'store_to_urec2_3D                       '
  name(1120) = 'renew_d_last_br_3D                      '
  name(1121) = 'scatter_chg_onto_d_3D                   '
  name(1122) = 'concentrate_d_to_chg_3D                 '
  name(1123) = 'simple_mix_large_Gc_3D                  '
  name(1124) = 'm_CD_mix_broyden2_3D                    '
  name(1125) = 'dF_F_d0_u_and_v_3D                      '
  name(1126) = 'm_CD_mix_DFP_3D                         '
  name(1127) = 'dF_F_d0_u_and_w_3D                      '
  name(1128) = 'renew_w_3D                              '
  name(1129) = 'renew_d_3D                              '
  name(1130) = 'renew_d_last_3D                         '
  name(1131) = 'm_CD_mix_pulay_3D                       '
  name(1132) = 'Resid_and_dd_into_urec                  '
  name(1133) = 'Ri_dot_Rj_3D                            '
  name(1134) = 'get_finv_3D                             '
  name(1135) = 'rdecomp                                 '
  name(1136) = 'rsolve                                  '
  name(1137) = 'mult1s10_3D                             '
  name(1138) = 'Rj_dot_d_3D                             '
  name(1139) = 'get_gmatrix_3D                          '
  name(1140) = 'renew_d_using_g_3D                      '
  name(1141) = 'm_OP_simple_mixing                      '
  name(1142) = 'm_OP_mix_broyden1                       '
  name(1143) = 'm_OP_mix_broyden2                       '
  name(1144) = 'm_CD_simple_mixing_hsr_3D               '
  name(1145) = 'm_CD_check_3D                           '
  name(1050) = 'm_ESlhxc_potential_3D                   '
  name(1051) = 'm_Dipole_calc_3D                        '
  name(1052) = 'electron_part_3D                        '
  name(1053) = 'ion_part                                '
  name(1054) = 'm_Dipole_potential_3D                   '
  name(1055) = 'edip_ion_part                           '
  name(1056) = 'm_ESlhxc_wd_vlhxc_3D                    '
  name(1057) = 'm_ESiVQ_integrate_VlhxcQlm_3D           '
  name(1058) = 'm_PP_find_maximum_l                     '
  name(1059) = 'm_pwBS_sphrp2_3D                        '
  name(1060) = 'substitute_vlhxcred                     '
  name(1061) = 'dp_Vlhxcq_exp_Q_div                     '
  name(1062) = 'veQ_dot_ylm_div                         '
  name(1063) = 'symmetrize_vlhxcq                       '
! -- KUKAN_11 / COMMUNICATION & DGEMM & DO / ID= 1146 -> 1197
  name(1146) = 'm_CD_prepare_precon_3D:allreduce_1      '
  name(1147) = 'm_CD_prepare_precon_3D:allreduce_2      '
  name(1148) = 'm_CD_simple_mixing_3D:DO_1              '
  name(1149) = 'alloc_chgqstore_recompose_chgq_3D       '
  name(1150) = 'precon_4_charge_mix_3D:DO_1             '
  name(1151) = 'precon_4_charge_mix_3D:DO_2             '
  name(1152) = 'precon_4_charge_mix_3D:DO_3             '
  name(1153) = 'compose_chgq_dealloc_chgqstore_3D:DO_1  '
  name(1154) = 'simple_mix1_3D:DO_1                     '
  name(1155) = 'simple_mix1_3D:DO_2                     '
  name(1156) = 'dF_F_d0_u_v_and_dd_3D:DO_1              '
  name(1157) = 'dF_F_d0_u_v_and_dd_3D:DO_2              '
  name(1158) = 'rotate_cmix_arrays:DO_1                 '
  name(1159) = 'rotate_cmix_arrays:DO_2                 '
  name(1160) = 'mult1s_3D:DO_1                          '
  name(1161) = 'mult1s_3D:allreduce_1                   '
  name(1162) = 'mult1s5_3D:DO_1                         '
  name(1163) = 'mult1s5_3D:allreduce_1                  '
  name(1164) = 'subtr_j_th_term_3D:DO_1                 '
  name(1165) = 'store_to_urec2_3D:DO_1                  '
  name(1166) = 'renew_d_last_br_3D:DO_1                 '
  name(1167) = 'renew_d_last_br_3D:DO_2                 '
  name(1168) = 'scatter_chg_onto_d_3D:DO_1              '
  name(1169) = 'scatter_chg_onto_d_3D:bcast_1           '
  name(1170) = 'scatter_chg_onto_d_3D:DO_2              '
  name(1171) = 'concentrate_d_to_chg_3D:DO_1            '
  name(1172) = 'concentrate_d_to_chg_3D:send/recv_1     '
  name(1173) = 'concentrate_d_to_chg_3D:DO_2            '
  name(1174) = 'simple_mix_large_Gc_3D:DO_1             '
  name(1175) = 'dF_F_d0_u_and_v_3D:DO_1                 '
  name(1176) = 'dF_F_d0_u_and_v_3D:DO_2                 '
  name(1177) = 'dF_F_d0_u_and_w_3D:DO_1                 '
  name(1178) = 'renew_d_last_3D:DO_1                    '
  name(1179) = 'Ri_dot_Rj_3D:DO_1                       '
  name(1180) = 'get_finv_3D:DO_1                        '
  name(1181) = 'get_finv_3D:DO_2                        '
  name(1182) = 'rdecomp:DO_1                            '
  name(1183) = 'rdecomp:DO_2                            '
  name(1184) = 'rsolve:DO_1                             '
  name(1185) = 'rsolve:DO_2                             '
  name(1186) = 'mult1s10_3D:DO_1                        '
  name(1187) = 'mult1s10_3D:allreduce_1                 '
  name(1188) = 'get_gmatrix_3D:DO_1                     '
  name(1189) = 'renew_d_using_g_3D:DO_1                 '
  name(1190) = 'renew_d_using_g_3D:DO_2                 '
  name(1191) = 'm_OP_mix_broyden1:DO_1                  '
  name(1192) = 'm_CD_simple_mixing_hsr_3D:DO_1          '
  name(1193) = 'm_CD_check_3D:DO_1                      '
  name(1194) = 'm_CD_check_3D:DO_2                      '
  name(1195) = 'm_CD_check_3D:DO_3                      '
  name(1196) = 'm_CD_check_3D:bcast_1                   '
  name(1197) = 'm_CD_check_3D:bcast_2                   '
  name(1064) = 'm_ESlhxc_potential_3D:DO_1              '
  name(1065) = 'm_ESlhxc_potential_3D:DO_2              '
  name(1066) = 'm_ESlhxc_potential_3D:DO_3              '
  name(1067) = 'electron_part_3D:DO_1                   '
  name(1068) = 'electron_part_3D:bcast_1                '
  name(1069) = 'electron_part_3D:allreduce_1            '
  name(1070) = 'ion_part:DO_1                           '
  name(1071) = 'ion_part:DO_2                           '
  name(1072) = 'm_Dipole_potential_3D:DO_1              '
  name(1073) = 'm_Dipole_potential_3D:DO_2              '
  name(1074) = 'edip_ion_part:DO_1                      '
  name(1075) = 'edip_ion_part:DO_2                      '
  name(1076) = 'm_ESiVQ_integrate_VlhxcQlm_3D:DO_1      '
  name(1077) = 'm_ESiVQ_integrate_VlhxcQlm_3D:allreduce1'
  name(1078) = 'm_ESiVQ_integrate_VlhxcQlm_3D:DO_2      '
  name(1079) = 'm_ESiVQ_integrate_VlhxcQlm_3D:allreduce2'
  name(1080) = 'm_PP_find_maximum_l:DO_1                '
  name(1081) = 'm_pwBS_sphrp2_3D:DO_1                   '
  name(1082) = 'm_pwBS_sphrp2_3D:DO_2                   '
  name(1083) = 'substitute_vlhxcred:DO_1                '
  name(1084) = 'dp_Vlhxcq_exp_Q_div:DO_1                '
  name(1085) = 'veQ_dot_ylm_div:DO_1                    '
  name(1086) = 'symmetrize_vlhxcq:DO_1                  '
  name(1087) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP0      '
  name(1088) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP1      '
  name(1089) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP2      '
  name(1090) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP3      '
  name(1091) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP4      '
  name(1092) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP5      '
  name(1093) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP6      '
  name(1094) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP7      '
  name(1095) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP8      '
  name(1096) = 'm_ESiVQ_integrate_VlhxcQlm_3D:OMP9      '
! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- KUKAN_FORCE
  name(1501) =  'm_Force_initialize                      '
  name(1502) =  'm_Force_term_drv_of_flmt_3D             '
  name(1503) =  'm_Force_term_drv_of_VlhxcQ_3D           '
  name(1504) =  'm_Force_term_Elocal_and_Epc_3D          '
  name(1505) =  'm_Force_term_dipole                     '
  name(1506) =  'm_Force_sumup_and_symmetrize            '
  name(1507) =  'm_Force_cal_forcmx                      '
  name(1508) =  'm_Force:comm1                           '
  name(1509) =  'm_Force:comm2                           '
  name(1510) =  'm_Force:comm3                           '
  name(1511) =  'm_Force:comm4                           '
  name(1512) =  'm_Force:comm5                           '
  name(1513) =  'm_Force:comm6                           '
  name(1514) =  'm_Force:comm7                           '
  name(1515) =  'm_Force:comm8                           '
  name(1516) =  'm_Force:comm9                           '
  name(1517) =  'gather_f_3d_to_2d                       '
  name(1518) =  'term_related_to_drv_of_flmt_NC          '
  name(1519) =  'term_related_to_drv_of_flmt_VDB         '
  name(1520) =  'substitute_ngabcred                     '
  name(1521) =  'calc_phase_div                          '
  name(1522) =  'pre_drv_of_betar_dot_WFs_div            '
  name(1523) =  'drv_of_betar_dot_WFs_div                '
  name(1524) =  'calc_phase_div                          '
  name(1525) =  'sum_hsr_dot_gauntc                      '
  name(1526) =  'hsr_dot_gaunt                           '
  name(1527) =  'add_hardpart_to_flhxcq_l_core           '
  name(1528) =  'add_hardpart_to_flhxcq_l_div            '
  name(1529) =  'force_from_Epc_for_each_atom            '
  name(1530) =  'force_from_Elocal_4_each_atom           '
  name(1531) =  'm_Force_term_Elocal_and_Epc_3D:DO_1     '
  name(1532) =  'm_Force_term_Elocal_and_Epc_3D:DO_2     '
  name(1533) =  'force_from_Epc_for_each_atom:DO_1       '
  name(1534) =  'force_from_Epc_for_each_atom:DO_2       '
  name(1535) =  'force_from_Elocal_4_each_atom:DO_1      '
  name(1536) =  'm_Force_term_drv_of_vlhxcQ_3D:DO_1      '
  name(1537) =  'm_Force_term_drv_of_vlhxcQ_3D:DO_2      '
  name(1538) =  'calc_phase_div:DO_1                     '
  name(1539) =  'sum_hsr_dot_gauntc:DO_1                 '
  name(1540) =  'add_hardpart_to_flhxcq_l_core:DO_1      '
  name(1541) =  'add_hardpart_to_flhxcq_l_core:DO_2      '
  name(1542) =  'add_hardpart_to_flhxcq_l_div:DO_1       '
  name(1543) =  'add_hardpart_to_flhxcq_l_div:DO_2       '
  name(1544) =  'm_Force_term_drv_of_flmt_3D:DO_1        '
  name(1545) =  'gather_f_3d_to_2d:DO_1                  '
  name(1546) =  'term_related_to_drv_of_flmt_NC:DO_1     '
  name(1547) =  'term_related_to_drv_of_flmt_VDB:DO_1    '
  name(1548) =  'substitute_ngabcred:DO_1                '
  name(1549) =  'calc_phase_div:DO_1                     '
  name(1550) =  'pre_drv_of_betar_dot_WFs_div:DO_1       '
  name(1551) =  'pre_drv_of_betar_dot_WFs_div:DO_2       '
  name(1552) =  'drv_of_betar_dot_WFs_div:DO_1           '
  name(1553) =  'drv_of_betar_dot_WFs_div:DO_2           '

! name(000) = '---------1---------2---------3---------4---------5---------6----'
! -- RMM SOLVERE
  name(1601) =  'RMM_Initialize                          '
  name(1602) =  'RMM_zajold2zaj_phi2zaj_old_3D           '
  name(1603) =  'RMM_Vnonlocal_W_RMM_3D                  '
  name(1604) =  'RMM_temporary_FFT                       '
  name(1605) =  'RMM_rmm1_3D                             '
  name(1606) =  'RMM_rmm_n_uda_3D                        '
  name(1607) =  'RMM_zajold2zaj_phi2zaj_old_all_3D       '
  name(1608) =  'RMM_Vnonlocal_W_RMMn_3D                 '
  name(1609) =  'RMM_temporary_FFT                       '
  name(1610) =  'RMM_rmm1_3D                             '
  name(1611) =  'RMM_rmm_n_uda_3D                        '
  name(1612) =  'RMM_Finalize                            '
  name(1613) =  'RMM2_calc_phasek_b_3D                   '
  name(1614) =  'RMM2_sumset_rmm_all4                    '
  name(1615) =  'RMM2_Vnonlocal_W_part_sum_ovr_lmt4      '
  name(1616) =  'RMM2_add_vnlph_l_part4                  '
  name(1617) =  'RMM2_sumset_rmm_all3                    '
  name(1618) =  'RMM2_Vnonlocal_W_part_sum_over_lmt1b    '
  name(1619) =  'RMM2_add_vnlph_l_without_eko_part3      '
  name(1620) =  'RMM2_add_vnlph_l_with_eko_part3         '
  name(1621) =  'RMM2_rr_avr_final_3D:allreduce          '
  name(1622) =  'RMM2_rmm1_3D:allreduce                  '
  name(1623) =  'RMM2_sumset_rmm:allreduce               '
  name(1624) =  'RMM2_sumset_rmm_all4:allreduce          '
  name(1625) =  'RMM2_sumset_rmm_all3:allreduce          '


  return
end subroutine timer_init

!-----------------------------------------------------------------------
subroutine timer_fin
  use m_IterationNumbers,   only  : iteration
  use m_Parallelization ,   only  : nrank_e,nrank_g,nrank_k
  use m_timer
  implicit none

#ifdef USE_WTIME
  include 'mpif.h'
#endif

  integer :: ierr
  integer :: i, it, m, n_err, n_err_g
  real(4) :: rdum, rarray(2)
!  real(4), external :: etime
  real(4) :: etime
  logical  :: infoflg

#ifdef USE_KTRACE
      integer ktrace_target(n_tim)
      common /c_ktrace_target/ ktrace_target

!$omp parallel
      call ktrace_exit
!$omp end parallel
#endif

  ! user and system cpu time
  rdum = etime(rarray)
#ifdef USE_WTIME
  ttot(1) = MPI_WTIME() - ttot(1)  ! real
#else
  call getclkreg(t_now)
  ttot(1) = t_now       - ttot(1)  ! real
#endif
  ttot(2) = rarray(1)   - ttot(2)  ! user
  ttot(3) = rarray(2)   - ttot(3)  ! sys

  call getenv('FJ_TIMER_DIR', dname)
  if (len_trim(dname) .eq. 0) then
    dname = '.'
  else
    call system( &
 &    'if [ ! -d '//dname//' ]; then mkdir -p '//dname//'; fi')
  endif
  call getenv('FJ_TIMER_RNK', crank)
  if (crank .ne. "ALL") then
     if (len_trim(crank) .eq. 0) then
        orank = 0
     else
        read(crank,*) orank
     endif
  endif
  if ((myrank .eq. orank) .or. (crank .eq. "ALL")) then

    ncount_all = 0
    timers_all = 0.0d0

    do it = 0, iteration
!
      if (it > n_ite) cycle
!
      if (crank .eq. "ALL") then
        write(fname,'(A,"/timer_ite_",I3.3,"_",I6.6,".out")') &
     &        dname(:len_trim(dname)), it, myrank
      else
        write(fname,'(A,"/timer_ite_",I3.3,".out")') &
     &        dname(:len_trim(dname)), it
      endif
      if (iunit1 .ne. 6) open(iunit1, file=fname)


      write(IUNIT1,'("cfj_info # nodes = ",i5,2x,"ne,ng,ng = ",i5,x,i5,x,i5)') &
     &                nrank_e*nrank_g*nrank_k,nrank_e,nrank_g,nrank_k
#ifdef USE_WTIME
      write(IUNIT1,'("cfj # rank = ",i6,2x,"real,user,sys = ",3f13.3)') &
     &                myrank, ttot
#else
      write(IUNIT1,'("cfj # rank = ",i6,2x,"real,user,sys = ",3f13.3)') &
     &                myrank, ttot*5e-10
#endif
      write(IUNIT1,'("cfj # rank = ", i6,2x," Iteration = ", i6)') &
     &                myrank, it
      do i = 1, N_TIM
        if ((name(i) .ne. '-').or.(ncount(i,it).ne.0)) then
#ifdef USE_WTIME
          write(IUNIT1,'("cfj ",i4,1x,a,f13.3,1x,i8)') &
         &  i, name(i), timers(i,it),       ncount(i,it)
#else
          write(IUNIT1,'("cfj ",i4,1x,a,f13.3,1x,i8)') &
         &  i, name(i), timers(i,it)*5e-10, ncount(i,it)
#endif
        endif
        timers_all(i) = timers_all(i) + timers(i,it)
        ncount_all(i) = ncount_all(i) + ncount(i,it)
      enddo
      write(IUNIT1,'()')
    enddo

    if (crank .eq. "ALL") then
      write(fname,'(A,"/timer_tot_",I6.6,".out")') &
   &        dname(:len_trim(dname)), myrank
    else
      write(fname,'(A,"/timer_tot.out")') &
   &        dname(:len_trim(dname))
    endif
    if (iunit2 .ne. 6) open(iunit2, file=fname)
    write(IUNIT2,'("cfj_info # nodes = ",i5,2x,"ne,ng,ng = ",i5,x,i5,x,i5)') &
   &                  nrank_e*nrank_g*nrank_k,nrank_e,nrank_g,nrank_k 
#ifdef USE_WTIME
    write(IUNIT2,'("cfj # rank = ",i6,2x,"real,user,sys = ",3f13.3)') &
   &                  myrank, ttot
#else
    write(IUNIT2,'("cfj # rank = ",i6,2x,"real,user,sys = ",3f13.3)') &
   &                  myrank, ttot*5e-10
#endif
    write(IUNIT2,'("cfj # rank = ", i6,2x," Total Iteration = ", i3)') &
   &                  myrank, iteration
    do i = 1, N_TIM
      if ((name(i) .ne. '-').or.(ncount_all(i).ne.0)) then
!      if (ncount_all(i).ne.0) then
#ifdef USE_WTIME
         write(IUNIT2,'("cfj ",i4,1x,a,f13.3,1x,i8)') &
        &  i, name(i), timers_all(i),       ncount_all(i)
#else
         write(IUNIT2,'("cfj ",i4,1x,a,f13.3,1x,i8)') &
        &  i, name(i), timers_all(i)*5e-10, ncount_all(i)
#endif
      endif
    enddo
    write(IUNIT2,'()')

#ifdef USE_MEM
    call flush(unitmw)
    close(unitmw)
#endif
    call flush(IUNIT1)
    call flush(IUNIT2)
    if (IUNIT1 .ne. 6) close(IUNIT1)
    if (IUNIT2 .ne. 6) close(IUNIT2)

  endif

  return
end subroutine timer_fin
!-----------------------------------------------------------------------

subroutine timer_sta(id)
  use m_IterationNumbers,   only  : iteration
  use m_timer
  implicit none
#ifdef USE_WTIME
  include 'mpif.h'
#endif
  integer :: id,idx
  character*5 cid
  character*3 citer
  character*11 ccc

! if (.not.tflag) return
!$omp master
  if(id < 0) then
!   idx = id * (-1) + 90000
    idx = id * (-1)
  else
    idx = id
  endif

  write(cid,'(I5.5)') idx
  write(citer,'(I3.3)') iteration
  ccc = "r_"//citer//"_"//cid

#ifdef USE_MEM
  if (myrank .eq. 0) then
    integer :: j
    character(1000) :: line
    character(15),dimension(24) :: str

    open(unitmr,file='/proc/self/stat', form="formatted", access="stream", status='old')
    read(unitmr,'(a)') line
    close(unitmr)
    read(line,*) (str(j),j=1,24)
#ifdef USE_WTIME
    write(unitmw,*) ccc, " sta", MPI_WTIME()-t_sta, str(23), str(24)
#else
    call getclkreg(t_now)
    write(unitmw,*) ccc, " sta", (t_now-t_sta)*5d-10, str(23), str(24)
#endif
  endif
#endif

#ifdef USE_TOFUPA
  call fj_tofupa_start(ccc,id,1)
#endif

#ifdef USE_KTRACE
  integer ktrace_target(n_tim)
  common /c_ktrace_target/ ktrace_target

  if (ktrace_target(id) .eq. 1) then
    call ktrace_start(id)
  endif
#endif

#ifdef USE_PROF
  call start_collection(ccc)
#endif

#ifdef USE_FAPP
  call FAPP_START(ccc,1,1)
#endif
!
  if(iteration <= n_ite) then
!
#ifdef USE_WTIME
  tsta(idx) = MPI_WTIME()
#else
  call getclkreg(tsta(idx))
#endif
  ncount(idx,iteration) = ncount(idx,iteration) + 1
!
  end if
!
!$omp end master
  return
end subroutine timer_sta

!-----------------------------------------------------------------------
subroutine timer_end(id)
  use m_IterationNumbers,   only  : iteration
  use m_timer
  implicit none
#ifdef USE_WTIME
  include 'mpif.h'
#endif
  integer :: id,idx
  character*5 cid
  character*3 citer
  character*11 ccc

! if (.not.tflag) return
!$omp master
  if(id < 0) then
!   idx = id * (-1) + 90000
    idx = id * (-1)
  else
    idx = id
  endif

  write(cid,'(I5.5)') idx
  write(citer,'(I3.3)') iteration
  ccc = "r_"//citer//"_"//cid

#ifdef USE_MEM
  if (myrank .eq. 0) then
    integer :: j
    character(1000) :: line
    character(15),dimension(24) :: str

    open(unitmr,file='/proc/self/stat', form="formatted", access="stream", status='old')
    read(unitmr,'(a)') line
    close(unitmr)
    read(line,*) (str(j),j=1,24)
#ifdef USE_WTIME
    write(unitmw,*) ccc, " end", MPI_WTIME()-t_sta, str(23), str(24)
#else
    call getclkreg(t_now)
    write(unitmw,*) ccc, " end", (t_now-t_sta)*5d-10, str(23),str(24)
#endif
  endif
#endif

#ifdef USE_KTRACE
      integer ktrace_target(n_tim)
      common /c_ktrace_target/ ktrace_target

      if (ktrace_target(id) .eq. 1) then
        call ktrace_stop(id)
      endif
#endif

#ifdef USE_PROF
  call stop_collection(ccc)
#endif

#ifdef USE_FAPP
  call FAPP_STOP(ccc,1,1)
#endif

#ifdef USE_TOFUPA
  call fj_tofupa_stop(ccc,id,1)
#endif
!
  if(iteration <= n_ite) then
!
#ifdef USE_WTIME
  timers(idx,iteration) = timers(idx,iteration) + (MPI_WTIME() - tsta(idx))
#else
  call getclkreg(t_now)
  timers(idx,iteration) = timers(idx,iteration) + (t_now       - tsta(idx))
#endif
  ncnt_e(idx,iteration) = ncnt_e(idx,iteration) + 1
!
  end if
!
!$omp end master
  return
end subroutine timer_end

!-----------------------------------------------------------------------
subroutine timer_bar(id)
  use m_IterationNumbers,   only  : iteration
  use m_timer
  implicit none
  include 'mpif.h'
  integer :: id
  integer :: ierr

!
  if(iteration > n_ite) return
!
!$omp master
#ifdef USE_WTIME
  tsta(id) = MPI_WTIME()
#else
  call getclkreg(tsta(id))
#endif
  ncount(id,iteration) = ncount(id,iteration) + 1
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#ifdef USE_WTIME
  timers(id,iteration) = timers(id,iteration) + (MPI_WTIME() - tsta(id))
#else
  call getclkreg(t_now)
  timers(id,iteration) = timers(id,iteration) + (t_now       - tsta(id))
#endif
  ncnt_e(id,iteration) = ncnt_e(id,iteration) + 1
!$omp end master

  return
end subroutine timer_bar

!-----------------------------------------------------------------------
subroutine timer_barrier(mpi_comm)
  use m_IterationNumbers,   only  : iteration
  use m_timer
  implicit none
  include 'mpif.h'
  integer :: mpi_comm
  integer :: ierr

!$omp master
  call MPI_BARRIER(MPI_COMM, ierr)
!$omp end master

  return
end subroutine timer_barrier

