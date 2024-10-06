subroutine Epsilon_Paramset
  use m_Const_Parameters    ,  only : MESH
  use m_Kpoints             ,  only : way_ksample
  use m_Epsilon_ek          ,  only : auto_mode, nppcorr, way_BZintegral, e_high

  auto_mode=1
  nppcorr = 2
  e_high=3.0d0
  if(way_ksample == MESH) then
     way_BZintegral= 2
  else
     way_BZintegral= 1
  end if  
end subroutine Epsilon_Paramset

