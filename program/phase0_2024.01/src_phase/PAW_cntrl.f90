logical function itr_first_paw()
    use m_IterationNumbers,   only : iteration
    use m_PseudoPotential,    only : flg_paw
    implicit none
        
    itr_first_paw=(iteration==1.and.flg_paw)
    return
end function itr_first_paw

logical function use_paw()
    use m_PseudoPotential,    only : flg_paw
    use_paw=flg_paw
    return
end function use_paw
    
    
    
    
