  !We build the transposed H here, so we take conjg part
  !to comply with the MPI decomposition of the matrix.
  do i=MpiIstart,MpiIend
     m  = Hsector%H(1)%map(i)
     ib  = bdecomp(m,2*Ns)
     !
     htmp = zero
     !
     do io=1,Ns
        do jo=1,Ns
           !UP electrons
           Jcondition = (Hij(1,io,jo)/=zero) .AND. (ib(jo)==1) .AND. (ib(io)==0)
           if (Jcondition) then
              call c(jo,m,k1,sg1)
              call cdg(io,k1,k2,sg2)
              j    = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(Hij(1,io,jo))*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0d,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0d,htmp,i,j)
              end select
           endif
           !DW electrons
           Jcondition = (Hij(Nspin,io,jo)/=zero) .AND. (ib(jo+Ns)==1) .AND. (ib(io+Ns)==0)
           if (Jcondition) then
              call c(jo+Ns,m,k1,sg1)
              call cdg(io+Ns,k1,k2,sg2)
              j    = binary_search(Hsector%H(1)%map,k2)
              htmp = conjg(Hij(Nspin,io,jo))*sg1*sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0d,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0d,htmp,i,j)
              end select
           endif
        enddo
     enddo

  enddo







  !MOVED OUT
  ! do io=1,iNs-1
  !    !UP impurity: C^+_i C_(i+1) i-->i+1
  !    Jcondition =  (t_imp/=0d0) .AND. (Nimp_up(io+1)==1) .AND. (Nimp_up(io)==0)
  !    if (Jcondition) then
  !       call c(ket_imp_index(io+1,1),m,k1,sg1)
  !       call cdg(ket_imp_index(io,1),k1,k2,sg2)
  !       j    = binary_search(Hsector%H(1)%map,k2)
  !       htmp = t_imp*sg1*sg2
  !       !
  !       select case(MpiStatus)
  !       case (.true.)
  !          call sp_insert_element(MpiComm,spH0d,htmp,i,j)
  !       case (.false.)
  !          call sp_insert_element(spH0d,htmp,i,j)
  !       end select
  !    endif
  !    !UP impurity: C^+_(i+1) C_i i<--i+1
  !    Jcondition =  (t_imp/=0d0) .AND. (Nimp_up(io)==1) .AND. (Nimp_up(io+1)==0)
  !    if (Jcondition) then
  !       call c(ket_imp_index(io,1),m,k1,sg1)
  !       call cdg(ket_imp_index(io+1,1),k1,k2,sg2)
  !       j    = binary_search(Hsector%H(1)%map,k2)
  !       htmp = t_imp*sg1*sg2
  !       !
  !       select case(MpiStatus)
  !       case (.true.)
  !          call sp_insert_element(MpiComm,spH0d,htmp,i,j)
  !       case (.false.)
  !          call sp_insert_element(spH0d,htmp,i,j)
  !       end select
  !    endif
  !    !
  !    !DW impurity C^+_i C_(i+1) i-->i+1
  !    Jcondition = (t_imp/=0d0) .AND. (Nimp_dw(io+1)==1) .AND. (Nimp_dw(io)==0)
  !    if (Jcondition) then
  !       call c(ket_imp_index(io+1,2),m,k1,sg1)
  !       call cdg(ket_imp_index(io,2),k1,k2,sg2)
  !       j    = binary_search(Hsector%H(1)%map,k2)
  !       htmp = t_imp*sg1*sg2
  !       !
  !       select case(MpiStatus)
  !       case (.true.)
  !          call sp_insert_element(MpiComm,spH0d,htmp,i,j)
  !       case (.false.)
  !          call sp_insert_element(spH0d,htmp,i,j)
  !       end select
  !    endif
  !    !DW impurity C^+_(i+1) C_i i<--i+1
  !    Jcondition = (t_imp/=0d0) .AND. (Nimp_dw(io)==1) .AND. (Nimp_dw(io+1)==0)
  !    if (Jcondition) then
  !       call c(ket_imp_index(io,2),m,k1,sg1)
  !       call cdg(ket_imp_index(io+1,2),k1,k2,sg2)
  !       j    = binary_search(Hsector%H(1)%map,k2)
  !       htmp = t_imp*sg1*sg2
  !       !
  !       select case(MpiStatus)
  !       case (.true.)
  !          call sp_insert_element(MpiComm,spH0d,htmp,i,j)
  !       case (.false.)
  !          call sp_insert_element(spH0d,htmp,i,j)
  !       end select
  !    endif
  ! enddo


