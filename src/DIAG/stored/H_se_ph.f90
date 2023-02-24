  if(jhflag)then
     !We build the transposed H so we take conjg. if not symmetric 
     !to comply with the MPI decomposition of the matrix.
     !A better MPI handling might be necessary here...
     do i=MpiIstart,MpiIend
        m  = Hsector%H(1)%map(i)
        ib  = bdecomp(m,2*Ns)
        !
        Nele_up = ib(1:eNs)
        Nele_dw = ib(eNs+1:2*eNs)
        Nimp_up = ib(2*eNs+1:2*eNs+iNs)
        Nimp_dw = ib(2*eNs+iNs+1:2*eNs+2*iNs)
        Nup     = [Nele_up,Nimp_up]
        Ndw     = [Nele_dw,Nimp_dw]
        !
        htmp = zero
        !
        ! SPIN-EXCHANGE (S-E) TERMS
        !    S-E: J Cdg_a.up Cdg_b.dw C_a.dw C_b.up
        if(Jx/=0d0)then
           do iorb=1,Norb
              do jorb=1,Norb
                 do isite=1,Nsites(iorb)
                    do jsite=1,Nsites(jorb)
                       if(isite/=jsite)cycle !local interaction only:
                       io = pack_indices(isite,iorb)
                       jo = pack_indices(isite,jorb)
                       Jcondition=(&
                            (io/=jo).AND.&
                            (Nele_Up(jo)==1).AND.&
                            (Nele_Dw(io)==1).AND.&
                            (Nele_Dw(jo)==0).AND.&
                            (Nele_Up(io)==0))
                       if(Jcondition)then
                          call c(jo,m,k1,sg1)
                          call c(io+eNs,k1,k2,sg2)
                          call cdg(jo+eNs,k2,k3,sg3)
                          call cdg(io,k3,k4,sg4)
                          j=binary_search(Hsector%H(1)%map,k4)
                          htmp = one*Jx*sg1*sg2*sg3*sg4
                          !
                          select case(MpiStatus)
                          case (.true.)
                             call sp_insert_element(MpiComm,spH0d,htmp,i,j)
                          case (.false.)
                             call sp_insert_element(spH0d,htmp,i,j)
                          end select
                          !
                       endif
                    enddo
                 enddo
              enddo
           enddo
        endif
        !
        ! PAIR-HOPPING (P-H) TERMS
        !    P-H: J Cdg_a.up Cdg_a.dw   C_b.dw   C_b.up
        if(Jp/=0d0)then
           do iorb=1,Norb
              do jorb=1,Norb
                 do isite=1,Nsites(iorb)
                    do jsite=1,Nsites(jorb)
                       if(isite/=jsite)cycle !local interaction only:
                       io = pack_indices(isite,iorb)
                       jo = pack_indices(isite,jorb)
                       Jcondition=(&
                            (Nele_Up(jo)==1).AND.&
                            (Nele_Dw(jo)==1).AND.&
                            (Nele_Dw(io)==0).AND.&
                            (Nele_Up(io)==0))
                       if(Jcondition)then
                          call c(jo,m,k1,sg1)
                          call c(jo+eNs,k1,k2,sg2)
                          call cdg(io+eNs,k2,k3,sg3)
                          call cdg(io,k3,k4,sg4)
                          j=binary_search(Hsector%H(1)%map,k4)
                          htmp = one*Jp*sg1*sg2*sg3*sg4
                          !
                          select case(MpiStatus)
                          case (.true.)
                             call sp_insert_element(MpiComm,spH0d,htmp,i,j)
                          case (.false.)
                             call sp_insert_element(spH0d,htmp,i,j)
                          end select
                          !
                       endif
                    enddo
                 enddo
              enddo
           enddo
        endif
        !
     enddo
  endif
