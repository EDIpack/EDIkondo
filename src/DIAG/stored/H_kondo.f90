  if(any([Jk_z,Jk_xy]/=0d0))then
     do i=MpiIstart,MpiIend
        m   = Hsector%H(1)%map(i)
        ib  = bdecomp(m,2*Ns)
        !
        Nele_up = ib(1:eNs)
        Nele_dw = ib(eNs+1:2*eNs)
        Nimp_up = ib(2*eNs+1:2*eNs+iNs)
        Nimp_dw = ib(2*eNs+iNs+1:2*eNs+2*iNs)
        Nup     = [Nele_up,Nimp_up]
        Ndw     = [Nele_dw,Nimp_dw]
        Sele_z  = 0.5d0*(Nele_up - Nele_dw)
        Simp_z  = 0.5d0*(Nimp_up - Nimp_dw)
        Sz      = 0.5d0*(Nup-Ndw)
        !
        htmp = zero
        !
        ! KONDO COUPLING
        do iimp=1,iNs
           if(Jkindx(iimp)==0)cycle
           isite=Jkindx(iimp)
           do iorb=1,Norb
              io = pack_indices(isite,iorb)
              htmp = htmp -2.d0*Jk_z*Sele_z(io)*Simp_z(iimp)
           enddo
        enddo
        select case(MpiStatus)
        case (.true.)
           call sp_insert_element(MpiComm,spH0d,htmp,i,i)
        case (.false.)
           call sp_insert_element(spH0d,htmp,i,i)
        end select
        !
        do iimp=1,iNs
           if(Jkindx(iimp)==0)cycle
           isite = Jkindx(iimp)
           do iorb=1,Norb
              io    = pack_indices(isite,iorb)
              io_up  = ket_site_index(io,1)
              io_dw  = ket_site_index(io,2)
              imp_up = ket_imp_index(iimp,1)
              imp_dw = ket_imp_index(iimp,2)
              !c^+_up d^+_dw c_dw  d_up
              Jcondition=(&
                   (Nimp_up(iimp)==1).AND.&
                   (Nele_dw(io)  ==1).AND.&
                   (Nimp_dw(iimp)==0).AND.&
                   (Nele_up(io)  ==0) )
              if(Jcondition)then
                 call c(imp_up,m,k1,sg1)    !d_up
                 call c(io_dw,k1,k2,sg2)    !c_dw
                 call cdg(imp_dw,k2,k3,sg3) !d^+_dw
                 call cdg(io_up,k3,k4,sg4)  !c^+_up
                 j=binary_search(Hsector%H(1)%map,k4)
                 htmp = one*Jk_xy*sg1*sg2*sg3*sg4
                 !
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0d,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0d,htmp,i,j)
                 end select
              endif
              !
              ! c^+_dw d^+_up c_up  d_dw                 
              Jcondition=(&
                   (Nimp_dw(iimp)==1).AND.&
                   (Nele_up(io)  ==1).AND.&
                   (Nele_dw(io)  ==0).AND.&
                   (Nimp_up(iimp)==0) )
              if(Jcondition)then
                 call c(imp_dw,m,k1,sg1)    !d_dw
                 call c(io_up,k1,k2,sg2)    !c_up
                 call cdg(imp_up,k2,k3,sg3) !d^+_up
                 call cdg(io_dw,k3,k4,sg4)  !c^+_dw
                 j=binary_search(Hsector%H(1)%map,k4)
                 htmp = one*Jk_xy*sg1*sg2*sg3*sg4
                 !
                 select case(MpiStatus)
                 case (.true.)
                    call sp_insert_element(MpiComm,spH0d,htmp,i,j)
                 case (.false.)
                    call sp_insert_element(spH0d,htmp,i,j)
                 end select
              endif
              !
           enddo
        enddo
        !
        !Kondo exchange.
        if(Vkondo/=0d0)then
           htmp=zero
           do iimp=1,iNs
              if(Jkindx(iimp)==0)cycle
              isite=Jkindx(iimp)
              do iorb=1,Norb
                 io = pack_indices(isite,iorb)
                 htmp = htmp + Vkondo*(&
                      Nele_Up(io)*Nimp_Up(iimp) + &
                      Nele_Up(io)*Nimp_Dw(iimp) + &
                      Nele_Dw(io)*Nimp_Up(iimp) + &
                      Nele_Dw(io)*Nimp_Dw(iimp) )
              enddo
           enddo
           select case(MpiStatus)
           case (.true.)
              call sp_insert_element(MpiComm,spH0d,htmp,i,i)
           case (.false.)
              call sp_insert_element(spH0d,htmp,i,i)
           end select
        endif
        !
     enddo
  endif
