  ! KONDO COUPLING
  if(any([Jk_z,Jk_xy]/=0d0))then     
     !
     htmp = zero
     !
     do iimp=1,iNs
        do iorb=1,Norb
           do isite=1,Nsites(iorb)
              if(isite/=Jkindx(iimp))cycle
              io = pack_indices(isite,iorb)
              htmp = -2.d0*Jk_z*Sele_z(io)*Simp_z(iimp)
           enddo
        enddo
     enddo
     i = j
     hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
     !
     !
     do iimp=1,iNs
        do iorb=1,Norb
           do isite=1,Nsites(iorb)
              if(isite/=Jkindx(iimp))cycle
              io = pack_indices(isite,iorb)
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
                 call c(imp_up,  m, k1,sg1)  !d_up
                 call c(io_dw,   k1,k2,sg2)  !c_dw
                 call cdg(imp_dw,k2,k3,sg3)  !d^+_dw
                 call cdg(io_up, k3,k4,sg4)  !c^+_up
                 i=binary_search(Hsector%H(1)%map,k4)
                 htmp = one*Jk_xy*sg1*sg2*sg3*sg4
                 !
                 if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
                 !
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
                 i=binary_search(Hsector%H(1)%map,k4)
                 htmp = one*Jk_xy*sg1*sg2*sg3*sg4
                 !
                 if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
                 !
              endif
              !
           enddo
        enddo
     enddo
     !
     !Kondo exchange Vkondo=(Vdir-Jk_z/4)
     if(Vkondo/=0d0)then
        htmp=zero
        do iimp=1,iNs
           do iorb=1,Norb
              do isite=1,Nsites(iorb)
                 if(isite/=Jkindx(iimp))cycle
                 io = pack_indices(isite,iorb)                 
                 htmp = htmp + Vkondo*(&
                      Nele_Up(io)*Nimp_Up(iimp) + &
                      Nele_Up(io)*Nimp_Dw(iimp) + &
                      Nele_Dw(io)*Nimp_Up(iimp) + &
                      Nele_Dw(io)*Nimp_Dw(iimp) )
              enddo
           enddo
        enddo
        i = j
        hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
     endif
  endif
