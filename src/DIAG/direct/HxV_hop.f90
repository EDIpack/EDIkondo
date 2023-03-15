  htmp = zero
  !
  do io=1,Ns
     do jo=1,Ns
        !UP electrons
        Jcondition = (Hij(1,io,jo)/=zero) .AND. (Nup(jo)==1) .AND. (Nup(io)==0)
        if (Jcondition) then
           call c(ket_index(jo,1),m,k1,sg1)
           call cdg(ket_index(io,1),k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = Hij(1,io,jo)*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
        !DW electrons
        Jcondition = (Hij(Nspin,io,jo)/=zero) .AND. (Ndw(jo)==1) .AND. (Ndw(io)==0)
        if (Jcondition) then
           call c(ket_index(jo,2),m,k1,sg1)
           call cdg(ket_index(io,2),k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = Hij(Nspin,io,jo)*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo
