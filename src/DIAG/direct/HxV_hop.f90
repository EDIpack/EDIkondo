  htmp = zero
  !
  do io=1,Ns
     do jo=1,Ns
        !UP electrons
        Jcondition = (Hij(1,io,jo)/=zero) .AND. (ib(jo)==1) .AND. (ib(io)==0)
        if (Jcondition) then
           call c(jo,m,k1,sg1)
           call cdg(io,k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = Hij(1,io,jo)*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
        !DW electrons
        Jcondition = (Hij(Nspin,io,jo)/=zero) .AND. (ib(jo+Ns)==1) .AND. (ib(io+Ns)==0)
        if (Jcondition) then
           call c(jo+Ns,m,k1,sg1)
           call cdg(io+Ns,k1,k2,sg2)
           i    = binary_search(Hsector%H(1)%map,k2)
           htmp = Hij(Nspin,io,jo)*sg1*sg2
           !
           if(i/=0)hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
        endif
     enddo
  enddo
