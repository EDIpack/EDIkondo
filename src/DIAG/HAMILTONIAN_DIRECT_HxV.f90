! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_DIRECT_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec
#ifdef _MPI
  public  :: directMatVec_MPI
#endif



contains


  subroutine directMatVec(Nloc,vin,Hv)
    integer                             :: Nloc
    complex(8),dimension(Nloc)          :: vin
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vt,Hvt
    complex(8),dimension(Nspin,eNs,eNs) :: Hij,Hloc
    real(8),dimension(Nspin,eNs)        :: Hdiag
    integer,dimension(2*Ns)             :: ib
    integer,dimension(Ns)               :: Nup,Ndw
    integer,dimension(eNs)              :: Nele_up,Nele_dw
    integer,dimension(iNs)              :: Nimp_up,Nimp_dw
    real(8),dimension(Ns)               :: Sz
    real(8),dimension(eNs)              :: Sele_z
    real(8),dimension(iNs)              :: Simp_z
    integer                             :: i,j,io_up,io_dw,imp_up,imp_dw
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = dreal(diagonal(Hloc(ispin,:,:)))
    enddo
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !
    do j=MpiIstart,MpiIend
       m   = Hsector%H(1)%map(j)
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
       !LOCAL HAMILTONIAN TERMS
       include "direct/HxV_diag.f90"
       !
       !NON-LOCAL INTERACTION HAMILTONIAN TERMS
       include "direct/HxV_se_ph.f90"
       !
       !KONDO COUPLING HAMILTONIAN TERMS
       include "direct/HxV_kondo.f90"
       !
       !HOPPING TERMS
       include "direct/HxV_hop.f90"
       !
    enddo
    !-----------------------------------------------!
    !
    !
  end subroutine directMatVec




#ifdef _MPI
  subroutine directMatVec_MPI(Nloc,v,Hv)
    integer                             :: Nloc,N
    complex(8),dimension(Nloc)          :: v
    complex(8),dimension(Nloc)          :: Hv
    complex(8),dimension(:),allocatable :: vin
    integer,dimension(2*Ns)             :: ib
    integer,dimension(Ns)               :: Nup,Ndw
    integer,dimension(eNs)              :: Nele_up,Nele_dw
    integer,dimension(iNs)              :: Nimp_up,Nimp_dw
    real(8),dimension(Ns)               :: Sz
    real(8),dimension(eNs)              :: Sele_z
    real(8),dimension(iNs)              :: Simp_z
    integer                             :: io_up,io_dw,imp_up,imp_dw
    complex(8),dimension(Nspin,eNs,eNs) :: Hij,Hloc
    real(8),dimension(Nspin,eNs)        :: Hdiag
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index    
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    call Hij_get(Hij)
    call Hij_get(Hloc)
    do ispin=1,Nspin
       Hdiag(ispin,:) = dreal(diagonal(Hloc(ispin,:,:)))
    enddo
    !
    !MPI part:
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    allocate(vin(N)) ; vin = zero
    call allgather_vector_MPI(MpiComm,v,vin)
    !
    Hv=zero
    !
    !-----------------------------------------------!
    !
    states: do j=MpiIstart,MpiIend
       m   = Hsector%H(1)%map(j)
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
       !LOCAL HAMILTONIAN TERMS
       include "direct/HxV_diag.f90"
       !
       !NON-LOCAL INTERACTION HAMILTONIAN TERMS
       include "direct/HxV_se_ph.f90"
       !
       !KONDO COUPLING HAMILTONIAN TERMS
       include "direct/HxV_kondo.f90"
       !
       !HOPPING TERMS
       include "direct/HxV_hop.f90"
       !
    enddo states
    !
    !-----------------------------------------------!
    deallocate(vin)
    return
  end subroutine directMatVec_MPI
#endif


END MODULE ED_HAMILTONIAN_DIRECT_HXV
