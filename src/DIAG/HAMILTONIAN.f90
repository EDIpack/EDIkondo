MODULE ED_HAMILTONIAN
  USE ED_HAMILTONIAN_COMMON
  USE ED_HAMILTONIAN_STORED_HxV
  USE ED_HAMILTONIAN_DIRECT_HxV
  !
  implicit none
  private


  !>Build sparse hamiltonian of the sector
  public  :: build_Hv_sector
  public  :: delete_Hv_sector
  public  :: vecDim_Hv_sector

  !> Tridiag sparse Hamiltonian of the sector
  public  :: tridiag_Hv_sector


  !>Sparse Mat-Vec product using stored sparse matrix
  public  :: spMatVec
#ifdef _MPI
  public  :: spMatVec_MPI
#endif


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec
#ifdef _MPI
  public  :: directMatVec_MPI
#endif




contains







  !####################################################################
  !####################################################################
  !                          KONDO
  !####################################################################
  !####################################################################
  subroutine build_Hv_sector(isector,Hmat)
    integer                            :: isector
    complex(8),dimension(:,:),optional :: Hmat
    integer                            :: irank
    integer                            :: i,j
    !
    !
    call build_sector(isector,Hsector)
    !
    Dim    = Hsector%Dim
    !
    !#################################
    !          MPI SETUP
    !#################################
    mpiAllThreads=.true.
    MpiQ = Dim/MpiSize
    MpiR = 0
    if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    !
    MpiIshift = MpiRank*mpiQ
    MpiIstart = MpiRank*mpiQ + 1
    MpiIend   = (MpiRank+1)*mpiQ + mpiR
    !
#ifdef _MPI
#ifdef _DEBUG
    if(MpiStatus.AND.ed_verbose>4)then
       write(LOGfile,*)&
            "         mpiRank,   mpi_Q,   mpi_R,   mpi_Istart,   mpi_Iend,   mpi_Iend-mpi_Istart"
       do irank=0,MpiSize-1
          call Barrier_MPI(MpiComm)
          if(MpiRank==irank)then
             write(LOGfile,*)MpiRank,MpiQ,MpiR,MpiIstart,MpiIend,MpiIend-MpiIstart+1
          endif
       enddo
       call Barrier_MPI(MpiComm)
    endif
#endif
#endif
    !
    !
    !#################################
    !          HxV SETUP
    !#################################
    if(present(Hmat))then
       spHtimesV_p => null()
       call ed_buildh(Hmat)          
       return
    endif
    !
    select case (ed_sparse_H)
    case (.true.)
       spHtimesV_p => spMatVec
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => spMatVec_MPI
#endif
       call ed_buildh()
    case (.false.)
       spHtimesV_p => directMatVec
#ifdef _MPI
       if(MpiStatus)spHtimesV_p => directMatVec_MPI
#endif
    end select
    !
  end subroutine build_Hv_sector





  

  subroutine delete_Hv_sector()
    !
    call delete_sector(Hsector)
    Dim    = 0
#ifdef _MPI
    if(MpiStatus)then
       call sp_delete_matrix(MpiComm,spH0d)
    else
       call sp_delete_matrix(spH0d)
    endif
#else
    call sp_delete_matrix(spH0d)
#endif
    !
    spHtimesV_p => null()
    !
#ifdef _MPI
    if(MpiStatus)then
       MpiComm = MpiComm_Global
       MpiSize = get_Size_MPI(MpiComm_Global)
       MpiRank = get_Rank_MPI(MpiComm_Global)
       mpiQ=0
       mpiR=0
       mpiIstart=0
       mpiIend=0
       mpiIshift=0
    endif
#endif
    !
  end subroutine delete_Hv_sector





  
  function vecDim_Hv_sector(isector) result(vecDim)
    integer :: isector
    integer :: vecDim
    integer :: Dim
    !
    Dim  = getdim(isector)
    !
#ifdef _MPI
    if(MpiStatus)then
       MpiQ = Dim/MpiSize
       MpiR = 0
       if(MpiRank==(MpiSize-1))MpiR=mod(Dim,MpiSize)
    else
       MpiQ = Dim
       MpiR = 0
    endif
#else
    MpiQ = Dim
    MpiR = 0
#endif
    !
    vecDim=MpiQ + MpiR
    !
  end function vecDim_Hv_sector





  
  subroutine tridiag_Hv_sector(isector,vvinit,alanc,blanc,norm2)
    integer                             :: isector
    complex(8),dimension(:)             :: vvinit
    real(8),dimension(:),allocatable    :: alanc,blanc
    real(8)                             :: norm2
    complex(8),dimension(:),allocatable :: vvloc
    integer                             :: vecDim
    !
    !
    if(MpiMaster)then
       norm2=dot_product(vvinit,vvinit)
       vvinit=vvinit/sqrt(norm2)
    endif
#ifdef _MPI
    if(MpiStatus)call bcast_MPI(MpiComm,norm2)
#endif
    call build_Hv_sector(isector)
    allocate(alanc(Hsector%Nlanc),blanc(Hsector%Nlanc))
    alanc=0d0 ; blanc=0d0
    if(norm2/=0d0)then
#ifdef _MPI
       if(MpiStatus)then
          vecDim = vecDim_Hv_sector(isector)
          allocate(vvloc(vecDim))
          call scatter_vector_MPI(MpiComm,vvinit,vvloc)
          call sp_lanc_tridiag(MpiComm,spHtimesV_p,vvloc,alanc,blanc)
       else
          call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
       endif
#else
       call sp_lanc_tridiag(spHtimesV_p,vvinit,alanc,blanc)
#endif
    endif
    call delete_Hv_sector()
  end subroutine tridiag_Hv_sector



end MODULE ED_HAMILTONIAN
