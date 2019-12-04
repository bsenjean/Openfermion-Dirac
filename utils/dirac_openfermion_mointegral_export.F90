!dirac_copyright_start
!      Copyright (c) 2015 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

module dirac_openfermion_mointegral_export

! Written by Bruno Senjean (October 2018, July 2019) based on the
! previous one by Lucas Visscher, VU University Amsterdam, December 2009, July 2010
! This program is aimed to interface DIRAC output (electronic integrals) with 
! the quantum computing package so-called OpenFermion.

  implicit none

  integer, parameter     :: filenumber_1e = 21
  integer, parameter     :: filenumber_2e = 22
  integer, parameter     :: filenumber_fcidump = 24
  integer, parameter     :: filenumber_mtable  = 25
  integer, parameter     :: filenumber_propint = 26
  integer, parameter     :: filenumber_ccamp = 27
  logical, parameter     :: generate_full_list = .true.
  logical, parameter     :: generate_lower_triangular = .false.
  character(10)          :: target
  character(10)          :: prop_name

  type SpinorInformation

      integer     :: irrep
      integer     :: abelian_irrep
      integer     :: occupation
      integer     :: index
      real(8)     :: energy

  end type SpinorInformation

  integer                 :: number_of_spinors, number_of_electrons, number_of_irreps, number_of_abelian_irreps
  integer                 :: inversion_symmetry, group_type, inop
  real(8)                 :: core_energy
  real(8)                 :: threshold = 1.0E-16
  integer, allocatable    :: kramer_to_spinor(:)
  integer, allocatable    :: multiplication_table(:,:)
  integer                 :: irrep_occupation(128)
  type(SpinorInformation), allocatable :: spinor(:)
  character(32)           :: ACHAR
  real(8), allocatable    ::  propr(:,:)
  real(8), allocatable    ::  propi(:,:)

  public initialize
  private process_1e, process_2e, make_index_to_occupied_first

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize ()

     integer       ::  number_of_kramerspairs, i, j
     character(10) ::  date_of_generation
     character(8)  ::  time_of_generation
     logical       ::  breit_included
     character(14) ::  representation_name

     open (filenumber_1e, file='MRCONEE', Form='UNFORMATTED')
     read (filenumber_1e) number_of_spinors, &
                          breit_included,    &
                          core_energy,       &
                          inversion_symmetry,&
                          group_type
     read (filenumber_1e) number_of_irreps,(representation_name,i=1,number_of_irreps), &
                          (irrep_occupation(i),i=1,number_of_irreps)
     read (filenumber_1e) number_of_abelian_irreps

!    count total number of electrons
     number_of_electrons = 0
     do i = 1, number_of_irreps
        number_of_electrons = number_of_electrons + irrep_occupation(i)
     end do

     allocate (multiplication_table(2*number_of_abelian_irreps,2*number_of_abelian_irreps))
     allocate (spinor(number_of_spinors))
     allocate (kramer_to_spinor(-number_of_spinors/2:number_of_spinors/2))

     read (filenumber_1e) ((multiplication_table(i,j),i=1,2*number_of_abelian_irreps),j=1,2*number_of_abelian_irreps)
     read (filenumber_1e) (spinor(i)%irrep,spinor(i)%abelian_irrep,spinor(i)%energy,i=1,number_of_spinors)

!    create indices for optional reordering according to occupation
!     call make_index_to_occupied_first
!    create indices for optional reordering according to energy
     call make_index_lowestenergy_first

     write (*,*) " Initialized reading from MRCONEE"
     write (*,*) " Core energy: ", core_energy
     write (*,*) " Breit interaction: ", breit_included
     write (*,*) " Group type (1:real, 2:complex, 4:quaternion) :",group_type

     open (filenumber_2e, file='MDCINT', Form='UNFORMATTED')
     read (filenumber_2e) date_of_generation,  &
                          time_of_generation,  &
                          number_of_kramerspairs, &
                          (kramer_to_spinor(i),kramer_to_spinor(-i),i=1,number_of_spinors/2) 

     if (2 * number_of_kramerspairs /= number_of_spinors ) then
        write (*,*) number_of_spinors, number_of_kramerspairs
        stop 'inconsistent MRCONEE and MDCINT files'
     end if 
     write (*,*) " Initialized reading from MDCINT"
     write (*,*) " MDCINT created: ",date_of_generation,"  ",time_of_generation
    
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine process_1e

     integer                :: i, j, rcw
     real(8), allocatable   :: integral(:,:,:)

     allocate (integral(number_of_spinors,number_of_spinors,2))
     read (filenumber_1e) ((integral(i,j,1),integral(i,j,2),i=1,number_of_spinors),j=1,number_of_spinors)

     select case (target)
     case ('fcidump')
        rcw = 1
        if (group_type .ne. 1) rcw = 2
        call print_1e_integral(filenumber_fcidump,integral,rcw)
     end select

     deallocate (integral)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_2e
  
     integer                :: nonzero, ikr, jkr, inz
     integer                :: i, rcw
     integer                :: ii, jj, kk, ll
     integer                :: ij, kl
     integer, allocatable   :: indk(:), indl(:)
     real(8), allocatable   :: integral(:)
     logical                :: select_integral

     allocate (integral(number_of_spinors**2))
     allocate (indk(number_of_spinors**2))
     allocate (indl(number_of_spinors**2))

     rcw = 1
     select case (target)
        case ('fcidump')

           if (group_type .ne. 1) rcw = 2
           do
              read (filenumber_2e) ikr, jkr, nonzero, (indk(inz), indl(inz), inz=1, nonzero), (integral(inz), inz=1, nonzero*rcw)
              if (ikr == 0) exit

              do inz = 1, nonzero

                 select_integral = .false.
                 if (generate_lower_triangular) then
                    ii = kramer_to_spinor(ikr)
                    jj = kramer_to_spinor(jkr)
                    kk = kramer_to_spinor(indk(inz))
                    ll = kramer_to_spinor(indl(inz))

                    ij = ii*(ii-1)/2 + jj
                    kl = kk*(kk-1)/2 + ll
                    
                    select_integral = (ii .ge. jj).and.(kk .ge. ll).and.(ij .ge. kl)
                 else
                    select_integral = .true.
                 endif

                 if ( select_integral ) then
                    call print_2e_integral(filenumber_fcidump,ikr,jkr,indk(inz),indl(inz),inz,integral,rcw)

!                   make also kramers-related integral if desired
!                   note that we assume real integrals, in which case the kr integral is identical
                    if (generate_full_list) then
                       call print_2e_integral(filenumber_fcidump,-ikr,-jkr,-indk(inz),-indl(inz),inz,integral,rcw)
                    end if   
                 end if
              end do
           end do
    end select
  
    deallocate(integral)
    deallocate(indk)
    deallocate(indl)
    
  
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_1e_integral(filenumber,integral,rcw)
     integer :: filenumber
     real(8) :: integral(:,:,:)
     integer :: i, j, rcw, end_index

     do i = 1, number_of_spinors
        if ( generate_lower_triangular ) then
           end_index = i
        else
           end_index = number_of_spinors
        end if
        do j = 1, end_index
           if (abs(integral(i,j,1)) > threshold .or. abs(integral(i,j,2)) > threshold) then
              if (rcw .ne. 1) then
                 write (filenumber_fcidump,'(1P,2E20.12,7i4)') &
                    integral(i,j,1),                &
                    integral(i,j,2),                &
!                    i,                              &
!                    j,                              &
                    spinor(i)%index,                              &
                    spinor(j)%index,                              &
                    0,                              &
                    0
              else
                 write (filenumber_fcidump,'(1P,E20.12,7i4)') &
                    integral(i,j,1),                &
!                    i,                              &
!                    j,                              &
                    spinor(i)%index,                              &
                    spinor(j)%index,                              &
                    0,                              &
                    0
              end if
           end if
        end do
     end do
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_2e_integral(filenumber,ikr,jkr,kkr,lkr,inz,integral,g_type)
     integer :: filenumber 
     integer :: i, g_type 
     integer :: ikr, jkr, kkr, lkr, inz
     real(8) :: integral(:)

     if ( g_type .ne. 1 ) then
        write (filenumber,'(1P,2E20.12,7i4)') (integral(g_type*(inz-1)+i), i=1,g_type), &
              spinor(kramer_to_spinor(ikr))%index,                  &
              spinor(kramer_to_spinor(jkr))%index,                  &
              spinor(kramer_to_spinor(kkr))%index,                  &
              spinor(kramer_to_spinor(lkr))%index
     else
        write (filenumber,'(1P,E20.12,7i4)') integral(inz), &
              spinor(kramer_to_spinor(ikr))%index,                  &
              spinor(kramer_to_spinor(jkr))%index,                  &
              spinor(kramer_to_spinor(kkr))%index,                  &
              spinor(kramer_to_spinor(lkr))%index
     end if

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_index_to_occupied_first

     integer :: index, i, j, n

!    We may want to reorder spinors such that we have all occupied first, create index for that
     do i = 1, number_of_irreps
        n = 0
        do j = 1, number_of_spinors
           if (spinor(j)%irrep == i) then
              n = n + 1
              if (n <= irrep_occupation(i)) then
                 index = index + 1
                 spinor(j)%index = index
                 spinor(j)%occupation = 1
              end if
           end if
        end do
     end do
     do i = 1, number_of_irreps
        n = 0
        do j = 1, number_of_spinors
           if (spinor(j)%irrep == i) then
              n = n + 1
              if (n > irrep_occupation(i)) then
                 index = index + 1
                 spinor(j)%index = index
                 spinor(j)%occupation = 0
              end if
           end if
        end do
     end do

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_index_lowestenergy_first

     integer :: index, i, j, k

!     do i = 1, number_of_irreps
!        do j = 1, number_of_spinors
!           if (spinor(j)%irrep == i) then
!              spinor(j)%index = i
!              do k = 1, number_of_spinors
!                 if ((spinor(k)%irrep == i) .and. (spinor(k)%energy < spinor(j)%energy)) then
!                    spinor(j)%index = spinor(j)%index + number_of_irreps
!                 end if
!              end do
!              do k = 1, j
!                 if ((spinor(k)%index == spinor(j)%index) .and. (j .ne. k)) then
!                    spinor(j)%index = spinor(j)%index + number_of_irreps
!                 end if
!              end do
!           end if
!        end do
!     end do

     do i = 1, number_of_spinors
        spinor(i)%index = 1
        do j = 1, number_of_spinors
           if (spinor(i)%energy > spinor(j)%energy) then
              spinor(i)%index = spinor(i)%index + 1
           end if
        end do
     end do

     do i = 1, number_of_spinors
        do j = i+1, number_of_spinors
           if ((spinor(i)%energy == spinor(j)%energy) .and. (spinor(i)%irrep == spinor(j)%irrep)) then
              spinor(j)%index = spinor(j)%index + 2
           end if
        end do
     end do

     do i = 1, number_of_spinors
        do j = i+1, number_of_spinors
           if ((spinor(i)%index == spinor(j)%index) .and. (spinor(i)%irrep .ne. spinor(j)%irrep)) then
              spinor(j)%index = spinor(j)%index + 1
           end if
        end do
     end do

 

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_fcidump_file

  integer               :: i, j

  open  (filenumber_fcidump, file='FCIDUMP', Form='FORMATTED')
  write (filenumber_fcidump,'(A,I5,A)') "&FCI NORB=",number_of_spinors,","
  write (filenumber_fcidump,'(A,I5,A)') "    NELEC=",number_of_electrons,","
  write (filenumber_fcidump,'(A,30(I2,A))') "    ORBSYM=", &
  & (spinor(i)%abelian_irrep,",",i=1,number_of_spinors) 
  write (filenumber_fcidump,'(A,I5,A)') "    ISYM=",(2 * number_of_abelian_irreps),","
  write (filenumber_fcidump,'(A)') "    IUHF=1,"
  write (filenumber_fcidump,'(A)') "&END"

  open  (filenumber_mtable, file='FCITABLE', Form='FORMATTED')
  write (filenumber_mtable,'(2i6)') 2 * number_of_abelian_irreps
  do j = 1, 2 * number_of_abelian_irreps
     write (filenumber_mtable,'(8i6)') (multiplication_table(i,j),i=1,2*number_of_abelian_irreps)
  end do

  call process_2e
  call process_1e

  do j = 1, number_of_spinors
    if (group_type .ne. 1) then
        write (filenumber_fcidump,'(1P,2E20.12,7i4)') (spinor(j)%energy),0.0d0,spinor(j)%index,0,0,0
    else
        write (filenumber_fcidump,'(1P,E20.12,7i4)') (spinor(j)%energy),spinor(j)%index,0,0,0
    end if
  end do

  if (group_type .ne. 1) then
     write (filenumber_fcidump,'(1P,2E20.12,7i4)') core_energy,0.0d0,0,0,0,0
  else
     write (filenumber_fcidump,'(1P,E20.12,7i4)') core_energy,0,0,0,0
  end if

  deallocate (multiplication_table)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_propint_file

      integer               :: i, j

      allocate (propr(number_of_spinors,number_of_spinors))
      allocate (propi(number_of_spinors,number_of_spinors))

      open (filenumber_propint,FILE='MDPROP',FORM='UNFORMATTED')
      inop = 0
    1 read (filenumber_propint,ERR=10,END=11) ACHAR
      if (ACHAR(1:8).NE.'********'.OR.ACHAR(25:32).NE.prop_name) GOTO 1 !Once ACHAR = prop_name, it continues to read the integrals. 
      write (*,1000) prop_name,ACHAR(9:16),ACHAR(17:24)
      read (filenumber_propint) ((propr(i,j),propi(i,j),i=1,number_of_spinors),j=1,number_of_spinors)
      close (filenumber_propint,STATUS='KEEP')
      GOTO 12
   10 INOP = 1
      GOTO 12
   11 INOP = 2
   12 CONTINUE

!
!     Error exit if the integrals could not be read
!
      if (inop.EQ.1) then
       write (*,*) ' Error reading property ',prop_name,' on file MDPROP'
       stop ' Error reading property integrals'
      end if
      if (inop.EQ.2) then
       write (*,*) ' Property ',prop_name,' not found on file MDPROP'
       stop ' Property integrals missing'
      end if
 
 1000 format (/' Read integral type ',A8,' created ',  &
             A8,' storage info : ',A8)

!  deallocate(propr) ; deallocate(propi)
  
      open (filenumber_propint, file='PROPINT', Form = 'FORMATTED')
      write (filenumber_propint,'(A,I5,A)') "&FCI NORB=",number_of_spinors,","
      write (filenumber_propint,'(A,I5,A)') "    NELEC=",number_of_electrons,","
      write (filenumber_propint,'(A,30(I2,A))') "    ORBSYM=", &
      & (spinor(i)%abelian_irrep,",",i=1,number_of_spinors)
      write (filenumber_propint,'(A,I5,A)') "    ISYM=",(2 * number_of_abelian_irreps),","
      write (filenumber_propint,'(A)') "    IUHF=1,"
      write (filenumber_propint,'(A,A)') "   PROP_NAME = ", prop_name
      write (filenumber_propint,'(A)') "&END"

      do i = 1, number_of_spinors
         do j = 1, number_of_spinors
            if (abs(propr(i,j)) > threshold .or. abs(propi(i,j)) > threshold) then
              write (filenumber_propint,'(1P,2E20.12,7i4)') &
                 propr(i,j),                &
                 propi(i,j),                &
                 spinor(i)%index,                              &
                 spinor(j)%index,                              &
                 0,                              &
                 0
            end if
         end do
      end do
    
      deallocate(propr) ; deallocate(propi)

  end subroutine

  subroutine write_ccamp_file

      ! NOT YET WORKING
      ! This subroutine should write CC amplitudes into a readable format that can be interpreted in the OpenFermion-Dirac interface.

      real(8), allocatable :: T1(:),T2(:), BUFF(:)
      integer              :: occ_spinors, vir_spinors, ini, N, NDIMT1, NDIMT2

      occ_spinors = number_of_electrons
      vir_spinors = number_of_spinors - occ_spinors
      NDIMT1 = occ_spinors*vir_spinors
      NDIMT2 = (occ_spinors*(occ_spinors-1)/2)*(vir_spinors*(vir_spinors-1)/2)

      print*,"occ_spinors",occ_spinors
      print*,"vir_spinors",vir_spinors
      print*,"NDIMT1",NDIMT1
      print*,"NDIMT2",NDIMT2

      allocate(T1(NDIMT1))
      allocate(T2(NDIMT2))

      open(filenumber_ccamp,FILE='MCCRES.0',ACCESS='DIRECT',STATUS='OLD')
        N = NDIMT1
        ini = 0
        do while (N .gt. 0)
           if (N .lt. 1024) then 
              allocate(BUFF(N))
              read(filenumber_ccamp) BUFF
              print*,"hey,N",N,BUFF
              T1(ini+1:ini+N) = BUFF
           else
              allocate(BUFF(1024))
              read(filenumber_ccamp) BUFF
              print*,"hey,N,",N,BUFF
              T1(ini+1:ini+1024) = BUFF
           end if
           deallocate(BUFF)
           ini = ini + 1024
           N = N - 1024
        end do
        N = NDIMT2
        do while (N .gt. 0)
           if (N .lt. 1024) then
              allocate(BUFF(N))
              read(filenumber_ccamp) BUFF
              T2(ini+1:ini+N) = BUFF
           else
              allocate(BUFF(1024))
              read(filenumber_ccamp) BUFF
              T2(ini+1:ini+1024) = BUFF
           end if
           deallocate(BUFF)
           ini = ini + 1024
           N = N - 1024
        end do
      close (filenumber_ccamp,STATUS='KEEP')

      open(filenumber_ccamp,FILE='CCAMP',FORM='FORMATTED')
        write(filenumber_ccamp,*) T1
        write(filenumber_ccamp,*) T2
      close (filenumber_ccamp,STATUS='KEEP')

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module

  program make_interface_files
  
  use dirac_openfermion_mointegral_export

  call get_command_argument(1, target)
  call get_command_argument(2, prop_name)
  write(*,*) target,prop_name
  call initialize 
  select case (target)
     case ('fcidump')
       write (*,*) ' Writing fcidump interface file ...'
       call write_fcidump_file
       write (*,*) ' fcidump file ready'
     case ('propint')
       write (*,*) ' Writing property integrals file ...'
       call write_propint_file
       write (*,*) ' propint file ready'
     case ('ccamp')
       write (*,*) ' Writing coupled-cluster amplitudes into file ... (NOT YET WORKING)'
       call write_ccamp_file
       write (*,*) ' ccamp file ready'
  write(*,*)
  end select

  end
