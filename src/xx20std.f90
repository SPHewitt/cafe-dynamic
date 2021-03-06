!$Id: xx20std.f90 2240 2017-04-18 16:36:34Z mexas $
PROGRAM xx20std
!-----------------------------------------------------------------------
! Sam Hewitt (University of Manchester)
! Anton Shterenlikht (University of Bristol)
!
! Program xx20std- linking ParaFEM with CGPACK, specifically
! modifying p129 from 5th edition to link with the cgca module.
!
! 12.9 --
!
! This version must conform to F2008 standard.
! This mainly means none of the routines using Cray extensions
! can be used.
!
! This version has been tested on Intel 17
!
! The CA CS is aligned with the FE CS, but not the same origin.
! Anyway, all cells are within the FE model, so mapping is easy.
!-----------------------------------------------------------------------
!USE mpi_wrapper  !remove comment for serial compilation

!*** CGPACK part *****************************************************72
! The CGPACK module must be used
use cgca
!*** end CGPACK part *************************************************72

USE precision
USE global_variables
USE mp_interface
USE input
USE output
USE loading
USE timing
USE maths
USE gather_scatter
USE steering
USE new_library

IMPLICIT NONE

! neq, ntot are now global variables - must not be declared
INTEGER,PARAMETER       ::  nodof=3,nst=6

INTEGER                 ::  nn,nr,nip,nod,i,j,k,l,m,iel,ndim=3,nstep
INTEGER                 ::  npri,iters,limit,ndof,nels,npes_pp
INTEGER                 ::  node_end,node_start,nodes_pp,loaded_nodes
INTEGER                 ::  nlen,nres,meshgen,partitioner,it,is,outProc
INTEGER                 ::  statu(MPI_STATUS_SIZE)

REAL(iwp),PARAMETER     ::  zero=0.0_iwp    
REAL(iwp)               ::  e,v,det,rho,alpha1,beta1,omega,theta
REAL(iwp)               ::  period,pi,dtim,volume,c1,c2,c3,c4,q
REAL(iwp)               ::  real_time,tol,up,alpha,beta
REAL(iwp)               :: outDisp(3),send,recv,tLoad
     
CHARACTER(LEN=15)       ::  element
CHARACTER(LEN=50)       ::  argv
CHARACTER(LEN=6)        ::  ch

LOGICAL                 ::consistent=.TRUE.,converged

!---------------------------- dynamic arrays -------------------------72
REAL(iwp),ALLOCATABLE ::  loads_pp(:),points(:,:),dee(:,:),fun(:)
REAL(iwp),ALLOCATABLE ::  jac(:,:),der(:,:),deriv(:,:),weights(:)
REAL(iwp),ALLOCATABLE ::  bee(:,:),g_coord_pp(:,:,:),fext_pp(:)
REAL(iwp),ALLOCATABLE ::  x1_pp(:),d1x1_pp(:),d2x1_pp(:),emm(:,:)
REAL(iwp),ALLOCATABLE ::  ecm(:,:),x0_pp(:),d1x0_pp(:),d2x0_pp(:)
REAL(iwp),ALLOCATABLE ::  store_km_pp(:,:,:),vu_pp(:),u_pp(:)
REAL(iwp),ALLOCATABLE ::  store_mm_pp(:,:,:),p_pp(:),d_pp(:),x_pp(:)
REAL(iwp),ALLOCATABLE ::  xnew_pp(:),pmul_pp(:,:),utemp_pp(:,:)
REAL(iwp),ALLOCATABLE ::  temp(:),diag_precon_pp(:)
REAL(iwp),ALLOCATABLE ::  diag_precon_tmp(:,:),temp_pp(:,:,:)
REAL(iwp),ALLOCATABLE ::  disp_pp(:),val(:,:),timest(:),eld_pp(:,:)
REAL(iwp),ALLOCATABLE ::  r_pp(:),tot_r_pp(:),eps(:),sigma(:)

INTEGER,ALLOCATABLE   ::  rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)    

!*** CGPACK part *****************************************************72
! CGPACK parameters
integer, parameter :: cgca_linum=5 ! number of loading iterations
logical( kind=ldef ), parameter :: cgca_yesdebug = .false.,             &
 cgca_nodebug = .false.
real( kind=rdef ), parameter :: cgca_zero = 0.0_rdef,                  &
 cgca_one = 1.0_rdef,                                                  &
 ! cleavage stress on 100, 110, 111 planes for BCC,
 ! see the manual for derivation, GPa.
 ! Material = Iron
 !cgca_scrit(3) = (/ 1.05e1_rdef, 1.25e1_rdef, 4.90e1_rdef /)
 cgca_scrit(3) = (/ 1.05e10_rdef, 1.25e10_rdef, 4.90e10_rdef /)


! CGPACK variables
integer( kind=idef ) ::                                                &
   cgca_ir(3),            & ! coarray codimensions
   cgca_img,              &
   cgca_nimgs,            &
   cgca_ng,               & ! number of grains in the whole model
   cgca_clvg_iter,        & ! number of cleavage iterations
!   cgca_lc(3),            & ! local coordinates of a cell with its image
!   cgca_lowr(3),          & ! local coordinates of the lower box corner
!   cgca_uppr(3),          & ! local coordinates of the upper box corner
!   cgca_iflag,            & ! 1 box in, 2 box out, 3 neither 
   cgca_liter               ! load iteration number
integer( kind=iarr ) :: cgca_c(3) ! coarray dimensions
integer( kind=iarr ), allocatable :: cgca_space(:,:,:,:) [:,:,:]

real( kind=rdef ) ::                                                   &
  cgca_qual,             & ! quality
  cgca_bsz(3),           & ! the given and the updated "box" size
  cgca_origin(3),        & ! origin of the "box" cs, in FE cloads
  cgca_rot(3,3),         & ! rotation tensor *from* FE cs *to* CA cs
  cgca_dm,               & ! mean grain size, linear dim, phys units
  cgca_res,              & ! resolutions, cells per grain
  cgca_bcol(3),          & ! lower phys. coords of the coarray on image
  cgca_bcou(3),          & ! upper phys. coords of the coarray on 
  cgca_stress(3,3),      & ! stress tensor
  cgca_length,           & ! fracture length scale
  cgca_time_inc,         & ! time increment
  cgca_lres,             & ! linear resolution, cells per unit of length
  cgca_charlen,          & ! characteristic element length
  cgca_fracvol             ! volume (number) of fractured cells per img
real( kind=rdef ), allocatable :: cgca_grt(:,:,:)[:,:,:]
logical( kind=ldef ) :: cgca_solid
! logical( kind=ldef ) :: cgca_lflag
character( len=6 ) :: cgca_citer
!*** end CGPACK part *************************************************72


!------------------------ input and initialisation ---------------------
ALLOCATE(timest(20))
timest  	=  zero
timest(1)	=  elap_time()

!*    Get the rank of the processes and the total number of processes
!* intent( out ):
!*          numpe - integer, process number (rank)
!*           npes - integer, total number of processes (size)
!CALL find_pe_procs( numpe, npes )
! CAnnot use MPI_INIT in a coarray program with ifort 16!
CALL find_pe_procs(numpe,npes)

CALL getname(argv,nlen)

CALL read_p129(argv,numpe,alpha1,beta1,e,element,limit,loaded_nodes,    &
   meshgen,nels,nip,nn,nod,npri,nr,nres,nstep,omega,partitioner,rho,     &
   theta,tol,v)

CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)

ndof  =  nod*nodof
ntot  =  ndof

ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
ALLOCATE(g_num_pp(nod,nels_pp))
ALLOCATE(rest(nr,nodof+1))

g_num_pp  	=  0
g_coord_pp	=  zero
rest		    =  0

CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)

IF(meshgen==2) CALL abaqus2sg(element,g_num_pp)
 
CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
CALL read_rest(argv,numpe,rest)

ALLOCATE(points(nip,ndim),bee(nst,ntot),dee(nst,nst),jac(ndim,ndim),     &
   weights(nip),der(ndim,nod),deriv(ndim,nod),g_g_pp(ntot,nels_pp),      &
   store_km_pp(ntot,ntot,nels_pp),utemp_pp(ntot,nels_pp),ecm(ntot,ntot), &
   pmul_pp(ntot,nels_pp),store_mm_pp(ntot,ntot,nels_pp),emm(ntot,ntot),  &
   temp_pp(ntot,ntot,nels_pp),diag_precon_tmp(ntot,nels_pp),fun(nod),    &
   eps(nst),sigma(nst))

!*** end of ParaFEM input and initialisation *************************72


!*** CGPACK part *****************************************************72
! *** CGPACK first executable statement ***
! In this test set the number of images via the env var,
! or simply as an argument to aprun.
! the code must be able to cope with any value >= 1.
cgca_img = this_image()
cgca_nimgs = num_images()

! Need to init RND before cgca_pfem_cenc, where there order of comms
! is chosen at *random* to even the remote access pattern.

! Initialise random number seed. Choose either a routine with
! reproducible seeds (cgca_ins), or random seeds (cgca_irs).
! Multiple runs with cgca_ins *on the same number of cores (images)*
! on the same platform should produce reproducible results.
!
! Argument:
! .false. - no debug output
!  .true. - with debug output
call cgca_ins( .false. )
!call cgca_irs( .true. )

! dump CGPACK parameters and some ParaFEM settings
if ( cgca_img .eq. 1 ) then
  call cgca_pdmp
  write (*,*) "Young's mod:", e, "Poisson's ratio", v
end if

! Try to separate stdout output
sync all

! The physical dimensions of the box, must be the same
! units as in the ParaFEM.
! Must be fully within the FE model, which for xx14
! is a cube with lower bound at (0,0,-10), and the
! upper bound at (10,10,0)

! xx20 The PARAFEM Model is (0,0,-1) to (1,5,0)

!cgca_bsz = (/ 10.0, 10.0, 10.0 /)
!cgca_bsz = (/ 1.0, 5.0, 1.0 /)
cgca_bsz = (/ 1.0, 1.0, 1.0 /)

! Origin of the box cs, in the same units.
! This gives the upper limits of the box at 0+10=10, 0+10=10, -10+10=0
! all within the FE model.

! xx20 Origin of box is (0,0,-1)
!cgca_origin = (/ 0.0, 0.0, -10.0 /)
!cgca_origin = (/ 0.0, 0.0, -1.0 /)
cgca_origin = (/ 0.0, 2.0, -1.0 /)

! Rotation tensor *from* FE cs *to* CA cs.
! The box cs is aligned with the box.
cgca_rot         = cgca_zero
cgca_rot( 1, 1 ) = 1.0
cgca_rot( 2, 2 ) = 1.0
cgca_rot( 3, 3 ) = 1.0

! mean grain size, also mm
! Beam is of size 1 * 5 * 1
! From Kato2015 - Iron has grain dimater 0.3 and 0.5 mm
!cgca_dm = 3.0e0_rdef
cgca_dm = 0.2e0_rdef ! 625 Grains
!cgca_dm = 0.5e0_rdef ! 40 Grains

! resolution, cells per grain
!cgca_res = 1.0e5_rdef
cgca_res = 1.0e5_rdef

! cgpack length scale, also in mm
! Equivalent to crack propagation distance per unit of time,
! i.e. per second. Let's say 1 km/s = 1.0e3 m/s = 1.0e6 mm/s. 
cgca_length = 1.0e6_rdef

! In p129_tiny, each finite element is 0.125 x 0.125 x 0.125 mm, so
! the charlen must be bigger than that.

cgca_charlen = 1.0

! each image calculates the coarray grid dimensions
call cgca_gdim( cgca_nimgs, cgca_ir, cgca_qual )

! calculate the resolution and the actual phys dimensions of the box
! subroutine cgca_cadim( bsz, res, dm, ir, c, lres, ng )
call cgca_cadim( cgca_bsz, cgca_res, cgca_dm, cgca_ir, cgca_c,         &
                 cgca_lres, cgca_ng )

! dump some stats from img 1
if (cgca_img .eq. 1 ) then
  write ( *, "(9(a,i0),tr1,g0,tr1,g0,3(a,g0),a)" )                     &
    "img: ", cgca_img  , " nimgs: ", cgca_nimgs,                       &
     " ("  , cgca_c (1), ","       , cgca_c (2), ",", cgca_c (3),      &
     ")["  , cgca_ir(1), ","       , cgca_ir(2), ",", cgca_ir(3),      &
     "] "  , cgca_ng   ,                                               &
    cgca_qual, cgca_lres,                                              &
         " (", cgca_bsz(1), ",", cgca_bsz(2), ",", cgca_bsz(3), ")"
  write (*,*) "dataset sizes for ParaView", cgca_c*cgca_ir
  write (*,"(a, es10.2, a, i0)") "Total cells in the model (real): ",  &
    product( real(cgca_c) * real(cgca_ir) ), " (int): ",               &
    product( int(cgca_c, kind=ilrg) * int(cgca_ir, kind=ilrg) )
  write(*,"(3(a,g0))")"Critical Stress: ",cgca_scrit(1),", ",            &
          cgca_scrit(2),", ",cgca_scrit(3)
end if

! Allocate space coarray with 2 layers, implicit SYNC ALL inside.
!subroutine cgca_as( l1, u1, l2, u2, l3, u3,                           &
!             col1, cou1, col2, cou2, col3, props, coarray )
call cgca_as( 1, cgca_c(1),  1, cgca_c(2),  1, cgca_c(3),              &
              1, cgca_ir(1), 1, cgca_ir(2), 1, 2, cgca_space )

! Calculate the phys. dim. of the coarray on each image
!subroutine cgca_imco( space, lres, bcol, bcou )
call cgca_imco( cgca_space, cgca_lres, cgca_bcol, cgca_bcou )

! dump box lower and upper corners from every image
write ( *,"(a,i0,2(a,3(es9.2,tr1)))" ) "img ", cgca_img,               &
       " bcol: ", cgca_bcol, "bcou: ", cgca_bcou

CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

! and now in FE cs:
write ( *,"(a,i0,2(a,3(es9.2,tr1)),a)" ) "img: ", cgca_img,            &
   " FE bcol: (",                                                      &
    matmul( transpose( cgca_rot ),cgca_bcol ) + cgca_origin,           &
  ") FE bcou: (",                                                      &
    matmul( transpose( cgca_rot ),cgca_bcou ) + cgca_origin, ")"

CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

! confirm that image number .eq. MPI process number
write (*,*) "img",cgca_img," <-> MPI proc", numpe

CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

! Allocate the tmp centroids array: cgca_pfem_centroid_tmp%r ,
! an allocatable array component of a coarray variable of derived type.
call cgca_pfem_ctalloc( ndim, nels_pp )

! Set the centroids array component on this image, no remote comms.
! first dim - coord, 1,2,3
! second dim - element number, always starting from 1
! g_coord_pp is allocated as g_coord_pp( nod, ndim, nels_pp )
cgca_pfem_centroid_tmp%r = sum( g_coord_pp(:,:,:), dim=1 ) / nod

         ! set cgca_pfem_centroid_tmp[*]%r
sync all ! must add execution segment
         ! use cgca_pfem_centroid_tmp[*]%r

! Set lcentr private arrays on every image. Choose one of the two
! routines that do this:
! - cgca_pfem_cenc - uses all-to-all algorithm.
! - cgca_pfem_map  - uses CO_SUM, CO_MAX and *large* tmp arrays
! Both routines have identical sets of input arguments.
call cgca_pfem_cenc( cgca_origin, cgca_rot, cgca_bcol, cgca_bcou )
!call cgca_pfem_map( cgca_origin, cgca_rot, cgca_bcol, cgca_bcou )

! Dump lcentr for debug
! *** a lot *** of data
!call cgca_pfem_lcentr_dump

! Allocate cgca_pfem_integrity%i(:), array component of a coarray of
! derived type. Allocating *local* array.
! i is set to 1.0 on allocation.
call cgca_pfem_integalloc( nels_pp )
 
! Allocate the Young's modulus 2D array
call cgca_pfem_ealloc( nip, nels_pp )
  
! initially set the Young's modulus to "e" everywhere
cgca_pfem_enew = e

! Dump the FE centroids in CA cs to stdout.
! Obviously, this is an optional step, just for debug.
! Can be called from any or all images.
! call cgca_pfem_cendmp ! NO DUMP IN THIS VERSION - TAKES TOO LONG!!!

! Generate microstructure

! Allocate rotation tensors, implicit SYNC ALL inside.
call cgca_art( 1, cgca_ng, 1, cgca_ir(1), 1, cgca_ir(2), 1, cgca_grt )

! Set initial values to all layers of the space array.
cgca_space( :, :, :, cgca_state_type_grain ) = cgca_liquid_state
cgca_space( :, :, :, cgca_state_type_frac  ) = cgca_intact_state

! Make sure all images set their space arrays, before calling
! the grain nucleation routine, which will update the state of
! the space array.
sync all

! Set grain nuclei, SYNC ALL inside.
! last argument:
! .false. - no debug output
!  .true. - with debug output
call cgca_nr( cgca_space, cgca_ng, .false. )

! Assign rotation tensors, SYNC ALL inside
call cgca_rt( cgca_grt )

! Solidify, SYNC ALL inside.
!subroutine cgca_sld( coarray, periodicbc, iter, heartbeat, solid )
! second argument:
!  .true. - periodic BC
! .false. - no periodic BC
call cgca_sld( cgca_space, .false., 0, 10, cgca_solid )

! initiate grain boundaries
call cgca_igb( cgca_space )
sync all
call cgca_hxi( cgca_space )
sync all

! Smoothen the GB, several iterations. cgca_gbs has not remote comms.
! cgca_hxi has remote comms, so need to sync before and after it.
call cgca_gbs( cgca_space )
sync all
call cgca_hxi( cgca_space )
sync all
call cgca_gbs( cgca_space )
sync all
call cgca_hxi( cgca_space )
sync all

! update grain connectivity, local routine, no sync needed
call cgca_gcu( cgca_space )


! IF CRACK
! Set a single crack nucleus in the centre of the x1=x2=0 edge
! 
if ( cgca_img .eq. 1 ) then
  !cgca_space( 1, 1, cgca_c(3)/2, cgca_state_type_frac )                &
  !            [ 1, 1, cgca_ir(3)/2 ] = cgca_clvg_state_100_edge
  cgca_space( 1, 50, 1, cgca_state_type_frac )                &
              [ 1, 1, 1 ] = cgca_clvg_state_100_edge


  WRITE(*,*)"Setting crack Nucleus"
  ! SET coarray position (1,1,)
  !! cgca_c = # Cells in space coarray
  !! cgca_ir= Rearranged Coarray dimensions
end if


! Must be slightly bigger than the characteristic length of FE.
! This Value is set earlier in the Code
!cgca_charlen = 0.4

! Each image will update its own cgca_space coarray.
!subroutine cgca_pfem_partin( coarray, cadim, lres, bcol, charlen, &
! debug )
call cgca_pfem_partin( cgca_space, cgca_c, cgca_lres, cgca_bcol,       &
  cgca_charlen, cgca_yesdebug )

         ! cgca_space changed locally on every image
sync all !
         ! cgca_space used

! Now can deallocate the temp array cgca_pfem_centroid_tmp%r.
! Could've done this earlier, but best to wait until sync all is
! required, to avoid extra sync.
call cgca_pfem_ctdalloc

! Img 1 dumps space arrays to files.
! Remote comms, no sync inside, so most likely want to sync afterwards
if ( cgca_img .eq. 1 ) write (*,*) "dumping model to files"
!call cgca_fwci( cgca_space, cgca_state_type_grain, "zg0text.raw" )
call cgca_fwci( cgca_space, cgca_state_type_frac,  "zf0text.raw" )
! HEWITT
call cgca_pswci( cgca_space, cgca_state_type_grain, "zg0.raw" )
call cgca_pswci( cgca_space, cgca_state_type_frac,  "zf0.raw" )
if ( cgca_img .eq. 1 ) write (*,*) "finished dumping model to files"

! Allocate the stress array component of cgca_pfem_stress coarray
! subroutine cgca_pfem_salloc( nels_pp, intp, comp )
call cgca_pfem_salloc( nels_pp, nip, nst )

! Need a sync after file write, because not sure what is coming next...
sync all

!*** end CGPACK part *************************************************72


!----------  find the steering array and equations per process ---------
CALL rearrange(rest)
g_g_pp  = 0
neq 	  = 0

! map equations to elements
elements_0: DO iel=1,nels_pp
  CALL find_g3( g_num_pp(:,iel), g_g_pp(:,iel), rest )
END DO elements_0

neq = MAXVAL(g_g_pp)
neq = max_p(neq)
CALL calc_neq_pp
CALL calc_npes_pp( npes, npes_pp )
CALL make_ggl( npes_pp, npes, g_g_pp )

DO i=1,neq_pp
  IF(nres==ieq_start+i-1)THEN
    it=numpe
    is=i
  END IF
END DO

CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)

ALLOCATE(x0_pp(neq_pp),d1x0_pp(neq_pp),x1_pp(neq_pp),vu_pp(neq_pp),     &
   diag_precon_pp(neq_pp),u_pp(neq_pp),d2x0_pp(neq_pp),loads_pp(neq_pp), &
   d1x1_pp(neq_pp),d2x1_pp(neq_pp),d_pp(neq_pp),p_pp(neq_pp),            &
   x_pp(neq_pp),xnew_pp(neq_pp),fext_pp(neq_pp),r_pp(neq_pp),            &
   tot_r_pp(neq_pp),disp_pp(nodes_pp*ndim),temp(nodes_pp),               &
   eld_pp(ntot,nels_pp))

   x0_pp=zero; d1x0_pp=zero; x1_pp=zero; vu_pp=zero; diag_precon_pp=zero
   u_pp=zero; d2x0_pp=zero; loads_pp=zero; d1x1_pp=zero; d2x1_pp=zero
   d_pp=zero; p_pp=zero; x_pp=zero; xnew_pp=zero; fext_pp=zero 
   r_pp=zero; tot_r_pp=zero ; disp_pp	= zero ; eld_pp = zero
   temp = 0

!---------------------------- Calculate Mass Matrix -------------------------
store_mm_pp=zero

!------ NEW nstep and dtim
! Omega is in Rad/s

pi=ACOS(-1._iwp); period=2._iwp*pi/omega; !dtim=period/20._iwp
dtim=0.001

c1=(1._iwp-theta)*dtim; c2=beta1-c1; c3=alpha1+1._iwp/(theta * dtim)
c4=beta1+theta*dtim

CALL sample(element,points,weights)
elements_1: DO iel=1,nels_pp;   
  volume=zero; emm=zero; ecm=zero
  gauss_points_1: DO i=1,nip     
    CALL shape_der(der,points,i)
    jac=MATMUL(der,g_coord_pp(:,:,iel))
    det=determinant(jac)
    CALL invert(jac)   
    volume=volume+det*weights(i)
    CALL shape_fun(fun,points,i)
    IF(consistent)THEN; CALL ecmat(ecm,fun,ntot,nodof)
       ecm=ecm*det*weights(i)*rho; emm=emm+ecm
    END IF
   END DO gauss_points_1   
   IF(.NOT.consistent)THEN
     DO i=1,ntot; emm(i,i)=volume*rho/13._iwp; END DO
     DO i=1,19,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=2,20,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=3,21,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=37,55,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=38,56,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
     DO i=39,57,6; emm(i,i)=emm(4,4)*.125_iwp; END DO
   END IF
   store_mm_pp(:,:,iel)=emm
 END DO elements_1
 

! Read Loads
IF ( loaded_nodes > 0 ) THEN
     ALLOCATE( node(loaded_nodes), source=0 )
     ALLOCATE( val(ndim, loaded_nodes), source=zero )
     CALL read_loads( argv, numpe, node, val )
     CALL load( g_g_pp, g_num_pp, node, val, tot_r_pp )
     DEALLOCATE( node )
     DEALLOCATE( val )
  end if
  DEALLOCATE(g_g_pp)

IF(numpe==it) THEN
 	OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
  	WRITE(11,'(A,I6,A)') "This job ran on ",npes," processes"
   	WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes ",nr,                &
     " restrained and ", neq," equations"
    WRITE(11,'(3(A,F10.2))') " Youngs Modulus",e," Poissons Ratio",v,                &
     " Density", rho
    WRITE(11,'(4(A,F10.2))') "Alpha1",alpha1," Beta1",beta1,                &
    " Omega", omega," Theta",theta
    write( 11,'(A,F10.4)') "Time after setup is:", elap_time()-timest(1)
  END IF 

!---------------------------- initial conditions -------------------------
 x0_pp=zero; d1x0_pp=zero; d2x0_pp=zero; real_time=zero

  ! Scaling Factor for Force
  tot_r_pp=tot_r_pp * 100
  
  tLoad=SUM_P(tot_r_pp)
  IF(numpe==it) THEN
   WRITE(11,'(A,E12.4)') "Total Load: ",tLoad
   WRITE(11,'(A)') "   Time t  cos(omega*t) Displacement Iterations"
   WRITE(11,'(3E12.4,I10)') real_time,cos(omega*real_time),x1_pp(is),iters
  END IF

!-------------------------------------------------------------------------
!---------------------------- time stepping loop ------------------------- 
!-------------------------------------------------------------------------
cgca_liter=1
timesteps: DO j=1,nstep
  cgca_liter = cgca_liter + 1  !counter for cafe output

!------ element stiffness integration and build the preconditioner ---72
     
  !HEWITT
  dee = zero

  CALL sample( element, points, weights )
  store_km_pp = zero ;

  ! Caluclate Updated Stiffness Matrix
  elements_2a: DO iel=1,nels_pp
    gauss_pts_1: DO i=1,nip
      ! from cgca_m3pfem.f90:
      ! allocate( cgca_pfem_enew( nip, nels_pp ) )
      ! HEWITT
      CALL deemat( dee, cgca_pfem_enew(i,iel), v )
      !CALL deemat( dee, e, v )
      CALL shape_der( der, points, i )
      jac = MATMUL( der, g_coord_pp(:,:,iel) )
      det = determinant(jac)
      CALL invert(jac)
      deriv = MATMUL( jac, der )
      CALL beemat( bee, deriv )
      store_km_pp(:,:,iel) = store_km_pp(:,:,iel) +                             &
      MATMUL( MATMUL( TRANSPOSE(bee), dee ), bee ) * det * weights(i)
    END DO gauss_pts_1
  END DO elements_2a
  
  ! Calculate Updated Diagonal Preconditioner
  diag_precon_tmp = zero

  elements_2b: DO iel=1,nels_pp
    DO k=1,ntot
      diag_precon_tmp(k,iel) = diag_precon_tmp(k,iel)			& 
	  	+ store_mm_pp(k,k,iel)*c3+store_km_pp(k,k,iel)*c4
    END DO
  END DO elements_2b

  diag_precon_pp = zero
  CALL scatter( diag_precon_pp, diag_precon_tmp )

!---------------------------- Displacement part ------------------------------
  elements_3: DO iel=1,nels_pp
    temp_pp(:,:,iel)=store_km_pp(:,:,iel)*c2+store_mm_pp(:,:,iel)*c3
  END DO elements_3 
 
  CALL gather(x0_pp,pmul_pp)

  DO iel=1,nels_pp
    utemp_pp(:,iel)=MATMUL(temp_pp(:,:,iel),pmul_pp(:,iel))
  END DO
  u_pp = zero
  CALL scatter(u_pp,utemp_pp)
!---------------------------- velocity part ------------------------------
  temp_pp=store_mm_pp/theta
  CALL gather(d1x0_pp,pmul_pp)

  DO iel=1,nels_pp
    utemp_pp(:,iel)=MATMUL(temp_pp(:,:,iel),pmul_pp(:,iel))
  END DO
  vu_pp=zero
  CALL scatter(vu_pp,utemp_pp)   ! doesn't add to last u_pp

!----------------------------- RHS and loading --------------------------
  real_time=real_time+dtim;
  r_pp = zero
  r_pp = tot_r_pp * (theta*dtim*cos(omega*real_time)+c1*                 &
                           cos(omega*(real_time-dtim)))
 
  IF ( loaded_nodes > 0 ) THEN
    q = SUM_P( r_pp )
    tLoad=q
    IF ( numpe==it ) then
      write (*, '(A,I12,A,F10.3)') "Time Step:", j," Real Time:",real_time
      write (*,'(A,E12.4)') "The total load is:", q
    END IF
  END IF

  r_pp = u_pp+vu_pp+r_pp
  diag_precon_pp = 1._iwp/diag_precon_pp
  d_pp = diag_precon_pp*r_pp
  p_pp = d_pp
  x_pp = zero

!--------------------- preconditioned cg iterations --------------------
  iters=0
  timest(3) = elap_time()
  temp_pp = store_mm_pp*c3+store_km_pp*c4
  iterations: DO
    iters = iters + 1
    u_pp = zero
    vu_pp = zero
    pmul_pp = zero
    utemp_pp = zero
    CALL gather(p_pp,pmul_pp)
    elements_4: DO iel=1,nels_pp
      utemp_pp(:,iel) = MATMUL(temp_pp(:,:,iel),pmul_pp(:,iel))
    END DO elements_4
    u_pp=zero
    CALL scatter(u_pp,utemp_pp)
!-------------------------- pcg equation solution ----------------------
    up = DOT_PRODUCT_P(r_pp,d_pp)
    alpha = up/DOT_PRODUCT_P(p_pp,u_pp)
    xnew_pp = x_pp+p_pp*alpha
    r_pp = r_pp-u_pp*alpha
    d_pp = diag_precon_pp*r_pp
    beta = DOT_PRODUCT_P(r_pp,d_pp)/up
    p_pp = d_pp+p_pp*beta
    CALL checon_par(xnew_pp,tol,converged,x_pp)
    IF ( converged .OR. iters == limit ) EXIT
  END DO iterations

  x1_pp=xnew_pp 
  d1x1_pp=(x1_pp-x0_pp)/(theta*dtim)-d1x0_pp*(1._iwp-theta)/theta
  d2x1_pp=(d1x1_pp-d1x0_pp)/(theta*dtim)-d2x0_pp*(1._iwp-theta)/theta
  x0_pp=x1_pp; d1x0_pp=d1x1_pp; d2x0_pp=d2x1_pp; utemp_pp=zero

  IF ( numpe==it ) THEN
    write(*,'(A,I6)')"The number of iterations to convergence was ",iters
    write(*,'(A,F10.4)')"Time to solve equations was  :",                &
      elap_time()-timest(3)
    write(*,'(A,E12.4)')"The central nodal displacement is :",xnew_pp(1)
  END IF

!--------------- recover stresses at centroidal gauss point ------------
! IF CAFE
if(.TRUE. .AND. j .GE. 1500) THEN
  CALL gather(xnew_pp(1:),eld_pp)

  elmnts: DO iel = 1, nels_pp
    intpts: DO i = 1, nip
      !    call deemat(dee, cgca_pfem_enew( i, iel ), v)
      call deemat(dee, e, v)
      ! Compute the derivatives of the shape functions at a Gauss point.
      ! http://parafem.googlecode.com/svn/trunk/parafem/src/modules/shared/new_library.f90
      CALL shape_der(der,points,i)
      jac = MATMUL(der,g_coord_pp(:,:,iel))
      CALL invert(jac)
      deriv = MATMUL(jac,der)
      CALL beemat(bee,deriv)
      eps = MATMUL(bee,eld_pp(:,iel))
      sigma = MATMUL(dee,eps)
      !    write (*,*) "MPI rank", numpe, "el", iel,  &
      !                "int. point", i, "stress", sigma

      ! set the stress array on this image
      ! from cgca_m3pfem.f90:
      ! allocate( cgca_pfem_stress%stress( nels_pp, intp, comp ),
      cgca_pfem_stress%stress( iel, i, : ) = sigma

    END DO intpts
  end do elmnts
!*** end ParaFEM part ************************************************72

!*** CGPACK part *****************************************************72
! debug: dump stresses to stdout
!call cgca_pfem_sdmp
! dump stresses from last image for element 1

!if ( cgca_img .eq. cgca_nimgs ) then
!  do i = 1, nip
!   write (*,"(2(a,i0),a,6es10.2)") "img: ", cgca_nimgs,               &
!         " FE 1 int. p. ", i, " stress: ",                            &
!          cgca_pfem_stress % stress( 1 , i , : )
!  end do
!end if

! all images sync here
sync all

! Calculate the mean stress tensor per image
!   subroutine cgca_pfem_simg( simg )
!   real( kind=rdef ), intent(out) :: simg(3,3)
call cgca_pfem_simg( cgca_stress )
write (*,*) "img:", cgca_img, " mean s tensor:", cgca_stress

! all images wait for each other, to make sure the stress arrays
! are not modified until all images calculate their mean values
sync all

! no real time increments in this problem
! I use the inverse of the length scale,
! which gives 1mm of crack propagation per increment maximum.
! I then can multiply it by a factor, e.g. a factor of 3 will mean
! that I get 3mm max ( 1/3 of the model ) per load increment.

! HEWITT
! Time Increment for this problem needs researching
cgca_time_inc = 0.2 * (1.0_rdef / cgca_length)
!cgca_time_inc=dtim

! run cleavage for a correct number of iterations, which is a function
! of the characteristic length and the time increment
cgca_clvg_iter = nint( cgca_length * cgca_lres * cgca_time_inc )
if ( cgca_img .eq. 1 ) write (*,*) "load inc:", cgca_liter,            &
                                   "clvg iter:", cgca_clvg_iter

! ===>>> sync all inside <<<===
! lower the crit stresses by a factor of 100.
! On Intel 16: no support for CO_SUM yet, so use _nocosum.
! On Intel 16: subroutine cgca_clvgp_nocosum( coarray, rt, t, scrit, &
!                     sub, gcus, periodicbc, iter, heartbeat, debug )


! rt - rotation tensor
! t  - stress tensor
! scrit - critical values of cleavage stress
! sub - cleavage state calculation routine,
!      cgca_clvgsd (deterministic) or cgca_clvgsp (probablistic).
! gcus - a subrotine to use to update the grain connectivity array
!        gcupd_a (all-to-all) or gcupd_n (nearest neighbour)
! periodicbc - if .true. periodic boundary conditions are used
! iter - # of clevage iterations
! heartbeat - if >0 then dump a message every heartbeat iterations

call cgca_clvgp_nocosum( cgca_space, cgca_grt, cgca_stress,            &
     0.01 * cgca_scrit, cgca_clvgsd, cgca_gcupdn, .false.,                   &
     cgca_clvg_iter, 10, cgca_yesdebug )

! dump the model out, no sync inside
if ( cgca_img .eq. 1 ) write (*,*) "dumping model to file"
!HEWITT
write ( cgca_citer, "(i0)" ) cgca_liter
call cgca_pswci( cgca_space, cgca_state_type_frac,                     &
                "zf"//trim( cgca_citer )//".raw" )
!call cgca_fwci( cgca_space, cgca_state_type_frac,  "zf"//trim( cgca_citer )//"text.raw" )
if ( cgca_img .eq. 1 ) write (*,*) "finished dumping model to file"

sync all
     
! Calculate number (volume) of fractured cells on each image.
! cgca_fracvol is a local, non-coarray, array, so no sync needed.
! Note : Cells of states "cgca_frac_states" are considered failed.
! This includes grain boundaries (1, -1 ... -6)
call cgca_fv( cgca_space, cgca_fracvol )

write (*,*) "img:", cgca_img, "fracvol:", cgca_fracvol

! calculate integrity, remote write
call cgca_pfem_intcalc1( cgca_c, cgca_fracvol )

! dump min integrity for all FE stored on this image
write (*,*) "img:", cgca_img, "min. integrity:",                       &
  minval( cgca_pfem_integrity % i )

! wait for integrity on all images to be calculated
sync all

! Young's modulus need to be updated on each image, local arrays only.
! The original Young's modulus value must be given as input.
! For each FE stored on this image.
!HEWITT
call cgca_pfem_uym( e, nels_pp )

!HEWITT
! dump updated Young's modulus
write (*,*) "img:", cgca_img, "*min* Young's mod:",                    &
            minval( cgca_pfem_enew )

sync all

endif ! CAFE
!*** end CGPACK part *************************************************72


!*** ParaFEM part ****************************************************72

!Find what processor node 12456
IF(j.EQ.1)THEN
  DO i = 1,nels_pp
    DO k =1,nod
	    if (g_num_pp(k,i) .EQ. 12465)THEN 
  	  outProc=numpe
      endif
    ENDDO
  ENDDO
ENDIF 

! Write Displacements out
IF(.true.)THEN
IF(j/npri*npri==j) THEN
     IF(numpe==it) THEN
       WRITE(11,'(3E12.4,I10)') real_time,cos(omega*real_time),x1_pp(is),iters
     ENDIF
     IF(numpe==1) THEN
      WRITE(ch,'(I6.6)') j
       OPEN(12+j,file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',   &
            action='write'); WRITE(12+j,'(A)')                             &
       "Alya Ensight Gold --- Vector per-node variable file"
       WRITE(12+j,'(A/A/A)') "part", "     1","coordinates"
     END IF
     CALL gather(x1_pp(1:),eld_pp)
     disp_pp=zero
     CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,        &
                         node_start,node_end,eld_pp,disp_pp,1)
     DO i=1,ndim
      temp=zero

      DO l=1,nodes_pp
        k=i+(ndim*(l-1))
        temp(l)=disp_pp(k)
      END DO

      ! Send Displacement of node to Master
	    IF(numpe==outProc) THEN
	      m=12456-node_start
	      send=temp(m)
        CALL MPI_SEND(send,1,MPI_REAL8,0,0,MPI_COMM_WORLD,ier)
      END IF

      IF(numpe==1) THEN
	      ! NOTE: recieve processor entered manually
        CALL MPI_RECV(recv,1,MPI_REAL8,7,0,MPI_COMM_WORLD,statu,ier)
	      outDisp(i)=recv
      END IF

      CALL dismsh_ensi_p(12+j,l,nodes_pp,npes,numpe,1,temp)
     
     END DO; IF(numpe==1) CLOSE(12+j)
   END IF 
ENDIF

IF (numpe==1)THEN
  IF(j.EQ.1)THEN
    OPEN(10,FILE='Displacement.dat',STATUS='replace',ACTION='write')
  ENDIF
  WRITE(10,'(E12.4,E12.4,3E12.4)') real_time,tLoad,outDisp
ENDIF

END DO timesteps

! deallocate all arrays, moved from inside the loop
DEALLOCATE( p_pp )
deallocate( r_pp )
deallocate( x_pp )
deallocate( u_pp )
deallocate( d_pp )
deallocate( diag_precon_pp )
deallocate( store_km_pp )
deallocate( pmul_pp )
DEALLOCATE( xnew_pp )
DEALLOCATE( g_coord_pp )

!------------------------ write out displacements ----------------------

  CALL calc_nodes_pp(nn,npes,numpe,node_end,node_start,nodes_pp)

  !IF ( numpe==1 ) THEN
  !   write(ch,'(I6.6)') numpe
  !   open( 12, file=argv(1:nlen)//".ensi.DISPL-"//ch,status='replace',    &
  !        action='write')
  !  write(12,'(A)') "Alya Ensight Gold --- Vector per-node variable file"
  !   write(12,'(A/A/A)') "part", "     1","coordinates"
  !END IF

 ! ALLOCATE( disp_pp(nodes_pp*ndim), source=zero )
 ! allocate( temp(nodes_pp), source=zero )
  !CALL scatter_nodes( npes, nn, nels_pp, g_num_pp, nod, ndim, nodes_pp,  &
  !     node_start, node_end, eld_pp, disp_pp, 1 )
  !DO i=1,ndim
  !   temp=zero
  !   DO j=1,nodes_pp
  !      k = i+(ndim*(j-1))
  !      temp(j) = disp_pp(k)
  !   END DO
  !   CALL dismsh_ensi_p(12,1,nodes_pp,npes,numpe,1,temp)
  !END DO

  IF ( numpe==it ) then
     write( 11, '(A,F10.4)') "This analysis took: ", elap_time()-timest(1)
     close( 11 )
     close( 12 )
  end if
!*** end ParaFEM part ************************************************72


!*** CGPACK part *****************************************************72
! deallocate all CGPACK arrays
call cgca_ds(  cgca_space )
call cgca_drt( cgca_grt )
call cgca_pfem_sdalloc
call cgca_pfem_edalloc
call cgca_pfem_integdalloc
!*** end CGPACK part *************************************************72


!*** ParaFEM part ****************************************************72
!CALL SHUTDOWN() ! cannot call MPI_FINALIZE with coarrays with Intel 16.
END PROGRAM xx20std
