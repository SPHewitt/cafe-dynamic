PROGRAM xx20_noio
!-----------------------------------------------------------------------
! Sam Hewitt (University of Manchester)
! Anton Shterenlikht (University of Bristol)
! Lee Margetts (University of Manchester)
!
! Program xx20- linking ParaFEM with CASUP, specifically
! modifying p129 from 5th edition to link with the cgca modules.
!
!
! This version does not contains any serious I/O, it simply prints
! the beam tip displacements and the timing information.
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
INTEGER                 ::  nlen,nres,meshgen,partitioner,it,is
INTEGER                 ::  statu(MPI_STATUS_SIZE),tstep,ierr
INTEGER                 ::  ii,jj,kk,ll,mm,idx1,idx2,open_flag
INTEGER                 ::  I_NUMBER
LOGICAL                 ::  I_OPEN,I_EXIST

REAL(iwp),PARAMETER     ::  zero=0.0_iwp    
REAL(iwp)               ::  e,v,det,rho,alpha1,beta1,omega,theta
REAL(iwp)               ::  period,pi,dtim,volume,c1,c2,c3,c4,q
REAL(iwp)               ::  real_time,tol,up,alpha,beta,ld_scale
REAL(iwp)               ::  outDisp(3),send(3),recv(3),tLoad
REAL(iwp)               ::  scrit_scale,crack_start,exit_disp
REAL(iwp)               ::  exit_flag,dw,ca_res,ca_charlen
REAL(iwp)               ::  tot_fracvol,fe_fracvol
REAL(iwp)               ::  updatek,precon,builde,solve,getstress,cleave

CHARACTER(LEN=15)       ::  element
CHARACTER(LEN=50)       ::  argv
CHARACTER(LEN=6)        ::  ch
CHARACTER(LEN=50)       ::  fname,job_name,label

LOGICAL                 ::  consistent=.TRUE.,converged

!---------------------------- dynamic arrays -------------------------72
REAL(iwp),ALLOCATABLE ::  loads_pp(:),points(:,:),dee(:,:),fun(:)
REAL(iwp),ALLOCATABLE ::  jac(:,:),der(:,:),deriv(:,:),weights(:)
REAL(iwp),ALLOCATABLE ::  bee(:,:),g_coord_pp(:,:,:),fext_pp(:)
REAL(iwp),ALLOCATABLE ::  x1_pp(:),d1x1_pp(:),d2x1_pp(:),emm(:,:)
REAL(iwp),ALLOCATABLE ::  ecm(:,:),x0_pp(:),d1x0_pp(:),d2x0_pp(:)
REAL(iwp),ALLOCATABLE ::  store_km_pp(:,:,:),vu_pp(:),u_pp(:)
REAL(iwp),ALLOCATABLE ::  store_mm_pp(:,:,:),p_pp(:),d_pp(:),x_pp(:)
REAL(iwp),ALLOCATABLE ::  xnew_pp(:),pmul_pp(:,:),utemp_pp(:,:)
REAL(iwp),ALLOCATABLE ::  temp(:),diag_precon_pp(:),value_shape(:)
REAL(iwp),ALLOCATABLE ::  diag_precon_tmp(:,:),temp_pp(:,:,:)
REAL(iwp),ALLOCATABLE ::  disp_pp(:),val(:,:),timest(:),eld_pp(:,:)
REAL(iwp),ALLOCATABLE ::  r_pp(:),tot_r_pp(:),eps(:),sigma(:)

REAL(iwp),ALLOCATABLE :: stress_integral_pp(:,:),stressnodes_pp(:)
REAL(iwp),ALLOCATABLE :: shape_integral_pp(:,:)
!REAL(iwp),ALLOCATABLE :: principal_integral_pp(:,:),princinodes_pp(:)
!REAL(iwp),ALLOCATABLE :: principal(:),reacnodes_pp(:)
!REAL(iwp),ALLOCATABLE :: strain_integral_pp(:,:)

INTEGER,ALLOCATABLE   :: rest(:,:),g_num_pp(:,:),g_g_pp(:,:),node(:)    

!*** CGPACK part *****************************************************72

! CGPACK parameters

integer, parameter :: cgca_linum=5 ! number of loading iterations
logical( kind=ldef ), parameter :: cgca_yesdebug = .false.,            &
 cgca_nodebug = .false.
real( kind=rdef ), parameter :: cgca_zero = 0.0_rdef,                  &
 cgca_one = 1.0_rdef,                                                  &
 ! cleavage stress on 100, 110, 111 planes for BCC,
 ! see the manual for derivation, GPa.
 ! Material = Iron
 cgca_scrit(3) = (/ 1.05e4_rdef, 1.25e4_rdef, 4.90e4_rdef /)


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
  cgca_fracvol,           & ! volume (number) of fractured cells per img
  cgca_scrit_scale         ! Scale for the Critical Stress 
real( kind=rdef ), allocatable :: cgca_grt(:,:,:)[:,:,:]
logical( kind=ldef ) :: cgca_solid
! logical( kind=ldef ) :: cgca_lflag
character( len=6 ) :: cgca_citer
!*** end CGPACK part *************************************************72

updatek = 0.0
precon = 0.0
builde = 0.0
solve = 0.0
getstress = 0.0
cleave = 0.0
!------------------------ input and initialisation ---------------------
ALLOCATE(timest(20))
timest     =  zero
timest(1)  =  elap_time()

!*    Get the rank of the processes and the total number of processes
!* intent( out ):
!*          numpe - integer, process number (rank)
!*           npes - integer, total number of processes (size)
!CALL find_pe_procs( numpe, npes )
! CAnnot use MPI_INIT in a coarray program with ifort 16!
CALL find_pe_procs(numpe,npes)

CALL getname(argv,nlen)

! Input file .def replaces read_p129
CALL READ_DEF(argv,alpha1,beta1,e,element,limit,loaded_nodes,           &
              meshgen,nels,nip,nn,nod,npri,nr,nres,nstep,omega,         &
              partitioner,rho,theta,tol,v,dtim,ld_scale,scrit_scale,       &
              crack_start,exit_disp,ca_res,ca_charlen)

cgca_scrit_scale=scrit_scale

!CALL read_p129(argv,numpe,alpha1,beta1,e,element,limit,loaded_nodes,    &
!   meshgen,nels,nip,nn,nod,npri,nr,nres,nstep,omega,partitioner,rho,    &
!   theta,tol,v)

CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

CALL calc_nels_pp(argv,nels,npes,numpe,partitioner,nels_pp)

ndof  =  nod*nodof
ntot  =  ndof

ALLOCATE(g_coord_pp(nod,ndim,nels_pp))
ALLOCATE(g_num_pp(nod,nels_pp))
ALLOCATE(rest(nr,nodof+1))

g_num_pp    =  0
g_coord_pp  =  zero
rest        =  0
open_flag   =  0

CALL read_g_num_pp(argv,iel_start,nn,npes,numpe,g_num_pp)

IF(meshgen==2) CALL abaqus2sg(element,g_num_pp)
 
CALL read_g_coord_pp(argv,g_num_pp,nn,npes,numpe,g_coord_pp)
CALL read_rest(argv,numpe,rest)

ALLOCATE(points(nip,ndim),bee(nst,ntot),dee(nst,nst),jac(ndim,ndim),     &
   weights(nip),der(ndim,nod),deriv(ndim,nod),g_g_pp(ntot,nels_pp),      &
   store_km_pp(ntot,ntot,nels_pp),utemp_pp(ntot,nels_pp),ecm(ntot,ntot), &
   pmul_pp(ntot,nels_pp),store_mm_pp(ntot,ntot,nels_pp),emm(ntot,ntot),  &
   temp_pp(ntot,ntot,nels_pp),diag_precon_tmp(ntot,nels_pp),fun(nod),    &
   eps(nst),sigma(nst),value_shape(nod),shape_integral_pp(nod,nels_pp),  &
   stress_integral_pp(nod*nst,nels_pp))


  shape_integral_pp     = zero
  stress_integral_pp    = zero
  stressnodes_pp        = zero

timest(2)  =  elap_time()
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

! xx20 The PARAFEM Model is (0,0,-0.1) to (0.1,0.5,0)
cgca_bsz = (/ 0.004, 0.004, 0.004 /)

! xx20 Origin of box is (0,0,-0.1)
cgca_origin = (/ 0.096, 0.248, -0.1 /)

! Rotation tensor *from* FE cs *to* CA cs.
! The box cs is aligned with the box.
cgca_rot         = cgca_zero
cgca_rot( 1, 1 ) = 1.0
cgca_rot( 2, 2 ) = 1.0
cgca_rot( 3, 3 ) = 1.0

! mean grain size, also mm
! Beam is of size 0.1 * 0.5 * 0.1
! From Kato2015 - Iron has grain dimater 0.3 and 0.5 mm
 cgca_dm = 0.0005e0_rdef

! resolution, cells per grain
! Coarse resolution for testing ** INCREASE FOR PRODUCTION RUNS
!  cgca_res = 1.0e5_rdef
  cgca_res = ca_res

! cgpack length scale, also in mm
! Equivalent to crack propagation distance per unit of time,
! i.e. per second. Let's say 1 km/s = 1.0e3 m/s = 1.0e6 mm/s. 
 cgca_length = 1.0e3_rdef!1.0e6_rdef

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

! Allocate cgca_pfem_integrity%i(:), array component of a coarray of
! derived type. Allocating *local* array.
! i is set to 1.0 on allocation.
call cgca_pfem_integalloc( nels_pp )
 
! Allocate the Young's modulus 2D array
call cgca_pfem_ealloc( nip, nels_pp )
  
! initially set the Young's modulus to "e" everywhere
 cgca_pfem_enew = e

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


! ** HEWITT
if ( cgca_img .eq. 1 ) then
  ! Setting Crack in the center of the domain
  WRITE(*,'(A)')"Setting crack Nucleus"
  cgca_space( cgca_c(1)/2, cgca_c(2)/2, cgca_c(3)/2, cgca_state_type_frac )  &
       [ int(cgca_ir(1)/2), int(cgca_ir(2)/2), int(cgca_ir(3)/2) ] = cgca_clvg_state_100_edge
end if

sync all


! Must be slightly bigger than the characteristic length of FE.
! Coarse = 0.004  m (4mm)
! Medium = 0.002  m (2mm)
! fine   = 0.001  m (1mm)
! vfine  = 0.0005 m (0.5mm) Iron grain size
! ** NEEDS TO BE CHANGED DEPENDING ON MESH
! cgca_charlen = 0.005
 cgca_charlen = ca_charlen

 IF(numpe==1)THEN
   WRITE(*,'(A)')" "
   WRITE(*,"(2(A,F10.5))")"Resolution: ", cgca_res, " Charlen: ", cgca_charlen
 ENDIF
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

! Allocate the stress array component of cgca_pfem_stress coarray
! subroutine cgca_pfem_salloc( nels_pp, intp, comp )
call cgca_pfem_salloc( nels_pp, nip, nst )

! Need a sync after file write, because not sure what is coming next...
sync all
timest(3)  =  elap_time()
!*** end CGPACK part *************************************************72

!----------  find the steering array and equations per process ---------
CALL rearrange(rest)
g_g_pp  = 0
neq     = 0

! map equations to elements
elements_0: DO iel=1,nels_pp
  CALL find_g3( g_num_pp(:,iel), g_g_pp(:,iel), rest )
END DO elements_0

neq = MAXVAL(g_g_pp)
neq = max_p(neq)
CALL calc_neq_pp
CALL calc_npes_pp( npes, npes_pp )
CALL make_ggl( npes_pp, npes, g_g_pp )

! Set to the last equation (z)
nres=neq

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
   eld_pp(ntot,nels_pp),stressnodes_pp(nodes_pp*nst))

   x0_pp=zero; d1x0_pp=zero; x1_pp=zero; vu_pp=zero; diag_precon_pp=zero
   u_pp=zero; d2x0_pp=zero; loads_pp=zero; d1x1_pp=zero; d2x1_pp=zero
   d_pp=zero; p_pp=zero; x_pp=zero; xnew_pp=zero; fext_pp=zero 
   r_pp=zero; tot_r_pp=zero ; disp_pp=zero ; eld_pp=zero
   temp = 0

timest(4)  =  elap_time()
!---------------------------- Calculate Mass Matrix -------------------------
store_mm_pp=zero

IF(numpe==1)THEN
  WRITE(*,"(2(A,F10.5))")"Time Step: ", dtim, " Load scale: ", ld_scale 
ENDIF

! Omega is in Rad/s
pi=ACOS(-1._iwp); period=2._iwp*pi/omega;

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
 timest(5)  =  elap_time()

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


  ! Scaling Factor for Force
  tot_r_pp=tot_r_pp * ld_scale
  tLoad=SUM_P(tot_r_pp)

 timest(6)  =  elap_time()
IF(numpe==it) THEN
     OPEN(11,FILE=argv(1:nlen)//".res",STATUS='REPLACE',ACTION='WRITE')
       WRITE(11,'(A,I6,A)') "This job ran on ",npes," processes"
       WRITE(11,'(A,3(I12,A))') "There are ",nn," nodes ",nr,                &
             " restrained and ", neq," equations"
       WRITE(11,'(A)')" "
       WRITE(11,'(A,E12.4,2(A,F10.2))') "Youngs Modulus",e," Poissons Ratio",v,                &
                                        " Density", rho
       WRITE(11,'(A)')" "
       WRITE(11,'(4(A,F10.2))') "Alpha1",alpha1," Beta1",beta1,                &
                                " Omega", omega," Theta",theta
       WRITE(11,'(A)')" "
       WRITE(11,'(2(A,E12.4))') "Total Load: ",tLoad," Load factor: ",ld_scale
       WRITE(11,'(A)')" "
       WRITE(11,'(A)')"Intilisation Timings: "
       WRITE(11,'(A)') "---------------------"
       WRITE(11,'(A,F10.4)') "ParaFEM File I/O:           ", timest(2)-timest(1)
       WRITE(11,'(A,F10.4)') "CASUP Initialise:           ", timest(3)-timest(2)
       WRITE(11,'(A,F10.4)') "ParaFEM Make ggl:           ", timest(4)-timest(3)
       WRITE(11,'(A,F10.4)') "ParaFEM Create Mass Matrix: ", timest(5)-timest(4)
       WRITE(11,'(A,F10.4)') "ParaFEM Read loads:         ", timest(6)-timest(5)
       WRITE(11,'(A,F10.4)') "Total Setup Time:           ", elap_time()-timest(1)
       WRITE(11,'(A)')" "
       WRITE(11,'(A)') "  Time        Update K    Precon      Build Eqns  Solve       Iters Sigma       Cleavage  "
       WRITE(11,'(A)') "------------------------------------------------------------------------------------------"
  END IF 


!---------------------------- initial conditions -------------------------
 x0_pp=zero; d1x0_pp=zero; d2x0_pp=zero; real_time=zero

!-------------------------------------------------------------------------
!---------------------------- time stepping loop ------------------------- 
!-------------------------------------------------------------------------
 timest(7)  =  elap_time()

cgca_liter=1
timesteps: DO tstep=1,nstep
  cgca_liter = cgca_liter + 1  !counter for cafe output

 timest(8)  =  elap_time()
!------ element stiffness integration and build the preconditioner ---72
  dee = zero

  CALL sample( element, points, weights )
  store_km_pp = zero ;

  ! Calculate Updated Stiffness Matrix
  elements_2a: DO iel=1,nels_pp
    gauss_pts_1: DO i=1,nip
      CALL deemat( dee, cgca_pfem_enew(i,iel), v )
      CALL shape_der( der, points, i )
      jac = MATMUL( der, g_coord_pp(:,:,iel) )
      det = determinant(jac)
      CALL invert(jac)
      deriv = MATMUL( jac, der )
      CALL beemat( bee, deriv )
      store_km_pp(:,:,iel) = store_km_pp(:,:,iel) +                    &
      MATMUL( MATMUL( TRANSPOSE(bee), dee ), bee ) * det * weights(i)
    END DO gauss_pts_1
  END DO elements_2a

   timest(9)  =  elap_time()
   updatek = updatek + (timest(9)-timest(8))
  ! Calculate Updated Diagonal Preconditioner
  diag_precon_tmp = zero

  elements_2b: DO iel=1,nels_pp
    DO k=1,ntot
      diag_precon_tmp(k,iel) = diag_precon_tmp(k,iel)                  & 
               + store_mm_pp(k,k,iel)*c3+store_km_pp(k,k,iel)*c4
    END DO
  END DO elements_2b

  diag_precon_pp = zero
  CALL scatter( diag_precon_pp, diag_precon_tmp )

   timest(10)  =  elap_time()
   precon = precon + timest(10)-timest(9)
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
  CALL scatter(vu_pp,utemp_pp)

!----------------------------- RHS and loading --------------------------
  real_time=real_time+dtim;
  r_pp = zero
  r_pp = tot_r_pp*(theta*dtim*cos(omega*real_time)+c1*                 &
                   cos(omega*(real_time-dtim)))
  ! Impulse Force
  !IF(tstep .eq. 1)THEN
  !  r_pp = tot_r_pp
  !ELSE
  ! r_pp = zero
  !ENDIF

  IF ( loaded_nodes > 0 ) THEN
    q = SUM_P( r_pp )
    tLoad=q
    IF ( numpe==it ) then
      write (*,*)" "
      write (*, '(A,I12,A,F10.5)') "Time Step:", tstep," Real Time:",real_time
      write (*,'(A,E12.4)') "The total load is:", q
    END IF
  END IF

  r_pp = u_pp+vu_pp+r_pp
  diag_precon_pp = 1._iwp/diag_precon_pp
  d_pp = diag_precon_pp*r_pp
  p_pp = d_pp
  x_pp = zero

  timest(11)  =  elap_time()
  builde = builde + timest(11)- timest(10)
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

  timest(12)  =  elap_time()
  solve = solve + timest(12)-timest(11)

  IF ( numpe==it ) THEN
    write(*,'(A,I6)')"The number of iterations to convergence was ",iters
    write(*,'(A,F10.4)')"Time to solve equations was  :",                &
      timest(12)-timest(11)
    write(*,'(A,E12.4)')"The central nodal displacement is :",xnew_pp(is)
  END IF

!--------------- recover stresses at centroidal gauss point ------------
! IF CAFE
! Manual Input

IF((.TRUE.) .AND. (tstep .EQ. crack_start)) THEN

  timest(13)  =  elap_time()
  CALL gather(xnew_pp(1:),eld_pp)
  elmnts: DO iel = 1, nels_pp
    intpts: DO i = 1, nip
      call deemat(dee, cgca_pfem_enew( i, iel ), v)
      CALL shape_der(der,points,i)
      jac = MATMUL(der,g_coord_pp(:,:,iel))
      CALL invert(jac)
      deriv = MATMUL(jac,der)
      CALL beemat(bee,deriv)
      eps = MATMUL(bee,eld_pp(:,iel))
      sigma = MATMUL(dee,eps)
      !    write (*,*) "MPI rank", numpe, "el", iel,  &
      !                "int. point", i, "stress", sigma
      cgca_pfem_stress%stress( iel, i, : ) = sigma

    END DO intpts
  end do elmnts

  timest(14)  =  elap_time()
  getstress = getstress + timest(14)-timest(13)
!*** end ParaFEM part ************************************************72

!*** CGPACK part *****************************************************72

if ( cgca_img .eq. cgca_nimgs ) then
  do i = 1, nip
   write (*,"(2(a,i0),a,6es10.2)") "img: ", cgca_nimgs,               &
         " FE 1 int. p. ", i, " stress: ",                            &
          cgca_pfem_stress % stress( 1 , i , : )
  end do

  PRINT*,cgca_stress
end if



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

! ** CAUSE FOR DEBATE
! Time Increment for this problem needs researching
! 0.02 mm/increment approxmiatley (1/5) of the CA domain
 cgca_time_inc = dtim * 1.0_rdef !*(1.0_rdef / cgca_length)

! run cleavage for a correct number of iterations, which is a function
! of the characteristic length and the time increment
 cgca_clvg_iter = nint( cgca_length * cgca_lres * cgca_time_inc )

if (cgca_clvg_iter .LE. 0) cgca_clvg_iter=1

!cgca_clvg_iter = 200
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
     cgca_scrit_scale * cgca_scrit, cgca_clvgsp, cgca_gcupdn,          &
     .false., cgca_clvg_iter, 100, cgca_yesdebug )

sync all
     
! Calculate number (volume) of fractured cells on each image.
! cgca_fracvol is a local, non-coarray, array, so no sync needed.
! Note : Cells of states "cgca_frac_states" are considered failed.
! This includes grain boundaries (1, -1 ... -6)
call cgca_fv( cgca_space, cgca_fracvol )

fe_fracvol = cgca_fracvol
tot_fracvol = 0

CALL MPI_REDUCE(cgca_fracvol,tot_fracvol,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

if (numpe .eq. 1 ) PRINT*,"Total Fracture Volume", tot_fracvol

! calculate integrity, remote write
call cgca_pfem_intcalc1( cgca_c, cgca_fracvol )

! dump min integrity for all FE stored on this image
!write (*,*) "img:", cgca_img, "min. integrity:",                       &
!  minval( cgca_pfem_integrity % i )

! wait for integrity on all images to be calculated
sync all

! Young's modulus need to be updated on each image, local arrays only.
! The original Young's modulus value must be given as input.
! For each FE stored on this image.
call cgca_pfem_uym( e, nels_pp )

! dump updated Young's modulus
!write (*,*) "img:", cgca_img, "*min* Young's mod:",                    &
!            minval( cgca_pfem_enew )
sync all


  timest(15)  =  elap_time()
  cleave = cleave  + timest(15) - timest(14)
endif ! CAFE
!*** end CGPACK part *************************************************72

!*** ParaFEM part ****************************************************72
IF((tstep .EQ. crack_start))THEN
  IF(numpe==it) THEN
    WRITE(11,'(5E12.4,I6,2E12.4)') real_time,timest(9)-timest(8), &
    timest(10)-timest(9),timest(11)-timest(10),timest(12)-timest(11), &
    iters, timest(14)-timest(13),timest(15)-timest(14)
  ENDIF
ELSE
  IF(numpe==it) THEN
    WRITE(11,'(5E12.4,I6)') real_time,timest(9)-timest(8), &
    timest(10)-timest(9),timest(11)-timest(10),timest(12)-timest(11), &
    iters
  ENDIF
ENDIF

! Send Displacement of node to Master
IF(numpe==it) THEN
  send=x1_pp(is-2:is)
  CALL MPI_SEND(send,3,MPI_REAL8,0,0,MPI_COMM_WORLD,ier)
END IF

IF(numpe==1) THEN
  ! NOTE: recieve processor entered manually
  CALL MPI_RECV(recv,3,MPI_REAL8,npes-1,0,MPI_COMM_WORLD,statu,ier)
  outDisp=recv
END IF

IF (numpe==1)THEN
  IF(tstep.EQ.1)THEN
    OPEN(10,FILE='Displacement.dat',STATUS='replace',ACTION='write')
  ENDIF
  WRITE(10,'(E12.4,E12.4,3E12.4,E12.4)') real_time,tLoad,outDisp,tot_fracvol
ENDIF

exit_flag = ABS(outDisp(3))
CALL MPI_BCAST(exit_flag,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

! Check for program exit
! Note outDisp(3) -> Z
IF( exit_flag .GE. exit_disp) THEN
  IF(numpe == 1)THEN    
    WRITE(*,*) "Displacement Exceeded specified Value"
    WRITE(*,*) "Ending Time Loop"
  ENDIF
  EXIT
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

DEALLOCATE(stress_integral_pp)
DEALLOCATE(stressnodes_pp)                        
!DEALLOCATE(principal_integral_pp)                
!DEALLOCATE(princinodes_pp)
!DEALLOCATE(reacnodes_pp)
DEALLOCATE(shape_integral_pp)

  IF ( numpe==it ) then
     WRITE(11, '(A)') " "
     WRITE(11,'(A)')"Summary of Timings:   "
     WRITE(11,'(A)') "--------------------------------------"
     WRITE(11,'(A,F10.4)') "Updating stiffness matrix: ", updatek
     WRITE(11,'(A,F10.4)') "Building preconditioner:   ", precon
     WRITE(11,'(A,F10.4)') "Building eqns:             ", builde
     WRITE(11,'(A,F10.4)') "Solving Eqns:              ", solve
     WRITE(11,'(A,F10.4)') "Getting Stress tensor:     ", getstress
     WRITE(11,'(A,F10.4)') "Cleavage:                  ", cleave
     WRITE(11,'(A)') " "
     WRITE(11,'(A,F10.4)') "Total analysis took:       ", elap_time()-timest(1)
     CLOSE(11 )
     CLOSE(12 )
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
END PROGRAM xx20_noio


SUBROUTINE READ_DEF(job_name,alpha1,beta1,e,element,           &
                    limit,loaded_nodes,mesh,nels,nip,nn,nod,   &
                    npri,nr,nres,nstep,omega,partitioner,rho,  &
                    theta,tol,v,dtim,load_scale,scrit_scale,   & 
                    crack_start,exit_disp,resolution,charlen)

!USE mpi_wrapper
USE precision
USE mp_interface

IMPLICIT none

CHARACTER(LEN=50), INTENT(IN)    :: job_name
CHARACTER(LEN=15), INTENT(INOUT) :: element
CHARACTER(LEN=50)                :: fname,nullname

REAL(iwp), INTENT(INOUT)         :: rho,e,v,alpha1,beta1,theta,omega
REAL(iwp), INTENT(INOUT)         :: tol,dtim,load_scale,scrit_scale
REAL(iwp), INTENT(INOUT)         :: crack_start,exit_disp
REAL(iwp), INTENT(INOUT)         :: resolution,charlen

INTEGER, INTENT(INOUT)           :: nels,nn,nr,nres,nod,nip,loaded_nodes
INTEGER, INTENT(INOUT)           :: nstep,npri,limit,mesh,partitioner 

!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------

INTEGER                          :: integer_store(12)
REAL(iwp)                        :: real_store(15)

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------

IF (numpe==1) THEN
  
  PRINT*,"Reading Definitions File"

  fname = job_name(1:INDEX(job_name, " ") -1)//".def"
  OPEN(8,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(8,*)
    READ(8,*)
    READ(8,*) element, nullname
    READ(8,*) mesh, nullname
    READ(8,*) partitioner, nullname
    READ(8,*) nels, nullname
    READ(8,*) nn, nullname
    READ(8,*) nr, nullname
    READ(8,*) nip, nullname
    READ(8,*) nod, nullname
    READ(8,*) loaded_nodes, nullname
    READ(8,*) nres, nullname
    READ(8,*)
    READ(8,*) rho, nullname
    READ(8,*) e, nullname
    READ(8,*) v, nullname
    READ(8,*) alpha1, nullname
    READ(8,*) beta1, nullname
    READ(8,*) nstep, nullname
    READ(8,*) npri, nullname
    READ(8,*) theta, nullname
    READ(8,*) omega, nullname
    READ(8,*) tol, nullname
    READ(8,*) limit, nullname
    READ(8,*) dtim, nullname
    READ(8,*)
    READ(8,*) charlen, nullname
    READ(8,*) resolution, nullname
    READ(8,*) load_scale, nullname
    READ(8,*) scrit_scale, nullname
    READ(8,*) crack_start, nullname
    READ(8,*) exit_disp, nullname
  CLOSE(8)

  integer_store      = 0

  integer_store(1)   = mesh
  integer_store(2)   = nels
  integer_store(3)   = nn
  integer_store(4)   = nr 
  integer_store(5)   = nip
  integer_store(6)   = nod
  integer_store(7)   = loaded_nodes
  integer_store(8)   = nstep
  integer_store(9)   = npri
  integer_store(10)  = limit
  integer_store(11)  = partitioner
  integer_store(12)  = nres

  real_store         = 0.0_iwp

  real_store(1)      = rho  
  real_store(2)      = e  
  real_store(3)      = v  
  real_store(4)      = alpha1  
  real_store(5)      = beta1  
  real_store(6)      = theta  
  real_store(7)      = omega  
  real_store(8)      = tol 
  real_store(9)      = dtim  
  real_store(10)     = load_scale
  real_store(11)     = scrit_scale
  real_store(12)     = crack_start
  real_store(13)     = exit_disp
  real_store(14)     = charlen
  real_store(15)     = resolution
  

ENDIF
!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------
bufsize = 12
CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

bufsize = 15
CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)

bufsize = 15
CALL MPI_BCAST(element,bufsize,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------

IF (numpe/=1) THEN

  mesh         = integer_store(1)
  nels         = integer_store(2)
  nn           = integer_store(3)
  nr           = integer_store(4)
  nip          = integer_store(5)
  nod          = integer_store(6)
  loaded_nodes = integer_store(7)
  nstep        = integer_store(8)
  npri         = integer_store(9)
  limit        = integer_store(10)
  partitioner  = integer_store(11)
  nres         = integer_store(12)

  rho          = real_store(1)
  e            = real_store(2)
  v            = real_store(3)
  alpha1       = real_store(4)
  beta1        = real_store(5)
  theta        = real_store(6)
  omega        = real_store(7)
  tol          = real_store(8)
  dtim         = real_store(9)
  load_scale   = real_store(10)
  scrit_scale  = real_store(11)
  crack_start  = real_store(12)
  exit_disp    = real_store(13)
  charlen      = real_store(14)
  resolution   = real_store(15)

END IF

RETURN

ENDSUBROUTINE READ_DEF

