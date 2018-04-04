IF(.NOT. ALLOCATED(shape_integral_pp))THEN
  ALLOCATE(shape_integral_pp(nod,nels_pp))
  ALLOCATE(stress_integral_pp(nod*nst,nels_pp))
  ALLOCATE(stressnodes_pp(nodes_pp*nst))
  !ALLOCATE(principal_integral_pp(nod*nodof,nels_pp))
  !ALLOCATE(princinodes_pp(nodes_pp*nodof))
  !ALLOCATE(reacnodes_pp(nodes_pp*nodof))
 !ALLOCATE(principal(ndim))
  !ALLOCATE(strain_integral_pp(nod*nst,nels_pp))
ENDIF

!------------------------------------------------------------------------------
! 15. Print out displacements, stress, principal stress and reactions
!------------------------------------------------------------------------------
  job_name = argv  

  WRITE(ch,'(I6.6)') tstep

  IF(numpe==1 .AND. open_flag==0) THEN
    open_flag=1
    fname = job_name(1:INDEX(job_name, " ")-1) //".dis"
    OPEN(24, file=fname, status='new', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) //".str"
    OPEN(25, file=fname, status='new', action='write')
    !fname = job_name(1:INDEX(job_name, " ")-1) //".pri"
    !OPEN(26, file=fname, status='UNKNOWN', action='write')
    !fname = job_name(1:INDEX(job_name, " ")-1) //".vms"
    !OPEN(27, file=fname, status='UNKNOWN', action='write')
    !fname = job_name(1:INDEX(job_name, " ")-1) //".rea"
    !OPEN(28, file=fname, status='UNKNOWN', action='write')
    !fname = job_name(1:INDEX(job_name, " ")-1) //".est"
    !OPEN(29, file=fname, status='UNKNOWN', action='write')
  END IF

 IF(numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ")-1) //".dis"
    OPEN(24, file=fname, status='old', action='write',position='append')
    fname = job_name(1:INDEX(job_name, " ")-1) //".str"
    OPEN(25, file=fname, status='old', action='write',position='append')
    !fname = job_name(1:INDEX(job_name, " ")-1) //".pri"
    !OPEN(26, file=fname, status='UNKNOWN', action='write')
    !fname = job_name(1:INDEX(job_name, " ")-1) //".vms"
    !OPEN(27, file=fname, status='UNKNOWN', action='write')
    !fname = job_name(1:INDEX(job_name, " ")-1) //".rea"
    !OPEN(28, file=fname, status='UNKNOWN', action='write')
    !fname = job_name(1:INDEX(job_name, " ")-1) //".est"
    !OPEN(29, file=fname, status='UNKNOWN', action='write')
  END IF



!------------------------------------------------------------------------------
! 16a. Displacements
!------------------------------------------------------------------------------
  eld_pp = zero

  CALL gather(xnew_pp(1:),eld_pp)

  disp_pp = zero

  label   = "*DISPLACEMENT"

  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,ndim,nodes_pp,              &
                     node_start,node_end,eld_pp,disp_pp,1)
  CALL write_nodal_variable(label,24,tstep,nodes_pp,npes,numpe,ndim,disp_pp)


  IF(numpe==1) CLOSE(24)

!------------------------------------------------------------------------------
! 16b. Stresses
!------------------------------------------------------------------------------ 
  !principal             = zero
  shape_integral_pp     = zero
  stress_integral_pp    = zero
  stressnodes_pp        = zero
  !principal_integral_pp = zero  
  !princinodes_pp        = zero
  !reacnodes_pp          = zero
  !utemp_pp              = zero

CALL sample(element,points,weights)

  DO iel = 1,nels_pp

    DO ii = 1,nip

      CALL deemat(dee, cgca_pfem_enew( ii, iel ), v)
      CALL shape_der(der,points,ii)
      jac   = MATMUL(der,g_coord_pp(:,:,iel))
      det   = DETERMINANT(jac) 
      CALL invert(jac)
      deriv = MATMUL(jac,der)
      CALL beemat(deriv,bee)
      eps   = MATMUL(bee,eld_pp(:,iel))
      sigma = MATMUL(dee,eps)
      !CALL PRINCIPALSTRESS3D(sigma,principal)
      !utemp_pp(:,iel) = utemp_pp(:,iel) +                                    &
                        !MATMUL(TRANSPOSE(bee),sigma)*det*weights(ii)

      CALL shape_fun(fun,points,ii)
      dw = det * weights(ii)
      shape_integral_pp(:,iel) = shape_integral_pp(:,iel) + fun(:) * dw

      DO jj = 1,nod
        idx1 = (jj-1)*nst
        idx2 = (jj-1)*nodof
        DO kk = 1,nst
          stress_integral_pp(idx1+kk,iel) = stress_integral_pp(idx1+kk,iel) +   &
                                           fun(jj)*sigma(kk)*det*weights(ii)

          !strain_integral_pp(idx1+j,iel) = strain_integral_pp(idx1+j,iel) +&
          !              fun(jj)*eps(kk)*det*weights(ii)
        END DO !nst
        !DO kk = 1,nodof
        !  principal_integral_pp(idx2+kk,iel) = principal_integral_pp(idx2+kk,iel) + &
        !                                      fun(jj)*principal(kk)*det*weights(ii)
        !END DO !nodof
      END DO !node

    END DO !gauss

  END DO !elements

!------------------------------------------------------------------------------
! 16c. Stress
!------------------------------------------------------------------------------

  label = "*STRESS"
  CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nst,nodes_pp,            &
                        node_start,node_end,shape_integral_pp,                &
                        stress_integral_pp,stressnodes_pp)

  CALL write_nodal_variable(label,25,tstep,nodes_pp,npes,numpe,nst,               &
                            stressnodes_pp)
                            
 IF(numpe==1) CLOSE(25)

!------------------------------------------------------------------------------
! 16d. Principal stress
!------------------------------------------------------------------------------

!  label = "*PRINCIPAL STRESS"
  
!  CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,          &
!                        node_start,node_end,shape_integral_pp,                &
!                        principal_integral_pp,princinodes_pp)

!  CALL write_nodal_variable(label,26,tstep,nodes_pp,npes,numpe,nodof,             &
!                            princinodes_pp)

  !IF(numpe==1) CLOSE(26)

!------------------------------------------------------------------------------
! 16e. Von Mises stress (rho_v)
!      rho_v = sqrt( ( (rho1-rho2)^2 + (rho2-rho3)^2 + (rho1-rho3)^2 ) / 2 )
!------------------------------------------------------------------------------
  
!  label = "*MISES STRESS"
  
!  DO ii = 1,nodes_pp
!    jj = ((ii-1)*nodof)+1
!    kk = jj + 1
!    ll = jj + 2
!    princinodes_pp(j) = SQRT(((princinodes_pp(jj)-princinodes_pp(kk)) **2 +     &
!                              (princinodes_pp(kk)-princinodes_pp(ll)) **2 +     &
!                              (princinodes_pp(jj)-princinodes_pp(ll)) **2)      &
!                              * 0.5_iwp)
!    princinodes_pp(kk:ll) = zero
!  END DO

!  CALL write_nodal_variable(label,27,tstep,nodes_pp,npes,numpe,nodof,             &
!                            princinodes_pp)

  !IF(numpe==1) CLOSE(27)

!------------------------------------------------------------------------------
! 16f. Reactions
!------------------------------------------------------------------------------

!  label = "*NODAL REACTIONS"
!  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,             &
!                     node_start,node_end,utemp_pp,reacnodes_pp,0)
!  CALL write_nodal_variable(label,28,tstep,nodes_pp,npes,numpe,nodof,             &
!                            reacnodes_pp)

  !IF(numpe==1) CLOSE(28)

!------------------------------------------------------------------------------
! 16c. Strains
!------------------------------------------------------------------------------
  ! Reuse stressnodes_pp for strains
  !stressnodes_pp = zero

  !label = "*STRAIN"
  !CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nst,nodes_pp,            &
  !                      node_start,node_end,shape_integral_pp,                &
  !                      strain_integral_pp,stressnodes_pp)

!  CALL write_nodal_variable(label,29,tstep,nodes_pp,npes,numpe,nst,               &
!                            stressnodes_pp)
                            
 !IF(numpe==1) CLOSE(29)




