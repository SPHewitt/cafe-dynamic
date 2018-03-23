!------------------------------------------------------------------------------
! 15. Print out displacements, stress, principal stress and reactions
!------------------------------------------------------------------------------
  job_name = argv  

  IF(numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ")-1)//".dis"
    OPEN(100+tstep, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".str"
    OPEN(200+tstep, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".pri"
    OPEN(300+tstep, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".vms"
    OPEN(400+tstep, file=fname, status='replace', action='write')
    fname = job_name(1:INDEX(job_name, " ")-1) // ".rea"
    OPEN(500+tstep, file=fname, status='replace', action='write')
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
  CALL write_nodal_variable(label,100+tstep,tstep,nodes_pp,npes,numpe,ndim,disp_pp)


  IF(numpe==1) CLOSE(100+tstep)

!------------------------------------------------------------------------------
! 16b. Stresses
!------------------------------------------------------------------------------ 
  principal             = zero
  shape_integral_pp     = zero
  stress_integral_pp    = zero
  stressnodes_pp        = zero
  principal_integral_pp = zero  
  princinodes_pp        = zero
  reacnodes_pp          = zero
  utemp_pp              = zero

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
      CALL PRINCIPALSTRESS3D(sigma,principal)
      utemp_pp(:,iel) = utemp_pp(:,iel) +                                    &
                        MATMUL(TRANSPOSE(bee),sigma)*det*weights(ii)

      CALL shape_fun(fun,points,ii)

      DO jj = 1,nod
        idx1 = (jj-1)*nst
        idx2 = (jj-1)*nodof
        shape_integral_pp(jj,iel) = shape_integral_pp(jj,iel) +                 &
                                   fun(jj)*det*weights(ii)
        DO kk = 1,nst
          stress_integral_pp(idx1+kk,iel) = stress_integral_pp(idx1+kk,iel) +   &
                                           fun(jj)*sigma(kk)*det*weights(ii)
        END DO
        DO kk = 1,nodof
          principal_integral_pp(idx2+kk,iel) = principal_integral_pp(idx2+kk,iel) + &
                                              fun(jj)*principal(kk)*det*weights(ii)
        END DO

      END DO !node

    END DO !gauss

  END DO !elements
!------------------------------------------------------------------------------
! 16c. Stress
!------------------------------------------------------------------------------
  PRINT*,"Writing Stress"

  label = "*STRESS"
  CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nst,nodes_pp,            &
                        node_start,node_end,shape_integral_pp,                &
                        stress_integral_pp,stressnodes_pp)

  PRINT*,"Nodal Projection"
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  CALL write_nodal_variable(label,200+tstep,tstep,nodes_pp,npes,numpe,nst,               &
                            stressnodes_pp)
                            

  IF(numpe==1) CLOSE(200+tstep)

!------------------------------------------------------------------------------
! 16d. Principal stress
!------------------------------------------------------------------------------
  IF(numpe==1)PRINT*,"Principal Stress"
  label = "*PRINCIPAL STRESS"
  
  CALL nodal_projection(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,          &
                        node_start,node_end,shape_integral_pp,                &
                        principal_integral_pp,princinodes_pp)

  CALL write_nodal_variable(label,300+tstep,tstep,nodes_pp,npes,numpe,nodof,             &
                            princinodes_pp)

  IF(numpe==1) CLOSE(300+tstep)

!------------------------------------------------------------------------------
! 16e. Von Mises stress (rho_v)
!      rho_v = sqrt( ( (rho1-rho2)^2 + (rho2-rho3)^2 + (rho1-rho3)^2 ) / 2 )
!------------------------------------------------------------------------------
  
  label = "*MISES STRESS"
  
  DO ii = 1,nodes_pp
    jj = ((ii-1)*nodof)+1
    kk = jj + 1
    ll = jj + 2
    princinodes_pp(j) = SQRT(((princinodes_pp(jj)-princinodes_pp(kk)) **2 +     &
                              (princinodes_pp(kk)-princinodes_pp(ll)) **2 +     &
                              (princinodes_pp(jj)-princinodes_pp(ll)) **2)      &
                              * 0.5_iwp)
    princinodes_pp(kk:ll) = zero
  END DO

  CALL write_nodal_variable(label,400+tstep,tstep,nodes_pp,npes,numpe,nodof,             &
                            princinodes_pp)

  IF(numpe==1) CLOSE(400+tstep)

!------------------------------------------------------------------------------
! 16f. Reactions
!------------------------------------------------------------------------------

  label = "*NODAL REACTIONS"
  CALL scatter_nodes(npes,nn,nels_pp,g_num_pp,nod,nodof,nodes_pp,             &
                     node_start,node_end,utemp_pp,reacnodes_pp,0)
  CALL write_nodal_variable(label,500+tstep,tstep,nodes_pp,npes,numpe,nodof,             &
                            reacnodes_pp)

  IF(numpe==1) CLOSE(500+tstep)


