IF(numpe==it) THEN
    WRITE(11,'(3E12.4,I10)') real_time,cos(omega*real_time),x1_pp(is),iters
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
    WRITE(10,'(E12.4,E12.4,3E12.4)') real_time,tLoad,outDisp
  ENDIF
