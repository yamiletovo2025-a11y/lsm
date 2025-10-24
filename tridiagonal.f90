      SUBROUTINE tridiagonal(ist,ifn,a1,a2,a3,b,delp,ierror)

      INTEGER MXD, ist, ifn, i, k
      INTEGER :: ierror
      PARAMETER (MXD = 50)
      REAL*8  a1(MXD),a2(MXD),a3(MXD),b(MXD)
      REAL*8  c1(MXD),d1(MXD)
      REAL*8  delp(MXD)
      REAL*8  temp
!
      a1(ist) = 0.d0
      a3(ifn) = 0.d0
      delp(:) = 0.d0

      c1(ist)=a3(ist)/a2(ist)
      d1(ist)=b(ist)/a2(ist)

      IF ( ierror == 1 ) THEN
        write(65, *) "=== tria ==="
        write(65,'(a5,4a17)')"i","temp","a2(i+1)","a1(i+1)","c1(i)"
      ENDIF

      DO i=ist,ifn-1
         temp=a2(i+1)-a1(i+1)*c1(i)

         IF ( ierror == 1 ) THEN
           write(65,'(i3,2x,4(d15.5,2x))')i,temp,a2(i+1),a1(i+1),c1(i)
         ENDIF
         if(temp.eq.0.) then
           stop 'tridiagonal: error in thomas solution'
         endif

         c1(i+1)=a3(i+1)/temp
         d1(i+1)=(b(i+1)-a1(i+1)*d1(i))/temp
      ENDDO
      delp(ifn)=d1(ifn)
     
      IF ( ierror == 1 ) THEN
        write(65,'(a5,5a17)')"i","delp(i)","d1(i)","c1(i)*delp(i+1)","c1(i)","delp(i+1)"
        write(65, *) "check the last line"
      ENDIF
 
      DO i=ifn-1,ist,-1
        delp(i)=d1(i)-c1(i)*delp(i+1)
        IF ( ierror == 1 ) THEN 
        write(65,'(i3,2x,5(d15.5,2x))')i,delp(i),d1(i),c1(i)*delp(i+1),c1(i),delp(i+1) 
        ENDIF
      ENDDO
        IF ( ierror == 1 ) THEN
        write(65,*) "delp(1) ~= d1(1)=b(1)/a2(1)", ist, d1(ist), b(ist), a2(ist)
        write(66,*) "delp(1) ~= d1(1)=b(1)/a2(1)", ist, d1(ist), b(ist), a2(ist)
        ENDIF 
      END SUBROUTINE tridiagonal
