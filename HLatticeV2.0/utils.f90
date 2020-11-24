Module Utils
  use Parameters
  implicit none

  Type file_pointer
     integer(IB) unit
     integer(IB) pt
  End type file_pointer

  Interface Num2str
     Module Procedure Int2Str,float2str
  End Interface


  Interface Dverk
     Module Procedure Dverk_S,Dverk_V,Dverk_Easy
  End Interface

  Interface SPlints
     Module Procedure Splint_S,Splint_V
  End Interface
  
  Interface Write_Binary_File
     Module Procedure Int_to_Binary_File, float_To_Binary_File
  End Interface

  Interface Read_Binary_File
     Module Procedure Int_From_Binary_File, float_From_Binary_File
  ENd Interface


contains

  SUBROUTINE FINDGEN(X,XSTART,XEND)	  
    REAL(DL),DIMENSION(:),INTENT(OUT)::X
    REAL(DL) XSTART,XEND
    INTEGER(IB) I,N
    N=SIZE(X)
    X(1)=XSTART
    X(N)=XEND
    DO I=2,N-1
       X(I)=XSTART+(I-1)*(XEND-XSTART)/(N-1)
    ENDDO
  END SUBROUTINE FINDGEN

  SUBROUTINE ReportErr(SubRtName,Action)
    CHARACTER(LEN=*)SubRtName
    LOGICAL,OPTIONAL::ACTION
    WRITE(*,*) "Execution Error in "//TRIM(SUBRTNAME)
    IF(PRESENT(ACTION))THEN
       IF(ACTION)THEN
          stop
       ELSE
          return
       ENDIF
    ELSE
       stop
    ENDIF
  END subroutine ReportErr

  FUNCTION GETDIM(SUBRTNAME,N1,N2,N3,N4,N5,N6)
    CHARACTER(LEN=*) SUBRTNAME
    INTEGER(IB),INTENT(IN)::N1,N2
    INTEGER(IB),OPTIONAL::N3,N4,N5,N6
    INTEGER(IB) GETDIM
    GETDIM=N1
    IF(N1.NE.N2)THEN
       CALL REPORTERR(SUBRTNAME,.TRUE.)
       return
    ELSE
       IF(PRESENT(N3))THEN
          IF(N1.NE.N3)THEN
             CALL REPORTERR(SUBRTNAME,.TRUE.)
             return
          Endif
          IF(PRESENT(N4))THEN
             IF(N1.NE.N4)then
                CALL REPORTERR(SUBRTNAME,.TRUE.)
                return
             endif
             if(Present(N5))then
                if(N5.ne.N1)then
                   call reportErr(SUBRTName,.true.)
                   return
                endif
                if(present(N6))then
                   if(N1.ne.N6)then
                      call reportErr(SUBRtName,.true.)
                      return
                   endif
                ENDIF
             ENDIF
          ENDIF
       endif
    endif
  END FUNCTION GETDIM


#ifdef __GFORTRAN__
  function iargc ()
    integer(IB) iargc
    iargc=command_argument_count()
  end function iargc
  
  subroutine getarg(num, res)
    integer(IB), intent(in) :: num
    character(len=*), intent(out) :: res
    integer(IB) l, err
    call get_command_argument(num,res,l,err)
  end subroutine getarg
#endif

  function GetParam(i)
    character(LEN=512) GetParam
    integer(IB), intent(in) :: i
    if (iargc() < i) then
       GetParam = ''
    else
       call getarg(i,GetParam)
    end if
  end function GetParam

  FUNCTION INT2STR(I)
    INTEGER(IB),INTENT(IN)::I
    CHARACTER(LEN=16) INT2STR
    WRITE (INT2STR,*) I
    INT2STR= Trim(ADJUSTL(INT2STR))
  END FUNCTION INT2STR

  Function Float2str(x,fmt) !!This is a smater algorithm that convert number to string
    Character(LEN=128) Float2str, str
    real(dl) x, ax
    integer(IB) ind, i, n
    character(LEN=*),optional::fmt
    if(present(fmt))then
       if(trim(fmt).ne."f" .and. trim(fmt).ne."")then
          write(Float2str,"("//trim(fmt)//")") x
          float2str=trim(adjustl(float2str))
          return
       endif
    endif
    ax = abs(x)
    if(ax.eq.0.)then
       Float2str="0"
       return
    endif
    ind=0
    do while(ax.lt.1.d0)
       ind = ind -1
       ax = ax*10.d0
    enddo
    do while(ax.ge.10.)
       ind = ind + 1
       ax = ax/10.d0
    enddo
    if(ind.gt.5 .or. ind.lt.-5)then
       write(Str,*) nint(ax*1.d5)
       Str=trim(adjustl(Str))
       n=verify(trim(str), "0", .true.)
       if(n.gt.1)then
          Float2str=str(1:1)//"."//str(2:n)
       else
          Float2str=Str(1:1)
       endif
       write(str,*) ind
       Float2str=Trim(Float2str)//"e"//Trim(adjustl(Str))
    else
       ax=abs(x)
       write(Str,'(F16.8)') ax
       str=trim(adjustl(str))
       n=scan(str,".")
       if(n>0)then
          n=max(n-1,verify(trim(str),"0",.true.))
       else
          n=Len_Trim(str)
       endif
       if(Str(n:n).eq.'.')then
          Float2str =  Str(1:n-1)
       else
          Float2str = Str(1:n)
       endif
    endif
    if(x.lt.0.d0) Float2str='-'//trim(Float2str)

  End Function Float2str



  SUBROUTINE LinearInterpolate(X,Y,xs,ys) 
    !!before calling this routine, the array x must be sorted, otherwise unexpected error will occur
    Real(dl),dimension(:),intent(in)::X,Y
    REAL(DL) XS,YS
    INTEGER(IB) N,J,L,R
    INTEGER(IB),SAVE::PRER=100000
    REAL(DL) A
    N=size(x)
    if(size(y).ne.N) stop "LinearInterpolate: sizes of array do not agree"
    if(x(1).lt.x(N))then
       IF(XS.LE.X(1))THEN
          YS=Y(1)
          RETURN
       ELSEIF(XS.GE.X(N))THEN
          YS=Y(N)
          RETURN
       ENDIF
       IF(PRER.LE.N)THEN
          IF(X(PRER).GE.XS .AND. X(PRER-1).LE.XS)THEN
             R=PRER
             GOTO 100
          ENDIF
       ENDIF
       L=1
       R=N
       DO WHILE(R-L.GT.1)
          J=ISHFT(R+L,-1)
          IF(X(J).GT.XS)THEN
             R=J
          ELSE
             L=J
          ENDIF
       ENDDO
    else
       IF(XS.GE.X(1))THEN
          YS=Y(1)
          RETURN
       ELSEIF(XS.LE.X(N))THEN
          YS=Y(N)
          RETURN
       ENDIF
       IF(PRER.LE.N)THEN
          IF(X(PRER).GE.XS .AND. X(PRER-1).LE.XS)THEN
             R=PRER
             GOTO 100
          ENDIF
       ENDIF
       L=1
       R=N
       DO WHILE(R-L.GT.1)
          J=ISHFT(R+L,-1)
          IF(X(J).LT.XS)THEN
             R=J
          ELSE
             L=J
          ENDIF
       ENDDO
    endif
    PRER=R
100 A=(xs-x(R-1))/(x(R)-x(R-1))
    ys=y(R)*A+Y(R-1)*(1._dl -A)
  end SUBROUTINE LinearInterpolate



  SUBROUTINE DVERK_Easy (FCN,X,Y,xend,TOL)
    REAL(DL)::X,xend
    REAL(DL),DIMENSION(:),INTENT(INOUT)::Y
    REAL(DL),OPTIONAL::TOL
    EXTERNAL FCN
    REAL(DL) W(SIZE(Y),9),TOLL,C(24)
    INTEGER(IB) IND
    IF(PRESENT(TOL))THEN
       TOLL=TOL
    ELSE
       TOLL=1.D-6
    ENDIF
    C=0._dl
    w=0._dl
    IND=1
    call Dverk_s(fcn,x,y,xend,ind,c,w,toll)
    IF(IND.NE.3) stop "Dverk Error"
  END	SUBROUTINE DVERK_Easy


  SUBROUTINE DVERK_V (FCN,X,Y,TOL)
    REAL(DL),DIMENSION(:),INTENT(IN)::X
    REAL(DL),DIMENSION(:,:),INTENT(INOUT)::Y
    REAL(DL),OPTIONAL::TOL
    EXTERNAL FCN
    REAL(DL) W(SIZE(Y,1),9),XSTART,TOLL,C(24)
    INTEGER(IB) N,M,I,IND
    IF(PRESENT(TOL))THEN
       TOLL=TOL
    ELSE
       TOLL=1.D-6
    ENDIF
    C=0._dl
    N=Size(x)
    if(N.ne.Size(Y,2)) stop "Dverk_V: array sizes do not agree"
    M=SIZE(Y,DIM=1)
    IND=1
    XSTART=X(1)
    DO I=2,N
       Y(1:M,I)=Y(1:M,I-1)
       CALL DVERK_S(FCN,XSTART,Y(1:M,I),X(I),IND,C,W,TOLL)
       IF(IND.NE.3) stop "Dverk Error"
    ENDDO
  END	SUBROUTINE DVERK_V

  SUBROUTINE DVERK_S (fcn,x,y,xend,ind,c,w,TOLER)
    !!subroutine fcn(n,x,y,yprime)	(where n=size(Y)=size(yprime))
    INTEGER(IB) ind, k,N
    REAL(DL),DIMENSION(:),INTENT(INOUT)::Y,C  !!C(24)
    REAL(DL) W(SIZE(Y),9)
    REAL(DL) x, xend
    REAL(DL),OPTIONAL::TOLER
    REAL(DL) tol,temp
    EXTERNAL FCN
    !     ******************************************************************
    !     * begin initialization, parameter checking, interrupt re-entries *
    !     ******************************************************************
    N=SIZE(Y)
    if(N.ne.SIZE(W,1)) stop "wrong dimension in DVERK"
    IF(PRESENT(TOLER))THEN
       TOL=MAX(ABS(TOLER),1.D-10)
    ELSE
       TOL=1.D-5
    ENDIF
    !  ......abort if ind out of range 1 to 6
    if (ind.lt.1 .or. ind.gt.6) GOTO 500
    !        cases - initial entry, normal re-entry, interrupt re-entries
    if (IND.NE.3) THEN
       if (ind==4) goto 1111
       if (ind==5 .or. ind==6) goto 2222
       !        case 1 - initial entry (ind .eq. 1 or 2)
       !  .........abort if n.gt.nw or tol.le.0
       if (ind.NE.2) THEN
          !              initial entry without options (ind .eq. 1)
          !              set c(1) to c(9) equal to 0
          c(1:9) = 0._dl
       ELSE
          !              initial entry with options (ind .eq. 2)
          !              make c(1) to c(9) non-negative
          c(1:9) = dabs(c(1:9))
          !              make floor values non-negative if they are to be used
          if (c(1).EQ.4._dl .OR. c(1).EQ.5._dl) c(31:30+N) = dabs(c(31:30+N)) 
       ENDIF
       !           initialize rreb, dwarf, prev xend, flag, counts
       c(10) = 1.3877787807814456755295395851135D-17	!!2^{-56}
       c(11) = 1.D-35
       !           set previous xend initially to initial value of x
       c(20) = x
       c(21:24) = 0._dl
       !        case 2 - normal re-entry (ind .eq. 3)
       !  .........abort if xend reached, and either x changed or xend not
    ELSE
       if (c(21).ne.0._dl .and. (x.ne.c(20) .or. xend.eq.c(20))) GOTO 500
       !           re-initialize flag
       c(21) = 0._dl
       !        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
       !           transfer control to the appropriate re-entry point..........
       !           this has already been handled by the computed GOTO        .
       !        end cases                                                     v
    ENDIF
50  continue
    !
    !     end initialization, etc.
    !
    !     ******************************************************************
    !     * loop through the following 4 stages, once for each trial  step *
    !     * until the occurrence of one of the following                   *
    !     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
    !     *        stage 4                                                 *
    !     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
    !     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
    !     *        requested, in stage 1 or stage 4                        *
    !     ******************************************************************
    !
    !
    !        ***************************************************************
    !        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
    !        * and some parameter  checking,  and  end  up  with  suitable *
    !        * values of hmag, xtrial and htrial in preparation for taking *
    !        * an integration step.                                        *
    !        ***************************************************************
    !
    !***********error return (with ind=-1) if no of fcn evals too great
    IF (c(7).NE.0._dl .AND. c(24).GE.c(7)) THEN
       ind = -1
       return
    ENDIF
    !           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
    IF (ind .NE. 6) THEN
       CALL FCN(n, x, y, w(1,1))
       c(24) = c(24) + 1._dl
    ENDIF
    !           calculate hmin - use default unless value prescribed
    c(13) = c(3)
    IF (c(3) .NE. 0._dl) GOTO 165
    !              calculate default value of hmin
    !              first calculate weighted norm y - c(12) - as specified
    !              by the error control indicator c(1)
    temp = 0._dl
    IF (c(1) .EQ. 1._dl) THEN
       !                 absolute error control - weights are 1                  
       temp = MAXVAL(dabs(y(1:N)))
       c(12) = temp
       GOTO 160
    ENDIF
    IF (c(1) .EQ. 2._dl) THEN
       !                 relative error control - weights are 1/dabs(y(k)) so
       !                 weighted norm y is 1
       c(12) = 1._DL 
       GOTO 160
    ENDIF
    IF (c(1) .EQ. 3._DL ) THEN
       !                 weights are 1/max(c(2),abs(y(k)))
       temp = dmax1(temp, MAXVAL(dabs(y(1:N))/c(2)))
       c(12) = dmin1(temp, 1._DL )
       GOTO 160
    ENDIF
    IF (c(1) .EQ. 4._DL ) THEN
       !                 weights are 1/max(c(k+30),abs(y(k)))
       temp = dmax1(temp, MAXVAL(dabs(y(k))/c(31:30+N)))
       c(12) = dmin1(temp, 1._DL )
       GOTO 160
    ENDIF
    IF (c(1) .EQ. 5._DL ) THEN
       !                 weights are 1/c(k+30)
       temp = dmax1(temp, MAXVAL(dabs(y(1:N))/c(31:30+N)))
       c(12) = temp
       GOTO 160
    ENDIF
    !                 default case - weights are 1/max(1,abs(y(k)))
    temp = dmax1(temp, MAXVAL(dabs(y(1:N))))
    c(12) = dmin1(temp, 1._DL )
160 continue
    c(13) = 10._DL *dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
165 continue
    !           calculate scale - use default unless value prescribed
    c(15) = c(5)
    if (c(5) .eq. 0._DL ) c(15) = 1._DL 
    !           calculate hmax - consider 4 cases
    !           case 1 both hmax and scale prescribed
    if (c(6).ne.0._DL  .and. c(5).ne.0._DL ) c(16) = dmin1(c(6), 2._DL /c(5))
    !           case 2 - hmax prescribed, but scale not
    if (c(6).ne.0._DL  .and. c(5).eq.0._DL ) c(16) = c(6)
    !           case 3 - hmax not prescribed, but scale is
    if (c(6).eq.0._DL  .and. c(5).ne.0._DL ) c(16) = 2._DL /c(5)
    !           case 4 - neither hmax nor scale is provided
    if (c(6).eq.0._DL  .and. c(5).eq.0._DL ) c(16) = 2._DL 
    !***********error return (with ind=-2) if hmin .gt. hmax
    IF (c(13) .GT. c(16))THEN
       ind = -2
       return
    ENDIF
    !           calculate preliminary hmag - consider 3 cases
    IF (ind .LE. 2) THEN
       !           case 1 - initial entry - use prescribed value of hstart, if
       !              any, else default
       c(14) = c(4)
       if (c(4) .eq. 0._DL ) c(14) = c(16)*tol**(1._DL /6._DL )
       GOTO 185
    ENDIF
    IF (c(23) .LE. 1._DL ) THEN
       !           case 2 - after a successful step, or at most  one  failure,
       !              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
       !              overflow. then avoid reduction by more than half.
       temp = 2._DL *c(14)
       if (tol .lt. (2._DL /.9_dl )**6*c(19)) &
            temp = .9_dl *(tol/c(19))**(1._DL /6._DL )*c(14)
       c(14) = dmax1(temp, .5_dl *c(14))
       GOTO 185
    ENDIF
    !           case 3 - after two or more successive failures
    c(14) = .5_dl *c(14)
185 continue
    !
    !           check against hmax
    c(14) = dmin1(c(14), c(16))
    !
    !           check against hmin
    c(14) = dmax1(c(14), c(13))
    !
    !***********interrupt no 1 (with ind=4) if requested
    if (c(8) .NE. 0._DL ) THEN
       ind = 4
       return
    ENDIF
    !           resume here on re-entry with ind .eq. 4   ........re-entry..
1111 continue
    !
    !           calculate hmag, xtrial - depending on preliminary hmag, xend
    if (c(14) .LT. dabs(xend - x)) THEN
       !    do not step more than half way to xend
       c(14) = dmin1(c(14), .5_dl *dabs(xend - x))
       c(17) = x + dsign(c(14), xend - x)
       GOTO 195
    ENDIF
    !              hit xend exactly
    c(14) = dabs(xend - x)
    c(17) = xend
195 continue
    !
    !           calculate htrial
    c(18) = c(17) - x
    !
    !        end stage 1
    !
    !        ***************************************************************
    !        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
    !        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
    !        * stage 3. w(*,9) is temporary storage until finally it holds *
    !        * ytrial.                                                     *
    !        ***************************************************************
    !
    temp = c(18)/1398169080000._DL 
    !
    w(1:N,9) = y(1:N) + temp*w(1:N,1)*233028180000._DL 
    CALL FCN(n, x + c(18)/6._DL , w(1,9), w(1,2))
    !
    w(1:N,9) = y(1:N) + temp*(   w(1:N,1)*74569017600._DL  &
         + w(1:N,2)*298276070400._DL   )
    CALL FCN(n, x + c(18)*(4._DL /15._DL ), w(1,9), w(1,3))
    !
    w(1:N,9) = y(1:N) + temp*(w(1:N,1)*1165140900000._DL  &
         - w(1:N,2)*3728450880000._DL  &
         + w(1:N,3)*3495422700000._DL  )
    CALL FCN(n, x + c(18)*(2._DL /3._DL ), w(1,9), w(1,4))
    !
    w(1:N,9) = y(1:N) + temp*( - w(1:N,1)*3604654659375._DL  &
         + w(1:N,2)*12816549900000._DL  &
         - w(1:N,3)*9284716546875._DL  &
         + w(1:N,4)*1237962206250._DL  )
    CALL FCN(n, x + c(18)*(5._DL /6._DL ), w(1,9), w(1,5))
    !
    w(1:N,9) = y(1:N) + temp*(   w(1:N,1)*3355605792000._DL  &
         - w(1:N,2)*11185352640000._DL  &
         + w(1:N,3)*9172628850000._DL  &
         - w(1:N,4)*427218330000._DL  &
         + w(1:N,5)*482505408000._DL   )
    CALL FCN(n, x + c(18), w(1,9), w(1,6))
    !
    w(1:N,9) = y(1:N) + temp*( - w(1:N,1)*770204740536._DL  &
         + w(1:N,2)*2311639545600._DL  &
         - w(1:N,3)*1322092233000._DL  &
         - w(1:N,4)*453006781920._DL  &
         + w(1:N,5)*326875481856._DL   )
    CALL FCN(n, x + c(18)/15._DL , w(1,9), w(1,7))
    !
    w(1:N,9) = y(1:N) + temp*(   w(1:N,1)*2845924389000._DL  &
         - w(1:N,2)*9754668000000._DL  &
         + w(1:N,3)*7897110375000._DL  &
         - w(1:N,4)*192082660000._DL  &
         + w(1:N,5)*400298976000._DL  &
         + w(1:N,7)*201586000000._DL   )
    CALL FCN(n, x + c(18), w(1,9), w(1,8))
    !
    !           calculate ytrial, the extrapolated approximation and store
    !              in w(*,9)
    w(1:N,9) = y(1:N) + temp*(   w(1:N,1)*104862681000._DL  &
         + w(1:N,3)*545186250000._DL  &
         + w(1:N,4)*446637345000._DL  &
         + w(1:N,5)*188806464000._DL  &
         + w(1:N,7)*15076875000._DL  &
         + w(1:N,8)*97599465000._DL    )
    !
    !           add 7 to the no of fcn evals
    c(24) = c(24) + 7._DL 
    !
    !        end stage 2
    !
    !        ***************************************************************
    !        * stage 3 - calculate the error estimate est. first calculate *
    !        * the  unweighted  absolute  error  estimate vector (per unit *
    !        * step) for the unextrapolated approximation and store it  in *
    !        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
    !        * specified by the error  control  indicator  c(1).  finally, *
    !        * modify  this result to produce est, the error estimate (per *
    !        * unit step) for the extrapolated approximation ytrial.       *
    !        ***************************************************************
    !
    !           calculate the unweighted absolute error estimate vector
    w(1:N,2) = (w(1:N,1)*8738556750._DL  &
         + w(1:N,3)*9735468750._DL  &
         - w(1:N,4)*9709507500._DL  &
         + w(1:N,5)*8582112000._DL  &
         + w(1:N,6)*95329710000._DL  &
         - w(1:N,7)*15076875000._DL  &
         - w(1:N,8)*97599465000._DL )/1398169080000._DL 
    !
    !           calculate the weighted max norm of w(*,2) as specified by
    !           the error control indicator c(1)
    temp = 0._DL 
    IF (c(1) .EQ. 1._DL ) THEN
       !              absolute error control
       temp = MAXVAL(dabs(w(1:N,2)))
       GOTO 360
    ENDIF
    IF (c(1) .EQ. 2._DL ) THEN
       !              relative error control
       temp = dmax1(temp, MAXVAL(dabs(w(1:N,2)/y(1:N))))
       GOTO 360
    ENDIF
    if (c(1) .EQ. 3._DL ) THEN
       !             weights are 1/max(c(2),abs(y(k)))
       temp = dmax1(temp, MAXVAL(dabs(w(1:N,2))/dmax1(c(2), dabs(y(1:N)))) )
       GOTO 360
    ENDIF
    if (c(1) .EQ. 4._DL ) THEN
       !              weights are 1/max(c(k+30),abs(y(k)))
       temp = dmax1(temp, MAXVAL(dabs(w(1:N,2))/dmax1(c(31:30+N), dabs(y(1:N)))) )
       GOTO 360
    ENDIF
    IF (c(1) .EQ. 5._DL ) THEN
       !              weights are 1/c(k+30)
       temp = dmax1(temp, MAXVAL(dabs(w(1:N,2)/c(31:30+N))))
       GOTO 360
    ENDIF
    !              default case - weights are 1/max(1,abs(y(k)))
    DO k = 1, n
       temp = dmax1(temp, dabs(w(K,2))/dmax1(1._DL , dabs(y(K))) ) 
    ENDDO
360 continue
    !
    !           calculate est - (the weighted max norm of w(*,2))*hmag*scale
    !              - est is intended to be a measure of the error  per  unit
    !              step in ytrial
    c(19) = temp*c(14)*c(15)
    !
    !        end stage 3
    !
    !        ***************************************************************
    !        * stage 4 - make decisions.                                   *
    !        ***************************************************************
    !
    !           set ind=5 if step acceptable, else set ind=6
    ind = 5
    if (c(19) .gt. tol) ind = 6
    !
    !***********interrupt no 2 if requested
    if (c(9) .NE. 0._DL ) RETURN
    !           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
2222 continue
    !
    IF (ind .NE. 6) THEN
       !              step accepted (ind .eq. 5), so update x, y from xtrial,
       !                 ytrial, add 1 to the no of successful steps, and set
       !                 the no of successive failures to zero
       x = c(17)
       y(1:N) = w(1:N,9)
       c(22) = c(22) + 1._DL 
       c(23) = 0._DL 
       !**************return(with ind=3, xend saved, flag set) if x .eq. xend
       if (x .EQ. xend) THEN
          ind = 3
          c(20) = xend
          c(21) = 1._DL 
          return
       ENDIF
    ELSE
       !              step not accepted (ind .eq. 6), so add 1 to the no of
       !                 successive failures
       c(23) = c(23) + 1._DL 
       !**************error return (with ind=-3) if hmag .le. hmin
       if (c(14) .LE. c(13)) THEN
          ind = -3
          return
       ENDIF
    ENDIF
    !
    !        end stage 4
    !
    GOTO 50
    !     end loop
    !
    !  begin abort action
500 write (*,*) 'Error in dverk, x =',x, 'xend=', xend
    STOP
    !  end abort action
    !
  END SUBROUTINE DVERK_S





  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !!%%%%%%%%%%%%%%%%%%% RANSEED %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Function ThreadSafeRand(seed)
#ifdef _OPENMP
    use omp_lib
    Integer(IB),parameter::Max_Num_Threads = 64
#else
    Integer(IB),parameter::Max_Num_Threads =1
#endif
    Integer(IB),Parameter :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    Integer(IB),optional::seed
    Integer(IB),dimension(Max_Num_Threads),save::idum=0
    Real(dl) ThreadSafeRand
    Real(dl),dimension(Max_Num_Threads),save:: am
    Integer(IB),dimension(Max_Num_Threads),save :: ix=-1
    Integer(IB),dimension(Max_Num_Threads),save:: iy=-1
    Integer(IB) i,th,k
#ifdef _OPENMP
    th=omp_get_thread_num()+1
#else
    th=1
#endif
    If(th.gt.Max_Num_Threads) stop "too many threads; ThreadSafeRand overflow"
    If(Present(seed))then
       if(th.gt.1)then
          write(*,*) "Can not initialize random seeds in multi-threads section"
          stop
       endif
       do i=1,Max_Num_Threads
          idum(i)=-i-mod(i*seed,IA)
       enddo
    endif
    if (idum(th) <= 0 .or. iy(th) < 0) then
       am(th)=nearest(1.0,-1.0)/IM
       iy(th)=ior(ieor(888889999,abs(idum(th))),1)
       ix(th)=ieor(777755555,abs(idum(th)))
       idum(th)=abs(idum(th))+1
    end if

    ix(th)=ieor(ix(th),ishft(ix(th),13))
    ix(th)=ieor(ix(th),ishft(ix(th),-17))
    ix(th)=ieor(ix(th),ishft(ix(th),5))
    k=iy(th)/IQ
    iy(th)=IA*(iy(th)-k*IQ)-IR*k
    if (iy(th) < 0) iy(th)=iy(th)+IM
    ThreadSafeRand = am(th)*ior(iand(IM,ieor(ix(th),iy(th))),1)
  End Function ThreadSafeRand


  Function Gaussian_random() !!thread safe
    Real(dl) Gaussian_Random,x(2),R2
    x(1)=ThreadSafeRand()
    x(2)=ThreadSafeRand()
    x=2._dl *x-1._dl 
    R2=X(1)*X(1)+X(2)*X(2)
    DO WHILE(R2.GT.1._DL )
       x(1)=ThreadSafeRand()
       x(2)=ThreadSafeRand()
       x=2._dl *x-1._dl        
       R2=X(1)*X(1)+X(2)*X(2)    
    ENDDO
    R2=SQRT(-2._DL *DLOG(R2)/R2)
    Gaussian_random=X(1)*R2
  End Function Gaussian_random

  Function Gaussian_cmplx() !!thread safe 
    complex(dlc) Gaussian_cmplx
    Real(dl) x(2),R2
    x(1)=ThreadSafeRand()
    x(2)=ThreadSafeRand()
    x=2._dl *x-1._dl 
    R2=X(1)*X(1)+X(2)*X(2)
    DO WHILE(R2.GT.1._DL )
       x(1)=ThreadSafeRand()
       x(2)=ThreadSafeRand()
       x=2._dl *x-1._dl        
       R2=X(1)*X(1)+X(2)*X(2)    
    ENDDO
    R2=SQRT(-DLOG(R2)/R2)
    Gaussian_cmplx=cmplx(X(1)*R2,X(2)*R2)
  End Function Gaussian_cmplx

  SUBROUTINE init_random_seed(rank)
    Real(dl) tmp
    Integer(IB) rank
    tmp=ThreadsafeRand(rank)
  END SUBROUTINE init_random_seed
  

  SUBROUTINE LOWERCASE(C)
    CHARACTER,INTENT(INOUT)::C
    INTEGER(IB) K
    K=ICHAR(C)
    IF(K.GE.65.AND.K.LE.90)THEN
       C=CHAR(K+32)
    ENDIF
  END SUBROUTINE LOWERCASE




  FUNCTION LOWER_CASE(C)
    CHARACTER C,LOWER_CASE
    LOWER_CASE=C
    CALL LOWERCASE(LOWER_CASE)
  END FUNCTION LOWER_CASE



  FUNCTION NOCASE_EQ(STR1,STR2)
    CHARACTER(LEN=*)STR1,STR2
    LOGICAL NOCASE_EQ
    INTEGER(IB) N1,N2,I
    N1=LEN_TRIM(STR1)
    N2=LEN_TRIM(STR2)
    IF(N1.NE.N2)THEN
       NOCASE_EQ=.FALSE.
       RETURN
    ENDIF
    DO I=1,N1
       IF(LOWER_CASE(STR1(I:I)).NE.LOWER_CASE(STR2(I:I)))THEN
          NOCASE_EQ=.FALSE.
          RETURN
       ENDIF
    ENDDO
    NOCASE_EQ=.TRUE.
  END FUNCTION NOCASE_EQ
  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE PRTSYSTIME(RESET)
    LOGICAL,OPTIONAL::RESET
    REAL(DL),SAVE::PRETIME=0._DL 
    INTEGER(IB) NOWTIME,COUNTRATE
    CALL SYSTEM_CLOCK(NOWTIME,COUNTRATE)
    IF(PRESENT(RESET))THEN
       IF(RESET)THEN
          PRETIME=NOWTIME/REAL(COUNTRATE)
          PRINT*,"===== time label reset to zero. ===="
          RETURN
       ENDIF
    ENDIF
    WRITE(*,'(A18,F10.3,A10)') "===== time label: ",NOWTIME/REAL(COUNTRATE)-PRETIME," sec ====="
  END SUBROUTINE PRTSYSTIME


  SUBROUTINE splines(X,Y,Y2,YPL,YPR)
    REAL(DL),DIMENSION(:),INTENT(IN)::X,Y
    REAL(DL),DIMENSION(:),INTENT(OUT)::Y2
    REAL(DL),OPTIONAL,INTENT(IN)::YPL,YPR
    INTEGER(IB) I,N
    REAL(DL) YIL,YIR,BET,DXR,DXL
    REAL(DL),DIMENSION(:),ALLOCATABLE::GAM
    N=SIZE(X)
    ALLOCATE(GAM(N-1))
    DXR=X(2)-X(1)
    YIR=(Y(2)-Y(1))/DXR
    IF(PRESENT(YPL))THEN
       Y2(1)=(YIR-YPL)/DXR*3._DL 
       GAM(1)= 0.5_DL 
    ELSE
       Y2(1)=0._DL 
       GAM(1)=0._DL 
    ENDIF
    DXR=DXR/6._DL 
    DO I=2,N-1
       DXL=DXR
       DXR=X(I+1)-X(I)
       BET=(X(I+1)-X(I-1))/3._DL -DXL*GAM(I-1)
       IF(DABS(BET) .EQ. 1.D-20) STOP 'Error in SPLINE.'
       YIL=YIR
       YIR=(Y(I+1)-Y(I))/DXR
       Y2(I)=(YIR-YIL-DXL*Y2(I-1))/BET
       DXR=DXR/6._DL 
       GAM(I)=DXR/BET
    ENDDO
    IF(PRESENT(YPR))THEN
       BET=(X(N)-X(N-1))/3._DL -DXR*GAM(N-1)
       IF(DABS(BET) .EQ. 1.D-20) STOP 'Error in SPLINE.'
       Y2(N)=(YPR-YIR-DXR*Y2(N-1))/BET
    ELSE
       Y2(N)=0._DL 
    ENDIF
    DO I=N-1,1,-1
       Y2(I)=Y2(I)-GAM(I)*Y2(I+1)
    ENDDO
    DEALLOCATE(GAM)
  END SUBROUTINE SPLINEs

  SUBROUTINE SPLINT_S(X,Y,Y2,XS,YS)
    REAL(DL),DIMENSION(:),INTENT(IN)::X,Y,Y2
    REAL(DL) XS,YS
    INTEGER(IB) N,J,L,R
    INTEGER(IB),SAVE::PRER=100000
    REAL(DL)::A,B
    N=SIZE(X)
    IF(XS.LE.X(1))THEN
       YS=Y(1)
       RETURN
    ELSEIF(XS.GE.X(N))THEN
       YS=Y(N)
       RETURN
    ENDIF
    IF(PRER.LE.N)THEN
       IF(X(PRER).GE.XS .AND. X(PRER-1).LE.XS)THEN
          R=PRER
          GOTO 100
       ENDIF
    ENDIF
    L=1
    R=N
    DO WHILE(R-L.GT.1)
       J=ISHFT(R+L,-1)
       IF(X(J).GT.XS)THEN
          R=J
       ELSE
          L=J
       ENDIF
    ENDDO
    PRER=R
100 A=(X(R)-XS)/(X(R)-X(R-1))
    B=1._DL -A
    YS=Y(R-1)*A+Y(R)*B+  &
         (Y2(R-1)*(A*A-1._DL )*A+Y2(R)*(B*B-1._DL )*B)/6._DL *(X(R)-X(R-1))**2
  END SUBROUTINE SPLINT_S

  SUBROUTINE SPLINT_V(X,Y,Y2,XS,YS)
    !! xS should be sorted s.t. xS(1)<xS(2)<...
    REAL(DL),DIMENSION(:),INTENT(IN)::X,Y,Y2,XS
    REAL(DL),DIMENSION(:),INTENT(OUT)::YS
    INTEGER(IB) I,N
    N=SIZE(XS)
    if(N.ne.SIZE(YS)) stop "SPLINT: Array sizes do not agree."
    DO I=1,N
       CALL SPLINT_s(X,Y,Y2,XS(I),YS(I))
    ENDDO
  END SUBROUTINE SPLINT_V

  FUNCTION FileExists(FILENAME)
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    LOGICAL FileExists
    INQUIRE(FILE=FILENAME, EXIST = FileExists)
  END FUNCTION FileExists


  SUBROUTINE qromb(func,a,b,ss,notconverge)
    INTEGER(IB) :: JMAX,JMAXP,K,KM
    REAL(DL)  :: a,b,func,ss,EPS
    logical,optional::notconverge
    EXTERNAL func
    PARAMETER (EPS=1.D-5, JMAX=25, JMAXP=JMAX+1, K=5, KM=K-1)
    !USES polint,trapzd
    INTEGER(IB) :: j
    REAL(DL)  :: dss,h(JMAXP),s(JMAXP)
    if(present(notconverge))notconverge=.false.
    h(1)=1.
    DO j=1,JMAX
       CALL D_trapzd(func,a,b,s(j),j)
       IF (j.GE.K) THEN
          if(maxval(s(j-km:j))-minval(s(j-km:j)).le.eps*maxval(abs(s(j-km:j))) .or. maxval(abs(s(km:j))).lt.1.d-30)THEN
             ss=sum(s(j-km:j))/k
             return
          endif
          CALL d_polint(h(j-KM),s(j-KM),K,0._dl ,ss,dss)
          IF (isnan(ss).or.ABS(dss).LE.EPS*ABS(ss).or.ABS(dss).LT.1.D-30) RETURN
       ENDIF
       s(j+1)=s(j)
       h(j+1)=0.25_DL *h(j)
    ENDDO
    WRITE(*,*) 'too many steps in qromb'
    if(present(notconverge))notconverge=.true.
  END SUBROUTINE qromb

 SUBROUTINE D_trapzd(func,a,b,s,n)
    INTEGER(IB) :: n
    REAL(DL)  :: a,b,s,func
    EXTERNAL func
    INTEGER(IB) :: j,it
    REAL(DL)  :: del,sum,tnm,x
    IF (n.EQ.1) THEN
       s=0.5*(b-a)*(func(a)+func(b))
    ELSE
       it=ishft(1,n-2)
       tnm=it	
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       DO j=1,it
          sum=sum+func(x)
          x=x+del
       ENDDO
       s=0.5_DL *(s+(b-a)*sum/tnm)
    ENDIF
    RETURN
  END SUBROUTINE D_trapzd

  SUBROUTINE D_polint(xa,ya,n,x,y,dy)
    INTEGER(IB) :: n,NMAX
    REAL(DL)  :: dy,x,y,xa(n),ya(n)
    PARAMETER (NMAX=10)
    INTEGER(IB) :: i,m,ns
    REAL(DL)  :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    ns=1
    dif=ABS(x-xa(1))
    DO i=1,n
       dift=ABS(x-xa(i))
       IF (dift.LT.dif) THEN
          ns=i
          dif=dift
       ENDIF
       c(i)=ya(i)
       d(i)=ya(i)
    ENDDO
    y=ya(ns)
    ns=ns-1
    DO m=1,n-1
       DO i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          IF(den.EQ.0.) THEN
             WRITE(*,*) 'failure in polint'
             STOP
          ENDIF
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       ENDDO
       IF (2*ns.LT.n-m)THEN
          dy=c(ns+1)
       ELSE
          dy=d(ns)
          ns=ns-1
       ENDIF
       y=y+dy
    ENDDO
    RETURN
  END SUBROUTINE D_polint

  Subroutine Matrix_Inverse(A,n)
    integer(IB) n
    Real(dl) A(n,n)
    Real(dl) tmp
    select case(n)
    case(1)
       A=1._dl /A
    case(2) 
       A=ReShape( (/ A(2,2),-A(2,1), -A(1,2),A(1,1) /) / (A(1,1)*A(2,2)-A(1,2)*A(2,1)), (/ 2, 2 /) )
    case default
       call MatInv(A,n)
    end select
  End Subroutine Matrix_Inverse

  SUBROUTINE SWAP(X,Y)
    REAL(DL),DIMENSION(:),INTENT(INOUT)::X,Y
    REAL(DL) TEMP(SIZE(X))
    temp=x
    x=y
    y=temp
  END SUBROUTINE SWAP

  FUNCTION imaxloc(arr)
    REAL(DL), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(IB) imaxloc
    INTEGER(IB),DIMENSION(1)::imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
  END FUNCTION imaxloc

  FUNCTION OUTERPROD(A,B)
    REAL(DL), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DL), DIMENSION(size(a),size(b)) :: outerprod
    OUTERPROD = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION OUTERPROD

  SUBROUTINE LUDCMP(A,INDX,D,n)
    integer(IB) n
    Real(dl) a(n,n)
    INTEGER(IB)  indx(n)
    Real(dl),INTENT(OUT) :: d
    Real(dl) vv(n)
    Real(dl),PARAMETER::TINY=1.E-20
    INTEGER(IB)  j,imax
    D=1._DL 
    vv=maxval(abs(a),dim=2)
    if (any(vv .EQ. 0.0)) stop 'singular matrix in ludcmp'
    vv=1._dl /vv
    DO j=1,n
       imax=(j-1)+imaxloc(vv(j:n)*DABS(a(j:n,j)))
       IF (j .NE. imax) THEN
          CALL SWAP(a(imax,1:N),a(j,1:N))
          d=-d
          vv(imax)=vv(j)
       ENDIF
       indx(j)=imax
       if (a(j,j) .EQ. 0.0) a(j,j)=TINY
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
    ENDDO
  END SUBROUTINE LUDCMP

  SUBROUTINE LUBKSB(A,indx,b,n)
    integer(IB) n
    Real(dl) A(n,n)
    INTEGER(IB) indx(n)
    Real(dl) b(n)
    INTEGER(IB):: i,ii,ll
    Real(dl) :: summ
    ii=0
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii=i
       end if
       b(i)=summ
    end do
    do i=n,1,-1
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
  END SUBROUTINE LUBKSB

  SUBROUTINE MATSOLVE(A,B,n,m)
    integer(IB) n,m
    Real(dl) A(n,n)
    Real(dl) B(n,m)
    INTEGER(IB) INDX(n)
    Real(dl) D
    INTEGER(IB) I
    CALL LUDCMP(A,INDX,D,n)
    DO I=1,m
       CALL LUBKSB(A,INDX,B(:,I),n)
    ENDDO
  END SUBROUTINE MATSOLVE

  SUBROUTINE MatInv(A,n)
    integer(IB) n
    Real(dl) A(n,n)
    Real(dl) ACOPY(n,n)
    integer(IB) i
    Acopy=A
    A=0._dl 
    do i=1,size(A)
       A(i,i)=1._dl 
    enddo
    CALL MATSOLVE(ACOPY,A,n,n)
  END SUBROUTINE MatInv

  SUBROUTINE choldc(n,a,p) !!a=LL^T,
!! only the upper left of a is required
    INTEGER(IB) n
    REAL(DL) a(n,n)
    REAL(DL) p(n)
    integer(IB) i
    REAL(DL) summ
    do i=1,n
       summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
       if (summ <= 0.0) stop 'choldc failed'
       p(i)=sqrt(summ)
       a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
    end do
  END SUBROUTINE choldc

  SUBROUTINE cholsl(n,a,p,b,x) !!solve LL^T x = b
    integer(IB) n
    Real(dl) a(n,n)
    Real(dl) p(n),b(n)
    Real(dl) x(n)
    integer(IB) i
    do i=1,n
       x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/p(i) !!solve L y = b
    end do
    do i=n,1,-1
       x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/p(i) !! solve L^T x = y
    end do
  END SUBROUTINE cholsl

  subroutine cholinv(n,a,p)  !!invert L
    integer(IB) n
    Real(dl) a(n,n),p(n)
    integer(IB) i,j
    do  i=1,n
       a(i,i)=1./p(i)
       do  j=i+1,n
          a(j,i)=-dot_product(a(j,i:j-1),a(i:j-1,i))/p(j)
       enddo
    enddo
  end subroutine cholinv

  subroutine SymMatInv(A,n)
    integer(IB) n
    Real(dl) A(n,n),p(n)
    integer(IB) i,j
    select case(n)
    case(1)
       A=1._dl /A
    case(2)
       A=ReShape( (/ A(2,2),-A(2,1), -A(1,2),A(1,1) /) / (A(1,1)*A(2,2)-A(1,2)*A(2,1)), (/ 2, 2 /) )
    case default
       call choldc(n,a,p)
       call cholinv(n,a,p)
       do i=1,n
          do j=i+1,n
             A(i,j)=dot_product(A(j:n,i),A(j:n,j))
          enddo
       enddo
       do i=1,n-1
          A(i,i)=dot_product(A(i:n,i),A(i:n,i))
          A(i+1:n,i)=A(i,i+1:n)
       enddo
       A(n,n)=A(n,n)**2
    end select
  end subroutine SymMatInv


  function open_file(filename,mode) 
    !!supported modes:
    !!    'w' (text mode for read/write),
    !!    'b' (binary mode for read/write), 
    !!    'r' (text mode for read),
    !!    'rb' (binary mode for read)
    !!    'a' (text mode for append)
    type(file_pointer) open_file
    character(LEN=*) filename
    character(LEN=*),optional::mode
    character(LEN=2)::open_mode
    integer(IB) funit
    funit = 10
    do while (unit_opened(funit))
       funit = funit+1
       if (funit.gt.99)then
          stop "error in open_txt_file: all file units are occupied"
       endif
    enddo
    open_file%unit = funit
    open_file%pt = 1
    if(present(mode))then
       open_mode = mode
    else
       open_mode = 'w'
    endif
    select case(open_mode)
    case("w")
       open(funit,FILE=trim(filename),FORM="FORMATTED",STATUS="UNKNOWN",ACCESS='SEQUENTIAL',ACTION="READWRITE", Err=200)
    case("r")
       open(funit,FILE=trim(filename),FORM="FORMATTED",STATUS="UNKNOWN", ACCESS='SEQUENTIAL', ACTION="READ", Err=200)
    case("a")
       open(funit,FILE=trim(filename),FORM="FORMATTED",STATUS="UNKNOWN", ACCESS='APPEND',ACTION="READWRITE", Err=200)
    case("b")
       open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READWRITE",RECL=IB,ERR=200) 
       open_file%pt = 1
    case("rb")
       open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="DIRECT",ACTION="READ",RECL=IB,ERR=200) 
    case('u')
       open(UNIT=funit,FILE=trim(filename),FORM="UNFORMATTED",STATUS="UNKNOWN", ACCESS="SEQUENTIAL",ACTION="READWRITE",ERR=200)        
    case default
       write(*,*) "Error in open_file: unknow mode "//trim(mode)
       goto 200
    end select

    unit_opened(funit) = .true.
    return
200 write(*,*) "Error: can not open file "//trim(filename)
    stop
  end function open_file


  subroutine int_to_binary_file(fp, i)
    type(file_pointer) fp
    integer(IB) i
    write(fp%unit,REC=fp%pt) i
    fp%pt=fp%pt+1
  end subroutine int_to_binary_file


  subroutine float_to_binary_file(fp, x)
    type(file_pointer) fp
    real(dl) x
    real(IB) xapp
    if(DL.le.IB)then
       write(fp%unit,REC=fp%pt) x
       fp%pt=fp%pt +1
    else
       xapp = x
       write(fp%unit,REC=fp%pt) xapp
       xapp = x-xapp
       write(fp%unit,REC=fp%pt+1) xapp
       fp%pt=fp%pt+2
    endif
  end subroutine float_to_binary_file

  subroutine int_from_binary_file(fp, i)
    type(file_pointer) fp
    integer(IB) i
    read(fp%unit, REC=fp%pt) i
    fp%pt=fp%pt+1
  end subroutine int_from_binary_file
 

  subroutine float_from_binary_file(fp, x)
    type(file_pointer) fp
    real(dl) x
    real(IB) x1,x2
    if(dl .le. IB)then
       read(fp%unit, REC=fp%pt) x
       fp%pt =fp%pt +1
    else
       read(fp%unit, REC=fp%pt) x1
       read(fp%unit, REC=fp%pt+1) x2
       x = x1 + (x2 + 0._dl)
       fp%pt=fp%pt+2
    endif
  end subroutine float_from_binary_file

  function new_file_unit()
    integer(IB) new_file_unit
    new_file_unit = 11
    do while(unit_opened(new_file_unit))
       new_file_unit = new_file_unit + 1
       if(new_file_unit .gt. 99)stop "All file units have been used."
    enddo
    unit_opened(new_file_unit)=.true.
    return
  end function new_file_unit

  subroutine close_file_unit(funit)
    integer(IB) funit
    close(funit)
    unit_opened(funit) = .false.
  end subroutine close_file_unit


  subroutine close_file(fp)
    type(file_pointer) fp
    if( fp%unit.ge.7 .and. fp%unit.le.99)then
       close(fp%unit)
       unit_opened(fp%unit)=.false.
       fp%pt = 1
       fp%unit = 0
    else
       write(*,*) "can not close unknow file unit: "//trim(int2str(fp%unit))
       stop
    endif    
  end subroutine close_file

  subroutine print_bar()
    write(*,*) "***********************************************************"
  end subroutine print_bar

  function fast_exp(x) !!for |x|<<1
    real(dl),parameter::c0 = 1._dl
    real(dl),parameter::c1 = c0/1._dl/16._dl
    real(dl),parameter::c2 = c1/2._dl/16._dl
    real(dl),parameter::c3 = c2/3._dl/16._dl
    real(dl),parameter::c4 = c3/4._dl/16._dl
    real(dl),parameter::c5 = c4/5._dl/16._dl
    real(dl),parameter::c6 = c5/6._dl/16._dl
    real(dl) x,fast_exp
    fast_exp =((( (c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+c6*x)))))) **2)**2)**2)**2
  end function fast_exp

End Module Utils
