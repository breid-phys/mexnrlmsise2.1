#include "fintrf.h"

C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
c      implicit none
      use msis_init, only          : msisinit
      use msis_calc, only          : msiscalc
C     Declarations
C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
      
C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN
      integer*4 mexPrintf
      integer*4 mxGetString
      
C     Local Variables      
      mwPointer in_ptr, alt_ptr, D_ptr, T_ptr
      
C     Array information:
      mwPointer mrows, ncols, malts
      mwSize asize, dsize, tsize

      DIMENSION D(10),T(2)
      real*8 IDAY
      real*8 input(3)
      real*8 SEC
      real*8 F107A, F107
      real*8 AP(7)
      real*8 SW(25)

      INTEGER i
      logical OK

      mwSize  maxbuf
      parameter(maxbuf = 1000)
      character*1000 input_buf
      mwPointer strlen

      dsize = 10
      tsize = 2
      
      if(nrhs .ne. 8) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:nInput',
     &  'Incorrect number of input arguments.')
      elseif(nlhs .ne. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:nOutput',
     &  'Incorrect number of output arguments.')
      endif
      
c     parameter path
C     Get the length of the input string.
      strlen = mxGetM(prhs(8))*mxGetN(prhs(8))
      if (strlen .gt. maxbuf) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:maxbuf',
     +       'Max parameter path length 1000.')
      endif

C     Get the string contents (dereference the input pointer).
      status = mxGetString(prhs(8), input_buf, maxbuf)

C     Check if mxGetString is successful.
      if (status .ne. 0) then 
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:readError',
     +             'Error reading parameter path.')
      endif
      
      inquire( file=trim(input_buf)//'msis21.parm', exist=OK)
      if (.not.OK) then
        call mexErrMsgIdAndTxt('MATLAB:NRLMSISE2_1:parm',
     +  trim(input_buf)//'msis21.parm does not exist.')
      endif
      
c     DAY       Day of year (1.0 to 365.0 or 366.0)
      in_ptr = mxGetPr(prhs(1))
      mrows = mxGetM(prhs(1))
      ncols = mxGetN(prhs(1))
      asize = mrows*ncols
      call mxCopyPtrToReal8(in_ptr,IDAY,asize)
      
      if(mrows*ncols.ne.1) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:DAY',
     &                           'INPUT: Day must be scalar')
      elseif((IDAY.lt.1.0).or.(IDAY.gt.366.0)) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:DAY',
     &                           'INPUT: Day out of range.')      
      endif
      
c     UTSEC     Universal time (seconds)
      in_ptr = mxGetPr(prhs(2))
      mrows = mxGetM(prhs(2))
      ncols = mxGetN(prhs(2))
      asize = mrows*ncols
      call mxCopyPtrToReal8(in_ptr,SEC,asize)
      
      if(mrows*ncols.ne.1) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:UTSEC',
     &                           'INPUT: UTSEC must be scalar')
      elseif((SEC.lt.0.0).or.(SEC.gt.86400.0)) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:UTSEC',
     &                           'INPUT: UTSEC out of range.')      
      endif
      
c     F107A
      in_ptr = mxGetPr(prhs(3))
      mrows = mxGetM(prhs(3))
      ncols = mxGetN(prhs(3))
      asize = mrows*ncols
      call mxCopyPtrToReal8(in_ptr,F107A,asize)
      
      if(mrows*ncols.ne.1) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:F107A',
     &                           'INPUT: F107A must be scalar')
      elseif(F107A.lt.0.0) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:F107A',
     &                           'INPUT: F107A out of range.')      
      endif

c     F107
      in_ptr = mxGetPr(prhs(4))
      mrows = mxGetM(prhs(4))
      ncols = mxGetN(prhs(4))
      asize = mrows*ncols
      call mxCopyPtrToReal8(in_ptr,F107,asize)

      if(mrows*ncols.ne.1) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:F107',
     &                           'INPUT: F107 must be scalar')
      elseif(F107.lt.0.0) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:F107',
     &                           'INPUT: F107 out of range.')      
      endif

c     AP
      in_ptr = mxGetPr(prhs(5))
      mrows = mxGetM(prhs(5))
      ncols = mxGetN(prhs(5))
      asize = mrows*ncols
      call mxCopyPtrToReal8(in_ptr,AP,asize)     
 
       if((mrows*ncols.ne.1).and.(mrows*ncols.ne.7)) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:AP',
     &   'INPUT: AP must have size [1] or size [7]')  
      endif
 
 
c     SW
      in_ptr = mxGetPr(prhs(6))
      mrows = mxGetM(prhs(6))
      ncols = mxGetN(prhs(6))
      asize = mrows*ncols
      call mxCopyPtrToReal8(in_ptr,SW,asize) 
      
      if(mrows*ncols.ne.25) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:SW',
     &                           'INPUT: SW must have size [25]')  
      endif
      
      
      call msisinit(switch_legacy=sngl(SW), 
     & parmpath=trim(input_buf))
            
c     get inputs
      in_ptr = mxGetPr(prhs(7))
      mrows = mxGetM(prhs(7))
      ncols = mxGetN(prhs(7))
      asize = 1*mrows
      
      if(mrows.ne.3) then
         call mexErrMsgIdAndTxt ('MATLAB:NRLMSISE2_1:mrows',
     &          'ALT/LAT/LON Input must have size [3,n]')
      endif
      
c     get output
      plhs(1) = mxCreateDoubleMatrix(dsize,ncols,0)
      D_ptr = mxGetPr(plhs(1))   
      plhs(2) = mxCreateDoubleMatrix(tsize,ncols,0)
      T_ptr = mxGetPr(plhs(2)) 
       
      DO I=1,ncols
      call mxCopyPtrToReal8(in_ptr+(8*mrows*(i-1)),input,asize)
      CALL MSISCALC(sngl(IDAY), sngl(SEC), 
     & sngl(input(1)), sngl(input(2)), sngl(input(3)), 
     & sngl(F107A), sngl(F107), sngl(AP), T(1), D, T(2))
      
      call mxCopyReal8ToPtr(dble(D),D_ptr+dsize*(8*(i-1)),dsize)
      call mxCopyReal8ToPtr(dble(T),T_ptr+tsize*(8*(i-1)),tsize)
      
      ENDDO
 
      return
      end

