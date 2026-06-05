! module containing some utilitary little parsing subroutines
! I am using a lot in my codes.
! Parsing functions : They allow to insert value through the
!   command line.
!   example : ./domywork val=10.
!   I may catch this value and assign it to the variable x
!   in the code through this code line :
!           if (parse_arg('val',x) != PARSE_OK) stop
!   The function will be chosen according to the argument type.
!   A third optional variable may be use to insert a file id in
!   order to use these functions also on a file.
!######################################################################
module forparse
use generic
implicit none

  integer, parameter :: PARSE_OK=0
  integer, parameter :: PARSE_UNKNOWN_KEYWORD=1
  integer, parameter :: PARSE_TYPE_ERROR=2
  integer, parameter :: PARSE_NO_ARGUMENT_FOUND=3
  integer, parameter :: PARSE_WRONG_NIT_VALUE=4
  integer, parameter :: PARSE_ERROR_PARSING_MULTI_VALUES=5
  integer, parameter :: PARSE_STRING_OVERFLOW=6

  private parse_int,parse_real,parse_string,parse_real_pr

  interface parse_arg
    module procedure parse_int,parse_real,parse_string,parse_real_pr
  end interface parse_arg

contains

!######################################################################
function helpisneeded(keyword)

  character(len=*) :: keyword
  logical :: helpisneeded

  integer :: n,nl
  character(100) :: line     

  helpisneeded=.false.
  nl=len_trim(keyword)
  n=command_argument_count()
  if (n == 0 .and. nl == 0) then
     helpisneeded=.true.
     return
  endif
  if (n == 1 .and. nl /= 0) then
     call get_command_argument(n,line)
     helpisneeded=(keyword == line(1:nl))
  endif
  return
end function helpisneeded
!######################################################################
function parse_real_pr(keyword,value,nit)
  use generic                        ! brings pr into scope
  real(pr) :: value
  character(len=*) :: keyword
  integer :: parse_real_pr
  integer,optional :: nit

  integer :: i,k,n,l,ios
  character(100) :: line     
  logical :: string2val

  if (present(nit)) then
     if (nit == 0) then
        parse_real_pr=PARSE_WRONG_NIT_VALUE
        return
     endif
     if (nit < 0) then
        n=command_argument_count()  
     else
        n=0
        rewind(nit)
     endif
  else
     n=command_argument_count()
     if (n == 0) then
        parse_real_pr=PARSE_NO_ARGUMENT_FOUND
        return
     endif
  endif
  i=1
  ios=0
  do while ((i <= n.and.n /= 0).or.(ios == 0.and.n == 0))
     if (n /= 0) then
        call get_command_argument(i,line)
     else
        read(nit,'(a)',iostat=ios) line
        if (ios /= 0.or.line(1:1) == '#') cycle
     endif
     if (linematchkeyword(line,keyword,k)) then
        string2val=string2real_pr(line(k:len_trim(line)),value)
        if (string2val) then
           parse_real_pr=PARSE_OK
        else
           parse_real_pr=PARSE_TYPE_ERROR
        endif
        return
     endif
     i=i+1
  enddo
  parse_real_pr=PARSE_UNKNOWN_KEYWORD
  return
end function parse_real_pr
!######################################################################
function parse_real(keyword,value,nit)

  real :: value
  character(len=*) :: keyword
  integer :: parse_real
  integer,optional :: nit

  integer :: i,k,n,l,ios
  character(100) :: line     
  logical :: string2val

  if (present(nit)) then
     if (nit == 0) then
        parse_real=PARSE_WRONG_NIT_VALUE
        return
     endif
     if (nit < 0) then
        n=command_argument_count()  
     else
        n=0
        rewind(nit)
     endif
  else
     n=command_argument_count()
     if (n == 0) then
        parse_real=PARSE_NO_ARGUMENT_FOUND
        return
     endif
  endif
  i=1
  ios=0
  do while ((i <= n.and.n /= 0).or.(ios == 0.and.n == 0))
     if (n /= 0) then
        call get_command_argument(i,line)
     else
        read(nit,'(a)',iostat=ios) line
        if (ios /= 0.or.line(1:1) == '#') cycle
     endif
     if (linematchkeyword(line,keyword,k)) then
        string2val=string2real(line(k:len_trim(line)),value)
        if (string2val) then
           parse_real=PARSE_OK
        else
           parse_real=PARSE_TYPE_ERROR
        endif
        return
     endif
     i=i+1
  enddo
  parse_real=PARSE_UNKNOWN_KEYWORD
  return
end function parse_real

!######################################################################
function parse_int(keyword,value,nit)

  integer :: value
  character(len=*) :: keyword
  integer :: parse_int
  integer,optional :: nit

  integer :: i,k,n,l,ios
  character(100) :: line     
  logical :: string2val

  if (present(nit)) then
     if (nit == 0) then
        parse_int=PARSE_WRONG_NIT_VALUE
        return
     endif
     if (nit < 0) then
        n=command_argument_count()  
     else
        n=0
        rewind(nit)
     endif
  else
     n=command_argument_count()
     if (n == 0) then
        parse_int=PARSE_NO_ARGUMENT_FOUND
        return
     endif
  endif
  i=1
  ios=0
  do while ((i <= n.and.n /= 0).or.(ios == 0.and.n == 0))
     if (n /= 0) then
        call get_command_argument(i,line)
     else
        read(nit,'(a)',iostat=ios) line
        if (ios /= 0.or.line(1:1) == '#') cycle
     endif
     if (linematchkeyword(line,keyword,k)) then
        string2val=string2int(line(k:len_trim(line)),value)
        if (string2val) then
           parse_int=PARSE_OK
        else
           parse_int=PARSE_TYPE_ERROR
        endif
        return
     endif
     i=i+1
  enddo
  parse_int=PARSE_UNKNOWN_KEYWORD
  return
end function parse_int

!######################################################################
function parse_string(keyword,value,nit)

  character(*) :: value
  character(len=*) :: keyword
  integer :: parse_string
  integer,optional :: nit

  integer :: i,k,n,l,ios
  character(100) :: line     

  if (present(nit)) then
     if (nit == 0) then
        parse_string=PARSE_WRONG_NIT_VALUE
        return
     endif
     if (nit < 0) then
        n=command_argument_count()  
     else
        n=0
        rewind(nit)
     endif
  else
     n=command_argument_count()
     if (n == 0) then
        parse_string=PARSE_NO_ARGUMENT_FOUND
        return
     endif
  endif
  i=1
  ios=0
  do while ((i <= n.and.n /= 0).or.(ios == 0.and.n == 0))
     if (n /= 0) then
        call get_command_argument(i,line)
     else
        read(nit,'(a)',iostat=ios) line
        if (ios /= 0.or.line(1:1) == '#') cycle
     endif
     if (linematchkeyword(line,keyword,k)) then
        value=line(k:len_trim(line))
        parse_string=PARSE_OK
        return
     endif
     i=i+1
  enddo
  parse_string=PARSE_UNKNOWN_KEYWORD
  return
end function parse_string
!######################################################################
function linematchkeyword(line,keyword,addr)

  logical :: linematchkeyword
  integer :: addr
  character(len=*) :: line,keyword

  integer :: k
      
  linematchkeyword=.false.
  k=index(line,'=')
  if (k.eq.0) return
  linematchkeyword=(trim(line(1:k-1))==keyword)

  if (linematchkeyword) addr=k+1
  return
end function linematchkeyword

!######################################################################
function string2int(string,val)

  integer :: val
  logical :: string2int
  character(len=*) :: string

  integer :: ios

  string2int=.false.
  read(string(1:len_trim(string)),*,iostat=ios) val
  string2int=(ios.eq.0)
  return
end function string2int
!######################################################################
function string2real_pr(string,val)

  real(pr) :: val
  logical :: string2real_pr
  character*(*) :: string
      
  integer :: ios

  string2real_pr=.false.
  read(string(1:len_trim(string)),*,iostat=ios) val
  string2real_pr=(ios.eq.0)
  return
end function string2real_pr
!######################################################################
function string2real(string,val)

  real :: val
  logical :: string2real
  character*(*) :: string
      
  integer :: ios

  string2real=.false.
  read(string(1:len_trim(string)),*,iostat=ios) val
  string2real=(ios.eq.0)
  return
end function string2real
!######################################################################
function parse_struct(keyword,npar,param,lpar,symb,nit)

  character*(*) :: keyword
  character*1 :: symb 
  character*20, dimension(10) :: param
  integer :: parse_struct,npar
  integer, dimension(10) :: lpar(10)
  integer,optional :: nit

  integer :: i,k,n,l,ios
  character(100) :: line

  if (present(nit)) then
     if (nit == 0) then
        parse_struct=PARSE_WRONG_NIT_VALUE
        return
     endif
     if (nit < 0) then
        n=command_argument_count()
     else
        n=0
        rewind(nit)
     endif
  else
     n=command_argument_count()
     if (n == 0) then
        parse_struct=PARSE_NO_ARGUMENT_FOUND
        return
     endif
  endif
  i=1
  ios=0
  do while ((i <= n.and.n /= 0).or.(ios == 0.and.n == 0))
     if (n /= 0) then
        call get_command_argument(i,line)
     else
        read(nit,'(a)',iostat=ios) line
        if (ios /= 0) cycle
     endif
     if (linematchkeyword(line,keyword,k)) then
        parse_struct=parse_line(line(k:len_trim(line)),npar,param,lpar,symb)
        return
     endif
     i=i+1
  enddo
  parse_struct=PARSE_UNKNOWN_KEYWORD
  return
end function parse_struct
!######################################################################
function parse_line(line,npar,param,lpar,symb)

  character*(*) :: line
  character*1 :: symb
  character*20, dimension(10) :: param
  integer :: parse_line,npar
  integer, dimension(10) :: lpar(10)

  integer :: k,kv

  npar=0
  k=0
  parse_line=PARSE_OK
  do while (.true.)
     npar=npar+1
     kv=k+1
     do while (line(kv:kv) /= symb) 
        kv=kv+1
        if (kv > len_trim(line)) then
           param(npar) (1:kv-k-1)=line(k+1:kv-1)
           lpar(npar)=kv-k-1
           if (lpar(npar) > 20) parse_line=PARSE_STRING_OVERFLOW
           return
        endif
     enddo
     param(npar) (1:kv-k-1)=line(k+1:kv-1)
     lpar(npar)=kv-k-1
     if (lpar(npar) > 20) then 
        parse_line=PARSE_STRING_OVERFLOW
        return
     endif
     k=kv
  enddo
end function parse_line
!######################################################################
end module forparse
