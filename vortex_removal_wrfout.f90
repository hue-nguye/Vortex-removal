      program tcrmwrf
!
!============================================================================

      implicit none
      INCLUDE 'netcdf.inc'
      integer  :: jdim
      parameter (jdim=6)
      integer ncid, status
      integer ishape(jdim)
      integer ishape2(jdim)
      character cval*50
      character name*31
      character (len=31),allocatable, dimension(:)      :: dname
      integer,           allocatable, dimension(:)       :: dval, dval2
      real,              allocatable, dimension(:,:,:,:) :: data, data2
      real,              allocatable, dimension(:,:)     :: xlat, xlon
      double precision,  allocatable, dimension(:,:,:,:) :: ddata, ddata2
      integer,           allocatable, dimension(:,:,:,:) :: idata, idata2
      character,         allocatable, dimension(:,:,:,:) :: text
      character omit(10)*80
      integer             :: start_dims(4)
      integer             :: dims(4)
      integer             :: dims_in(4), dims_out(4), box_start(3), box_end(3)
      integer             :: firstS,firstE, secondS,secondE, thirdS,thirdE
      integer             :: idm, ndims, nvars, natt, ngatts, nunlimdimid, iratio
      integer             :: i, ii, j,jj, iweg, jweg, isng, jsng, ibtg, jbtg, ix, iy
      integer             :: i_shape_we, i_shape_sn, i_shape_bt
      integer             :: i_shape_we_stag, i_shape_sn_stag, i_shape_bt_stag
      integer             :: ilen, itype, ival, na
      integer             :: mcid
      real                :: dx, rval, StormID
      real                :: new_cen,box_icent,box_jcent,V_max
      real                :: okm
      character (len=80)  :: input_file, output_file
      character (len=10)  :: option
      logical             :: debug=.FALSE.
      logical             :: debug_tc=.TRUE.
      logical             :: x_ave=.FALSE.
      logical             :: y_ave=.FALSE.
      logical             :: bit64
      logical             :: landuse
      integer, parameter  :: box_rad=55
!     
      integer                ::icent,jcent
      real,              allocatable, dimension(:,:) :: psfc_in,psfc_out
      real,              allocatable, dimension(:,:,:) :: u_in,u_out
      real,              allocatable, dimension(:,:,:) :: v_in,v_out
      real,              allocatable, dimension(:,:,:) :: t_in,t_out
      real,              allocatable, dimension(:,:,:) :: qv_in,qv_out
      real,              allocatable, dimension(:,:,:) :: p_in,p_out
      real,              allocatable, dimension(:,:,:) :: ph_in,ph_out
      real,              allocatable, dimension(:,:) :: landmask
      integer             :: nx, ny, nz, nx_u, ny_v, nz_p
      character(len=10) :: ID
!
      call read_args(input_file,option,box_icent,box_jcent,V_max,StormID,debug)
      write(ID,'(f10.0)')StormID
      output_file  = trim(ID)//trim(input_file)

      write(6,*) "#########################################"
      write(6,*) "Running TC RM WRF "
      write(6,*) 
      write(6,*) "INPUT FILE:         ",trim(input_file)
      write(6,*) "OUTPUT FILE:        ",trim(output_file)
      write(6,*) "OPTION:             ",option   
      write(6,*) "BOX_CENTER:         ",box_icent, box_jcent
      write(6,*) "VMAX:         ",V_max



! OPEN INPUT AND OUTPUT FILE
! output_file is input_file_new
      status = nf_open(input_file, 0, ncid)
      if (status .ne. nf_noerr) call handle_err(status)
      status = nf_create(output_file, 0, mcid)
      if (status .ne. nf_noerr) call handle_err(status)

! GET BASIC INFORMTION ABOUT THE FILE
! most important 
!   ndims:  number of dimensions
!   nvars:  number of variables
!   ngatts: number of global attributes
      status = nf_inq(ncid, ndims, nvars, ngatts, nunlimdimid)
      if (status .ne. nf_noerr) call handle_err(status)
      IF (debug) THEN
        write(6,*) 
        write(6,'(4(A,i4))') ' ndims = ',ndims,'    nvars = ',nvars,'    ngatts = ',ngatts, &
                   '    nunlimdimid =',nunlimdimid
        write(6,*) 
      ENDIF

! ALLOCATE SOME VARIABLES
      allocate (dval(ndims))
      allocate(dval2(ndims))
      allocate(dname(ndims))

! GET SOME BASIC DIMS FROM INPUT_FILE
      dx = -99.
      status = nf_get_att_real (ncid, nf_global, 'DX', dx)
      status = nf_get_att_int (ncid, nf_global, 'WEST-EAST_GRID_DIMENSION', iweg)
      status = nf_get_att_int (ncid, nf_global, 'SOUTH-NORTH_GRID_DIMENSION', isng)
      status = nf_get_att_int (ncid, nf_global, 'BOTTOM-TOP_GRID_DIMENSION', ibtg)
      IF (debug) THEN
        write(6,*) "BASICS from input file:"
        write(6,*) "       DX= ", dx
        write(6,*) "        X= ", iweg
        write(6,*) "        Y= ", isng
        write(6,*) "        Z= ", ibtg
      ENDIF
      if (dx .lt. 0.) stop 'dx is bad'

! CALCULATE BOX START, BOX END

    allocate(xlat(iweg,isng),xlon(iweg,isng))

    status = nf_inq_varid(ncid,'XLAT',i)
    status = nf_get_var_real(ncid, i, xlat)
    IF (debug_tc) write(6,*) '   SAMPLE VALUE = ',xlat(1,1),xlat(1,2),xlat(1,3)
    IF (debug_tc) write(6,*) '   SAMPLE VALUE = ',xlat(1,1),xlat(2,1),xlat(3,1)

    status = nf_inq_varid(ncid,'XLONG',i)
    status = nf_get_var_real(ncid, i, xlon)
    IF (debug_tc) write(6,*) '   SAMPLE VALUE = ',xlon(1,1),xlon(1,2),xlon(1,3)
    IF (debug_tc) write(6,*) '   SAMPLE VALUE = ',xlon(1,1),xlon(2,1),xlon(3,1)
   
    icent = -9999
    jcent = -9999

    DO jj = 1, iweg-1
       IF (xlat(1,jj).le.box_jcent .and. xlat(1,jj+1).ge.box_jcent) THEN
           jcent = jj
       ENDIF 
    END DO

    DO ii = 1, isng-1
       IF (xlon(ii,1).le.box_icent .and. xlon(ii+1,1).ge.box_icent) THEN
           icent = ii
       ENDIF 
    END DO

    DO ii = 1, 2
       IF (icent .eq. -9999 .or. jcent .eq. -9999) THEN
          write(6,*) "ERROR : Cannot found the TC center"
          EXIT
       ENDIF
    ENDDO
      
    write(6,*) "TC CENTER: ",icent,jcent

    box_start(1)=max(icent-box_rad,1)
    box_end(1)=min(icent+box_rad,iweg)

    box_start(2)=max(jcent-box_rad,1)
    box_end(2)=min(jcent+box_rad,isng)

    write(6,*) "BOX INFO : ",box_start(1),box_end(1),box_start(2),box_end(2)


! CALCULATE DIMS FOR OUTPUT FILE
      IF ( option == '-box' ) THEN
        okm = dx
        jweg = iweg
        jsng = isng
        jbtg = ibtg
        if ( box_end(1) .ne. 0 ) jweg = int(box_end(1) - box_start(1)) + 1
        if ( box_end(2) .ne. 0 ) jsng = int(box_end(2) - box_start(2)) + 1
        if ( box_end(3) .ne. 0 ) jbtg = int(box_end(3) - box_start(3)) + 1
      ELSE
        okm = dx
        jweg = iweg
        jsng = isng
        jbtg = ibtg
      ENDIF
      IF (debug) THEN
        write(6,*) "BASICS for output file:"
        write(6,*) "       DX= ", okm
        write(6,*) "        X= ", jweg
        write(6,*) "        Y= ", jsng
        write(6,*) "        Z= ", jbtg
      ENDIF
      !! We also need to fix the CEN_LAT and CEN_LON later, so get
      !! the middle of the new domain
      ix = int((jweg-1)/2.)
      iy = int((jsng-1)/2.)
      if ( ix .eq. int(jweg/2.) ) x_ave = .TRUE.
      if ( iy .eq. int(jsng/2.) ) y_ave = .TRUE.
      ix = int(jweg/2.)
      iy = int(jsng/2.)

! READ ALL DIMS FROM INPUT FILE AND CREATE DIMS FOR OUTPUT FILE
      IF (debug) THEN
        write(6,*) 
        write(6,*) "FILE dimensions:"
      ENDIF
      i_shape_we      = 0
      i_shape_sn      = 0
      i_shape_bt      = 0
      i_shape_we_stag = 0
      i_shape_sn_stag = 0
      i_shape_bt_stag = 0

      do i = 1, ndims
        status = nf_inq_dim(ncid, i, dname(i), dval(i))
        dval2(i) = dval(i)
!       CAUTION -- this stuff is hard-wired
        if (dname(i) .eq. 'west_east_stag') then
          dval2(i) = jweg
          i_shape_we_stag = i
        else if (dname(i) .eq. 'west_east') then
          dval2(i) = jweg-1
          i_shape_we = i
        else if (dname(i) .eq. 'south_north_stag') then
          dval2(i) = jsng
          i_shape_sn_stag = i
        else if (dname(i) .eq. 'south_north') then
          dval2(i) = jsng-1
          i_shape_sn = i
        else if (dname(i) .eq. 'bottom_top_stag') then
          dval2(i) = jbtg
          i_shape_bt_stag = i
        else if (dname(i) .eq. 'bottom_top') then
          dval2(i) = jbtg-1
          i_shape_bt = i
        endif
        if ( dname(i) == "Time" ) then
          status = nf_def_dim(mcid, dname(i), NF_UNLIMITED, i)
        else
          status = nf_def_dim(mcid, dname(i), dval2(i), i)
        end if
        IF (debug) THEN
          write(6,'(i4," : ",A," in = ",i4," (out = ",i4,")")') &
                i,dname(i),dval(i),dval2(i)
        ENDIF
      enddo
      IF (.not. debug) THEN
        write(6,*)
        write(6,*) "Set up file DIMENSIONS"
      ENDIF

! DEALING WITH THE GLOBAL ATTRIBUTES
      IF (debug) THEN
        write(6,*) 
        write(6,*) "FILE attributes:"
      ENDIF
      do i = 1, ngatts
        status = nf_inq_attname(ncid, nf_global, i, name)
        status = nf_inq_atttype(ncid, nf_global, name, itype)
        status = nf_inq_attlen(ncid, nf_global, name, ilen)

        if ( itype .eq. 2 ) then        ! characters
          status = nf_get_att_text (ncid, nf_global, name, cval)
          IF (debug) THEN
            write(6,'(i4," : ",A," in = ",A," (out = ",$)') &
                  i,name,cval(1:ilen)
          ENDIF
          if(name(1:5) .eq. 'TITLE') then
             cval = cval(1:ilen)//" : tcwrf"//option
             ilen = len_trim(cval)
          endif
          IF (debug) write(6,'(A,")")') cval(1:ilen)
          status = nf_put_att_text(mcid, nf_global, name, ilen,&
                    cval(1:ilen))

        elseif ( itype .eq. 4 ) then     ! integers
          status = nf_get_att_int (ncid, nf_global, name, ival)
          IF (debug) THEN
            write(6,'(i4," : ",A," in = ",i7," (out = ",$)') &
                  i,name,ival        
          ENDIF
          if(name .eq. 'WEST-EAST_GRID_DIMENSION') ival = jweg
          if(name .eq. 'SOUTH-NORTH_GRID_DIMENSION') ival = jsng
          if(name .eq. 'BOTTOM-TOP_GRID_DIMENSION') ival = jbtg
          IF (debug) write(6,'(i7,")")') ival
          status = nf_put_att_int(mcid, nf_global, name, itype,&
                    ilen, ival)

        elseif ( itype .eq. 5 ) then    ! real
          status = nf_get_att_real (ncid, nf_global, name, rval)
          IF (debug) THEN
            write(6,'(i4," : ",A," in = ",G18.10E2," (out = ",$)') &
                  i,name,rval        
          ENDIF
          if(name(1:2) .eq. 'DX' .or. name(1:2) .eq. 'DY') rval = okm
          IF (debug) write(6,'(G18.10E2,")")') rval
          status = nf_put_att_real(mcid, nf_global, name, itype,&
                    ilen, rval)
        endif
      enddo
      IF ( .not. debug ) THEN
        write(6,*) "Write file ATTRIBUTES"
        write(6,*) 
      ENDIF

! TRAIN FILE
      do i = 1, nvars
        status = nf_inq_var(ncid, i, cval, itype, idm, ishape, natt)
        ishape2 = ishape

        status = nf_def_var(mcid, cval, itype, idm, ishape2, i)
        do na = 1, natt
          status = nf_inq_attname(ncid, i, na, name)
          status = nf_copy_att(ncid, i, name, mcid, i)
        enddo
      enddo
      status = nf_enddef(mcid)

! Get variables for TC remove
! Get dimensions information

      do i = 1, nvars
        status = nf_inq_var(ncid, i, cval, itype, idm, ishape, natt)
        dims  = 1
        do ii = 1,idm
          dims(ii)  = dval(ishape(ii))
        enddo

        if ( trim(cval) .eq. 'PSFC' ) then
           nx=dims(1)
           ny=dims(2)
        elseif ( trim(cval) .eq. 'U' ) then
           nx_u=dims(1)
           nz=dims(3)
        elseif ( trim(cval) .eq. 'V' ) then
           ny_v=dims(2)
        elseif ( trim(cval) .eq. 'PH' ) then
           nz_p=dims(3)
        endif

      enddo ! End of i = 1, nvars

! ========================================================== !
! Get PSFC and LANDMASK to find TC center

    ALLOCATE(psfc_in(nx,ny),psfc_out(nx,ny))
    ALLOCATE(u_in(nx_u,ny,nz),u_out(nx_u,ny,nz))
    ALLOCATE(v_in(nx,ny_v,nz),v_out(nx,ny_v,nz))
    ALLOCATE(t_in(nx,ny,nz),t_out(nx,ny,nz))
    ALLOCATE(qv_in(nx,ny,nz),qv_out(nx,ny,nz))
    ALLOCATE(p_in(nx,ny,nz),p_out(nx,ny,nz))
    ALLOCATE(ph_in(nx,ny,nz_p),ph_out(nx,ny,nz_p))
    ALLOCATE(landmask(nx,ny))
!
    status = nf_inq_varid(ncid,'PSFC',i)
    status = nf_get_var_real(ncid, i, psfc_in)
    IF (debug_tc) write(6,*) '   SAMPLE VALUE = ',psfc_in(1,1)
!
    status = nf_inq_varid(ncid,'LANDMASK',i)
    status = nf_get_var_real(ncid, i, landmask)
! 
    status = nf_inq_varid(ncid,'U',i)
    status = nf_get_var_real(ncid, i, u_in)
!
    status = nf_inq_varid(ncid,'V',i)
    status = nf_get_var_real(ncid, i, v_in)
!
    status = nf_inq_varid(ncid,'T',i)
    status = nf_get_var_real(ncid, i, t_in)
!
    status = nf_inq_varid(ncid,'QVAPOR',i)
    status = nf_get_var_real(ncid, i, qv_in)
!
    status = nf_inq_varid(ncid,'P',i)
    status = nf_get_var_real(ncid, i, p_in)
!
    status = nf_inq_varid(ncid,'PH',i)
    status = nf_get_var_real(ncid, i, ph_in)

! Now, we remove TC

! For PSFC
    psfc_out(:,:)=psfc_in(:,:)
    CALL Smooth(psfc_in,landmask,nx,ny,icent,jcent,V_max,.FALSE.)
    CALL check_land_mask(landmask,nx,ny,psfc_in,psfc_out)
    deallocate(psfc_in)

! For U
    do i = 1,nz
       u_out(:,:,i)=u_in(:,:,i)
       if (i .eq. 1) then
         CALL Smooth(u_in(:,:,i),landmask,nx_u,ny,icent,jcent,V_max,.FALSE.)
         CALL check_land_mask(landmask,nx_u,ny,u_in(:,:,i),u_out(:,:,i))
       else
          CALL Smooth(u_out(:,:,i),landmask,nx_u,ny,icent,jcent,V_max,.TRUE.)
       endif
    enddo
    deallocate(u_in)

! For V
    do i = 1,nz
       v_out(:,:,i)=v_in(:,:,i)
       if (i .eq. 1) then
         CALL Smooth(v_in(:,:,i),landmask,nx,ny_v,icent,jcent,V_max,.FALSE.)
         CALL check_land_mask(landmask,nx,ny_v,v_in(:,:,i),v_out(:,:,i))
       else
          CALL Smooth(v_out(:,:,i),landmask,nx,ny_v,icent,jcent,V_max,.TRUE.)
       endif
    enddo
    deallocate(v_in)

    IF (debug_tc) write(6,*) '   SAMPLE VALUE OUT = ',v_out(1,1,1)
    write(6,*) '  Finish calulate'

! For T
    do i = 1,nz
       t_out(:,:,i)=t_in(:,:,i)
       if (i .eq. 1) then
         CALL Smooth(t_in(:,:,i),landmask,nx,ny,icent,jcent,V_max,.FALSE.)
         CALL check_land_mask(landmask,nx,ny,t_in(:,:,i),t_out(:,:,i))
       else
          CALL Smooth(t_out(:,:,i),landmask,nx,ny,icent,jcent,V_max,.TRUE.)
       endif
    enddo
    deallocate(t_in)

! For QV
    do i = 1,nz
       qv_out(:,:,i)=qv_in(:,:,i)
          CALL Smooth(qv_out(:,:,i),landmask,nx,ny,icent,jcent,V_max,.TRUE.)
    enddo
    deallocate(qv_in)

! For P
    do i = 1,nz
       p_out(:,:,i)=p_in(:,:,i)
       CALL Smooth(p_out(:,:,i),landmask,nx,ny,icent,jcent,V_max,.TRUE.)
    enddo
    deallocate(p_in)

! For PH
    do i = 1,nz_p
       ph_out(:,:,i)=ph_in(:,:,i)
       CALL Smooth(ph_out(:,:,i),landmask,nx,ny,icent,jcent,V_max,.TRUE.)
    enddo
    deallocate(ph_in)

    write(6,*) '  Finish calulate'
!
! LOOP THROUGH THE DATA 

      IF (debug) THEN
        write(6,*) 
        write(6,*) 
      ENDIF
      write(6,*) "Write file VARIABLES:"
      start_dims = 1
      do i = 1, nvars
        status = nf_inq_var(ncid, i, cval, itype, idm, ishape, natt)
        ishape2 = ishape
        IF (debug) THEN
          write(6,*) 
        ENDIF
        write(6,*) 'VARIABLE: ',trim(cval)

! GET THE DIMS FOR INPUT AND OUTPUT FROM THE SHAPE

        dims_in  = 1
        dims_out = 1
        do ii = 1,idm
          dims_in(ii)  = dval(ishape(ii))
          dims_out(ii) = dval2(ishape2(ii))
        enddo
        IF (debug) THEN
          write(6,*) '   DIMS  IN: ',dims_in
          write(6,*) '   DIMS OUT: ',dims_out
        ENDIF

        IF ( option == '-box' ) THEN
           !! Get the start and end dimensions of the box in the input file
           firstS  = 1
           firstE  = dims_out(1)
           secondS = 1
           secondE = dims_out(2)
           thirdS  = 1
           thirdE  = dims_out(3)

          if (idm.eq.2 .and. dims_out(1).ge.jbtg-1 .and. box_end(3).ne.0) then
            firstS = box_start(3)
            firstE = box_end(3)
            if (dims_out(3) .eq. jbtg-1) firstE = firstE-1
          endif
          if (idm .ge. 3) then
            if (box_end(1) .ne. 0) then
              firstS = box_start(1)
              firstE = box_end(1)
              if (dims_out(1) .eq. jweg-1) firstE = firstE-1
            endif
            if (box_end(2) .ne. 0) then
              secondS = box_start(2)
              secondE = box_end(2)
              if (dims_out(2) .eq. jsng-1) secondE = secondE-1
            endif
            if (idm == 4 .and. box_end(3).ne.0) then
              thirdS = box_start(3)
              thirdE = box_end(3)
              if (dims_out(3) .eq. jbtg-1) thirdE = thirdE-1
            endif
          endif
        ENDIF

! ALLOCATE THE INPUT AND OUTPUT ARRAYS
! READ THE DATA FROM INPUT FILE
! THIN THE GRID IF NEEDED, OR GET THE CORRECT BOX

        IF     (itype .eq. 2) THEN          ! character
          allocate (text(dims_in(1), dims_in(2), dims_in(3), &
                         dims_in(4)))
          status = nf_get_var_text(ncid, i, text)
          IF (debug) write(6,*) '   SAMPLE VALUE = ',text(:,1,1,1)
          status = nf_put_vara_text (mcid, i, start_dims, dims_in, text)
          deallocate (text)

        ELSEIF (itype .eq. 4) THEN          ! integer
          allocate (idata(dims_in(1), dims_in(2), dims_in(3), &
                          dims_in(4)))
          allocate(idata2(dims_out(1),dims_out(2),dims_out(3),&
                          dims_out(4)))
          status = nf_get_var_int(ncid, i, idata)
          IF (debug) write(6,*) '   SAMPLE VALUE = ',idata(1,1,1,1)
          IF ( option == '-box') THEN
            IF (debug) write(6,*) '   a BOX is extracted from the input domain '  
            idata2 = idata(firstS:firstE,secondS:secondE,thirdS:thirdE,:)
          ENDIF
          status = nf_put_vara_int (mcid, i, start_dims, dims_out, idata2)
          deallocate (idata)
          deallocate (idata2)

        ELSEIF (itype .eq. 5) THEN          ! real



          allocate (data(dims_in(1), dims_in(2), dims_in(3), &
                         dims_in(4)))
          allocate(data2(dims_out(1),dims_out(2),dims_out(3), &
                         dims_out(4)))
          status = nf_get_var_real(ncid, i, data)
          IF (debug) write(6,*) '   SAMPLE VALUE = ',data(1,1,1,1)
          IF ( option == '-box') THEN
            IF (debug) write(6,*) '   a BOX is extracted from the input domain '  
            data2 = data(firstS:firstE,secondS:secondE,thirdS:thirdE,:)

            IF ( trim(cval) .eq. 'PSFC' ) THEN
              write(6,*) 'PSFC new'
              data2(:,:,1,1) = psfc_out(firstS:firstE,secondS:secondE)

            ELSEIF ( trim(cval) .eq. 'U' ) THEN
              write(6,*) 'U new'
              data2(:,:,:,1) = u_out(firstS:firstE,secondS:secondE,thirdS:thirdE)

            ELSEIF ( trim(cval) .eq. 'V' ) THEN
              write(6,*) 'V new'
              data2(:,:,:,1) = v_out(firstS:firstE,secondS:secondE,thirdS:thirdE)

            ELSEIF ( trim(cval) .eq. 'T' ) THEN
              write(6,*) 'T new'
              data2(:,:,:,1) = t_out(firstS:firstE,secondS:secondE,thirdS:thirdE)

            ELSEIF ( trim(cval) .eq. 'QVAPOR' ) THEN
              write(6,*) 'QV new'
              data2(:,:,:,1) = qv_out(firstS:firstE,secondS:secondE,thirdS:thirdE)

            ELSEIF ( trim(cval) .eq. 'PH' ) THEN
              write(6,*) 'PH new'
              data2(:,:,:,1) = ph_out(firstS:firstE,secondS:secondE,thirdS:thirdE)

            ELSEIF ( trim(cval) .eq. 'P' ) THEN
              write(6,*) 'P new'
              data2(:,:,:,1) = p_out(firstS:firstE,secondS:secondE,thirdS:thirdE)

            ENDIF
          ENDIF
          status = nf_put_vara_real (mcid, i, start_dims, dims_out, data2)



          IF ( cval == 'XLAT' .or. cval == 'XLONG' ) THEN
!            We need fix the box's center long and lat
             new_cen = data2(ix,iy,1,1)
                 if ( x_ave .and. y_ave ) then
               new_cen = (data2(ix,  iy,1,1)+data2(ix  ,iy+1,1,1)+  &
                          data2(ix+1,iy,1,1)+data2(ix+1,iy+1,1,1))/4.
             elseif ( x_ave .and. .not. y_ave ) then
               new_cen = (data2(ix,  iy,1,1)+data2(ix+1,iy  ,1,1))/2.
             elseif ( .not. x_ave .and. y_ave ) then
               new_cen = (data2(ix,  iy,1,1)+data2(ix  ,iy+1,1,1))/2.
             endif
          ENDIF
          IF ( cval == 'XLAT' ) THEN
             IF (debug) write(6,*) '   Fix global attribute CEN_LAT: now = ', new_cen
             status = nf_inq_atttype(ncid, nf_global, 'CEN_LAT', itype)
             status = nf_inq_attlen(ncid, nf_global, 'CEN_LAT', ilen)
             status = nf_put_att_real(mcid, nf_global, 'CEN_LAT', itype,&
                       ilen, new_cen)
          ELSEIF ( cval == 'XLONG' ) THEN
             IF (debug) write(6,*) '   Fix global attribute CEN_LON: now = ', new_cen
             status = nf_inq_atttype(ncid, nf_global, 'CEN_LON', itype)
             status = nf_inq_attlen(ncid, nf_global, 'CEN_LON', ilen)
             status = nf_put_att_real(mcid, nf_global, 'CEN_LON', itype,&
                       ilen, new_cen)
          ENDIF
          deallocate (data)
          deallocate (data2)

        ELSEIF (itype .eq. 6) THEN          ! double
          allocate (ddata(dims_in(1), dims_in(2), dims_in(3), &
                          dims_in(4)))
          allocate(ddata2(dims_out(1),dims_out(2),dims_out(3),&
                          dims_out(4)))
          status = nf_get_var_double(ncid, i, ddata)
          IF (debug) write(6,*) '   SAMPLE VALUE = ',ddata(1,1,1,1)
          IF ( option == '-box') THEN
            IF (debug) write(6,*) '   a BOX is extracted from the input domain '  
            ddata2 = ddata(firstS:firstE,secondS:secondE,thirdS:thirdE,:)
          ENDIF
          status = nf_put_vara_double (mcid, i, start_dims, dims_out, ddata2)
          deallocate (ddata)
          deallocate (ddata2)
        ELSE
            stop 'trouble - do not know the variable type'
        ENDIF

      ENDDO     ! END OF VARIABLE LOOP
      status = nf_close(mcid)

      write(6,*) 
      write(6,*) "SUCCESS - we are out of here"      
      write(6,*) "#########################################"

      end program tcrmwrf
!---------------------------------------------------------------------
      subroutine handle_err(status)
      integer status
      write(6,*) 'Error number ',status
      stop
      end subroutine
!---------------------------------------------------------------------

  subroutine read_args(input_file,option,box_icent,box_jcent,Vmax,stid,debug)

  implicit none
  character (len=80)    :: input_file
  character (len=10)    :: option

  real               :: box_icent, box_jcent,idummy1,stid
  logical               :: debug

  integer               :: numarg, i, idummy
  real                  :: rdummy, Vmax
  integer, external     :: iargc
  character (len=80)    :: dummy, dir

! set up some defaults first
  input_file = " "
  option = " "
  idummy1 = 0.0
  box_icent = 0.0
  box_jcent   = 0.0
  Vmax = 0.0
  stid = 0.0
  numarg = iargc()
  i = 1

  if (numarg .lt. 1) call help_info

  do while (i <= numarg)
    call getarg(i,dummy)

    if (dummy(1:1) == "-") then    ! We have an option, else it is the filename

      SELECTCASE (trim(dummy))
          CASE ("-help")
               call help_info
          CASE ("-h")
               call help_info
          CASE ("-debug")
               debug = .TRUE.     
          CASE ("-box")
               option = dummy
               DO 
                 i = i+1
                 call getarg(i,dir)
                 if (dir(1:1) == '-') then
                   i=i-1
                   exit
                 endif
                 if (dir.ne.'lon') then
                 if (dir.ne.'lat') then
                 if (dir.ne.'vmax') then
                 if (dir.ne.'stormid') then
                   i=i-1
                   exit
                 endif
                 endif
                 endif
                 endif
                 i = i+1
                 call getarg(i,dummy)
                 read(dummy,*)idummy1
                 if (dir.eq.' '.or.idummy1.eq.0.0) exit
                 if ( dir == 'lon' ) then
                   box_icent = idummy1
                 endif
                 if ( dir == 'lat' ) then
                   box_jcent = idummy1
                 endif
                 if ( dir == 'vmax' ) then
                   Vmax = idummy1
                 endif
                 if ( dir == 'stormid' ) then
                   stid = idummy1
                 endif
                 idummy1 = 0.0
               ENDDO
          CASE DEFAULT
               call help_info
      END SELECT
    else
      input_file = dummy
    endif

      i = i+1

  enddo

  if (input_file == " ") call help_info

  end subroutine read_args
!------------------------------------------------------------------------------

  subroutine help_info

  print*," "
  print*," tcrmwrf   wrf_data_file_name  [-options] "
  print*," "
  print*," Current options available are:"
  print*," -help     : Print this information"                           
  print*," -h        : Print this information"                           
  print*," "
  print*," -box [ ]  : Will extract a box out of the input grid"
  print*,"             The box can have values for x/y/z  "
  print*,"             Examples:"
  print*,"                -box x 10 30 y 20 40 z 5 15" 
  print*,"                -box x 10 30 y 20 40 " 
  print*,"                -box x 10 30 " 
  print*,"                -box y 20 40 " 
  print*," "
  end subroutine help_info

!------------------------------------------------------------------------------

  subroutine Smooth(b_in,l_mask,b1,b2,b_ic,b_jc,v_max,use_land)
  implicit none 
  integer b1, b2
  real,dimension(b1,b2) :: b_in, l_mask
  real v_max
  integer i, j,i_n,j_n, ii,jj,k
  integer i1,i2,j1,j2,i3,i4,j3,j4,ncount
  integer rad, box_rad, max_loop, check
  integer b_jc,b_ic
  real Var1,Var2,Var3,Var4,delta_Var, delta_Var2,delta_crt
  logical keep_smooth
  logical use_land
!!
  rad = 60 !radius from TCs center
  box_rad=40 ! 
! 
  if (v_max .lt. 30.) then
      max_loop = 20
  else if (v_max .ge. 30.) then
      max_loop = 30
  endif
!
  i1 = max(b_ic-rad,1+1)
  i2 = min(b_ic+rad,b1-1)
  j1 = max(b_jc-rad,1+1)
  j2 = min(b_jc+rad,b2-1)
! Calculate value of each point in a box (i1:i2,j1:j2) by average value of a box with box_rad=40 to remove vortex
  do  k = 1,max_loop
    do ii=i1,i2
      do jj=j1,j2
        call box_average(b_in,l_mask,b1,b2,ii,jj,box_rad,use_land)
      end do
    end do
  end do
! ==========================
!Smooth value of each point in a box by average value of 9 points around 
  do k =1,max_loop
    do ii=i1,i2
      do jj=j1,j2
        if (use_land) then
          b_in(ii,jj)=(b_in(ii,jj)+b_in(ii-1,jj)+b_in(ii+1,jj)+b_in(ii-1,jj-1)+b_in(ii+1,jj+1)+b_in(ii,jj+1)+b_in(ii,jj-1)+b_in(ii-1,jj+1)+b_in(ii+1,jj-1))/9.
        else
            if (l_mask(ii,jj) .eq. 0 ) then
              b_in(ii,jj)=(b_in(ii,jj)+b_in(ii-1,jj)+b_in(ii+1,jj)+b_in(ii-1,jj-1)+b_in(ii+1,jj+1)+b_in(ii,jj+1)+b_in(ii,jj-1)+b_in(ii-1,jj+1)+b_in(ii+1,jj-1))/9.
            endif
        endif
      end do
    end do
  end do

  end subroutine Smooth

!================================================================================

  subroutine box_average(d_in,ll_mask,d1,d2,il,jl,boxhalf,useland)
  
  implicit none 
  integer d1, d2, il, jl, boxhalf
  real,dimension(d1,d2) :: d_in, ll_mask
  real tmp_s
  integer i3,i4,j3,j4,ncount,i_n,j_n
  logical useland

  i3=max(il-boxhalf,1)
  i4=min(il+boxhalf,d1)
  j3=max(jl-boxhalf,1)
  j4=min(jl+boxhalf,d2)
           
  tmp_s=d_in(il,jl)
!
  ncount=0
  do i_n=i3,i4
     do j_n=j3,j4
       if (useland) then
           ncount=ncount+1
           tmp_s = tmp_s + d_in(i_n,j_n)
       else
          if (ll_mask(i_n,j_n) .eq. 0 ) then
             ncount=ncount+1
             tmp_s = tmp_s + d_in(i_n,j_n)
          endif
       endif
     end do
  end do
!
  d_in(il,jl)=tmp_s/float(ncount)
!
  end subroutine box_average
!
!=======================================================================

  subroutine check_land_mask(lmask,lm1,lm2,var_in,var_out)
  implicit none
  integer lm1, lm2
  real,dimension(lm1,lm2) :: lmask,var_in,var_out
  integer i, j
  !
  do i=1,lm1
     do j=1,lm2
          if ( lmask(i,j) .eq. 0 ) then
             var_out(i,j)=var_in(i,j)
          endif
       enddo
  enddo
  end subroutine check_land_mask

!=======================================================================!
