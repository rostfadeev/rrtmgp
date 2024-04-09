

      Module mod_rrtmgp
        Use netcdf
        
        Implicit None
        
        Real(Kind=8) :: PI
        
      CONTAINS
        
        Subroutine rrtmgp()
          Implicit None
          Integer(Kind=4) ncid, did, nlon,nlat,nlev, varid
          Character(300) in_fn, varname
          Logical stat
          Real(Kind=4), Allocatable :: lon(:),lat(:),lev(:)
          Real(Kind=4), Allocatable :: r4(:), r42(:,:), r43(:,:,:)
          
          in_fn = "zz_diag_eb-latm.nc"
          
          call check( "Open file "//trim(in_fn), nf90_open (path = trim(in_fn), mode = IOR(NF90_NOWRITE,NF90_NETCDF4), ncid = ncid) )
          
          varname = "column"
          call check( "looking for "//trim(varname), nf90_inq_dimid        (ncid = ncid, name = trim(varname), dimid = did) )
          call check( "get NLON"        , nf90_inquire_dimension(ncid = ncid, dimid = did , len = nlon  ) )
          
          varname = "latm"
          call check( "looking for "//trim(varname), nf90_inq_dimid        (ncid = ncid, name = trim(varname), dimid = did) )
          call check( "get NLAT"        , nf90_inquire_dimension(ncid = ncid, dimid = did , len = nlat  ) )
          
          varname = "level"
          call check( "looking for "//trim(varname), nf90_inq_dimid        (ncid = ncid, name = trim(varname), dimid = did) )
          call check( "get NLEV"        , nf90_inquire_dimension(ncid = ncid, dimid = did , len = nlev  ) )
          
          print *,"nlon = ",nlon
          print *,"nlat = ",nlat
          print *,"nlev = ",nlev
          
          Allocate(lon(nlon), lat(nlat), lev(nlev))
          
          varname = "lon"; print *,trim(varname)
          call check( "looking for "//trim(varname), nf90_inq_varid (ncid, trim(varname), varid))
          call check( "get "        //trim(varname), nf90_get_var   (ncid, varid, lon, (/ 1 /), (/ nlon /)) )
          write(*,'(1000(f6.1))') lon
          
          varname = "latm"; print *,trim(varname)
          call check( "looking for "//trim(varname), nf90_inq_varid (ncid, trim(varname), varid))
          call check( "get "        //trim(varname), nf90_get_var   (ncid, varid, lat, (/ 1 /), (/ nlat /)) )
          write(*,'(1000(f6.1))') lat
          
          varname = "lev"; print *,trim(varname)
          call check( "looking for "//trim(varname), nf90_inq_varid (ncid, trim(varname), varid))
          call check( "get "        //trim(varname), nf90_get_var   (ncid, varid, lev, (/ 1 /), (/ nlev /)) )
          write(*,'(1000(f6.3))') lev
          
          !==========================
          
          Allocate(r42(nlon,nlat)); r42 = 0.
          Allocate(r43(nlev,nlon,nlat)); r43 = 0.
          
          varname = "skin_temperature"
          call check( "looking for "//trim(varname), nf90_inq_varid (ncid, trim(varname), varid))
          call check( "get "        //trim(varname), nf90_get_var   (ncid, varid, r42, (/ 1, 1 /), (/ nlon, nlat /)) )
          print *,trim(varname),minval(r42),maxval(r42)
          
          
          varname = "pressure_hl"
          call check( "looking for "//trim(varname), nf90_inq_varid (ncid, trim(varname), varid))
          call check( "get "        //trim(varname), nf90_get_var   (ncid, varid, r43, (/ 1, 1, 1 /), (/ nlev, nlon, nlat /)) )
          print *,trim(varname),minval(r43),maxval(r43)
          
          
          
          call check("Close file ", nf90_close (ncid))
          
          if (Allocated(lon)) Deallocate(lon)
          if (Allocated(lat)) Deallocate(lat)
          if (Allocated(lev)) Deallocate(lev)
        End Subroutine
        
        Subroutine check(code, stat)
          Implicit None
          Character(*),    Intent(In) :: code
          Integer(Kind=4), Intent(In) :: stat
          if (stat /= NF90_NOERR) then
            write(*,"(3x,2x,a,a,2x,a,2x,a)") "var (code): ",code,"mess: ",trim(NF90_STRERROR(stat))
            STOP
          end if
        End Subroutine
        
      End Module




