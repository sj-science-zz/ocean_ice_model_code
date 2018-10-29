/*
*************************************************** Hernan G. Arango ***
**Copyright (c) 2002-2015 The ROMS/TOMS Group Benjamin K. Galton-Fenzi**
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines iceshelf model input parameters in output NetCDF files.   **
**  It is included in routine "def_info.F".                           **
**                                                                    **
************************************************************************
*/

!
!  Define iceshelf model parameters.
!
      Vinfo( 1)='zice'
      Vinfo( 2)='ice draft'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               2, t2dgrd, Aval, Vinfo, ncname,                    &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
