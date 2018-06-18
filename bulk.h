!-------------------------------------------------------------------------------
! BULK configuration and constants
!-------------------------------------------------------------------------------
      logical
     &  calc_swr        ! Calculate shortwave radiation
     &, use_coare       ! Use COARE algo to estimate heat fluxes

      integer(kind=1)
     &  lwrad_formula   ! Formula to calculate net longwave radiation

      common/bulk_flags/
     &  calc_swr, use_coare, lwrad_formula

!
! Named constants
!
      integer(kind=1) lwBIGNAMI, lwMAY, lwHERZFELD, lwBERLIAND
      parameter ( lwBIGNAMI  = 0,
     &            lwMAY      = 1,
     &            lwHERZFELD = 2,
     &            lwBERLIAND = 3 )