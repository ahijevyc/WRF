!WRF:MODEL_LAYER:PHYSICS
!
MODULE MODULE_LTNG_FED

! USE WRF PHYSICS CONSTANTS
      USE MODULE_MODEL_CONSTANTS, ONLY:                           &
                                        PICONST                   &
                                      , GRAV => G

! USE GAMMA FUNCTION FROM THOMPSON MICROPHYSICS
      USE MODULE_MP_THOMPSON,     ONLY: WGAMMA

CONTAINS
!===================================================================
!
      SUBROUTINE CALCFED(                                         &
                         QV,QC,QR,QI,QS,QG,QH                     &
                        ,W,Z,DZ8W,PI_PHY,TH_PHY,P_PHY,RHO_PHY     & ! DZ8W is thickness (meters) of each model layer
                        ,DX,DY                                    &
                        ,IDS,IDE, JDS,JDE, KDS,KDE                & ! domain
                        ,IMS,IME, JMS,JME, KMS,KME                & ! memory
                        ,ITS,ITE, JTS,JTE, KTS,KTE                & ! tile
                        ,FED_rate                                 &
                        ,FED1_rate                                & 
                        ,FED2_rate                                &
                        ,FED3_rate                                &
                        ,FED4_rate                                &
                        ,FED_SNOW                                 &
                        ,FED_ICE                                  &
                        ,FED_GRAUP                                &
                        ,FED_IWP                                  &
                        ,FED_WVOL                                 &
                        ,FED_WVOL2                                &
                        ,FED_ICE_MASS_FLUX                        &
                        ,FED_PRECIP_MASS_FLUX                     &
                                                                  )

      USE MODULE_WRF_ERROR
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
!-------------------------------------------------------------------

      INTEGER, INTENT(IN) ::          IDS,IDE, JDS,JDE, KDS,KDE , &
                                      IMS,IME, JMS,JME, KMS,KME , &
                                      ITS,ITE, JTS,JTE, KTS,KTE
      REAL, DIMENSION( IMS:IME , KMS:KME , JMS:JME ),             &
            INTENT(IN) ::                                         &
                                                             QV , &
                                                             QC , &
                                                             QR , &
                                                             QI , &
                                                             QS , &
                                                             QG , &
                                                             QH

      REAL, DIMENSION( IMS:IME, KMS:KME, JMS:JME ),               &
            INTENT(IN ) :: W,Z,DZ8W,PI_PHY,TH_PHY,P_PHY,RHO_PHY

      REAL, INTENT(IN ) :: DX, DY ! GRID SPACING IN METERS

      REAL, INTENT(OUT), DIMENSION(IMS:IME,JMS:JME)::  FED_rate   &
                        ,FED1_rate                                & 
                        ,FED2_rate                                &
                        ,FED3_rate                                &
                        ,FED4_rate                                &
                        ,FED_SNOW                                 &
                        ,FED_ICE                                  &
                        ,FED_GRAUP                                &
                        ,FED_IWP                                  &
                        ,FED_WVOL                                 &
                        ,FED_WVOL2                                &
                        ,FED_ICE_MASS_FLUX                        &
                        ,FED_PRECIP_MASS_FLUX                     


! DEFINE LOCAL DATA
      REAL, DIMENSION(KTS:KTE):: TEMPC
      REAL :: AIR_VOL, AIR_MASS, GRAUP_MASS, ICE_MASS, SNOW_MASS
      REAL :: TMP1, TMP2, VTG, VTS
      REAL :: ICE_MASS_FLUX, PRECIP_MASS_FLUX
      REAL :: IWP
      REAL :: FLUXPROD, LIGHTFLUX, LIGHTIWP, LIGHTPRECIP_ICE
      REAL :: LIGHT_GRAUPEL_VOL, LIGHT_WMAX
      REAL :: W_VOL, W_VOL2, LIGHTWVOL, LIGHTWVOL2

      INTEGER         ::            I,J,K

      REAL, PARAMETER :: TMAXC    =   0.    ! TEMPERATURE THRESHOLD IN CELSIUS
      REAL, PARAMETER :: WMIN     =   5.    ! VERTICAL VEL THRESHOLD (M/S)
      REAL, PARAMETER :: WMIN2    =   3.    ! VERTICAL VEL THRESHOLD (M/S)
      REAL, PARAMETER :: CELS2K   = 273.15
      REAL, PARAMETER :: XNOG     =   4.0E4 !LIN SNOW FALL SPEED
      REAL, PARAMETER :: CDRAG    =   0.6
      REAL, PARAMETER :: RHOGRAUL = 917. ! DENSITY OF GRAUPEL KG/M^3 ! ONLY 500 IN MODULE_MP_THOMPSON.F
      REAL, PARAMETER :: XNOS     =   3.0E6
      REAL, PARAMETER :: RHOSNOW  = 100.! KG/M3 DENSITY OF SNOW
      REAL, PARAMETER :: CONSTD   =   0.25
      REAL, PARAMETER :: EXPON    = 1-CONSTD
      REAL, PARAMETER :: CONSTC   = 152.93*0.01**EXPON
      REAL, PARAMETER :: GAMDP4   = GAMMA(CONSTD+4.)
      REAL, PARAMETER :: GAM4PT5  = GAMMA(4.5)

      CHARACTER*256 :: MESSAGE
!-------------------------------------------------------------------

      CALL WRF_DEBUG(150, ' LTNG_FED: ENTERED CALCFED')

      DO J = JTS,JTE
         DO I = ITS,ITE

            ! INITIALIZE VERTICAL TOTALS TO ZERO
            AIR_MASS         = 0.
            AIR_VOL          = 0.
            GRAUP_MASS       = 0.
            ICE_MASS         = 0.
            SNOW_MASS        = 0.
            PRECIP_MASS_FLUX = 0.
            ICE_MASS_FLUX    = 0. ! UPWARD
            W_VOL            = 0.
            W_VOL2           = 0.

            ! CALCULATE MASS AND FLUX WITH A SINGLE LOOP THRU THE VERTICAL
            DO K = KTS,KTE
               WRITE(MESSAGE,'(" LTNG_FED: I,J,K=",I8,I8,I8)') I,J,K
               CALL WRF_DEBUG(600, MESSAGE)

               ! TEMPERATURE OF LAYER (CELSIUS)
               TEMPC(K)     = TH_PHY(I,K,J)*PI_PHY(I,K,J)-CELS2K

               ! IF TEMPERATURE <= TMAXC
               ! ADD MASS OF AIR IN GRID VOLUME TO AIR_MASS
               ! ADD MASS OF GRAUPEL IN GRID VOLUME
               ! ADD MASS OF ICE IN GRID VOLUME
               ! ADD MASS OF SNOW IN GRID VOLUME
               ! ADD DOWNWARD PRECIP MASS FLUX
               ! ADD UPWARD ICE MASS FLUX 
               IF(TEMPC(K).LE.TMAXC)THEN 
                   
                  ! COULD FORCE POSITIVE MIXING RATIOS HERE, BUT NEGATIVE VALUES ARE SMALL?

                  AIR_VOL    =      DX * DY * DZ8W(I,K,J)
                  AIR_MASS   = AIR_MASS   + (AIR_VOL * RHO_PHY(I,K,J)) 
                  GRAUP_MASS = GRAUP_MASS + (QG(I,K,J) * AIR_MASS)
                  ICE_MASS   = ICE_MASS   + (QI(I,K,J) * AIR_MASS)
                  SNOW_MASS  = SNOW_MASS  + (QS(I,K,J) * AIR_MASS)

                  ! ADD TO PRECIP MASS FLUX IF QG>1E-8
                  IF(QG(I,K,J).GT.1E-8)THEN

                     TMP1 = 1./(PICONST*RHOGRAUL*XNOG/(RHO_PHY(I,K,J)*QG(I,K,J)))**0.25
                     TMP2 = SQRT(4.*0.33334*GRAV*RHOGRAUL/(CDRAG*RHO_PHY(I,K,J)))
                     VTG = (1./6.)*GAM4PT5*TMP2*SQRT(TMP1)

                     PRECIP_MASS_FLUX = PRECIP_MASS_FLUX + ( GRAUP_MASS * VTG ) ! IF YOU SUBTRACT W FROM VTG, YOU MAY GET NEGATIVE FLUX
                  END IF

                  ! IF QI>1E-8 AND W>0 ADD ICE TO ICE MASS FLUX
                  IF(QI(I,K,J).GT.1E-8 .AND. W(I,K,J).GT.0) THEN
                     ICE_MASS_FLUX = ICE_MASS_FLUX + ICE_MASS * W(I,K,J)
                  END IF

                  ! IF QS>1E-8 ADD SNOW TO EITHER ICE MASS FLUX OR PRECIP MASS FLUX
                  ! DEPENDING ON FALL SPEED AND VERTICAL MOTION
                  IF(QS(I,K,J).GT.1E-8)THEN

                     TMP1 = (PICONST*RHOSNOW*XNOS/(RHO_PHY(I,K,J)*QS(I,K,J)))**0.25
                     TMP2 = SQRT(1.29/RHO_PHY(I,K,J))
                     VTS = (1./6.) * CONSTC * GAMDP4 * TMP2 / (TMP1**CONSTD)

                     ! IF SNOW IS BEING LOFTED UPWARD, ADD SNOW TO ICE MASS FLUX
                     IF(W(I,K,J).GT.VTS) THEN
                        ICE_MASS_FLUX = ICE_MASS_FLUX + SNOW_MASS * (W(I,K,J)-VTS) 
                     ELSE
                        ! ELSE IF SNOW IS FALLING, ADD SNOW TO PRECIP MASS FLUX
                        PRECIP_MASS_FLUX = PRECIP_MASS_FLUX +  SNOW_MASS * VTS 
                     END IF
                   END IF

                   ! ADD TO UPDRAFT VOLUME IF W EXCEEDS W THRESHOLD
                   IF (W(I,K,J).GT.WMIN)  W_VOL  = W_VOL  + AIR_VOL
                   IF (W(I,K,J).GT.WMIN2) W_VOL2 = W_VOL2 + AIR_VOL

               END IF
            END DO

            IWP       = SNOW_MASS + ICE_MASS + GRAUP_MASS

            FLUXPROD  = ICE_MASS_FLUX * PRECIP_MASS_FLUX


            ! LIGHTNING FLASH EXTENT DENSITY ESTIMATES
            LIGHTFLUX       = FLUXPROD * 9E-17
            LIGHTIWP        =     IWP  * 0.20 / 2000**2 * DX * DY ! McCaul et al. 2009
            LIGHTPRECIP_ICE =     3.4E-8 * GRAUP_MASS - 18.1     ! flashes per minute; Deierling et al. 2008
            LIGHTWVOL       =    W_VOL * 6.75E-10                ! Deierling and Petersen 2008
            LIGHTWVOL2      =   W_VOL2 * 6.75E-10                ! Deierling and Petersen 2008


            FED_rate(I,J)             = LIGHTIWP
            FED1_rate(I,J)            = LIGHTWVOL
            FED2_rate(I,J)            = LIGHTWVOL2
            FED3_rate(I,J)            = LIGHTPRECIP_ICE
            FED4_rate(I,J)            = LIGHTFLUX
            FED_SNOW(I,J)             = SNOW_MASS
            FED_ICE(I,J)              = ICE_MASS
            FED_GRAUP(I,J)            = GRAUP_MASS
            FED_IWP(I,J)              = IWP
            FED_WVOL(I,J)             = W_VOL
            FED_WVOL2(I,J)            = W_VOL2
            FED_ICE_MASS_FLUX(I,J)    = ICE_MASS_FLUX
            FED_PRECIP_MASS_FLUX(I,J) = PRECIP_MASS_FLUX


            ! TODO: Convert other rates to flashes per second
            ! flash rates multiplied by lightning_dt (a time interval in seconds) in 
            ! module_lightning_driver.F.  Output needs to be in flashes per second.
            FED3_rate(I,J)  = FED3_rate(I,J) / 60. ! Convert from fl/min to fl/sec



         END DO
      END DO
      RETURN
      END SUBROUTINE CALCFED

  END MODULE MODULE_LTNG_FED
