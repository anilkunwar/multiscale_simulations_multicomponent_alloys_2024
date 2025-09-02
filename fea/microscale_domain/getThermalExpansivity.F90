    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Thermal Expansivity of solid AlSiMgZr fitted as a function of temperature
    ! (alpha_alloy)solid = As*(T)**2 +Bs*T+  Cs, where As = -3.7262790630140875e-10 K^{-3}, Bs = 4.2255653792479394e-07 K^{-2}  and Cs = -8.371435292934836e-05 K
    ! 298 < T < Tm where Tm = 842.710 K
    ! Reference: Jiang et al., Mater. Research Letters, 6 (2019) 125806.
    ! https://iopscience.iop.org/article/10.1088/2053-1591/ab5b5b
    ! Thermal Expansivity of liquid AlSiMgZr fitted as a function of temperature
    ! (alpha_alloy)liquid  = 0 
    ! The expansivity of liquid is not zero. However, in order to render zero thermal stress in molten liquid, this value is considered to be zero
    !-----------------------------------------------------
    FUNCTION getThermalExpansivity( model, n, temp ) RESULT(expansivity)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, expansivity, tscaler

    ! variables needed inside function
    !REAL(KIND=dp) :: refSolExpansivity, refLiqExpansivity, refTemp,  &
    !alphas, alphal, betas, betal
    REAL(KIND=dp) :: refSolExpansivity, refTemp,  &
    alphas, betas
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getThermalExpansivity', 'No material found')
    END IF

    ! read in reference conductivity at reference temperature
    refSolExpansivity = GetConstReal( material, 'Reference Thermal Expansivity Solid AlMgSiZr',GotIt)
    !refDenst = GetConstReal( material, 'Solid_ti_rho_constant',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'Reference Thermal Expansivity Solid AlMgSiZr not found')
    END IF

    ! read in Temperature Coefficient of Resistance
    alphas = GetConstReal( material, 'Thermal Expansivity Coeff As Solid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'T2 term of Thermal Expansivity-temperature curve solid not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance
    betas = GetConstReal( material, 'Thermal Expansivity Coeff Bs Solid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'Slope  term of solid AlMgSiZr not found')
    END IF
    
    ! read in  Ds in Ds*ln(T) term
    !deltas = GetConstReal( material, 'Thermal Expansivity Coeff Ds Liquid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getThermalExpansivity', 'Coefficient of logT term liquid AlMgSiZr not found')
    !END IF
    
    ! read in reference Thermal Expansivity at reference temperature
    !refLiqExpansivity = GetConstReal( material, 'Reference Thermal Expansivity Liquid AlMgSiZr',GotIt)
    !refDenst = GetConstReal( material, 'Solid_ti_rho_constant',GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getThermalExpansivity', 'Reference Thermal Expansivity Solid AlMgSiZr not found')
    !END IF
    
    ! read in pseudo reference conductivity at reference temperature of liquid
    !alphal = GetConstReal( material, 'Thermal Expansivity Coefficient Liquid AlMgSiZr',GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getThermalExpansivity', 'Thermal Expansivity Coefficient Al of Liquid AlMgSiZr not found')
    !END IF
    
    ! read in Temperature Coefficient of Resistance
    !betal = GetConstReal( material, 'Thermal Expansivity Coeff Bl Liquid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getThermalExpansivity', 'Slope  term of liquid AlMgSiZr not found')
    !END IF

    ! read in reference temperature
    refTemp = GetConstReal( material, 'Melting Point Temperature of AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'Reference Melting Temperature of AlMgSiZr not found')
    END IF
    
    ! read in the temperature scaling factor
    tscaler = GetConstReal( material, 'Tscaler', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Scaling Factor for T not found')
    END IF

    ! compute Thermal Expansivity conductivity
    IF (refTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getThermalExpansivity', 'The AlMgSiZr material is in liquid state.')
            !CALL Warn('getThermalExpansivity', 'Using Thermal Expansivity reference value')
    !denst = 1.11*(refDenst + alpha*(temp))
    !denst = refLiqExpansivity + alphal*((tscaler)*(temp-910))
    expansivity = 0.0
    ELSE
    ! To add the expansivity
    expansivity = refSolExpansivity + alphas*((tscaler)*(temp))**2 + betas*(tscaler)*(temp)
    END IF

    END FUNCTION getThermalExpansivity

