    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Density of solid Ti6Al4V fitted as a function of temperature
    ! (rho_alloy)solid = As*(T-298) +  Cs, where As = -0.1898 kg/m3K  and Cs = 2672.85 kg/m3
    ! 298 < T < Tm where Tm = 842.710 K
    ! Reference: Schmitz et al., Journal of Materials Science, 47 (2012) 3706-3712.
    ! https://link.springer.com/article/10.1007/s10853-011-6219-8
    ! Density of liquid Ti6Al4V fitted as a function of temperature
    !(rho_alloy)liquid  = Al*(T-910) + Cl, where Al = -0.3153 kg/m3K  and Cl = 2409.200 kg/m3 (900.0 K < T < 2500.0 K)
    ! Reference: Schmon et al, EPJ Web of Conferences, 151 (2017) 04003.
    ! Reference: Schmitz et al., Journal of Materials Science, 47 (2012) 3706-3712.
    ! https://link.springer.com/article/10.1007/s10853-011-6219-8
    !-----------------------------------------------------
    FUNCTION getDensity( model, n, temp ) RESULT(denst)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, denst, tscaler

    ! variables needed inside function
    REAL(KIND=dp) :: refSolDenst, refLiqDenst, refTemp,  &
    alphas, alphal
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getDensity', 'No material found')
    END IF

    ! read in reference conductivity at reference temperature
    refSolDenst = GetConstReal( material, 'Reference Density Solid AlMgSiZr',GotIt)
    !refDenst = GetConstReal( material, 'Solid_ti_rho_constant',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Reference Density Solid AlMgSiZr not found')
    END IF

    ! read in Temperature Coefficient of Resistance
    alphas = GetConstReal( material, 'Density Coeff Solid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'slope of Density-temperature curve solid not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance
    !betas = GetConstReal( material, 'Density Coeff Bs Solid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getDensity', 'Coefficientt of Bs*T2 term solid AlMgSiZr not found')
    !END IF
    
    ! read in  Ds in Ds*ln(T) term
    !deltas = GetConstReal( material, 'Density Coeff Ds Liquid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getDensity', 'Coefficient of logT term liquid AlMgSiZr not found')
    !END IF
    
    ! read in reference density at reference temperature
    refLiqDenst = GetConstReal( material, 'Reference Density Liquid AlMgSiZr',GotIt)
    !refDenst = GetConstReal( material, 'Solid_ti_rho_constant',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Reference Density Solid AlMgSiZr not found')
    END IF
    
    ! read in pseudo reference conductivity at reference temperature of liquid
    alphal = GetConstReal( material, 'Density Coefficient Liquid AlMgSiZr',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Density Coefficient Al of Liquid AlMgSiZr not found')
    END IF

    ! read in reference temperature
    refTemp = GetConstReal( material, 'Melting Point Temperature of AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getDensity', 'Reference Melting Temperature of AlMgSiZr not found')
    END IF
    
    ! read in the temperature scaling factor
    tscaler = GetConstReal( material, 'Tscaler', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Scaling Factor for T not found')
    END IF

    ! compute density conductivity
    IF (refTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getDensity', 'The AlMgSiZr material is in liquid state.')
            !CALL Warn('getDensity', 'Using density reference value')
    !denst = 1.11*(refDenst + alpha*(temp))
    denst = refLiqDenst + alphal*((tscaler)*(temp-910))
    ELSE
    denst = refSolDenst + alphas*((tscaler)*(temp-298)) 
    END IF

    END FUNCTION getDensity

