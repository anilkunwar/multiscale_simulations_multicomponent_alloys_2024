    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Thermal Expansivity of solid Si fitted as a function of temperature
    ! (alpha_alloy)solid = Cs + As*[exp(temp/theta)*(temp/theta)**2]/(exp(temp/theta)-1)**2 + Bs*(((temp/phi)-1)**2)/(1+0.316*(temp/phi))
    ! where As = 5.0E-6 K^{-1}, Bs =0.22e-06 K^{-1}  and Cs = -0.687e-06 K^{-1}, theta=685.0 K and phi = 395 K
    ! 90 < T < 850 
    ! Reference: Swenson, J. Phys. Chem. Ref. Data, Vol. 12, No.2, 1983.
    ! https://doi.org/10.1063/1.555681
    ! Tm = 1687 K ,  Mills and Courtney, Thermophysical Properties of Silicon, ISIJ International 2000
    ! Thermal Expansivity of liquid Si fitted as a function of temperature
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
    REAL(KIND=dp) :: temp, expansivity

    ! variables needed inside function
    !REAL(KIND=dp) :: refSolExpansivity, refLiqExpansivity, refTemp,  &
    !alphas, alphal, betas, betal
    REAL(KIND=dp) :: refSolExpansivity, theta, refphi,  &
    alphas, betas, refTemp
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getThermalExpansivity', 'No material found')
    END IF

    ! read in reference conductivity at reference temperature
    refSolExpansivity = GetConstReal( material, 'Cs of Silicon',GotIt)
    !refDenst = GetConstReal( material, 'Solid_ti_rho_constant',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'Cs silicon not found')
    END IF

    ! read in Temperature Coefficient of Resistance
    alphas = GetConstReal( material, 'Thermal Expansivity Coeff As Solid Si', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'Heat capacity term of Thermal Expansivity-temperature curve solid not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance
    betas = GetConstReal( material, 'Thermal Expansivity Coeff Bs Solid Si', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'Deviation  term of solid Si not found')
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
    theta = GetConstReal( material, 'Debye Theta Temperature of Si', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'Theta of  Si not found')
    END IF
    
    ! read in reference temperature
    refphi = GetConstReal( material, 'Phi Temperature of Si', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'Phi of  Si not found')
    END IF

    ! read in reference temperature
    refTemp = GetConstReal( material, 'Melting Point Temperature of Si', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalExpansivity', 'Reference Melting Temperature of AlMgSiZr not found')
    END IF
    
    ! read in the temperature scaling factor
    !tscaler = GetConstReal( material, 'Tscaler', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getThermalConductivity', 'Scaling Factor for T not found')
    !END IF

    ! compute Thermal Expansivity conductivity
    IF (refTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getThermalExpansivity', 'The AlMgSiZr material is in liquid state.')
            !CALL Warn('getThermalExpansivity', 'Using Thermal Expansivity reference value')
    !denst = 1.11*(refDenst + alpha*(temp))
    !denst = refLiqExpansivity + alphal*((tscaler)*(temp-910))
    expansivity = 0.0
    !expansivity = refSolExpansivity + alphas*((tscaler)*(temp))**2 + betas*(tscaler)*(temp)
    ELSE
    ! To add the expansivity
    expansivity = refSolExpansivity + &
          alphas*(exp(temp/theta)*(temp/theta)**2)/(exp(temp/theta)-1)**2 + &
          betas*(((temp/refphi)-1)**2)/(1+0.316*(temp/refphi))
    END IF

    END FUNCTION getThermalExpansivity

