    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Density of gaseous air fitted as a function of temperature
    ! (rho_air)gas = Agas*T^2 + Bgas*T  + Dgas*ln(Egas*T) +Cgas, where Agas = -7.65e-07 kg/m3K2  and  Bgas=  3.16e-03 kg/m3K and Dgas=-1.94 kg/m3 and Egas = 1.0 K-1 and Cgas = 1.14e+01  kg/m3
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Initial porosity (porosity) = 0.45
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
    ! Method for defining equivalent density in powder bed material
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Written By: Anil Kunwar (Original 2015-03-13) (Modification 2023-11-08)
    !-----------------------------------------------------
    FUNCTION getDensity( model, n, temp ) RESULT(effdenst)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, effdenst, rhobulk, rhopowder, rhogas, tscaler

    ! variables needed inside function
    REAL(KIND=dp) :: refSolidDenst, refLiquidDenst, refGasDenst,  &
    alphas, alphal, alphag, betag, deltag,  & 
    porosity,  refStTemp, refMeltTemp
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getMaterialDensity', 'No material found')
    END IF

    ! read in reference density at reference temperature
    refSolidDenst = GetConstReal( material, 'Reference Density Solid AlMgSiZr',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Reference Density Solid AlMgSiZr not found')
    END IF
    
    ! read in the coefficient of  Acu*(T-300)^2 term
    alphas = GetConstReal( material, 'Density Coeff Solid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Coefficientt of T3 term Solid AlMgSiZr not found')
    END IF
    
    ! read in the coefficient of  Bcu*(T-300)^2 term
    !betas = GetConstReal( material, 'Density Coeff Bcu Solid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getMaterialDensity', 'Coefficientt of T2 term Solid AlMgSiZr not found')
    !END IF
    
    ! read in the coefficient of Dcu*(T-300) term
    !deltas = GetConstReal( material, 'Density Coeff Dcu of Solid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getMaterialDensity', 'slope of Density-temperature curve solid not found')
    !END IF
    
    ! read in reference density of Liquid AlMgSiZr
    refLiquidDenst = GetConstReal( material, 'Reference Density Liquid AlMgSiZr',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Reference Density Liquid AlMgSiZr not found')
    END IF
    
    ! read in the coefficient of  Acu*(T-1330)^3 term
    alphal = GetConstReal( material, 'Density Coeff Liquid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Coefficientt of T3 term Solid AlMgSiZr not found')
    END IF
    
    ! read in the coefficient of BCu*(T-1330)^2 term
    !betal = GetConstReal( material, 'Density Coeff Bcu Liquid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getMaterialDensity', 'Coefficient of T2 term Liquid AlMgSiZr not found')
    !END IF
    
    ! read in the coefficient of Dcu*(T-1330) term
    !deltal = GetConstReal( material, 'Density Coeff Dcu of Liquid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getMaterialDensity', 'slope of Density-temperature curve liquid not found')
    !END IF
        
    ! read in pseudo reference density of gaseous air
    refGasDenst = GetConstReal( material, 'Reference Density Cgas of Air',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Reference Density for Gaseous Air not found')
    END IF
        
    ! read in the coefficient of  Agas*(T)^2 term
    alphag = GetConstReal( material, 'Density Coeff Agas of Air', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Coefficient of T^2 term for air density not found')
    END IF
    
    ! read in the coefficient of Bgas*T term for air density
    betag = GetConstReal( material, 'Density Coeff Bgas of Air',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Slope of rhogas-T not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    deltag = GetConstReal( material, 'Density Coeff Dgas of Air', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Coefficient of lnT term Gaseous Air not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    porosity = GetConstReal( material, 'Initial porosity of AlMgSiZr powder', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Porosity of AlMgSiZr powder not found')
    END IF
    
    ! read in reference sintering temperature of Ag powder
    refStTemp = GetConstReal( material, 'Sintering Temperature of AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Reference Sintering Temperature Cu not found')
    END IF

    ! read in reference melting temperature of Ag powder
    refMeltTemp = GetConstReal( material, 'Melting Point Temperature of AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Reference Melting Temperature Cu not found')
    END IF
    
    ! read in the temperature scaling factor
    tscaler = GetConstReal( material, 'Tscaler', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Scaling Factor for T not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    porosity = GetConstReal( material, 'Initial porosity of AlMgSiZr powder', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getMaterialDensity', 'Porosity of AlMgSiZr powder not found')
    END IF

     
    ! compute density conductivity
    ! https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap03/else-if.html
    IF (refMeltTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getMaterialDensity', 'The AlMgSiZr material is in liquid state.')
            !CALL Warn('getMaterialDensity', 'Using density reference value')
    !denst = 1.11*(refDenst + alpha*(temp))
    effdenst = refLiquidDenst + alphal*(tscaler*(temp-910.00))
    ELSE IF (refMeltTemp > temp .AND. refStTemp < temp) THEN
       CALL Warn('getMaterialDensity', 'The Cu material is being sintered.')
    rhobulk = refSolidDenst +  alphas*(tscaler*(temp-298.0)) 
    rhogas = refGasDenst +  alphag*(tscaler*temp)**2 + betag*tscaler*temp + deltag*LOG(tscaler*temp)
    rhopowder = (1-porosity)*rhobulk + porosity*rhogas
    effdenst = ((rhobulk - rhopowder)*(temp-refStTemp))/(refMeltTemp-refStTemp)+rhopowder
    ELSE
    effdenst = (1-porosity)*(refSolidDenst +  alphas*(tscaler*(temp-298.0)))  + &
    porosity*(refGasDenst +  alphag*(tscaler*temp)**2 + betag*tscaler*temp + deltag*LOG(tscaler*temp))
    END IF

    END FUNCTION getDensity

