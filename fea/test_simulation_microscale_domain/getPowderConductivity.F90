    !-----------------------------------------------------
    ! material property user defined function for ELMER:
    ! Conductivity of gaseous air fitted as a function of temperature
    ! Expression 1; (kth_air)gas = Agas*T^2 + Bgas*T    where Agas = 0.00012 W/mK3  and  Bgas=  -0.001557 W/mK2 (250.0 < T < 1358.0 K) 
    ! Reference for (kth_air)gas: C.Y. Ho, R. W. Powell and P. E. Liley, Thermal Conductivity of the Elements: A Comprehensive Review,  Journal of Physical and Chemical Reference Data
    ! Volume 3 (1974) 1-756
    ! The following function yields better result for thermal conductivity of air
    ! Expression 2; (kth)_air = Agas*T + Cgas , where Agas = 1.7082E-4 W/(m  K2) and Cgas = -7.488E-3 W/m K
    ! This function is adapted from the expression ((kth)_air = Agas*T + Bgas*T**2 + Dgas*T**3 + Egas*T**4 + Fgas*T**5+ Cgas ) provided in the following reference
    ! https://www.cambridge.org/core/books/abs/gas-turbines/equations-of-air-thermophysical-properties/9572106E068EFF1B7C0896124C17A196
    ! It means a value of 0 has been put in Bgas, Dgas, Egas and Fgas to enable thermal conducitivity to be well around the order of 10^{-2} for given temperature range
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Initial porosity (porosity) = 0.45
    ! Thermal conductivity of AlMgSiZr fitted as a function of temperature
    ! (kth_alloy)solid = As*(T-298)^2  + Bs*(T-298)+  Cs, where As = -1.17E-04 W /mK3, Bs = 9.29E-03 W/mK2 and    Cs = 199.7 W/mK (298.0 K < T < 840.0 K)
    ! References: (i) Dong et al., Journal of Alloys and Compounds, 970 (2024) 172682.
    ! https://www.sciencedirect.com/science/article/pii/S0925838823039853
    !(ii) Zhang et al., Themochimica Acta, 635 (2016) 8-16.
    ! https://www.sciencedirect.com/science/article/abs/pii/S0040603116300892
    ! (kth_alloy)liquid = Al*(T-900) + Cl, where Al = 3.69E-02  W/mK2 and Cl = 64.549 W/m K (900.0 K < T < 2000.0 K) 
    ! Reference: Sun et al., International Journal of Thermophysics, 40 (2019) 31.
    ! https://link.springer.com/article/10.1007/s10765-019-2497-1
    ! Sintering temperature of AlMgSiZr = 450 oC ! https://journals.sagepub.com/doi/full/10.1080/00325899.2020.1719688
    !https://www.tandfonline.com/doi/pdf/10.1179/174329007X223947
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Written By: Anil Kunwar (Original 2015-03-13) (Modification 2023-11-08)
    !-----------------------------------------------------
    FUNCTION getThermalConductivity( model, n, temp ) RESULT(effcondt)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, effcondt, kthbulk, kthpowder, kthgas, tscaler

    ! variables needed inside function
    REAL(KIND=dp) ::  refMeltTemp, refStTemp, porosity,  &
    alphas, betas, alphal,  alphag, betag,  &
    refSolThCond, refLiqThCond
     
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getPowderConductivity', 'No material found')
    END IF

    ! read in reference Conductivity at reference temperature
    refSolThCond = GetConstReal( material, 'Reference Thermal Conductivity Solid AlMgSiZr',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Reference Conductivity Solid AlMgSiZr not found')
    END IF
    
    ! read in the coefficient of  Acu*(T-300)^2 term
    alphas = GetConstReal( material, 'Cond Coeff As Solid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Coefficientt of T2 term Solid AlMgSiZr not found')
    END IF
    
    ! read in the coefficient of  Bcu*(T-298)^2 term
    betas = GetConstReal( material, 'Cond Coeff Bs Solid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Coefficientt of T term Solid AlMgSiZr not found')
    END IF
    
    ! read in the coefficient of Dcu*(T-300) term
    !deltas = GetConstReal( material, 'Conductivity Coeff Dcu of Solid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getPowderConductivity', 'slope of Conductivity-temperature curve solid not found')
    !END IF
    
    ! read in reference Conductivity of Liquid AlMgSiZr
    refLiqThCond = GetConstReal( material, 'Reference Thermal Conductivity Liquid AlMgSiZr',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Reference Conductivity Liquid AlMgSiZr not found')
    END IF
    
    ! read in the coefficient of  Acu*(T-1330)^3 term
    alphal = GetConstReal( material, 'Cond Coeff Liquid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Coefficientt of T term Solid AlMgSiZr not found')
    END IF
    
    ! read in the coefficient of BCu*(T-1330)^2 term
    !betal = GetConstReal( material, 'Conductivity Coeff Bcu Liquid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getPowderConductivity', 'Coefficient of T2 term Liquid AlMgSiZr not found')
    !END IF
    
    ! read in the coefficient of Dcu*(T-1330) term
    !deltal = GetConstReal( material, 'Conductivity Coeff Dcu of Liquid AlMgSiZr', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getPowderConductivity', 'slope of Conductivity-temperature curve liquid not found')
    !END IF
        
    ! read in pseudo reference Conductivity of gaseous air
    !refGasDenst = GetConstReal( material, 'Reference Conductivity Cgas of Air',GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getPowderConductivity', 'Reference Conductivity for Gaseous Air not found')
    !END IF
        
    ! read in the coefficient of  Agas*(T)^2 term
    alphag = GetConstReal( material, 'Conductivity Coeff Agas of Air', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Coefficient of T^2 term for air Conductivity not found')
    END IF
    
    ! read in the coefficient of Bgas*T term for air Conductivity
    betag = GetConstReal( material, 'Conductivity Coeff Bgas of Air',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Slope of rhogas-T not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    !deltag = GetConstReal( material, 'Conductivity Coeff Dgas of Air', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getPowderConductivity', 'Coefficient of lnT term Gaseous Air not found')
    !END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    porosity = GetConstReal( material, 'Initial porosity of AlMgSiZr powder', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Porosity of AlMgSiZr powder not found')
    END IF
    
    ! read in reference sintering temperature of Ag powder
    refStTemp = GetConstReal( material, 'Sintering Temperature of AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Reference Sintering Temperature AlMgSiZr not found')
    END IF

    ! read in reference melting temperature of Ag powder
    refMeltTemp = GetConstReal( material, 'Melting Point Temperature of AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Reference Melting Temperature AlMgSiZr not found')
    END IF
    
     ! read in the temperature scaling factor
    tscaler = GetConstReal( material, 'Tscaler', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getPowderConductivity', 'Scaling Factor for T not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance for Liquid Ti
    !porosity = GetConstReal( material, 'Initial porosity of Cu powder', GotIt)
    !IF(.NOT. GotIt) THEN
    !CALL Fatal('getPowderConductivity', 'Porosity of Cu powder not found')
    !END IF

     
    ! compute Conductivity conductivity
    ! https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap03/else-if.html
    IF (refMeltTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getPowderConductivity', 'The AlMgSiZr material is in liquid state.')
            !CALL Warn('getPowderConductivity', 'Using Conductivity reference value')
    !denst = 1.11*(refDenst + alpha*(temp))
    !effcondt = alphal*(tscaler*temp)**3 + betal*(tscaler*temp)**2 +deltal*(tscaler*temp)
    effcondt = refLiqThCond + alphal*((tscaler)*(temp-900))
    ELSE IF (refMeltTemp > temp .AND. refStTemp < temp) THEN
       CALL Warn('getPowderConductivity', 'The AlMgSiZr material is being sintered.')
    !kthbulk = alphas*(tscaler*temp)**3 + betas*(tscaler*temp)**2 + deltas*(tscaler*temp)
    kthbulk = refSolThCond + betas*((tscaler)*(temp-298)) + alphas*((tscaler)*(temp-298))**2 
    !kthgas = alphag*(tscaler*temp)**2 + betag*(tscaler*temp) !Expression 1
    kthgas = alphag*(tscaler*temp)  + betag !Expression 2
    kthpowder = (1-porosity)*kthbulk + porosity*kthgas
    effcondt = ((kthbulk - kthpowder)*(temp-refStTemp))/(refMeltTemp-refStTemp)+kthpowder
    ELSE
    effcondt = (1-porosity)*(refSolThCond + betas*((tscaler)*(temp-298)) + alphas*((tscaler)*(temp-298))**2 )+ &
    porosity*(alphag*(tscaler*temp)**2 + betag*tscaler*temp)
    END IF

    END FUNCTION getThermalConductivity

