    !-----------------------------------------------------
    ! Written By: Anil Kunwar (Original 2015-03-13) (Modification 2023-11-12)
    ! material property user defined function for ELMER:
    ! Thermal conductivity of AlMgSiZr fitted as a function of temperature
    ! (kth_alloy)solid = As*(T-298)^2  + Bs*(T-298)+  Cs, where As = -1.17E-04 W /mK3, Bs = 9.29E-03 W/mK2 and    Cs = 199.7 W/mK (298.0 K < T < 840.0 K)
    ! References: (i) Dong et al., Journal of Alloys and Compounds, 970 (2024) 172682.
    ! https://www.sciencedirect.com/science/article/pii/S0925838823039853
    !(ii) Zhang et al., Themochimica Acta, 635 (2016) 8-16.
    ! https://www.sciencedirect.com/science/article/abs/pii/S0040603116300892
    ! (kth_alloy)liquid = Al*(T-900) + Cl, where Al = 3.69E-02  W/mK2 and Cl = 64.549 W/m K (900.0 K < T < 2000.0 K) 
    ! Reference: Sun et al., International Journal of Thermophysics, 40 (2019) 31.
    ! https://link.springer.com/article/10.1007/s10765-019-2497-1
    !-----------------------------------------------------
    FUNCTION getThermalConductivity( model, n, temp ) RESULT(thcondt)
    ! modules needed
    USE DefUtils
    IMPLICIT None
    ! variables in function header
    TYPE(Model_t) :: model
    INTEGER :: n
    REAL(KIND=dp) :: temp, thcondt, tscaler

    ! variables needed inside function
    REAL(KIND=dp) :: refSolThCond, refLiqThCond,refTemp, &
    alphas, betas, alphal 
    Logical :: GotIt
    TYPE(ValueList_t), POINTER :: material

    ! get pointer on list for material
    material => GetMaterial()
    IF (.NOT. ASSOCIATED(material)) THEN
    CALL Fatal('getThermalConductivity', 'No material found')
    END IF

    ! read in reference conductivity at reference temperature
    refSolThCond = GetConstReal( material, 'Reference Thermal Conductivity Solid AlMgSiZr',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Reference Thermal Conductivity Solid AlMgSiZr not found')
    END IF

    ! read in Temperature Coefficient of Resistance
    alphas = GetConstReal( material, 'Cond Coeff As Solid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Coefficient of (T-298)^2 term of solid not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance
    betas = GetConstReal( material, 'Cond Coeff Bs Solid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'slope of thermal conductivity-temperature curve solid not found')
    END IF
    
    ! read in reference conductivity at reference temperature
    refLiqThCond = GetConstReal( material, 'Reference Thermal Conductivity Liquid AlMgSiZr',GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Reference Thermal Conductivity Liquid AlMgSiZr not found')
    END IF
    
    ! read in Temperature Coefficient of Resistance
    alphal = GetConstReal( material, 'Cond Coeff Liquid AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Coefficientt of T  term Liquid AlMgSiZr not found')
    END IF
    
    
    ! read in reference temperature
    refTemp = GetConstReal(material, 'Melting Point Temperature of AlMgSiZr', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Reference Temperature not found')
    END IF
    
    ! read in the temperature scaling factor
    tscaler = GetConstReal( material, 'Tscaler', GotIt)
    IF(.NOT. GotIt) THEN
    CALL Fatal('getThermalConductivity', 'Scaling Factor for T not found')
    END IF


    ! compute density conductivity
    IF (refTemp <= temp) THEN ! check for physical reasonable temperature
       CALL Warn('getThermalConductivity', 'The Ti material is in liquid state.')
            !CALL Warn('getThermalConductivity', 'Using density reference value')
    !thcondt = 1.11*(refThCond + alpha*(temp))
    thcondt = refLiqThCond + alphal*((tscaler)*(temp-900))
    ELSE
    thcondt = refSolThCond + betas*((tscaler)*(temp-298)) + alphas*((tscaler)*(temp-298))**2 
    END IF

    END FUNCTION getThermalConductivity

