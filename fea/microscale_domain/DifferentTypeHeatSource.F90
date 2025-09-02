!Gaussian 
FUNCTION TravellingHeatSource( Model, n, t ) RESULT(f)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f


  INTEGER :: timestep, prevtimestep = -1
  REAL(KIND=dp) :: Alpha, Coeff, xspeed, yspeed,  Dist, &
      Time, x, y, z, s1, s2, sper, r, xzero, yzero, Omega
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, NewTimestep
  
  SAVE Mesh, Params, prevtimestep, time, Alpha, Coeff, xspeed, yspeed,  &
       Dist, Omega
  
  timestep = GetTimestep()
  NewTimestep = ( timestep /= prevtimestep )

  IF( NewTimestep ) THEN
    Mesh => GetMesh()
    Params => Model % Simulation
    time = GetTime()
    Alpha = GetCReal(Params,'Heat source width')
    Coeff = GetCReal(Params,'Heat source coefficient')
    xspeed = GetCReal(Params,'Heat source speed x')
    yspeed = GetCReal(Params,'Heat source speed y')
    Dist = GetCReal(Params,'Heat source distance')
    xzero = GetCReal(Params,'Heat source initial position', Found)
    yzero = GetCReal(Params,'y coordinate initial position', Found)
    Omega = GetCReal(Params,'Absorptance of Surface Material')
    prevtimestep = timestep
  END IF

  x = Mesh % Nodes % x(n)   
  y = Mesh % Nodes % y(n)   
  z = Mesh % Nodes % z(n)   

  s1 = xzero + time * xspeed
  s2 = yzero + time * yspeed 
  r = SQRT((x-s1)**2 + (y-s2)**2)   
  !sper = MODULO( s, 2 * Dist ) 
  !IF( sper > Dist ) sper = 2 * Dist - sper

  !r = x-s!sper
  ! in 3D this could be the radius
   !r = SQRT((x-s)**2 + (y-yzero)**2)
  
  !f = Coeff * EXP( -2*r**2 / Alpha**2 )
  f = Coeff * EXP( -2*r**2 / Alpha**2 -Omega * ABS(z))
    
END FUNCTION TravellingHeatSource

!Gaussian
FUNCTION FixedHeatSource( Model, n, t ) RESULT(f)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f


  INTEGER :: timestep, prevtimestep = -1
  REAL(KIND=dp) :: Alpha, Coeff, Dist0, &
      Time, x, y, z, r
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, NewTimestep
  
  SAVE Mesh, Params, prevtimestep, time, Alpha, Coeff, Dist0
  
  timestep = GetTimestep()
  NewTimestep = ( timestep /= prevtimestep )

  IF( NewTimestep ) THEN
    Mesh => GetMesh()
    Params => Model % Simulation
    time = GetTime()
    Alpha = GetCReal(Params,'Heat source width')
    Coeff = GetCReal(Params,'Heat source coefficient')
    Dist0 = GetCReal(Params,'Heat source initial position', Found)
    prevtimestep = timestep
  END IF

  x = Mesh % Nodes % x(n)   
  y = Mesh % Nodes % y(n)   
  z = Mesh % Nodes % z(n)   

  ! not a function of time
  r = x-Dist0

  ! in 3D this could be the radius
  ! r = SQRT((x-s)**2 + y**2)
  
  f = Coeff * EXP( -2*r**2 / Alpha**2 )
    
END FUNCTION FixedHeatSource

!Slab movement
FUNCTION TravellingSlabVelo( Model, n, t ) RESULT(f)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f


  INTEGER :: timestep, prevtimestep = -1
  REAL(KIND=dp) :: Alpha, Coeff, Speed, Dist, Dist0, &
      Time, x, y, z, s, sper, r
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, NewTimestep
  
  SAVE Mesh, Params, prevtimestep, time, Alpha, Coeff, Speed, Dist, &
      Dist0
  
  timestep = GetTimestep()
  NewTimestep = ( timestep /= prevtimestep )

  IF( NewTimestep ) THEN
    Mesh => GetMesh()
    Params => Model % Simulation
    time = GetTime()
    Speed = GetCReal(Params,'Heat source speed')
    Dist = GetCReal( Params,'Heat source distance')
    prevtimestep = timestep
  END IF

  x = Mesh % Nodes % x(n)   
  y = Mesh % Nodes % y(n)   
  z = Mesh % Nodes % z(n)   

  s = Dist0 + time * Speed  
  sper = MODULO( s, 2 * Dist ) 

  IF( sper > Dist ) THEN
    sper = 2 * Dist - sper
    f = Speed
  ELSE
    f = -Speed
  END IF
    
END FUNCTION TravellingSlabVelo




!Flat top for travelling source
FUNCTION FlatTopHeatSource( Model, n, t ) RESULT(f)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f


  INTEGER :: timestep, prevtimestep = -1
  REAL(KIND=dp) :: Alpha, Coeff, xspeed, yspeed,  Dist, &
      Time, x, y, z, s1, s2,  sper, r, xzero, yzero, Omega, &
      sgo, m1, m2,  rsgo
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, NewTimestep
  
  !Save allows taking coordinate memory?
  SAVE Mesh, Params, prevtimestep, time, Alpha, Coeff, xspeed, Dist, &
        yspeed, xzero, yzero, sgo, m1, m2, rsgo !, Sigma, Omega, Dist0
  
  timestep = GetTimestep()
  NewTimestep = ( timestep /= prevtimestep )

  IF( NewTimestep ) THEN
    Mesh => GetMesh()
    Params => Model % Simulation
    time = GetTime()
    Alpha = GetCReal(Params,'Heat source width')
    Coeff = GetCReal(Params,'Heat source coefficient')
    xspeed = GetCReal(Params,'Heat source speed x')
    yspeed = GetCReal(Params,'Heat source speed y')
    Dist = GetCReal(Params,'Heat source distance')
    xzero = GetCReal(Params,'Heat source initial position x', Found)
    yzero = GetCReal(Params,'Heat source initial position y', Found)
    Omega = GetCReal(Params,'Absorptance of Surface Material')
    sgo = GetCReal(Params,'Super gaussian order n')
    rsgo = GetCReal(Params,'reciproccal of Super gaussian order 1/n')
    m1 = GetCReal(Params,'prefactor within amplitude term')
    m2 = GetCReal(Params,'prefactor within exponential term')
    prevtimestep = timestep
  END IF

  x = Mesh % Nodes % x(n)   
  y = Mesh % Nodes % y(n)   
  z = Mesh % Nodes % z(n)   

  s1 = xzero + time * xspeed
  s2 = yzero + time * yspeed 
  r = SQRT((x-s1)**2 + (y-s2)**2)   
  
  !Coeff = P*eta/(pi*alpha^2)
  f = m1**(rsgo)*sgo*Coeff * EXP( -m2*r**sgo / Alpha**sgo )/gamma(rsgo)
  !f = Coeff * EXP( -2*r**2 / Alpha**2 -Omega * ABS(z))
  
    
END FUNCTION FlatTopHeatSource


!Flat top for travelling source
FUNCTION GaussianHeatSource( Model, n, t ) RESULT(f)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f


  INTEGER :: timestep, prevtimestep = -1
  REAL(KIND=dp) :: Alpha, Coeff, xspeed, yspeed,  Dist, &
      Time, x, y, z, s1, s2,  sper, r, xzero, yzero, Omega
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, NewTimestep
  
  !Save allows taking coordinate memory?
  SAVE Mesh, Params, prevtimestep, time, Alpha, Coeff, xspeed, Dist, &
        yspeed, xzero, yzero !, Sigma, Omega, Dist0
  
  timestep = GetTimestep()
  NewTimestep = ( timestep /= prevtimestep )

  IF( NewTimestep ) THEN
    Mesh => GetMesh()
    Params => Model % Simulation
    time = GetTime()
    Alpha = GetCReal(Params,'Heat source width')
    Coeff = GetCReal(Params,'Heat source coefficient')
    xspeed = GetCReal(Params,'Heat source speed x')
    yspeed = GetCReal(Params,'Heat source speed y')
    Dist = GetCReal(Params,'Heat source distance')
    xzero = GetCReal(Params,'Heat source initial position x', Found)
    yzero = GetCReal(Params,'Heat source initial position y', Found)
    Omega = GetCReal(Params,'Absorptance of Surface Material')
    prevtimestep = timestep
  END IF

  x = Mesh % Nodes % x(n)   
  y = Mesh % Nodes % y(n)   
  z = Mesh % Nodes % z(n)   

  s1 = xzero - time * xspeed
  s2 = yzero + time * yspeed 
  r = SQRT((x-s1)**2 + (y-s2)**2)   
  
  f = Coeff * EXP( -2*r**2 / Alpha**2 )
  !f = Coeff * EXP( -2*r**2 / Alpha**2 -Omega * ABS(z))
  
    
END FUNCTION GaussianHeatSource
