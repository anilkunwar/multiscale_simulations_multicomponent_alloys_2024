!Super-Gaussian Heat Source
FUNCTION FlatTopHeatSource( Model, n, t ) RESULT(f)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Model_t) :: Model
  INTEGER :: n
  REAL(KIND=dp) :: t, f


  INTEGER :: timestep, prevtimestep = -1
  REAL(KIND=dp) :: Alpha, Coeff, xspeed, yspeed,  &
      Time, x, y, z, s1, s2,  sper, r, xzero, yzero, Omega, &
      sgo, m1, m2,  rsgo
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, NewTimestep
  
  !Save allows taking coordinate memory?
  SAVE Mesh, Params, prevtimestep, time, Alpha, Coeff, xspeed,  &
        yspeed, xzero, yzero, sgo, m1, m2, rsgo 
  
  timestep = GetTimestep()
  NewTimestep = ( timestep /= prevtimestep )

  IF( NewTimestep ) THEN
    Mesh => GetMesh()
    Params => Model % Simulation
    time = GetTime()
    Alpha = GetCReal(Params,'Beam spot radius')
    Coeff = GetCReal(Params,'Heat source coefficient term in the amplitude term')
    xspeed = GetCReal(Params,'Laser scan speed in x direction')
    yspeed = GetCReal(Params,'Laser speed in x direction')
    xzero = GetCReal(Params,'Heat source initial position x', Found)
    yzero = GetCReal(Params,'Heat source initial position y', Found)
    sgo = GetCReal(Params,'Super gaussian order k')
    rsgo = GetCReal(Params,'reciproccal of Super-gaussian order 1/k')
    m1 = GetCReal(Params,'prefactor within amplitude term')
    m2 = GetCReal(Params,'prefactor within exponential term')
    prevtimestep = timestep
  END IF

  x = Mesh % Nodes % x(n)   
  y = Mesh % Nodes % y(n)   
  z = Mesh % Nodes % z(n)   

  s1 = xzero - time * xspeed ! + for forward and - for reverse journey
  s2 = yzero + time * yspeed 
  r = SQRT((x-s1)**2 + (y-s2)**2)   
  
  !Coeff = P*eta/(pi*alpha^2)
  f = m1**(rsgo)*sgo*Coeff * EXP( -m2*r**sgo / Alpha**sgo )/gamma(rsgo)
 
  
    
END FUNCTION FlatTopHeatSource




