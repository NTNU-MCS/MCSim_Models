% Marine Guidance, Navigation and Control Toolbox
% Version 1.6 (R13) 24-Mar-2003
%
% Ref.: T. I. Fossen (2002). Marine Control Systems: Guidance, Navigation and Control
%       of Ships, Rigs and Underwater Vehicles (Marine Cybernetics AS, Trondheim).
%
%       Book Examples: >>help examples 
% -------------------------------------------------------------------------
% Marine GNC Simulink Library: >>marineGNC
%
% M-File Demonstrations:
%    GNCDemo    - Demonstration menu
%    KinDemo    - Euler angle and quaternion kinematics.
%    ManDemo    - Maneuvering trials.
%    SimDemo    - User editable scripts for simulation and control of m-file vessels models
%    StabDemo   - Straight-line, directional and positional motion stability
%    WaveDemo   - Wave spectra demonstrations
%
% GNC:
%    ucalloc    - Unconstrained control allocation for tau = T*K*u (T=K=constant)
%    lqtracker  - Computes the LQ tracker gain matrices for LTI systems 
%    order3     - Path generation using cubic (3rd-order) polynominals 
%    order5     - Path generation using 5th order polynominals
%
% Standard ship, semi-submersible and underwater vehicle models:
%
%    Ships:
%    mariner     - Mariner class vessel, L=160 m (nonlinear maneuvering model)
%    tanker      - Esso Osaka tanker, L=304 m (nonlinear course unstable maneuvering model)
%    container   - Container ship, L=175 m (nonlinear maneuvering model including the roll mode)
%    Lcontainer  - Container ship, L=175 m (LINEAR maneuvering model including the roll mode)
%
%    Underwater Vehicles:
%    dsrv        - Deep submergence rescue vehicle (DSRV), L = 5.0 m
%    npsauv      - Naval Postgraduate School autnomous underwater vehicle (AUV), L = 5.3 m
% 
% Maneuvering trials/Bode plots:
%    nomoto      - Bode plots of Nomoto's 1st- and 2nd-order models
%    TurnCircle  - performs a turning circle for a given ship model
%    ZigZag      - performs a zig-zag maneuver for a given ship model
%    Pullout     - performs a pullout maneuver for a given ship model
%
% Model transformations and conversion factors:
%    conversion - loads a set of useful conversion factors to workspace, D2R, R2D, MS2KNOTS etc.
%    Smtrx      - matrix vector cross product operator
%    Hmtrx      - system transformation matrix
%    Gmtrx      - Gravitational-buoyancy matrix (floating vessels)
%    gvect      - Gravitational-buoyancy vector (submerged vehicles)
%    m2c        - System inertia to Coriolis-centripetal transformation matrix
%
% Kinematics:
%    eulerang   - Euler angle transformation matrices J,J1 and J2
%    ecef2llh   - ECEF xyz-coordinates to longitude-lattitude-height
%    euler2q    - Euler angles (roll-pitch-yaw) to quaternions
%    quest      - QUEST algorithm for attitude determination
%    quest6DOF  - 6 DOF position/attitude vector from 3 position measurements
%    llh2ecef   - longitude-lattitude-height to ECEF xyz-coordinates.
%    quatern    - Quaternion transformation matrices J,J1 and J2
%    euler2q    - Computes the Euler angles from the unit quaternions.
%    rad2pipi   - Converts an angle in rad to the interval <-pi pi>
%    Rll        - Euler angle rotation matrix (longitude-lattitude)
%    Rzyx       - Euler angle rotation matrix (roll-pitch-yaw)
%    Rquat      - Quaternion rotation matrix
%    q2euler    - Quaternions to Euler angles (roll-pitch-yaw)
%    R2euler    - Rotation matrix to Euler angles (roll-pitch-yaw)
%
% Wave Spectra, encounter frequency and motion sickness incidence (MSI):
%    encounter  - Computes the encounter frequency
%    HMmsi      - Motion Sickness Incidence using the criterion of O'Hanlon and McCauley
%    ISOmsi     - Motion Sickness Incidence using the ISO 2631-1 criterion 
%    mpierson   - Modified Pierson-Moskowitz (PM) wave spectrum for fully developed sea
%    pierson    - The "classical" one parameter Pierson-Moskowitz (PM) wave spectrum
%    jonswap    - JONSWAP (Joint North Sea Wave Project) wave spectrum for developing sea
%    torset     - Torsethaugen two peaked spectrum for swell and developing sea
%
% Wind forces and moments:
%
%    windcoef   - Wind forces and moments using the data of Isherwood (1972)
%
% Numerical integration:
%    euler2     - 2nd-order Euler integration, fixed step
%    rk4        - 4th-order Runge Kutta method, fixed step
% 
% Copyright (c) 2001 by Marine Cybernetics AS <http://www.marinecybernetics.com/>