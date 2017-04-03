%% KINGSROCKsim // Finian P.G. Russell // Kingston Uni MEng Astro 2016/17
%% Final Script - Dev note: search %%%DEV:NB%%%
%{
               /\
              /  \
             /    \
            /      \
           |        |
           |        |
           |        |
           |        |
           |        |
           |        |
           |        |
          /|        |\
         / |        | \
        /  |        |  \
       /___|________|___\
            ********
           **********
          ************
         **************


  -Requires 8 complementary functions:

    *atmosDensity
    *atmosTemp
    *calcParasiteDrag.m
    *changeInput.m
    *extendSimData.m
    *liftCoef
    *makeGraph.m
    *oswaldFactor_mod.m

  -Make sure these are in the same directory as this program before running it.


  -DEV: Use ctrl+F to go to a section/subsection using its associated code (see below).

%% Section             |   Code

   Variables           |   var01
   Initialization      |   ini02
   Rocket Design       |   roc03
   Inputs              |   inp04
   Simulation Setup    |   set05
   Engine Import       |   eng06
   COP Calculations    |   cop07
   COM Calculations    |   com08
   Stability           |   sta09
   Control systems     |   con10
   Flight simulation   |   fli11
   Results             |   res12
   Microgravity        |   mic13
   Graphs              |   gra14
   Animations          |   ani15

%}

close all
clear
clc
commandwindow
    %imshow('kingsrock.png','Border','tight')
    alphaVer = 1.03;
    betaVer = 1.01;
    disp(['KINGSROCKsim - KINGSton Model ROCKet Trajectory SIMulator - beta ',num2str(betaVer),' - by Finian P.G. Russell'])
    fprintf('%s\n',[],'Credits:','Dr Malcolm Claus & Dr Adam Baker, Kingston University - supervision','Rick Newlands, Aspirespace - advice','Niklaus Kamm - peer review',...
        'Sampo Niskanen (OpenRocket) | Charles Rogers/David Cooper (RASAero)','Simon Box et al (Cambridge Rocketry Toolbox) - inspiration',...
        'Newton Launch Systems (Newton Trajectory sim) - inspiration','See script file for others...')
    
    %{
    Sources:
        http://openrocket.sourceforge.net/
        http://www.rasaero.com/
        http://cambridgerocket.sourceforge.net/
        https://www.apogeerockets.com/Rocket_Software/RockSim
        ...
        http://openrocket.sourceforge.net/techdoc.pdf
        http://cambridgerocket.sourceforge.net/AerodynamicCoefficients.pdf
        http://web.aeromech.usyd.edu.au/AERO2705/Resources/Research/Drag_Coefficient_Prediction.pdf
        https://web.archive.org/web/20110411143013/http://www.if.sc.usp.br/~projetosulfos/artigos/NoseCone_EQN2.PDF
        http://www.rocketmime.com/rockets/Barrowman.html
        https://www.apogeerockets.com/downloads/PDFs/barrowman_report.pdf
        https://apogeerockets.com/education/downloads/Newsletter238.pdf
        http://ftp.demec.ufpr.br/foguete/bibliografia/Labudde_1999_CP.pdf
        http://argoshpr.ch/joomla1/articles/pdf/sentinel39-galejs.pdf
        http://cnqzu.com/library/Anarchy%20Folder/Rocketry/Missiles%20and%20Warheads/Missile%20Guidance%20&%20Control%20Systems.pdf
        http://www.philsrockets.org.uk/forces.pdf
        http://www.philsrockets.org.uk/wind.pdf
        http://math.etsu.edu/multicalc/prealpha/Chap3/Chap3-1/printversion.pdf
        http://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
        http://www.metoffice.gov.uk/public/weather/forecast/gcpvj0v07
        http://drømstørre.dk/wp-content/wind/miller/windpower%20web/en/tour/wres/weibull.htm
        http://www.wind-power-program.com/wind_statistics.htm
        https://uk.mathworks.com/help/stats/prob.truncatabledistribution.pdf.html?requestedDomain=www.mathworks.com#bttw6d3-1
        http://wpage.unina.it/agodemar/DSV-DQV/Digital_Datcom_Users_Manual_1.2.pdf
        http://dtic.mil/dtic/tr/fulltext/u2/a548461.pdf
        http://www.engineeringtoolbox.com/wind-shear-d_1215.html
        http://www.aerospaceweb.org/question/airfoils/q0150b.shtml
        https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/rktcock.html
        https://www.grc.nasa.gov/www/k-12/airplane/induced.html
        https://uk.mathworks.com/matlabcentral/fileexchange/38800-oswald-efficiency-estimation-function?requestedDomain=www.mathworks.com
        https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19690005172.pdf
    %}
    
    if	exist('datetime','builtin')
	disp('The current date/time cannot be displayed. Sorry!')
    else
	fprintf('%s\n',[],datetime('now'),[])
    end
    
    if      exist('atmosDensity.m') == 0
            error('You do not have the ''atmosDensity.m'' function located within the same directory as this script.')
    elseif  exist('atmosTemp.m') == 0
            error('You do not have the ''atmosTemp.m'' function located within the same directory as this script.')
    elseif  exist('calcParasiteDrag.m') == 0
            error('You do not have the ''calcParasiteDrag.m'' function located within the same directory as this script.')
    elseif  exist('changeInput.m') == 0
            error('You do not have the ''changeInput.m'' function located within the same directory as this script.')
    elseif  exist('extendSimData.m') == 0
            error('You do not have the ''extendSimData.m'' function located within the same directory as this script.')
    elseif  exist('liftCoef.m') == 0
            error('You do not have the ''liftCoef.m'' function located within the same directory as this script.')
    elseif  exist('makeGraph.m') == 0
            error('You do not have the ''makeGraph.m'' function located within the same directory as this script.')
    elseif  exist('oswaldFactor_mod.m') == 0
            error('You do not have the ''oswaldFactor_mod.m'' function located within the same directory as this script.')
    end
    
    fprintf('%s\n','If you make a mistake in entering values for user input, you''ll need to CTRL+C and restart the program.',[])

    
%% VARIABLES [var01]
%%
%{
 animCounter            //	Counter for animated trajectory graph, for loop
 anim_line              //  Animated trajectory graph
 bodyDiameter           //  Rocket body diameter
 centreOfMass           //  Rocket centre of mass
 centreOfMass_noProp    //  Rocket centre of mass without propellant
 centreOfPressure       //  Rocket centre of pressure
 componentCount         //  Counter for numbering components
 component_count        //  Count vector for recording component data
 component_name         //  Cell array of component names
 component_mass         //  Component mass
 component_position     //  Component position along rocket longitudinal axis
 component_MOI          //  Moment of inertia per component
 MOInertia              //  MOI at launch
 MOInertia_noProp       //  MOI at propellant depletion
 MOI_method
 keep_adding            //  Used to request input from user for adding more components
 staticMargin           //  Static margin at launch
 stabCaliber            //  Rocket stability at launch, in calibers
 staticMargin_noProp    //  Static margin at propellant depletion
 stabCaliber_noProp     //  Rocket stability at propellant depletion, in calibers
 time_past              //  Time past during loop for animated trajectory
 timeStep               //  Time step for simulation chosen by user
 totalSteps             //  Total number of time steps taken in simulation
 simEndCount            //  Count vector used to find time of rocket landing
 simEndIndex            //  Index found using simEndCount loop

 stage_propMass         //  Mass of propellant consumed by the stage
 stage_ISP              //  Appropriate specific impulse for the propellant combination, 
                            chamber pressure, expansion ratio and ambient pressure.
 stage_thrust           //  Model assumes a CONSTANT THRUST throughout the stage
 stage_diameter         //  Cross-sectional diameter of the stage.
 stage_preCoast         //  Time lag between staging and ignition of next stage.
 stage_postCoast        //  Time lag between end of thrust and staging
 stage_impulse          //  Calculated stage impulse based on propellant mass and specific impulse.
 stage_burnTime         //  Calculated from impulse and thrust.
 stage_area             //  Calculated cross-sectional area.
 stage_inertMass        //  Formula to be edited to reflect the expected inert mass fraction
 stage_mass             //  Inert mass plus propellant mass for the stage.
 stage_thrustWeight     //  This ratio should be used to guide the choice of thrust magnitude.
 stage_massFlow         //  Calculated propellant mass flow rate for the stage.

 grossLiftoffMass       //  Gross lift-off mass
 numOfStages            //  Number of stages
 stage_count            //  Count vector for stages, for loop
 ignition               //  Ignition time event
 burnout                //  Burnout time event
 separation             //  Stage separation time event

 stageCounter           //  Counter for stages, for loop
 phaseCounter           //  Counter for phases, for loop
 phase                  //  
 phase_endTime          //  Time when current phase ends.
 phase_thrust           //  Phase thrust
 phase_area             //  Phase area
 phase_massFlow         //  Phase mass flow rate
 phase_mass             //  Phase mass
 phase_count            //  Count vector for phases, for loop

 sim_time               //  Time at the end of the current time step
 sim_thrust             //  Predicted rocket motor thrust corresponding to the current time
 sim_massFlow           //  Mass flow rate - calculated rate of proellant usage, relating to the thrust magnitude
 sim_area               //  Cross-sectional area of the launch vehicle, used to calculate drag.
 sim_mass1              //  Mass at start of time step.
 sim_mass2              //  Mass at end of time step.
 sim_gravity            //  gravitational acceleration calculated for current altitude
 sim_velocity1          //  Velocity magnitude (i.e. speed measured parallel to the velocity vector) at start of current step
 sim_angle1RAD          //  Inclination (pitch angle) measured from the horizontal in degrees at the start of the time step.
                            It is assumed that the rocket always points in the direction of flight, 
                            so the pitch angle = the direction of the velocity vector.
 sim_angle1DEG          //  Inclination angle converted from radians to degrees
 sim_range1             //  Down range distance at start of time step.
 sim_altitude1          //  Altitude at start of time step.
 sim_accel2             //  Acceleration parallel to velocity vector at end of time step.
 sim_velocity2          //  Speed parallel to velocity vector at end of time step.
 sim_machNum            //  Mach number at end of time step.
 sim_deltaAngle         //  Change in pitch angle in radians.
 sim_angle2RAD          //  Pitch angle in radians at the end of the time step.
 sim_angle2DEG          //  Pitch angle in degrees at the end of the time step.
 sim_angleAvg           //  Average inclination during the time increment in radians. 
                        //  This value is used to update the displacement vector
 sim_deltaDisplace      //  Change in magnitude of the displacement vector (measured parallel to the velocity vector).
 sim_range2             //  Down range distance at end of time step.
 sim_altitude2          //  Altitude at end of time step.
 sim_atmosDensity       //  Atmospheric density based on altitude at start of time step.
 sim_atmosTemp          //  Atmospheric temperature based on altitude at start of time step.
 sim_soundSpeed         //  Speed of sound based on altitude at start of time step.
 sim_dragCoef           //  Drag coefficient CD -- CALCULATIONS NEED VERIFICATION
 sim_dynPressure1       //  Dynamic pressure based on conditions at start of time step
 sim_dynPressure2       //  Dynamic pressure based on conditions at end of time step
 sim_dragForce          //  Drag force based on conditions at start of time step
 sim_accel1             //  Estimated acceleration parallel to flight based on initial conditions
 sim_convergence        //  Confirmation that the iterations have converged.
 sim_phaseEnd           //  Time corresponding to the end of the current phase.
 sim_flightPhase        //  Flight phase corresponding to the end of the current phase.
 sim_phaseMass          //  Steady mass during the particular phase of the flight.
 sim_potEnergy          //  Potential energy at a given time
 sim_kinEnergy          //  Kinetic energy at a given time
 sim_totalEnergy        //  Total energy (potential plus kinetic) at a given time
 sim_dynViscosity       //  Dynamic viscosity of air at a given point
 sim_Reynolds           //  Reynolds number at a given point
 dynViscosity           //  Array of data for dynamic viscosity wrt temperature

 Results.Apogee                 //  Apogee altitude and time to reach it
 Results.maxRange               //  Maximum down range distance and time it was reached
 Results.maxSpeed               //  Maximum speed and time it was reached
 Results.maxMachNum             //  Maximum Mach number and time it was reached
 Results.maxAccel               //  Maximum acceleration and time it was reached
 Results.maxDynPressure         //  Maximum dynamic pressure and time it was reached
 Results.maxDragForce           //  Maximum drag force experienced and time it was reached
 Results.max_Potential_Energy   //  Maximum potential energy and time reached
 Results.max_Potential_Energy   //  Maximum kinetic energy and time reached
 Results.max_Potential_Energy   //  Maximum total energy and time reached
 Results.Time_of_Flight         //  Time of flight

 graph_trajectory       //  Graph of the trajectory (altitude wrt range)
 graph_altitude         //  Graph of the altitude wrt time
 graph_speed            //  Graph of the speed wrt time
 graph_accel            //  Graph of the acceleration wrt time
 graph_mass             //  Graph of the mass wrt time
 graph_machNum          //  Graph of the Mach number wrt time
 graph_dragForce        //  Graph of the drag force experienced wrt time
 graph_DynPress         //  Graph of the dynamic pressure wrt time
 graph_atmosDens        //  Graph of the atmospheric density wrt time
 graph_angle            //  Graph of the inclination/pitch angle wrt time

 x_increment            //  X-axis increment used when calculating values for Y-axis in rocket design
 x_iArray               //  Array of X-axis values
 y_i                    //  Y-axis value for a given X-axis location, either above or below the X-axis
 y_i_projected          //  Equivalent projected  Y-axis value for fins at an angle
 y_i_total              //  Total Y-axis value for a given X-axis location above AND below the X-axis
 y_iArray               //  Array of y_i_total values
 COPcount               //  Count vector used in loop for calculating COP
 xAxisRange             //  Vector used to identify type of rocket part is present at each x-axis location
 xRangeCount            //  Used to remove zero values at end of xAxisRange array
 new_xPosition          //  Index used to identify start of a new length of rocket to assign a type of rocket part to in xAxisRange
 findTubeCount          //  Count vector used to indentify postion of last body tube, to assign fins to
 nonTube_startIndex     //  Index for use as above
 nonTube_endIndex       //  Index for use as above
 nonTube_indexRange     //  Difference between startIndex & endIndex, for use as above
 y_i_projOverlay        //  Array of values to plot for projected image of underside fin

 COP_options.Barrowman_Equations	//  Used for selecting Barrowman equations as the method for calculating COP
 COP_options.Centre_of_Area         //  Used for selecting centre of area as the method for calculating COP
 COP_selected                       //  Key value for the COP calculation method selected
 COP_numerator                      //  Numerator for quotient used in calculating COP changes with time
 COP_denominator                    //  Denominator, as above
 nose_normalForce                   //  Normal force on nose cone (for use in Barrowman equations)
 nose_COP                           //  COP of the nose cone (for use in Barrowman equations)
 trans_normalForce                  //  Normal force on transitions (for use in Barrowman equations)
 trans_COP                          //  COP of transitions (for use in Barrowman equations)
 fin_normalForce                    //  Normal force on fins (for use in Barrowman equations)
 fin_COP                            //  COP of the fins (for use in Barrowman equations)
 transPointLength                   //  Used in calculation to determine trans_COP
 finStartPoint                      //  Captures data on length of rocket from nose tip to fin leading edge

 bodyTubeNUM            //  Counter for number of body tubes
 transitionNUM          //  Counter for number of transitions
 finSetNUM              //  Counter for number of fin/canard sets
 bodyTubeDiameter       //  Diameter of a given body tube
 bodyTubeLength         //  Length of a given body tube
 trans_frontDiameter    //  Front diameter of a given transition/boat tail
 trans_rearDiameter     //  Rear diameter of a given transition/boat tail
 trans_length           //  Length of a given transition/boat tail
 numOfFins              //  Number of fins in the rocket design
 finAngle               //  Angle between each pair of fins
 finCount               //  Count vector used in loop(s) for calculating Y-axis values for fin edges
 fin_leadAngle          //  Inside angle between fin leading edge and root chord
 fin_span               //  Fin span length
 fin_sweepDist          //  Fin sweep (leading edge) length
 trans_frontLength      //  X-axis distance covered by fin leading edge
 trans_rearLength       //  X-axis distance covered by fin trailing edge
 fin_tipChord           //  Fin tip chord length
 fin_rootChord          //  Fin root chord length
 fin_trailAngle         //  Inside angle between fin trailing edge and root chord
 ellipseCounter         //  Counter used in determing elliptical fin shape points
 y_ellipse              //  Vector of Y-axis points for elliptical fins
 finShapeCount          //  Used to find a fin shape identifier in the part_identity cell array
 finShapeCounter        //  As above
 finXcoords             //  Custom fin design X-coordinate data
 finYcoords             //  Custom fin design Y-coordinate data
 finXpoint              //  Used to determine X-axis position to apply Y-coordinate data for fins
 finCoordsArray         //  Array of X/Y coordinates for custom fin designs
 finCoordCounter        //  Counter used in process of adding new fin coordinates
 coordsCount            //  Count vector used in loop for displaying current collection of fin coordinates
 xCoord                 //  As above - current X-coordinate being displayed
 yCoord                 //  As above - current Y-coordinate being displayed

 finShape.Trapezoidal               //  Used to select trapezoidal fins as the desired fin shape
 finShape.Elliptical                //  As above, but for elliptical fins
 finShape.Freeform_image_file       //  As above, but for importing an image of a custom fin design
 finShape.Freeform_coordinates      //  As above, but for entering custom fin design coordinates manually


 microG_array           //  Microgravity calculations array
 microG_count           //  Count vector for microgravity calculatuions, for loop
 microG_duration        //  Microgravity duration calculation
 microG_limit           //  Microgravity acceleration magnitude limit
 microG_limit_minus     //  Negative limit
 microG_limit_plus      //  Positive limit
 microG_6_time          //  Duration of time the rocket will experience microgravity of 1x10^-6g quality
 microG_5_time          //  As above (1x10^-5g quality)
 microG_4_time          //  As above (1x10^-4g quality)
 microG_3_time          //  As above (1x10^-3g quality)
 microG_2_time          //  As above (1x10^-2g quality)
 microG_1_time          //  As above (1x10^-1g quality)
 Quality                //  List of above qualities
 Duration               //  List of above durations
 microG_table           //  Table of above two lists
 microG_analysis        //  Used to determine if user wants to see analysis of rocket microgravity experience

 engineCounter      //  Used to count number of engine files added
 engineFile         //  Name of engine file
 engine             //  Cell array for engine data extracted from engine file
 thrustData         //  Cell data converted to numerical array
 thrustCount        //  Counter used to find appropriate thrust range to draw data from
 columnCount        //  Counter for column data in thrustData
 rowCount           //  Counter for row data in thrustData
 numOfEngines       //  Total number of engine files added
 firingCounter      //  Counter used to determine which engine file to use

 microgFig          //  Used to get the handle for a microgravity figure with multiple graphs
 maxResult_Names    //  Cell array of names for max output values
 maxResult_Cell     //  Cell array of names and data for max output values
 maxResult_Table    //  Table version of above cell
 burnLine           //  Animated line for section of trajectory when engine is on
 recoverLine        //  Animated line for section of trajectory beyond parachute

 COP_method         //  Used to determine whether COP will be entered manually or calculated
 COM_method         //  Used to determine whether COM will be entered manually or calculated
 staticStability    //  Used to determine whether static stability should be analysed.
 instabilityArray   //  Array of ones/zeros, used to determine when rocket is stable/unstable during flight [WIP]
 stabilityCount     //  Count vector used to determine the above [WIP]
 unstableTime       //  Total time during flight that the rocket will be unstable for [WIP]
 sim_accuracy       //  Estimate of accuracy of flight simulation based on per-cent time that rocket is stable for [WIP]

 numOfChutes            //  Number of parachutes
 chuteDiameter          //  Diameter of a given parachute
 chuteDragCoef          //  Drag coefficient of a given parachute
 chuteDeployMethod      //  Deployment method for a given parachute
 chuteTimeDelay         //  Time delay (seconds) for a given parachute
 chuteSetAltitude       //  Altitude during descent for a given parachute deployment
 chuteCounter           //  Used to determine which number parachute is deploying at a given time
 chuteCount             //  Count vector for all parachutes, used to gather data on parachute attributes
 chuteDrag              //  Drag on all deployed parachutes, calculated without dynamic pressure
 deployTime             //  Records times at which parachutes deploy
 parachute.Apogee       //  Key for selecting deployment at apogee

 parachute.Time_delay_post_launch               //  Key for selecting deployment after time delay post-launch
 parachute.Time_delay_post_previous_deployment  //  Key for selecting deployment after time delay post-previous chute deployment
 parachute.Set_altitude_during_descent          //  Key for selecting deployment at a set altitude during descent

 vertLineCounter        //  Counter used to count number of instances to draw a vertical line to separate rocket components
 rocketDraw_count       //  Count vector used for loop to draw above described vertical lines
 rocketDraw_vertLineX   //  X-axis coordinates for vertical lines
 rocketDraw_vertLineY   //  Y-axis coordinates for vertical lines
 skipDescentAnim        //  Used to determine whether to skip plotting the animated descent
 deployArray1           //  Used to find the index within the sim_time array that 
                            corresponds to the deployment time of the first parachute
 deployArray2           //  As above
 deployCount1           //  As above
 deployCount2           //  As above
 deployIndex            //  As above
 COMpoint1              //  Plot point for the centre of mass at launch
 COMpoint2              //  As above
 COMpoint3              //  As above
 COMpointText           //  Plot text for COM

 Transformed_COP        //  Centre of pressure location with coordinates transformed 
                        //  to show its position at a given time during rocket FPA visualization.
 TransformMatrix        //  Matrix used to transform coordinates
 Ty_i                   //  Tranformed rocket body coordinates (upper surface)
 Ty_i_projected         //  Tranformed rocket body coordinates (lower surface)
 TxyEnd                 //  Tranformed rocket body coordinates (rear end)
 COMpoint               //  Legend label for transformed COM
 COPpoint               //  Legend label for transformed COP
 xticklabel             //  Axis labels for transformed rocket body X-coordinates
 yticklabel             //  Axis labels for transformed rocket body Y-coordinates
 xyTick                 //  Used to size the axis tick marks
 simRunTime             //  Time elapsed whilst running a given process
 simZeroCount           //  Count vector used in loop for setting negative altitude/speed vealues to zero

 xEnd_firing            //  X-coordinate array used in animating motor plume
 yEnd_firing            //  Y-coordinate array used as above
 TxyEnd_firing          //  Transformed coordinates array used as above
 plumeSpreadFctr        //  Imaginary factor used to animate spread in plume
 fireCounter            //  Counter used to alternate between narrow/spread plume drawing
 xChute                 //  X-coordinate array used in drawing a parachute
 yChute                 //  Y-coordinate array used as above
 TxyChute               //  Transformed coordinates array for drawing first deployed chute
 TxyCords               //  Transformed coordinates array for drawing chute cords
 distToOriginX
 distToOriginY

 inputsFromFile         //  Gives user option to import simulation input information from a pre-made file
 inputsFileName         //  File name for the above
 fileNameCounter        //  Counter used to increment data set number when creating new inputs data file
 inputsSetNum           //  Set number from the above
 inputsData             //  Array of imported data from the above
 rocketFromFile         //  Gives user option to import rocket shape design data from file
 rocketName             //  File name for the above
 rocketData             //  Array of imported data from the above
 rocketDataSize         //  Size of array of the above
 rocketCount            //  Count vector used in loop for generating a rocket design from the imported data
 rocketIntFromFile      //  Gives user option to import rocket design components data from file
 rocketIntName          //  File name for the above
 rocketIntData          //  Array of imported data from the above
 recoveryCount          //  Count vector used in loop for generating 
                            recovery information from the imported parachute data

 pauseFix           //  Fix for pause issue with R2016a (hangs otherwise)
 simRunNum          //  Number of the current simulation run (accounts for repeat runs)
 expFlightDur       //  Expected (predicted) flight duration
 R_length           //  Total rocket length, from nose tip to base
 int_range          //  Initial range value (usually zero)
 int_altitude       //  Initial altitude value (usually zero)
 pitchAngle         //  Launch inclination angle
 int_velocity       //  Initial velocity at launch (usually zero)
 towerHeight        //  Height/length of the launch tower/rails
 payloadMass        //  Payload mass (usually set to zero)
 yawAngle           //  Angle (in the horizontal plane) between the rocket length and wind velocity vector
 wind_velocity      //  Mean/average wind velocity
 gust_velocity      //  Wind gust velocity
 angleDiff          //  Angle moved through per time increment due to action of wind
 angularVelo        //  Angular velocity of the rocket due to wind effects
 angularAccel       //  Angular acceleration of the rocket due to wind effects
 numOfLaunchLugs    //  Number of launch lugs used in rocket design
 
 metricLength       //  Length dimension selected for use throughout simulation
 metricMass         //  Mass dimension selected for use throughout simulation
 flight_analysis    //  Option to conduct a flight analysis (obselete)
 simDescriptor      //  Description of simulation input data file being constructed
 descriptFileID     //  ID associated with the above
 sameInputs         //  Option to use same inputs data on repeat simulation runs
 descentPhase       //  Used to determine if the rocket is descending after reaching apogee
 
 deployChute        //  Determines whether to deploy a given parachute at a given time
 frustrumAngle      //  Angle used in calculations for certain nosecone designs
 delta_noseLength   //  Partial length of nosecone, used as above
 skinRoughness      //  Coefficient of surface roughness
 finThickness       //  Thickness of fins
 finArea            //  Surface area of fins
 finSpanTotal       //  Span from fin-to-fin

 sim_windVelocity   //  Wind velocity over time
 sim_angleOfAttack  //  Angle of attack over time
 sim_liftCoef       //  Lift coefficient over time
 sim_windLift       //  Lift due to wind over time
 sim_COP            //  Current COP during flight simulation
 sim_COM            //  Current COM during flight simulation
 sim_MOI            //  Current MOI during flight simulation
 sim_stability      //  Current (dynamic) stability (calibers) during flight simulation
 sim_dragCoefO      //  Current parasite/zero-lift drag coefficient during flight simulation
 sim_dragCoefI      //  Current lift-induced drag coefficient during flight simulation

 massCount          //  Count vector used in calculating current COM & MOI
 freeformFins       //  Used for determining if the fis are freeform designs
 setDragCoef        //  Requests user to set a constant parasite drag coefficient value 
                    //  in scenarios where the function isn't applicable
 finFile            //  File ID for importing monochromatic bitmap fin images
 finPic             //  Variable containing data from the above
 sizeFinPic         //  Size of the array data from the above
 fin_area           //  Fin planform area (includes attached body tube)
 chuteDataStart     //  Used when importing recovery system data
 tubeCounter        //  Counter used for the above
 frustumCounter     //  Counter used when importing transition data
 finsCounter        //  Counter used when importing fin design data
 slantLength        //  
 partSurfaceArea    //  Surface area calculated for a particular rocket component
 bodybaseDiameter   //  Base diameter of the rocket
 R_effLength        //  'Effective' length of the rocket
 bodyMaxDiameter    //  Maximum rocket diameter
 R_aftLength        //  'Aft' length of the rocket
 aftLengthCount     //  Count vector used with the above
 R_frontLength      //  'Front' length of the rocket

 launchLugCount     //  Count vector used in gathering launch lug data
 protub_length      //  Protuberance (e.g. launch lugs) length
 protub_dist        //  Protuberance distance from a reference point
 protub_crossArea   //  Protuberance cross-sectional area
 protub_surfArea    //  Protuberance surface area

 windEffect_method      //  Used when choosing to randomize wind speeds or not
 kCell                  //  Used in determing parameters for Weibull wind speed distribution
 kCount                 //  As above
 xWind                  //  Linearly-spaced vector of wind speeds from 0 to gust speed
 shapeParameter         //  Shape parameter of the Weibull distribution
 scaleParameter         //  Scale parameter of the Weibull distribution
 windCount              //  Count vector used to calculate random wind speed values
 windVelocity_preAlt    //  Wind velocity vector before accounting for altitude effects (shear)
 dragCalcMethod         //  Used when choosing which drag calculation method to employ
 
 CDcalcWarnID       //  Warning ID for when parasite drag function shouldn't be used
 CDcalcWarnMSG      //  Warning message for the above
 critWarnID         //  Warning ID for when parasite drag function may be inaccurate
 critWarnMSG        //  Warning message for the above
 rocketDataWarnID   //  Warning ID for when parasite drag function cannot be used
 rocketDataWarnMSG  //  Warning message for the above
 
 failureChoice      //  Used when choosing whether to simulate failure mode(s)
 failureOption      //  Used to select option for how failure mode(s) will be simulated
 failDescriptor     //  Descriptions of failure modes
 failureMode        //  Actual failure mode type
 failureChance      //  Chance of a failrue mode occuring
 failTypeChance     //  Chance of a particular type of failure mode occuring
 motorFailTime      //  Randomized time when motor experiences a failure mode
 
 rocketCompsFromFile    //  Option used to select importing rocket interior components design file
 rocketCompsName        //  Name of the above file
 namesFileID            //  ID to to identify the above file
 rocketCompsCount       //  Count vector used to assign data from the above imported file
 
 moveableFinSet         //  Set number for fins that will be moveable (control system)
 control_system         //  Used to choose whether to simulate a control system or not
 controlSystem          //  Used to choose a type of control system to simulate
 finDeflection          //  Current deflection of fins during flight
 finRotationRate        //  Rate at which fins can be rotated
 maxDeflection          //  Maximum deflection angle for fins or gimballed nozzle
 gimbalingRate          //  Rate at which nozzle can be gimballed
 sim_thrustAngle        //  Angle of the thrust plume at a given time
 COPvsAOA               //  
 moveableFinAOA         //  Angle of attack of the fins
 moveableFinCOP         //  Specific COP of the fins alone
 moveableFinLift        //  Lift acting on moveable fins at a given time
 moveableFinLiftCoef    //  Lift coefficient of the above
 gimbalCount            //  Count vector used to draw rocket plume at a given angle

 aspectRatio            //  Aspect ratio of fins
 ozFactor               //  Ozwald factor for the above
 saveVelocity           //  Used to save current velocity temporarily
 defaultDragCoefO       //  Default value for parasite drag

 windSpeedAxis          //  Linearly-spaced vector of wind speed values
 probDens_Fnct          //  Weibull probability density function
 probDens_Data          //  Data calculated using Weibull probability 
 windSpeedAxis_ext      //  Extended version of windSpeedAxis (100 to 1000)
 cumuProbDens_Fnct      //  Cumulative Weibull probability density function
 windStandardDev        //  Estimated wind speed standard deviation
 sigma                  //  Used to calculate the above
 landingPosition        //  Final range value, i.e. landing position
 circleAngles           //  Used to create drawing of circular landing zones
 circlePoints           //  Used as above
 stdDevCount            //  Count vector used to calculate probabailities of
                        //  rocket falling within particular landing zones
 SDplusIndex            //  Used to calculate landing positions (considering
                        //  deviations) forward of expected landing spot
 SDminusIndex           //  As above, but backwards from expected position
 SDprob                 //  Probabilities of rocket landing within particular landing zones
 landingRadius          //  Radius of landing zones
 maxMachIndex           //  Index for array value for maximum Mach number
 maxReynoldsIndex       //  As above, but for Reynolds number
 sim_ballisticCoef      //  Ballistic coefficient calculated for ascent

 frame                  //  Retrieves frame during rocket animation
 frameImage             //  Saves the above as an image
 filename               //  Name for GIF file
 gifCount               //  Count vector used to create GIF file
 indexedImage           //  Indexed version of images used in GIF
 colourMap              //  Colour map data used for the above
 
%}


%% INITIALIZATION [ini02]
%%
warning('off','MATLAB:hg:AutoSoftwareOpenGL')
reRunSim = true;
simRunNum = 0;
%
%%%
while reRunSim == true % un-comment this line when integrating scripts %%%DEV:NB%%%
simRunNum = simRunNum+1;
%%%
%
if  (exist(['Inputs_and_Responses_(',num2str(simRunNum),').txt'],'file')==2)
    fopen(['Inputs_and_Responses_(',num2str(simRunNum),').txt'],'w');
end

diary(['Inputs_and_Responses_(',num2str(simRunNum),').txt'])
diary on
lastwarn('')
               
%%

keep_adding = true;
global yes
global no
yes = {'Y','y'};
no = {'N','n'};

%%
if simRunNum == 1
%%% %%%
[numOfStages,...
expFlightDur,...
timeStep,...
R_length,...
int_range,...
int_altitude,...
pitchAngle,...
int_velocity,...
towerHeight,...
payloadMass,...
numOfChutes,...
dragCalcMethod,...
yawAngle,...
wind_velocity,...
gust_velocity,...
angleDiff,...
angularVelo,...
angularAccel,...
numOfLaunchLugs ]=deal(0);


fprintf('%s\n','The options for length dimensions when describing your rocket design are as such:',...
            '0. Metres (m)','1. Decimetres (dm)','2. Centimetres (cm)','3. Millimetres (mm)')
        metricLength = nan;
    while isnan(metricLength) == true
        metricLength = input('Enter a number corresponding to a length dimension to use based on the above key: ');
        switch metricLength
            case {0,1,2,3}
                % placeholder
            otherwise
                disp('You failed to enter a value within the limits of the provided key. Try again.')
                metricLength = nan;
        end
    end
    
    
    fprintf('%s\n','The options for mass dimensions when describing your rocket design are as such:',...
            '0. Kilograms (kg)','3. Grams (g)')
        metricMass = nan;
    while isnan(metricMass) == true
        metricMass = input('Enter a number corresponding to a mass dimension to use based on the above key: ');
        switch metricMass
            case {0,3}
                % placeholder
            otherwise
                disp('You failed to enter a value within the limits of the provided key. Try again.')
                metricMass = nan;
        end
    end
    
%%% %%%
end

%%
    %flight_analysis =  max(strcmp((input('Do you want a flight simulation analysis for the rocket? [Y/N]: ','s')),yes));
    flight_analysis = true;
if  flight_analysis == true
    
	inputsFromFile =  max(strcmp((input('Do you want to import simulation input data from a previously generated file? [Y/N]: ','s')),yes));
if  inputsFromFile == false
    
    [numOfStages] = changeInput(simRunNum,'number of stages','How many stages are in your rocket? ',numOfStages);
    
    [expFlightDur] = changeInput(simRunNum,'predicted flight duration','How long (maximum) do you expect the flight to last (s)? ',expFlightDur);
    disp(['The recommended time step is ',num2str(expFlightDur/10000),' seconds, based on 10000 steps.']);
    
    [timeStep] = changeInput(simRunNum,'time step','What time step do you want to use? ',timeStep);
    
    
    inputsFileName = [];
    fileNameCounter = 1;
    
    while isempty(inputsFileName)
        if      exist(['Simulation Input Data (Set ',num2str(fileNameCounter),').txt']) == 2
                inputsFileName = [];
                fileNameCounter = fileNameCounter + 1;
            
        elseif  exist(['Simulation Input Data (Set ',num2str(fileNameCounter),').txt']) == 0
                inputsFileName = ['Simulation Input Data (Set ',num2str(fileNameCounter),')'];
                dlmwrite([inputsFileName,'.txt'],[numOfStages,expFlightDur,timeStep],'delimiter',',');
        end
    end
    
    simDescriptor = input('Type a description for these simulation settings: ','s');
    if      exist('Simulation Descriptors.txt','file') == 0
            descriptFileID = fopen(['Simulation Descriptors','.txt'],'wt');
    elseif  exist('Simulation Descriptors.txt','file') == 2
            descriptFileID = fopen(['Simulation Descriptors','.txt'],'at');
    end
    fprintf(descriptFileID,'%s\n',['Set ',num2str(fileNameCounter),': ',simDescriptor]);
    fclose(descriptFileID);

    
elseif inputsFromFile == true
    
        sameInputs = false;
    if  simRunNum > 1
        sameInputs =  max(strcmp((input('Use the same simulation input values? [Y/N]: ','s')),yes));
    end
    
    if  sameInputs == false
        inputsSetNum = input('Enter the set number (see file name) for the inputs data file to use: ');
        inputsFileName = ['Simulation Input Data (set ',num2str(inputsSetNum),')'];
    end
    
    inputsData = dlmread([inputsFileName,'.txt']);
    
    numOfStages     = inputsData(1,1);
    expFlightDur    = inputsData(1,2);
    timeStep        = inputsData(1,3);
    
end % {inputs from file}

    numOfSteps = expFlightDur/timeStep;

end % {flight analysis}
%%
    %staticStability =  max(strcmp((input('Do you want an analysis of the rocket''s static stability? [Y/N]: ','s')),yes));
    staticStability = true;
    
if  simRunNum == 1
%%% %%%
    %%
    if  staticStability == true
    
[   component_mass,...
    component_position,...
    component_isProp,...
    deployChute,...
    nose_COP,...
    nose_normalForce,...
    noseFactor,...
    frustrumAngle,...
    delta_noseLength,...
    skinRoughness   ] = deal(0);

    %%
    
[   bodyTubeDiameter,...
    bodyTubeLength,...
    ...
	trans_frontDiameter,...         
	trans_rearDiameter,...
	trans_length,...
    ...
	numOfFins,...
	fin_distFromBase,...
	fin_span,...
	fin_sweepDist,...
	fin_tipChord,...
	fin_rootChord,...
    finThickness,...
    finAngle,...
    finTubeAssign,...
    finArea,...
    finSpanTotal,...
    trans_normalForce,...
    trans_COP,...
	fin_normalForce,...
    fin_COP     ] = deal(zeros(1,10));

    end

    %%
    if  flight_analysis == true
    
[   stage_propMass,...
    stage_ISP,...
    stage_thrust,...
    stage_diameter,...
    stage_preCoast,...
    stage_postCoast,...
    stage_impulse,...
    stage_burnTime,...
    stage_area,...
    stage_inertMass,...
    stage_mass,...
    stage_thrustWeight,...
    stage_massFlow,...
    ...
	ignition,...
	burnout,...
	separation	] = deal(zeros(1,numOfStages));


[   phase,...
    phase_endTime,...
    phase_thrust,...
    phase_area,...
    phase_massFlow,...
    phase_mass  ] = deal(zeros(1,(3*numOfStages)));

    %%

[   sim_time,...
    sim_thrust,...
    sim_massFlow,...
    sim_area,...
    sim_mass1,...
    sim_mass2,...
    sim_gravity,...
    sim_velocity1,...
    sim_angle1RAD,...
    sim_angle1DEG,...
    sim_range1,...
    sim_altitude1,...
    sim_accel2,...
    sim_velocity2,...
    sim_machNum,...
    sim_deltaAngle,...
    sim_angle2RAD,...
    sim_angle2DEG,...
    sim_angleAvg,...
    sim_deltaDisplace,...
    sim_range2,...
    sim_altitude2,...
    sim_atmosDensity,...
    sim_atmosTemp,...
    sim_soundSpeed,...
    sim_dragCoefO,...
    sim_dragCoefI,...
    sim_dynPressure1,...
    sim_dynPressure2,...         
    sim_dragForce,...
    sim_accel1,...
    sim_convergence,...
    sim_phaseEnd,...
    sim_flightPhase,...
    sim_phaseMass,...
    sim_potEnergy,...
    sim_kinEnergy,...
    sim_totalEnergy,...
    sim_dynViscosity,...
    sim_Reynolds,...
    sim_windVelocity,...
    sim_angleOfAttack,...
    sim_liftCoef,...
    sim_windLift,...
    sim_COP,...
    sim_COM,...
    sim_MOI,...
    sim_stability   ] = deal(zeros(1,(expFlightDur/timeStep)));

    %%
    
[   chuteDiameter,...
    chuteDragCoef,...
    chuteDeployMethod,...
    chuteTimeDelay,...
    chuteSetAltitude    ] = deal(zeros(1,10));
    
    end
    
    %%
    
elseif simRunNum > 1

    if  flight_analysis == true
        
[   sim_time,sim_thrust,sim_massFlow,sim_area,...
    sim_mass1,sim_mass2,sim_gravity,sim_velocity1,...
    sim_angle1RAD,sim_angle1DEG,sim_range1,sim_altitude1,...
    sim_accel2,sim_velocity2,sim_machNum,sim_deltaAngle,...
    sim_angle2RAD,sim_angle2DEG,sim_angleAvg,sim_deltaDisplace,...
    sim_range2,sim_altitude2,sim_atmosDensity,sim_atmosTemp,...
    sim_soundSpeed,sim_dragCoefO,sim_dynPressure1,sim_dynPressure2,...         
    sim_dragForce,sim_accel1,sim_convergence,sim_phaseEnd,...
    sim_flightPhase,sim_phaseMass,sim_potEnergy,sim_kinEnergy,...
    sim_totalEnergy, sim_dynViscosity, sim_Reynolds, sim_windVelocity,...
    sim_angleOfAttack, sim_liftCoef, sim_windLift, sim_dragCoefI,...
    sim_COP, sim_COM, sim_MOI, sim_stability ...
    ...
    ] = extendSimData(expFlightDur, timeStep, simRunNum, numOfSteps, ...
    ...
    sim_time,sim_thrust,sim_massFlow,sim_area,...
    sim_mass1,sim_mass2,sim_gravity,sim_velocity1,...
    sim_angle1RAD,sim_angle1DEG,sim_range1,sim_altitude1,...
    sim_accel2,sim_velocity2,sim_machNum,sim_deltaAngle,...
    sim_angle2RAD,sim_angle2DEG,sim_angleAvg,sim_deltaDisplace,...
    sim_range2,sim_altitude2,sim_atmosDensity,sim_atmosTemp,...
    sim_soundSpeed,sim_dragCoefO,sim_dynPressure1,sim_dynPressure2,...         
    sim_dragForce,sim_accel1,sim_convergence,sim_phaseEnd,...
    sim_flightPhase,sim_phaseMass,sim_potEnergy,sim_kinEnergy,...
    sim_totalEnergy,sim_dynViscosity,sim_Reynolds,sim_windVelocity,...
    sim_angleOfAttack,sim_liftCoef,sim_windLift,sim_dragCoefI,...
    sim_COP,sim_COM,sim_MOI,sim_stability);

    end

%%% %%%
end
%%

[   microG_6_time,...
    microG_5_time,...
    microG_4_time,...
    microG_3_time,...
    microG_2_time,...
    microG_1_time,...
    realTime,...
    part_count,...
    vertLineCounter,...
    time_past,...
    freeformFins,...
    skipDescentAnim,...
    bodyTubeNUM,...
    transitionNUM,...
    finSetNUM,...
    chuteDataStart,...
    componentCount,...
    y_iArray,...
    setDragCoef,...
    sim_dragCoefO(1)     ]=deal(0);

[animCounter,...
 fpaCounter]= deal(1);
errorMargin = 0.01;
tCounter    = zeros(1,2);
deployTime  = zeros(1,10);
chuteCounter = ones(1,(numOfSteps+1));


    rocket.Nose = 1;
    rocket.Body_Tube = 2;
    rocket.Canards = 'See Fins';
    rocket.Transition = 3;
    rocket.Fins = 4;
    rocket.Boat_Tail = 'See Transition';
    
	noseShape.Conical = 1;  noseShape.Bi_conic = 2;
	noseShape.Tangent_ogive = 3;    noseShape.Secant_ogive = 4;
	noseShape.Ellipsoid = 5;    noseShape.Power_series = 6;
	noseShape.Parabolic_series = 7; noseShape.Haack_series = 8;

	finShape.Trapezoidal = 1; finShape.Elliptical = 2;
	finShape.Freeform_image_file = 3; finShape.Freeform_coordinates = 4;
    
    parachute.Apogee = 1;
    parachute.Time_delay_post_launch = 2;
    parachute.Time_delay_post_previous_deployment = 3;
    parachute.Set_altitude_during_descent = 4;
                

%% Rocket Design [roc03]
%%
%
if  staticStability == true
%%%
fprintf('%s\n','You have the option to either describe the parameters of your rocket design,',...
        'or import data from a file generated after previously going through this process.')
    
    rocketFromFile =  max(strcmp((input('Do you want to import rocket design data from a previously generated file? [Y/N]: ','s')),yes));
if  rocketFromFile == false
    %%
    
    [R_length] = (changeInput(simRunNum,'rocket length','Enter a value for the length of your rocket: ',R_length))*10^-metricLength;
    
    fprintf('%s\n','Surface Finish Options','1: Totally smooth surface','2: Polished metal or wood','3: Natural sheet metal',...
            '4: Smooth matte paint finish','5: Standard camouflage paint finish')
	[skinRoughness] = changeInput(simRunNum,'surface finish','Enter a value based on the above key for the rocket''s surface finish:  ',skinRoughness,[1,2,3,4,5]);

        rocketName = [];
        while isempty(rocketName) == true
            try
                rocketName = input('Type a name for your rocket design: ','s');
                dlmwrite([rocketName,'.txt'],[R_length,skinRoughness],'delimiter',',');
            catch ME
                if  strcmp(ME.identifier,'MATLAB:dlmwrite:FileOpenFailure')
                    warning('Illegal characters used in file name. Try again.')
                    rocketName = [];
                end
            end
        end
        
        switch  skinRoughness
            case 1
                skinRoughness = 0.00001;
            case 2
                skinRoughness = (0.00002+0.00008)/2;
            case 3
                skinRoughness = 0.00016;
            case 4
                skinRoughness = 0.00025;
            case 5
                skinRoughness = (0.0004+0.0012)/2;
        end
        
    x_TotalSteps = 10000;
    x_increment = R_length/x_TotalSteps;

    xAxisRange = zeros(1,x_TotalSteps);
    new_xPosition = 1;


    disp('From nose to nozzle, describe the structure of your rocket using the key below: ')

    while keep_adding == true
   
        disp(rocket)
    
        part_count = part_count+1;
    
        part_identity{1,part_count} = input('Add a rocket part using the above key, entering a number from 1 to 4: ');
    
        
        if  (new_xPosition >= length(xAxisRange)) && (part_identity{1,part_count} ~= 4)
            error('You attempted to add a non-fin/canard component beyond the length of the rocket.')
        end
        
    
        switch part_identity{1,part_count}
        
            case 1  % Nose cone
                
                disp(noseShape)
                part_identity{2,part_count} = input(['Add a nose cone shape using the above key,'...
                                                     'entering a number from 1 to 8 corresponding to the desired shape: ']);

                switch part_identity{2,part_count}
                    case 4
                        fprintf('%s\n','The nose shape parameter (ogive radius) can be almost any value for a Secant Ogive nose,', ...
                                'but must at least be greater than half the nose length.');
                    case 6
                        fprintf('%s\n','The nose shape parameter must be one of the following for a Power Series nose:',...
                                'Cone = 1','3/4 Power = 0.75','1/2 Power (Parabola) = 0.5','Cylinder = 0');
                    case 7
                        fprintf('%s\n','The nose shape parameter must be one of the following for a Parabolic Series nose:',...
                                'Cone = 0','1/2 Parabola = 0.5','3/4 Parabola = 0.75','Parabola = 1');
                    case 8
                        fprintf('%s\n','The nose shape parameter must be one of the following for a Haack Series nose:',...
                                'LV-Haack = 1/3','LD-Haack (Von Karman (Ogive)) = 0');
                end
                
                switch part_identity{2,part_count}
                    case {4,6,7,8}
                        nose_shapePrmtr = input('Enter a value for the nose shape parameter based on the above key: ');
                    otherwise
                        nose_shapePrmtr = nan;
                end
                
                if (part_identity{2,part_count} == noseShape.Secant_ogive) && (nose_shapePrmtr > (nose_length/2)) == 0
                    error('You entered a value for the ogive radius outside of the allowable range.')
                end
                
                if  part_identity{2,part_count} == noseShape.Bi_conic
                    nose_baseRadius{1}  =(input('Enter a value for the nose cone base radius: '))*10^-metricLength;
                    nose_baseRadius{2}  =(input('Enter a value for the frustrum  base radius: '))*10^-metricLength;
                    nose_length{1}      =(input('Enter a value for the nose cone length: '))*10^-metricLength;
                    nose_length{2}      =(input('Enter a value for the frustrum  length: '))*10^-metricLength;
                else
                    nose_baseRadius     =(input('Enter a value for the nose base radius: '))*10^-metricLength;
                    nose_length         =(input('Enter a value for the nose length: '))*10^-metricLength;
                end
                    nose_tipRadius      =(input('Enter a value for the nose tip radius: '))*10^-metricLength;
            
                xAxisRange( new_xPosition : floor(nose_length/x_increment) ) = rocket.Nose;
                new_xPosition = floor(nose_length/x_increment) + 1;
                
                dlmwrite([rocketName,'.txt'],...
                    [part_identity{1,part_count},part_identity{2,part_count},nose_shapePrmtr,...
                    nose_baseRadius,nose_length,nose_tipRadius],...
                    '-append','delimiter',',');
            
            
            case 2  % Body tube
            
                bodyTubeNUM = bodyTubeNUM + 1;
                disp(['This is for body tube No.',num2str(bodyTubeNUM),'.'])
            
                bodyTubeDiameter(bodyTubeNUM)    =(input('Enter a value for the body tube diameter: '))*10^-metricLength;
                bodyTubeLength(bodyTubeNUM)      =(input('Enter a value for the body tube length: '))*10^-metricLength;
            
                xAxisRange( new_xPosition : new_xPosition + floor(bodyTubeLength(bodyTubeNUM)/x_increment) ) = rocket.Body_Tube;
                new_xPosition = new_xPosition + floor( bodyTubeLength(bodyTubeNUM)/x_increment ) + 1;
                
                dlmwrite([rocketName,'.txt'],...
                    [part_identity{1,part_count},bodyTubeDiameter(bodyTubeNUM),bodyTubeLength(bodyTubeNUM)],...
                    '-append','delimiter',',');
            
        
            case 3  % Transition (incl. boat tails)
            
                transitionNUM = transitionNUM + 1;
                disp(['This is for transition No.',num2str(transitionNUM),'.'])
            
                trans_frontDiameter(transitionNUM) =(input('Enter a value for the transition front diameter: '))*10^-metricLength;            
                trans_rearDiameter(transitionNUM)  =(input('Enter a value for the transition rear diameter: '))*10^-metricLength;
                trans_length(transitionNUM)        =(input('Enter a value for the transition length: '))*10^-metricLength;
            
                xAxisRange( new_xPosition : new_xPosition + floor(trans_length(transitionNUM)/x_increment) ) = rocket.Transition;
                new_xPosition = new_xPosition + floor( trans_length(transitionNUM)/x_increment ) + 1;
                
                dlmwrite([rocketName,'.txt'],...
                    [part_identity{1,part_count},trans_frontDiameter(transitionNUM),...
                    trans_rearDiameter(transitionNUM),trans_length(transitionNUM)],...
                    '-append','delimiter',',');
            
        
            case 4  % Fins (incl. canards)
            
                finSetNUM = finSetNUM + 1;
                finTubeAssign(finSetNUM) = bodyTubeNUM;
                disp(['This is for fin/canard set No.',num2str(finSetNUM),'.'])
                
                disp(finShape)
                disp('Note: Freeform fins can be used by either importing an image of the fin design (option 3) or entering coordinates manually (option 4)')
                part_identity{2,part_count} = input(['Add a fin shape using the above key, '...
                                                     'entering a number from 1 to 4 corresponding to the desired shape: ']);

                                                 
                numOfFins(finSetNUM)           = input('Enter a value for the number of fins or canards around the tube: ');
                fin_distFromBase(finSetNUM)    =(input('Enter a value for the fin/canard set''s distance from the base of the tube: '))*10^-metricLength;
                finThickness(finSetNUM)        =(input('Enter a value for the thickness of the fins: '))*10^-metricLength;
                fin_rootChord(finSetNUM)       =(input('Enter a value for the fin/canard root chord: '))*10^-metricLength;
                
                
                switch part_identity{2,part_count}
                    
                    case {finShape.Trapezoidal,finShape.Elliptical}
                        
                        fin_span(finSetNUM)    =(input('Enter a value for the fin/canard span: '))*10^-metricLength;
                
                        if (part_identity{2,part_count} == finShape.Trapezoidal)
                            fin_sweepDist(finSetNUM) =(input('Enter a value for the fin/canard sweep distance: '))*10^-metricLength;
                            fin_tipChord(finSetNUM)  =(input('Enter a value for the fin/canard tip chord: '))*10^-metricLength;
                            fin_midChord(finSetNUM)  =(input('Enter a value for the fin/canard mid-chord line length: '))*10^-metricLength;
                        else
                            [fin_sweepDist(finSetNUM),fin_tipChord(finSetNUM),fin_midChord(finSetNUM)] = deal(nan);
                        end
                        
                        [finXcoords,finYcoords] = deal(nan);
                        
                        
                    case finShape.Freeform_image_file
                        
                        [fin_span(finSetNUM),fin_sweepDist(finSetNUM),...
                         fin_tipChord(finSetNUM),fin_midChord(finSetNUM)] = deal(nan);
                     
                        clear finCoordsArray
        
                        disp(['The imported image needs to be a monochrome (fin black, background white) bitmap file,'... 
                              ' with the root chord flush with the image base.'])
                        finFile = input(['Type the file name for your set (',num2str(finSetNUM),...
                                         ') fin monochrome bitmap image file (with the extension): '],'s');
                        finPic = ~(imread(finFile));
                        sizeFinPic = size(finPic);
                        finCoordCounter = 0;
        
                        for     finCoordCount1 = 1:sizeFinPic(2) % X-coordinate data
                            for finCoordCount2 = 1:sizeFinPic(1) % Y-coordinate data
                
                                if  finPic(finCoordCount2,finCoordCount1) == 1
                                    % logical 1 indicates colour black for fin image data
                                    finCoordCounter = finCoordCounter + 1;
                                    finCoordsArray{finCoordCounter,1} = finCoordCount1;
                                    finCoordsArray{finCoordCounter,2} = (sizeFinPic(1)-finCoordCount2);
                                    break
                                end
                
                            end
                        end
        
                        finCoordsArray = cell2mat(finCoordsArray);
                        finXcoords = finCoordsArray(:,1) - finCoordsArray(1,1);
                        finYcoords = finCoordsArray(:,2);

                        finYcoords = finYcoords*(fin_rootChord(finSetNUM)/finXcoords(end));
                        finXcoords = finXcoords*(fin_rootChord(finSetNUM)/finXcoords(end));
        
        
                    case finShape.Freeform_coordinates
                        
                        [fin_span(finSetNUM),fin_sweepDist(finSetNUM),...
                         fin_tipChord(finSetNUM),fin_midChord(finSetNUM)] = deal(nan);
                     
                        clear finCoordsArray
        
                        keep_adding = true;
                        finCoordCounter = 0;
                        [finCoordsArray{1,1},finCoordsArray{2,1}]=deal(0);
        
                        while   keep_adding == true
            
                                finCoordCounter = finCoordCounter + 1;
                                finCoordsArray{1,finCoordCounter+1} = input('Enter a value for a fin/canard X coordinate: ');
                                finCoordsArray{2,finCoordCounter+1} = input('Enter a value for a the associated fin/canard Y coordinate: ');
        
                                keep_adding = max(strcmp((input('Do you want to add another set of fin/canard coordinates? [Y/N]: ','s')),yes));
                                
                                for coordsCount = 1:length(finCoordsArray)
                                    xCoord = finCoordsArray{1,coordsCount};
                                    yCoord = finCoordsArray{2,coordsCount};
                                    disp(['X',num2str(coordsCount-1),' = ',num2str(xCoord),'; Y',num2str(coordsCount-1),' = ',num2str(yCoord)])
                                end
                        end
                        
                            if  finCoordsArray{2,end} ~= 0
                                finCoordsArray{1,end+1} = finCoordsArray{1,end}+1;
                                finCoordsArray{2,end} = 0;
                                disp(['X',num2str(coordsCount),' = ',num2str(finCoordsArray{1,end}),...
                                    '; Y',num2str(coordsCount),' = ',num2str(finCoordsArray{2,end})])
                            end

                        finCoordsArray = cell2mat(finCoordsArray);
                        finXcoords = finCoordsArray(1,:) * (fin_rootChord(finSetNUM)/finCoordsArray(1,end));
                        finYcoords = finCoordsArray(2,:) * (fin_rootChord(finSetNUM)/finCoordsArray(1,end));
                end
                
                if (part_identity{2,part_count} == 3) || (part_identity{2,part_count} == 4)
                    freeformFins = true; % used later for warning regarding use of Barrowman equations
                    for     coordsCount = 2:length(finXcoords)
                        if  finXcoords(coordsCount) == finXcoords(coordsCount-1)
                            if  coordsCount == length(finXcoords)
                                finXcoords(coordsCount-1) = finXcoords(coordsCount-1)-(1*10^-6);
                            else
                                finXcoords(coordsCount) = finXcoords(coordsCount)+(1*10^-6);
                            end
                        end
                        if  finYcoords(coordsCount) == finYcoords(coordsCount-1)
                            if  coordsCount == length(finYcoords)
                                finYcoords(coordsCount-1) = finYcoords(coordsCount-1)-(1*10^-6);
                            else
                                finYcoords(coordsCount) = finYcoords(coordsCount)+(1*10^-6);
                            end
                        end
                    end     % fix for 'interp1' function's inability to use vectors with same points
                end
                
                        nonTube_endIndex = 0;
                for     findTubeCount = length(xAxisRange):-1:1
                    if  nonTube_endIndex == 0 && xAxisRange(findTubeCount) > 0
                        nonTube_endIndex = findTubeCount;
                    end
                    if  xAxisRange(findTubeCount) == rocket.Body_Tube
                        nonTube_startIndex = findTubeCount;
                        nonTube_indexRange = nonTube_endIndex - nonTube_startIndex;
                        break
                    end
                end

                xAxisRange(   new_xPosition - floor(fin_distFromBase(finSetNUM)/x_increment) ...
                            - nonTube_indexRange - floor(fin_rootChord(finSetNUM)/x_increment) ...
                            : new_xPosition - floor(fin_distFromBase(finSetNUM)/x_increment) ...
                            - nonTube_indexRange ) = rocket.Fins;
                            %  special rules for fins/canards - overrides section of assigned body tube
                            
                dlmwrite([rocketName,'.txt'],...
                    [part_identity{1,part_count},part_identity{2,part_count},finTubeAssign(finSetNUM),...
                    numOfFins(finSetNUM),fin_distFromBase(finSetNUM),finThickness(finSetNUM),fin_rootChord(finSetNUM),...
                    fin_span(finSetNUM),fin_sweepDist(finSetNUM),fin_tipChord(finSetNUM),fin_midChord(finSetNUM),...
                    nan,finXcoords,nan,finYcoords,nan,nonTube_indexRange],...
                    '-append','delimiter',',');
                
                
                disp(['This fin/canard set has been assigned to body tube No.',num2str(finTubeAssign(finSetNUM)),'.'])
            

            otherwise
                error('You entered a rocket part number outside the bounds of the key provided.')
        end


        if  new_xPosition >= length(xAxisRange)
            disp(['The rocket has been fully defined up to its total length.'...
                  ' You may now only add (more) fins/canards.'])
        end
    
        keep_adding = max(strcmp((input('Do you want to add another rocket part? [Y/N]: ','s')),yes));
    
        %if keep_adding == true; part_identity{2,length(part_identity)+1}=[]; end
    
    end

        if  new_xPosition < length(xAxisRange)
            error('The total length of the rocket parts you entered does not match the total rocket length you specified prior.')
        end

    if  part_identity{1,end} == 3 ...
        || ( part_identity{1,end} == 4 && part_identity{1,end-1} == 3 )
        disp('The final transition has been identified as a boat tail.')
    end

    rocket.Boat_Tail = 3;
    rocket.Canards = 4;

    
elseif rocketFromFile == true
    %%
        sameRocket = false;
    if  simRunNum > 1
        sameRocket =  max(strcmp((input('Use the same rocket design? [Y/N]: ','s')),yes));
    end
    if  sameRocket == false
        rocketName = [];
        while isempty(rocketName) == true
            try
                rocketName = input('Type the name of the rocket design data file to import: ','s');
                rocketData = dlmread([rocketName,'.txt']);
            catch ME
                if  strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                    warning('A file with that name doesn''t exist. Try again.')
                    rocketName = [];
                end
            end
        end
    end
    
    R_length = rocketData(1,1);
    skinRoughness = rocketData(1,2);
    rocketDataSize = size(rocketData);
    
	switch  skinRoughness
            case 1
                skinRoughness = 0.00001;
            case 2
                skinRoughness = (0.00002+0.00008)/2;
            case 3
                skinRoughness = 0.00016;
            case 4
                skinRoughness = 0.00025;
            case 5
                skinRoughness = (0.0004+0.0012)/2;
	end
    
	x_TotalSteps = 10000;
    x_increment = R_length/x_TotalSteps;

    xAxisRange = zeros(1,x_TotalSteps);
    new_xPosition = 1;
    
    for rocketCount = 2:rocketDataSize(1)
        switch rocketData(rocketCount,1)
            case 1
                part_identity{1,rocketCount-1}  = rocketData(rocketCount,1);
                part_identity{2,rocketCount-1}  = rocketData(rocketCount,2);
                nose_shapePrmtr                 = rocketData(rocketCount,3);
                nose_baseRadius                 = rocketData(rocketCount,4);
                nose_length                     = rocketData(rocketCount,5);
                nose_tipRadius                  = rocketData(rocketCount,6);
                
                xAxisRange( new_xPosition : floor(nose_length/x_increment) ) = rocket.Nose;
                new_xPosition = floor(nose_length/x_increment) + 1;
                
            case 2
                bodyTubeNUM = bodyTubeNUM + 1;
                part_identity{1,rocketCount-1}  = rocketData(rocketCount,1);
                bodyTubeDiameter(bodyTubeNUM)   = rocketData(rocketCount,2);
                bodyTubeLength(bodyTubeNUM)     = rocketData(rocketCount,3);
                
                xAxisRange( new_xPosition : new_xPosition + floor(bodyTubeLength(bodyTubeNUM)/x_increment) ) = rocket.Body_Tube;
                new_xPosition = new_xPosition + floor( bodyTubeLength(bodyTubeNUM)/x_increment ) + 1;
                
            case 3
                transitionNUM = transitionNUM + 1;
                part_identity{1,rocketCount-1}      = rocketData(rocketCount,1);
                trans_frontDiameter(transitionNUM)  = rocketData(rocketCount,2);
                trans_rearDiameter(transitionNUM)	= rocketData(rocketCount,3);
                trans_length(transitionNUM)         = rocketData(rocketCount,4);
                
                xAxisRange( new_xPosition : new_xPosition + floor(trans_length(transitionNUM)/x_increment) ) = rocket.Transition;
                new_xPosition = new_xPosition + floor( trans_length(transitionNUM)/x_increment ) + 1;
                
            case 4
                finSetNUM = finSetNUM + 1;
                part_identity{1,rocketCount-1}  = rocketData(rocketCount,1);
                part_identity{2,rocketCount-1}  = rocketData(rocketCount,2);
                finTubeAssign(finSetNUM)        = rocketData(rocketCount,3);
                numOfFins(finSetNUM)            = rocketData(rocketCount,4);
                fin_distFromBase(finSetNUM)     = rocketData(rocketCount,5);
                finThickness(finSetNUM)         = rocketData(rocketCount,6);
                fin_rootChord(finSetNUM)        = rocketData(rocketCount,7);
                fin_span(finSetNUM)             = rocketData(rocketCount,8);
                fin_sweepDist(finSetNUM)        = rocketData(rocketCount,9);
                fin_tipChord(finSetNUM)         = rocketData(rocketCount,10);
                fin_midChord(finSetNUM)         = rocketData(rocketCount,11);
                
                if ~(isnan(rocketData(rocketCount,13))) % if coordinate data for freeform fins exists
                    [finXcoords,finYcoords] = deal([]);
                    coordImportCounter = 13;
                    while ~(isnan(rocketData(rocketCount,coordImportCounter)))
                        finXcoords = [finXcoords,rocketData(rocketCount,coordImportCounter)];
                        coordImportCounter = coordImportCounter + 1;
                    end
                    coordImportCounter = coordImportCounter + 1;
                    while ~(isnan(rocketData(rocketCount,coordImportCounter)))
                        finYcoords = [finYcoords,rocketData(rocketCount,coordImportCounter)];
                        coordImportCounter = coordImportCounter + 1;
                    end
                    
                    for     coordsCount = 2:length(finXcoords)
                        if  finXcoords(coordsCount) == finXcoords(coordsCount-1)
                            finXcoords(coordsCount) = finXcoords(coordsCount)+(1*10^-6);
                        end
                        if  finYcoords(coordsCount) == finYcoords(coordsCount-1)
                            finYcoords(coordsCount) = finYcoords(coordsCount)+(1*10^-6);
                        end
                    end
                end
                
                nonTube_indexRange              = rocketData(rocketCount,end);
                
                xAxisRange(   new_xPosition - floor(fin_distFromBase(finSetNUM)/x_increment) ...
                            - nonTube_indexRange - floor(fin_rootChord(finSetNUM)/x_increment) ...
                            : new_xPosition - floor(fin_distFromBase(finSetNUM)/x_increment) ...
                            - nonTube_indexRange ) = rocket.Fins;
                        
            otherwise % {nan}
                chuteDataStart = rocketCount;
                break
        end
    end
    
    
end % {import rocket design}

    for xRangeCount = length(xAxisRange):-1:1
        if xAxisRange(xRangeCount) > 0
            xAxisRange = xAxisRange(1:xRangeCount);
            break
        end
    end
    
    
    for baseCount = length(part_identity):-1:1
        if      (cell2mat(part_identity(1,baseCount)) == 2)
                bodybaseDiameter = bodyTubeDiameter(max(bodyTubeNUM));
        elseif  (cell2mat(part_identity(1,baseCount)) == 3)
                bodybaseDiameter = trans_rearDiameter(max(transitionNUM));
        end
    end % gets diameter of rocket body base/rear end
    
                bodyMaxDiameter = max(bodyTubeDiameter);
     

disp('You have the option to enter the rocket''s wetted surface area manually, or the program will estimate it from your design.')
    surfaceArea_method =  max(strcmp((input('Enter a value for the surface area directly? [Y/N]: ','s')),yes));
if  surfaceArea_method == true
    
                surfaceArea = 0;
        [surfaceArea] = changeInput(simRunNum,'total wetted surface area',...
        'Enter a value for the rocket''s total wetted surface area (m^2):  ',surfaceArea);

elseif surfaceArea_method == false
    
    [surfaceArea,tubeCounter,frustumCounter,finsCounter] = deal(0);
	for surfaceCount = 1:length(part_identity)
        switch cell2mat(part_identity(1,surfaceCount))
            case 1 % nose
                slantLength = sqrt(nose_baseRadius^2 + nose_length^2);
                partSurfaceArea = pi*nose_baseRadius*slantLength;
            case 2 % tube
                tubeCounter = tubeCounter + 1;
                partSurfaceArea = pi*bodyTubeDiameter(tubeCounter)*bodyTubeLength(tubeCounter);
            case 3 % frustum
                frustumCounter = frustumCounter + 1;
                slantLength = sqrt(trans_length(frustumCounter)^2 + ...
                             (abs(trans_frontDiameter(frustumCounter)-trans_rearDiameter(frustumCounter)))^2);
                partSurfaceArea = pi*slantLength*(trans_frontDiameter(frustumCounter)/2 + trans_rearDiameter(frustumCounter)/2);
            case 4 % fins
                finsCounter = finsCounter + 1;
                switch cell2mat(part_identity(2,surfaceCount))
                    case 1
                        partSurfaceArea = fin_span(finsCounter)*(fin_rootChord(finsCounter)+fin_tipChord(finsCounter));
                        finTrailingEdge = sqrt(  fin_span(finsCounter)^2 + ( fin_rootChord(finsCounter)-fin_tipChord(finsCounter)-sqrt(fin_sweepDist(finsCounter)^2-fin_span(finsCounter)^2) )  );
                        partSurfaceArea = partSurfaceArea + finThickness(finsCounter)*( fin_sweepDist(finsCounter) + fin_tipChord(finsCounter) + finTrailingEdge );
                    case 2
                        partSurfaceArea = pi*fin_rootChord(finsCounter)*fin_span(finsCounter);
                        a = fin_rootChord(finsCounter)/2; b = fin_span(finsCounter);
                        partSurfaceArea = partSurfaceArea + finThickness(finsCounter)*0.5*pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));
                        
                    case {3,4}
                        approx_finSpan = max(finYcoords(finsCounter))*( fin_rootChord(finsCounter)/max(finXcoords(finsCounter)) );
                        partSurfaceArea = pi*fin_rootChord(finsCounter)*approx_finSpan;
                        for finPerimeterCount = 2:length(finXcoords)
                            finCoordPairLength{finPerimeterCount-1} = sqrt( (finXcoords(finPerimeterCount)-finXcoords(finPerimeterCount-1))^2 ...
                                                                          + (finYcoords(finPerimeterCount)-finYcoords(finPerimeterCount-1))^2 ) ...
                                                                    * ( fin_rootChord(finsCounter)/max(finXcoords(finsCounter)) );
                        end
                        partSurfaceArea = partSurfaceArea + finThickness(finsCounter)*sum(cell2mat(finCoordPairLength));
                end
                        partSurfaceArea = partSurfaceArea*numOfFins(finsCounter);
        end
        surfaceArea = surfaceArea + partSurfaceArea;
	end
        surfaceArea = surfaceArea + (pi*bodybaseDiameter^2)/4;
end     %%%DEV:NB%%% assumes fin thickness/sides have square shape/cross-section (no round or airfoil shape)
    
    
    if      bodybaseDiameter < bodyTubeDiameter(bodyTubeNUM)
            R_effLength = R_length; % rocket effective length
            
    elseif  bodybaseDiameter == bodyTubeDiameter(bodyTubeNUM)
            R_effLength = R_length - bodyTubeLength(bodyTubeNUM);
    end
    
    
    if      bodyMaxDiameter == (2*nose_baseRadius)
            R_aftLength = R_effLength - nose_length; % rocket length aft of max diameter point
        if  R_aftLength < 0
            R_aftLength = 0;
        end
            
    elseif  bodyMaxDiameter > (2*nose_baseRadius)
            [R_frontLength,tubeCounter,frustumCounter,finsCounter] = deal(0);
            for aftLengthCount = 1:length(part_identity)
                switch  cell2mat(part_identity(1,aftLengthCount))
                    case 1
                        R_frontLength = R_frontLength + nose_length;
                    case 2
                        tubeCounter = tubeCounter + 1;
                        if      bodyTubeDiameter(tubeCounter) < bodyMaxDiameter
                                R_frontLength = R_frontLength + bodyTubeLength(tubeCounter);
                        elseif  bodyTubeDiameter(tubeCounter) == bodyMaxDiameter
                                break
                        end
                    case 3
                        frustumCounter = frustumCounter + 1;
                        if      max(trans_frontDiameter(frustumCounter),trans_rearDiameter(frustumCounter)) < bodyMaxDiameter
                                R_frontLength = R_frontLength + trans_length(frustumCounter);
                        elseif  max(trans_frontDiameter(frustumCounter),trans_rearDiameter(frustumCounter)) < bodyMaxDiameter
                                break
                        end
                    case 4
                        % placeholder
                end
            end
            R_aftLength = R_effLength - R_frontLength;
    end
    
    
disp(['As of alpha ',num2str(alphaVer),', you can now add launch lugs to your rocket design.'])
[numOfLaunchLugs] = changeInput(simRunNum,'number of launch lugs',...
        'Enter a value for the number of launch lugs on your rocket: ',numOfLaunchLugs);
if  numOfLaunchLugs ~= 0
    for launchLugCount = 1:numOfLaunchLugs
        disp(['For launch lug no.',num2str(launchLugCount),','])
        protub_length{launchLugCount}       = input('-enter a value for its length: ')*10^-metricLength;
        protub_dist{launchLugCount}         = input('-enter a value for its distance from the rocket nose tip: ')*10^-metricLength;
        protub_crossArea{launchLugCount}    = input('-enter a value for its cross-sectional area: ')/10^(metricLength^2);
        protub_surfArea{launchLugCount}     = input('-enter a value for its total surface area: ')/10^(metricLength^2);
    end
else
    [protub_length,protub_dist,protub_crossArea,protub_surfArea] = deal([]);
end
    
    
%%%
end % {static stability analysis}


%% INPUTS [inp04]
%%
%
if  flight_analysis == true
%%%

if  inputsFromFile == false
    
%% Initial Conditions

disp('Initial Conditions:')

[int_range] = (changeInput(simRunNum,'initial down range distance', ...
    'Enter a value for the down range distance from your reference point at launch: ',int_range))*10^-metricLength;

[int_altitude] = (changeInput(simRunNum,'launch altitude', ...
    'Enter a value for the launch altitude above sea level: ',int_altitude))*10^-metricLength;

[pitchAngle] = changeInput(simRunNum,'launch inclination', ...
    'Enter a value for the inclination (measured from horizontal, in degrees) at launch: ',pitchAngle);

[int_velocity] = (changeInput(simRunNum,'initial velocity', ...
    'Enter a value for the velocity magnitude at launch: ',int_velocity))*10^-metricLength;

[towerHeight] = (changeInput(simRunNum,'launch tower height', ...
    'Enter a value for the height of your launch tower: ',towerHeight))*10^-metricLength;


[wind_velocity] = (changeInput(simRunNum,'average wind velocity', ...
    'Enter a value for the average wind velocity magnitude at the launch site: ',wind_velocity))*10^-metricLength;

if  wind_velocity ~= 0
    [gust_velocity] = (changeInput(simRunNum,'wind gust velocity', ...
    'Enter a value for the wind gust velocity (maximum wind velocity) magnitude at the launch site: ',gust_velocity))*10^-metricLength;
    gust_velocity = gust_velocity*cos(yawAngle*(pi/180));
 
    
    if  simRunNum == 1
        disp('The next question determines whether to launch the rocket into or against the flow of the wind.')
        yawAngle = [];
    end
    %while isempty(find([0,180,360] == yawAngle)) == true
    while isempty(yawAngle) == true
    [yawAngle] = changeInput(simRunNum,'wind angle', ...
    'Enter a value for the angle that the wind vector makes with your rocket at launch (0 or 180), in degrees: ',yawAngle);
    %%%DEV:NB%%% assume = 0 or 180 for now, for simplicity (2D plane, 2 or 3 DOF only)
        switch yawAngle
            case {0,180}
                % placeholder
            otherwise
                disp('You failed to enter an appropriate value for the wind direction. Try again.')
                yawAngle = [];
        end
    end
else
    [gust_velocity,yawAngle] = deal(0);
end
    wind_velocity = wind_velocity*cos(yawAngle*(pi/180));

 
%% Stage Information

disp('Stage Information:')

[payloadMass] = (changeInput(simRunNum,'payload mass', ...
    'What is your payload mass? This can be zero: ',payloadMass))*10^-metricMass;

dlmwrite([inputsFileName,'.txt'],...
                 [int_range,int_altitude,pitchAngle,int_velocity,towerHeight,...
                  yawAngle,wind_velocity,gust_velocity,payloadMass],...
                 '-append','delimiter',',');
             

for stage_count = 1:numOfStages

    [stage_propMass(stage_count)] = (changeInput(simRunNum,['stage ',num2str(stage_count),' propellant mass'], ...
        ['Enter the propellant mass of stage ',num2str(stage_count),': '],stage_propMass(stage_count)))*10^-metricMass;
    
    [stage_ISP(stage_count)] = changeInput(simRunNum,['stage ',num2str(stage_count), ...
        ' specific impulse'],['Enter the specific impulse (s) of stage ',num2str(stage_count),': '],stage_ISP(stage_count));
    
    [stage_thrust(stage_count)] = changeInput(simRunNum,['stage ',num2str(stage_count),' thrust'], ...
        ['Enter the thrust (kN) of stage ',num2str(stage_count),': '],stage_thrust(stage_count));
    
    [stage_diameter(stage_count)] = (changeInput(simRunNum,['stage ',num2str(stage_count),' diameter'], ...
        ['Enter the diameter of stage ',num2str(stage_count),': '],stage_diameter(stage_count)))*10^-metricLength;
    
    [stage_preCoast(stage_count)] = changeInput(simRunNum,['stage ',num2str(stage_count),' pre-ignition cost time'], ...
        ['Enter the pre-ignition coast time (s) of stage ',num2str(stage_count),': '],stage_preCoast(stage_count));
    
    [stage_postCoast(stage_count)] = changeInput(simRunNum,['stage ',num2str(stage_count),' post-burn coast time'], ...
        ['Enter the post-burn coast time (s) of stage ',num2str(stage_count),': '],stage_postCoast(stage_count));
    
    
	[stage_inertMass(stage_count)] = (changeInput(simRunNum,['stage ',num2str(stage_count),' inert mass'], ...
        ['Enter the inert mass of stage ',num2str(stage_count),': '],stage_inertMass(stage_count)))*10^-metricMass;
    
    
    stage_impulse(stage_count)      = stage_propMass(stage_count)*9.807*stage_ISP(stage_count)/1000;
    stage_burnTime(stage_count)     = stage_impulse(stage_count)/stage_thrust(stage_count);
    stage_area(stage_count)         = pi*(stage_diameter(stage_count)^2)/4;
    
    stage_mass(stage_count)         = stage_propMass(stage_count)+stage_inertMass(stage_count);
    stage_thrustWeight(stage_count) = (stage_thrust(stage_count)*1000)/(stage_mass(stage_count)*9.807);
    stage_massFlow(stage_count)     = stage_propMass(stage_count)/stage_burnTime(stage_count);
    
    
    dlmwrite([inputsFileName,'.txt'],...
                 [stage_propMass(stage_count),stage_ISP(stage_count),stage_thrust(stage_count),...
                 stage_diameter(stage_count),stage_preCoast(stage_count),stage_postCoast(stage_count),...
                 stage_inertMass(stage_count),...
                 stage_impulse(stage_count),stage_burnTime(stage_count),stage_area(stage_count),...
                 stage_mass(stage_count),stage_thrustWeight(stage_count),stage_massFlow(stage_count)],...
                 '-append','delimiter',',');
    
end


elseif  inputsFromFile == true
        
        inputsCellArray = num2cell(inputsData);
        
        [int_range,int_altitude,pitchAngle,int_velocity,towerHeight,...
            yawAngle,wind_velocity,gust_velocity,payloadMass]...
        = deal(inputsCellArray{2,1:9});
        
        stageInfoCounter = 1;
        
        for stage_count = 1:numOfStages
            
            [stage_propMass(stage_count),stage_ISP(stage_count),stage_thrust(stage_count),...
             stage_diameter(stage_count),stage_preCoast(stage_count),stage_postCoast(stage_count),...
             stage_inertMass(stage_count),...
             stage_impulse(stage_count),stage_burnTime(stage_count),stage_area(stage_count),...
             stage_mass(stage_count),stage_thrustWeight(stage_count),stage_massFlow(stage_count)] ...
             ...
             = deal(inputsCellArray{3,stageInfoCounter:stageInfoCounter+12});
         
            stageInfoCounter = stageInfoCounter + 13;
        end
end

grossLiftoffMass = payloadMass + sum(stage_mass);


if  wind_velocity ~= 0
    fprintf('%s\n','You have the option to either compute the effect of wind using the average wind velocity alone,',...
        'or to use randomized distributed data that incorporates the effect of gusts (wind speed variation).',...
        'The latter method will produce different results on each run due to randomization of wind speed.')
        windEffect_method =  max(strcmp((input('Use the latter method (randomize wind speeds)? [Y/N]: ','s')),yes));

        clear kCell
        xWind = linspace(0,gust_velocity);
        k0 = 1;
        for kCount = 1:100
            k1 = sum((xWind(2:end).^k0).*log(xWind(2:end)));
            k2 = sum((xWind(2:end).^k0));
            k3 = (1/(length(xWind(2:end))-1))*sum(log(xWind(2:end)));
            k  = 1/( (k1/k2) - k3 );
            k0 = k;
            kCell{kCount} = k;
            if (kCount > 1) && (kCell{kCount} == kCell{kCount-1})
                break
            end
        end
            shapeParameter = k;
            scaleParameter = ( (1/(length(xWind)-1))*sum(xWind(2:end).^k) )^(1/k);
        
    if  windEffect_method == true
            for windCount = 1:numOfSteps
                sim_windVelocity(windCount) = ( scaleParameter*(-log(rand(1)))^(1/shapeParameter) )...
                                            * cos(yawAngle*(pi/180));
            end
    else
                sim_windVelocity(1:end) = wind_velocity;
    end
else
                sim_windVelocity(1:end) = wind_velocity;
end
                windVelocity_preAlt = sim_windVelocity; % wind velocity irrespective of altitude


%% Recovery

if  (staticStability == true) && (rocketFromFile == false)
    
[numOfChutes] = changeInput(simRunNum,'number of parachutes', ...
    'Enter a value for the number of parachutes being used: ',numOfChutes);

dlmwrite([rocketName,'.txt'],[nan,numOfChutes],...
    '-append','delimiter',',');

if numOfChutes > 0
    
    disp('Recovery Options')
    
    for chuteCount = 1:numOfChutes
        
        [chuteDiameter(chuteCount)] = (changeInput(simRunNum,['parachute number ',num2str(chuteCount),' diameter'], ...
            ['Enter a value for the diameter of chute No.',num2str(chuteCount),': '],chuteDiameter(chuteCount)))*10^-metricLength;
                                        
        [chuteDragCoef(chuteCount)] = changeInput(simRunNum,['parachute number ',num2str(chuteCount),' drag coefficient'], ...
            ['Enter a value for the drag coefficient of chute No.',num2str(chuteCount),': '],chuteDragCoef(chuteCount));
        
        disp(parachute)
        
        [chuteDeployMethod(chuteCount)] = changeInput(simRunNum,['parachute number ',num2str(chuteCount),' deployment method'], ...
            ['Enter a value, based on the above key, for the deployment method of chute No.',num2str(chuteCount),': '],chuteDeployMethod(chuteCount),[1,2,3,4]);
                                         
        switch chuteDeployMethod(chuteCount)
            case 1
                chuteTimeDelay(chuteCount) = nan;
                chuteSetAltitude(chuteCount) = nan;
                
            case {2,3}
                [chuteTimeDelay(chuteCount)] = changeInput(simRunNum,['parachute number ',num2str(chuteCount),' deployment time delay'], ...
                    ['Enter a value for the deployment time (s) delay for chute No.',num2str(chuteCount),': '],chuteTimeDelay(chuteCount));
                 chuteTimeDelay(chuteCount)  = sum(chuteTimeDelay(chuteCount:-1:1));
                
                 chuteSetAltitude(chuteCount) = nan;
                
            case 4
                [chuteSetAltitude(chuteCount)] = (changeInput(simRunNum,['parachute number ',num2str(chuteCount),' deployment altitude'], ...
                    ['Enter a value for the deployment altitude for chute No.',num2str(chuteCount),': '],chuteSetAltitude(chuteCount)))*10^-metricLength;
                 chuteTimeDelay(chuteCount) = nan;
            
            otherwise
                error('You entered a value outside the range of the key.')
        end
        
        
        dlmwrite([rocketName,'.txt'],[chuteDiameter(chuteCount),chuteDragCoef(chuteCount),...
            chuteDeployMethod(chuteCount),chuteTimeDelay(chuteCount),chuteSetAltitude(chuteCount)],...
            '-append','delimiter',',');
        
    end
end


elseif (staticStability == true) && (rocketFromFile == true) && (chuteDataStart > 0)
    
        numOfChutes = rocketData(chuteDataStart,2);
    
    for recoveryCount = chuteDataStart+1:rocketDataSize(1)
        chuteDiameter(recoveryCount-chuteDataStart)        = rocketData(recoveryCount,1);
        chuteDragCoef(recoveryCount-chuteDataStart)        = rocketData(recoveryCount,2);
        chuteDeployMethod(recoveryCount-chuteDataStart)    = rocketData(recoveryCount,3);
        chuteTimeDelay(recoveryCount-chuteDataStart)       = rocketData(recoveryCount,4);
        chuteSetAltitude(recoveryCount-chuteDataStart)     = rocketData(recoveryCount,5);
    end
    

end % {import rocket design, incl. recovery}


if  numOfChutes > 0
    fprintf('%s\n','There are two methods that can be employed to calculate the overall drag force following parachute deployment:',...
        '1) Calculate using the more significant drag component, ignoring the lesser ones (e.g. deploy drogue, ignore rocket drag)',...
        '2) Calculate air flow speed reduction going from part to part (e.g. rocket to drogue), finally obtaining an overall drag force')
    [dragCalcMethod] = changeInput(simRunNum,'descent drag calculation method', ...
        'Enter a value based on the above key for the method you want to use: ',dragCalcMethod,[1,2]);
end

%%%
end % {flight analysis}


%% SIMULATION SETUP [set05]
%%
if  flight_analysis == true
%%%

%% Time Events - Stage


for stage_count = 1:numOfStages
    
    if stage_count == 1
        ignition(stage_count) = 0;
        burnout(stage_count) = ignition(stage_count) + stage_burnTime(stage_count);
        separation(stage_count) = stage_burnTime(stage_count) + stage_postCoast(stage_count);
    else
        ignition(stage_count) = separation(stage_count-1) + stage_preCoast(stage_count);
        burnout(stage_count) = ignition(stage_count) + stage_burnTime(stage_count);
        if stage_count == length(numOfStages)
            separation(stage_count) = expFlightDur*10;
        else
            separation(stage_count) = burnout(stage_count) + stage_postCoast(stage_count);
        end
    end
    
end


%% Time Events - Phase

stageCounter = 1;
phaseCounter = 1;

for phase_count = 1:3*numOfStages
    
    phase(phase_count) = phase_count;
    
    switch phaseCounter
        
        case 1
            phase_endTime(phase_count)   = ignition(floor(stageCounter));
            phase_thrust(phase_count)    = 0;
            phase_massFlow(phase_count)  = 0;
            
            if      floor(stageCounter) == 1
                    phase_mass(phase_count) = grossLiftoffMass;
            elseif  floor(stageCounter) == 2
                    phase_mass(phase_count) = grossLiftoffMass - stage_mass(1);
            else
                    phase_mass(phase_count) = phase_mass(phase_count-3) - stage_mass(floor(stageCounter)-1);
            end
            
        case 2
            phase_endTime(phase_count)   = burnout(floor(stageCounter));
            phase_thrust(phase_count)    = stage_thrust(floor(stageCounter))*1000;
            phase_massFlow(phase_count)  = stage_massFlow(floor(stageCounter));
            phase_mass(phase_count)      = -1;
            
        case 3
            phase_endTime(phase_count)   = separation(floor(stageCounter));
            phase_thrust(phase_count)    = 0;
            phase_massFlow(phase_count)  = 0;
            phase_mass(phase_count)      = phase_mass(phase_count-2) - stage_propMass(floor(stageCounter));
    end   
            phase_area(phase_count)      = stage_area(floor(stageCounter));
    
            
    if      phaseCounter < 3
            phaseCounter = phaseCounter+1;
    elseif  phaseCounter == 3
            phaseCounter = phaseCounter-2;
    end
    
    stageCounter = stageCounter + (1/3);
end


%% Warning Messages

        CDcalcWarnID  =  'KINGSROCKsim:CDcalc';
        CDcalcWarnMSG = ['Your rocket has a nose length to body length ratio outside the range ',...
                           'of 0.2 - 0.6, which means parasite drag calculations ',...
                           'shouldn''t be utilized. You can set the zero-lift drag coefficient ',...
                           'to a particular value instead.'];
                       
        critWarnID  =  'KINGSROCKsim:criticalRatio';
        critWarnMSG = ['Your rocket has a nose length to effective body length ratio ',...
                       'of more than 0.6, which means the transonic wave drag calculations ',...
                       '(and subsequently supersonic wave drag) may be inaccurate.'];
                       
        rocketDataWarnID  =  'KINGSROCKsim:rocketData';
        rocketDataWarnMSG = ['Since you have chosen not to use rocket design data for stability analysis, ',...
                             'the function for calculating parasite drag cannot be utilized. ',...
                             'You will instead set the zero-lift drag coefficient to a particular value.'];

                         
%% Data

dynViscosity=[
175	1.182;
200	1.329;
225	1.467;
250	1.599;
275	1.725;
300	1.846;
325	1.962;
350	2.075;
375	2.181;
400	2.286;
450	2.485;
500	2.67;
550	2.849;
600	3.017;
650	3.178;
700	3.332;
750	3.482;
800	3.624;
850	3.763;
900	3.897;
950	4.026;
1000	4.153;
1050	4.276;
1100	4.396;
1150	4.511;
1200	4.626;
1250	4.736;
1300	4.846;
1350	4.952;
1400	5.057;
1500	5.264;
1600	5.457;
1700	5.646;
1800	5.829;
1900	6.008
];


%% Failure modes

disp('As of alpha 1.04, you can also account for rocket failure modes.')
failureChoice = max(strcmp(input('Do you want to simulate the possibility of failure mode(s) occuring? [Y/N]: ','s'),yes));
if failureChoice == true
%
fprintf('%s\n','There are two methods for achieving failure mode simulation:',...
    '-Select a particular failure mode to simulate, or...',...
    '-Include a random factor in the simulation based on the likelihood of any given failure occuring.')
failureOption = max(strcmp(input('Use the latter method (randomize failure mode occurances)? [Y/N]: ','s'),yes));

if	exist('failDescriptor','var') == false
	failDescriptor = true;
elseif exist('failDescriptor','var') == true
    failDescriptor = false;
end
    
if  failDescriptor == true
    fprintf('%s\n','The failure modes that can be simulated are as such:')
    tic; while toc < 1; end; fprintf('%s\n',...
        '1. Unstable - If the rocket points nose-down during the motor burn.',...
        '              Will be apparent during normal simulation, depending on your rocket design and launch conditions.')
    tic; while toc < 1; end; fprintf('%s\n',...
        '2. Lawn dart - Recovery failure, but nosecone remains attached.',...
        '               The result is a ballistic flight (no parachutes).')
    tic; while toc < 1; end; fprintf('%s\n',...
        '3. Separation - Any component separation (without recovery device attached) event.',...
        '                NOT accounted for due to unpredictability and complexity.')
    tic; while toc < 1; end; fprintf('%s\n',...
        '4. Motor CATO - Motor failure at ignition or during firing.',...
        '                Accounted for by negating thrust at a random time prior to expected burnout.')
    tic; while toc < 1; end; fprintf('%s\n',...
        '5. Core sample - Ballistic flight like lawn dart, but nosecone ejects whilst remaining attached by cord.',...
        '                 Accounted for by removing parachutes, and increasing descent drag area.')
    tic; while toc < 1; end; fprintf('%s\n',...
        '6. Motor unrestrained - Motor separates from system prior to burnout.',...
        '                        Same method as ''Motor CATO'' used here.')
    tic; while toc < 1; end; fprintf('%s\n',...
        '7. Shred - Rocket complete separation/destruction during ascent.',...
        '           Simulation terminates at apogee, since a random time prior to apogee cannot be determined.')
    tic; while toc < 1; end; fprintf('%s\n',...
        '8. No chute - Non-ballistic trajectory with no parachutes deployed.',...
        '              Due to difficulty in predicting non-ballistic descent, the simulation terminates at apogee.')
    tic; while toc < 1; end;
end

switch failureOption
    case 0 % select failure mode
        [failureMode] = changeInput(simRunNum,'type of failure mode to simulate',...
                        'Enter a value from the above key for the type of failure mode to simulate [1-8]: ',failureMode,[1,2,3,4,5,6,7,8]);
    case 1 % randomize failure mode
        if  failDescriptor == true
            fprintf('%s\n','Based on historical data, failure modes typically occur 8.5% (523 of 6169) of the time.',...
            ['Therefore, a failure will occur in approximately every 1 in ',num2str(100/8.5),' runs.'],...
            'The chances of each failure mode occuring is listed below:',...
            ['1. Unstable - 19% of failures, so ',num2str(0.19*(523/6169)),'% of flights'],...
            ['2. Lawn dart - 23% of failures, so ',num2str(0.23*(523/6169)),'% of flights'],...
            ['3. Separation - 28% of failures, so ',num2str(0.28*(523/6169)),'% of flights'],...
            ['4. Motor CATO - 6% of failures, so ',num2str(0.06*(523/6169)),'% of flights'],...
            ['5. Core sample - 5% of failures, so ',num2str(0.05*(523/6169)),'% of flights'],...
            ['6. Motor unrestrained - 1% of failures, so ',num2str(0.01*(523/6169)),'% of flights'],...
            ['7. Shred - 2% of failures, so ',num2str(0.02*(523/6169)),'% of flights'],...
            ['8. No chute - 16% of failures, so ',num2str(0.16*(523/6169)),'% of flights'])
        end
        tic; while toc < 1.5; end;
            failureChance = rand(1);
        if  failureChance > (523/6169)
            fprintf('%s\n','Based on random variable calculation, this particular',...
                    'simulation run will NOT suffer a failure mode.')
                    failureMode = 0;
        else
                    failTypeChance = rand(1);
            if      failTypeChance <= 0.19
                    failureMode = 1;
                    fprintf('%s\n','Based on random variable calculation, this particular',...
                            'simulation run would suffer the ''Unstable'' failure mode.',...
                            'However, dynamic stability calculations will determine whether your rocket is actually unstable or not.')
            elseif  failTypeChance <= (0.19+0.23)
                    failureMode = 2;
                    fprintf('%s\n','Based on random variable calculation, this particular',...
                            'simulation run will suffer the ''Lawn dart'' failure mode.')
            elseif  failTypeChance <= (0.19+0.23+0.28)
                    failureMode = 3;
                    fprintf('%s\n','Based on random variable calculation, this particular',...
                            'simulation run will suffer the ''Separation'' failure mode.')
            elseif  failTypeChance <= (0.19+0.23+0.28+0.06)
                    failureMode = 4;
                    fprintf('%s\n','Based on random variable calculation, this particular',...
                            'simulation run will suffer the ''Motor CATO'' failure mode.')
            elseif  failTypeChance <= (0.19+0.23+0.28+0.06+0.05)
                    failureMode = 5;
                    fprintf('%s\n','Based on random variable calculation, this particular',...
                            'simulation run will suffer the ''Core sample'' failure mode.')
            elseif  failTypeChance <= (0.19+0.23+0.28+0.06+0.05+0.01)
                    failureMode = 6;
                    fprintf('%s\n','Based on random variable calculation, this particular',...
                            'simulation run will suffer the ''Motor unrestrained'' failure mode.')
            elseif  failTypeChance <= (0.19+0.23+0.28+0.06+0.05+0.01+0.02)
                    failureMode = 7;
                    fprintf('%s\n','Based on random variable calculation, this particular',...
                            'simulation run will suffer the ''Shred'' failure mode.')
            else
                    failureMode = 8;
                    fprintf('%s\n','Based on random variable calculation, this particular',...
                            'simulation run will suffer the ''No chute'' failure mode.')
            end
        end
end
%
else
    failureMode = 0;
end

%%%
end % {flight analysis}


%% IMPORT ENGINE (motor/thrust curve files) [eng06]
%%
%
if  flight_analysis == true
%%%

%
fprintf('%s\n','You now have the option to replace the average thrust data with more',...
    'accurate thrust wrt time data from an engine file.')

engineCounter = 0;
stageCounter = 1;

keep_adding = true;


while engineCounter < numOfStages
    
    keep_adding = max(strcmp((input(['Add engine data to replace constant thrust for stage No.',num2str(stageCounter),'? [Y/N]: '],'s')),yes));
    
    if  keep_adding == true
        %
        engineCounter = engineCounter + 1;
        %
            sameEngine = false;
        if  simRunNum > 1
            sameEngine =  max(strcmp((input(['Use the same engine file for stage ',num2str(stageCounter),'? [Y/N]: '],'s')),yes));
        end
        if  sameEngine == false
        engineFile = [];
        while isempty(engineFile) == true
            try
                engineFile = input(['Type the file name for your stage ',num2str(stageCounter),' engine file (with the extension): '],'s');
                disp('Please wait...')
                %
                engine = importdata(engineFile);
            catch ME
                if  strcmp(ME.identifier,'MATLAB:importdata:FileNotFound')
                    warning('File cannot be located. Try again.')
                    engineFile = [];
                end
            end
        end
        end
        %
        thrustData{engineCounter} = zeros(size(engine.textdata,1)-3,2); 
        % removes first 4 rows (not thrust/time data), and adds a row for t=0
        %
        for columnCount = 1:2
            for rowCount = 5:(size(engine.textdata,1))
                thrustData{engineCounter}(rowCount-3,columnCount) = str2num(cell2mat(engine.textdata(rowCount,columnCount)));
            end
        end
        %
        burnout(stageCounter) = ignition(stageCounter) + engine.data; % updates burnout accordingly
        if  stage_burnTime(stageCounter) ~= burnout(stageCounter)
            stage_burnTime(stageCounter)  = burnout(stageCounter);
        end
        %
        % thrustData contains all the thrust data (in the second column) for a range
        % of times specified by the first column
        %
        stageCounter = stageCounter + 1;
        
        
    elseif keep_adding == false
            break
    end

end

numOfEngines = engineCounter;

%%%
end % {flight analysis}


%% COP calculation [cop07]
%%
%
% Alt. method to Barrowman: calculate width of rocket (y_i_total) wrt length (x_i),
% then multiply to get areas of narrow strips, then use trapezium rule to get total area,
% then: sum((2*y_i*x_increment)*x_i)/sum(2*y_i*x_increment) = x_i Average = COP
%

if staticStability == true
%%%

fprintf('%s\n','Now we will calculate the COP (Centre of Pressure) of the rocket.')%,...
       %'Alternatively, if you know the COP location you can enter it directly.')

%COP_method = max(strcmp((input('Enter COP location directly? [Y/N]: ','s')),yes));
COP_method = 0;

[bodyTubeNUM,transitionNUM,finSetNUM]=deal(0);
finShapeCounter = 1;
y_i_projOverlay = [];
fin_area = 0;

[y_i,y_i_projected] = deal(zeros(1,length(xAxisRange)));


%%
for COPcount = 1:length(xAxisRange)
    
	switch xAxisRange(COPcount)
    
        %%
        case rocket.Nose
            
            switch part_identity{2,1} % assumes first rocket part is the nose, and only one nose in design
                
                case noseShape.Conical
                y_i(COPcount) = (COPcount*x_increment * nose_baseRadius) / nose_length;
                noseFactor = 2/3;
            
                case noseShape.Bi_conic
                    if      COPcount*x_increment <= nose_length{1}
                            y_i(COPcount) = (COPcount*x_increment * nose_baseRadius{1}) / nose_length{1};
                    elseif  COPcount*x_increment <= nose_length{2}
                            y_i(COPcount) = nose_baseRadius{1} + ( (COPcount*x_increment-nose_length{1})*...
                            (nose_baseRadius{2}-nose_baseRadius{1}) )/nose_length{2};
                    end
                frustrumAngle = atan(nose_length{2}/(nose_baseRadius{2}-nose_baseRadius{1}));
                delta_noseLength = tan(frustrumAngle)*nose_baseRadius{2} - (nose_length{1}+nose_length{2});
                noseFactor = 2/3;
            
                case noseShape.Tangent_ogive
                rho_i = ( nose_baseRadius^2 + nose_length^2 ) / (2*nose_baseRadius);
                y_i(COPcount) = sqrt( rho_i^2 - (COPcount*x_increment - nose_length)^2 ) + (nose_baseRadius - rho_i);
                noseFactor = 0.466;
                
                case noseShape.Secant_ogive
                alpha_i = atan(nose_baseRadius/nose_length) - ...
                    acos( sqrt(nose_length^2 + nose_baseRadius^2) / (2*nose_shapePrmtr) );
                y_i(COPcount) = sqrt( nose_shapePrmtr^2 - (COPcount*x_increment - ...
                    nose_shapePrmtr*cos(alpha_i))^2 ) + (nose_shapePrmtr*sin(alpha_i));
                noseFactor = 0.466;
            
                case noseShape.Ellipsoid
                y_i(COPcount) = nose_baseRadius*sqrt( 1 - ((COPcount*x_increment)^2)/(nose_length^2) );
                noseFactor = 1/3;
            
                case noseShape.Power_series
                y_i(COPcount) = nose_baseRadius*((COPcount*x_increment)/nose_length)^nose_shapePrmtr;
                noseFactor = 1/2;
            
                case noseShape.Parabolic_series
                y_i(COPcount) = nose_baseRadius*( (2*((COPcount*x_increment)/nose_length) - ...
                    nose_shapePrmtr*((COPcount*x_increment)/nose_length)^2) / (2-nose_shapePrmtr) );
                noseFactor = 1/2;
            
                case noseShape.Haack_series
                theta_i = acos(( 1 - (2*(COPcount*x_increment) / nose_length) ));
                y_i(COPcount) = (  nose_baseRadius*sqrt( theta_i - (sin(2*theta_i)/2) + ...
                    nose_shapePrmtr*(sin(theta_i))^3 )  )/sqrt(pi);
                switch nose_shapePrmtr
                    case 1/3 % LV_haack
                        noseFactor = 0.437;
                    case 0 % LD-Haack (Von Karman)
                        noseFactor = 0.5;
                end
            end

            
            y_i_projected(COPcount) = y_i(COPcount);

            
            % BARROWMAN EQUATIONS
            nose_normalForce = 2;
            nose_COP = noseFactor*nose_length;
            if  part_identity{2,1} == noseShape.Bi_conic
                nose_COP = nose_COP - (1-noseFactor)*delta_noseLength;
                % from: Xn = 2/3(L + dL) - dL = 2/3(L) + 2/3(dL) - dL
                %          = 2/3(L) - 1/3(dL)
                % approximation for other nose shapes e.g. bi-conic
            end
            
            
        %%    
        case rocket.Body_Tube
            
            if  ( xAxisRange(COPcount) ~= xAxisRange(COPcount-1) ) && ...
                ( xAxisRange(COPcount-1) ~= rocket.Fins )
                bodyTubeNUM = bodyTubeNUM + 1;
            end
            
            y_i(COPcount) = bodyTubeDiameter(bodyTubeNUM)/2;
            y_i_projected(COPcount) = y_i(COPcount);
            
            
        %%    
        case {rocket.Transition,rocket.Boat_Tail}
            
            if xAxisRange(COPcount) ~= xAxisRange(COPcount-1)
                transitionNUM = transitionNUM + 1;
                transPoint = COPcount;
                transPointLength = transPoint *x_increment;
            end

            x1  = x_increment*(transPoint);
            x2  = x1 + trans_length(transitionNUM);
            x_i = x_increment*COPcount;
            % y_i = y1 + ( (y2-y1)/(x2-x1) )*(x_i-x1) {linear equation}
            y_i(COPcount) = trans_frontDiameter(transitionNUM)/2 + ...
                          ((trans_rearDiameter(transitionNUM)/2 - ...
                            trans_frontDiameter(transitionNUM)/2)/(x2-x1))*(x_i-x1);
            y_i_projected(COPcount) = y_i(COPcount);
            
            
            % BARROWMAN EQUATIONS {conical transitions ONLY}
            trans_normalForce(transitionNUM) = 2*(  (trans_rearDiameter(transitionNUM) / 2*nose_baseRadius)^2 - ...
                                    (trans_frontDiameter(transitionNUM) / 2*nose_baseRadius)^2 );
                    
            trans_COP(transitionNUM)  =   (trans_length(transitionNUM)/3)*(  1 + ...
                                    ( 1- (trans_frontDiameter(transitionNUM)/trans_rearDiameter(transitionNUM)) ) / ...
                                    ( 1- (trans_frontDiameter(transitionNUM)/trans_rearDiameter(transitionNUM))^2 )  ) ...
                                    + transPointLength;
            
        %%    
        case {rocket.Fins,rocket.Canards}
            
            if xAxisRange(COPcount) ~= xAxisRange(COPcount-1)
                finSetNUM = finSetNUM + 1;
                transPoint = COPcount;
                finStartPoint = x_increment*transPoint;
                        finAngleCount = 1;
                while   isnan(finAngleCount) == false
                        finAngle(finSetNUM)= (360/numOfFins(finSetNUM))*finAngleCount;
                    if  finAngle(finSetNUM) > 180
                        finAngle(finSetNUM)= (360/numOfFins(finSetNUM))*(finAngleCount-1) - 90;
                        finAngleCount = nan;
                    else
                        finAngleCount = finAngleCount+1;
                    end
                end
                finAngle(finSetNUM) = (finAngle(finSetNUM))*(pi/180);
                
                for finShapeCount = (finShapeCounter+1):length(part_identity)
                    finShapeCounter = finShapeCount;
                    if ~(isempty(part_identity{2,finShapeCounter}))
                        break % determines the shape of the fins for a particular set
                    end
                end
                    
                switch  part_identity{2,finShapeCounter}
                    
                    case finShape.Trapezoidal
                    
                    fin_leadAngle       = asin(fin_span(finSetNUM)/fin_sweepDist(finSetNUM));
                    trans_frontLength   = fin_sweepDist(finSetNUM)*cos(fin_leadAngle);
                    trans_rearLength    = fin_rootChord(finSetNUM)-(trans_frontLength+fin_tipChord(finSetNUM));
                    fin_trailAngle      = atan(fin_span(finSetNUM)/trans_rearLength);
            
                    finLeadVertDist     = fin_sweepDist(finSetNUM)*asin((pi/2)-fin_leadAngle);

                    finLead     = finStartPoint + trans_frontLength;
                    finMid      = finLead + fin_tipChord(finSetNUM);
                    finTrail    = finMid + trans_rearLength;
                    
                    
                    case finShape.Elliptical
                        
                    ellipseCounter = 1;
 
                end
            end
            
            
            switch part_identity{2,finShapeCounter}
                
                case finShape.Trapezoidal
                    
                    % find inside angle: sin(angle) = opp/hyp = span/sweep
                    % use angle: tan(angle) = opp/adj = y_i / x_increment*count
                    % therefore: y_i = x_increment*COPcount*tan(angle)
        
                    % (1): leading edge {treated like a transition}
                    if COPcount*x_increment <= finLead
                        y_i(COPcount) = (trans_frontLength-(finLead - x_increment*COPcount))*tan(fin_leadAngle) ...
                                        + bodyTubeDiameter(finTubeAssign(finSetNUM))/2;

                    % (2): tip chord {treated like a body tube in calculations}
                    elseif COPcount*x_increment <= finMid
                        y_i(COPcount) = fin_span(finSetNUM) + bodyTubeDiameter(finTubeAssign(finSetNUM))/2;

                    % (3): trailing edge {treated like a transition}
                    elseif COPcount*x_increment <= finTrail
                        y_i(COPcount) = ( trans_rearLength - (x_increment*COPcount - finMid) )*tan(fin_trailAngle) ...
                                        + bodyTubeDiameter(finTubeAssign(finSetNUM))/2;   
                    end
            
                
                case finShape.Elliptical
                    
                    y_ellipse{ellipseCounter} = fin_span(finSetNUM)*...
                        sin(pi*((ellipseCounter)/(fin_rootChord(finSetNUM)/x_increment)));
                
                    y_i(COPcount) = y_ellipse{ellipseCounter} + bodyTubeDiameter(finTubeAssign(finSetNUM))/2;
                
                    ellipseCounter = ellipseCounter+1;
                    
                    
                case {finShape.Freeform_image_file,finShape.Freeform_coordinates}
                    
                        finXpoint = ( (COPcount-transPoint+1)*x_increment );
                    if  finXpoint > finXcoords(end)
                        finXpoint = finXcoords(end); % fix for fin trailing edge not connecting with body
                    end

                    y_i(COPcount) = interp1(finXcoords,finYcoords,finXpoint,'linear')...
                                  + bodyTubeDiameter(finTubeAssign(finSetNUM))/2;
            end
            
            
                y_i_projected(COPcount) = y_i(COPcount)*sin(finAngle(finSetNUM));
            if  y_i_projected(COPcount) < bodyTubeDiameter(finTubeAssign(finSetNUM))/2
                y_i_projOverlay = [y_i_projOverlay; COPcount*x_increment , y_i_projected(COPcount)];
                y_i_projected(COPcount) = bodyTubeDiameter(finTubeAssign(finSetNUM))/2;
                % accounts for when y_i_projected is less than body tube radius
            end
            
            
                finArea(finSetNUM) = finArea(finSetNUM) + x_increment*(y_i(COPcount)+y_i_projected(COPcount));
                fin_area = fin_area + x_increment*( y_i(COPcount) + y_i_projected(COPcount) - bodyTubeDiameter(finTubeAssign(finSetNUM)) );
            if  (y_i(COPcount)+y_i_projected(COPcount)) > (y_i(COPcount-1)+y_i_projected(COPcount-1))
                finSpanTotal(finSetNUM) = y_i(COPcount) + y_i_projected(COPcount);
            end
            
            
            % BARROWMAN EQUATIONS
            % may become inaccurate for a large number of fins (>4)
            % and cannot apply to freeform fins
            switch part_identity{2,finShapeCounter}
                
                case finShape.Trapezoidal
                    
                fin_normalForce(finSetNUM) = ( 1 + (bodyTubeDiameter(bodyTubeNUM)/2)/(fin_span(finSetNUM)+bodyTubeDiameter(bodyTubeNUM)) ) * ...
                            ( 4*numOfFins(finSetNUM)*(fin_span(finSetNUM)/(2*nose_baseRadius))^2 ) / ...
                            ( 1+sqrt(1+((2*fin_midChord(finSetNUM))/(fin_rootChord(finSetNUM)+fin_tipChord(finSetNUM)))^2) );
                    
                fin_COP(finSetNUM) = ( finLeadVertDist*(fin_rootChord(finSetNUM)+2*fin_tipChord(finSetNUM)) ) / ...
                        ( 3*(fin_rootChord(finSetNUM)+fin_tipChord(finSetNUM)) ) + ...
                        (1/6)* ( fin_rootChord(finSetNUM) + fin_tipChord(finSetNUM) - ...
                        ((fin_rootChord(finSetNUM)*fin_tipChord(finSetNUM))/(fin_rootChord(finSetNUM)+fin_tipChord(finSetNUM))) ) ...
                        + (finLead - trans_frontLength); % like transPointLength for transitions
                    
                    
                case finShape.Elliptical
                    
                fin_normalForce(finSetNUM) = ( 1 + (bodyTubeDiameter(bodyTubeNUM)/2)/(fin_span(finSetNUM)+bodyTubeDiameter(bodyTubeNUM)) ) * ...
                            ( 4*numOfFins(finSetNUM)*(fin_span(finSetNUM)/(2*nose_baseRadius))^2 ) / ...
                            ( 1+sqrt(1+1.623*(fin_span(finSetNUM)/fin_rootChord(finSetNUM))^2) );
                
                fin_COP(finSetNUM) = 0.288*fin_rootChord(finSetNUM) + finStartPoint;
                
                
                case {finShape.Freeform_image_file,finShape.Freeform_coordinates}
                    
                    fin_span(finSetNUM) = max(finYcoords);
                    
                fin_normalForce(finSetNUM) = ( 1 + (bodyTubeDiameter(bodyTubeNUM)/2)/(fin_span(finSetNUM)+bodyTubeDiameter(bodyTubeNUM)) ) * ...
                            ( 4*numOfFins(finSetNUM)*(fin_span(finSetNUM)/(2*nose_baseRadius))^2 ) / ...
                            ( 1+sqrt(1+1.623*(fin_span(finSetNUM)/fin_rootChord(finSetNUM))^2) );
                        
                fin_COP(finSetNUM) = 0.288*fin_rootChord(finSetNUM) + finStartPoint;
                    
            end
            
        %%
        otherwise
        error('You entered a rocket part number outside the bounds of the key provided.')
	end
            %%

            y_i_total       = y_i(COPcount) + y_i_projected(COPcount);            
            y_iArray        = [y_iArray,y_i_total];

end
            x_iArray        = x_increment*(1:length(y_iArray));
            y_i = [0,y_i];
            y_i_projected = [0,y_i_projected];

%

centreOfPlanformArea = sum( (y_iArray*x_increment).*x_iArray ) ...
                     / sum( y_iArray*x_increment );
            
COP_options.Barrowman_Equations = 1;
COP_options.Centre_of_Area = 2;

fprintf('%s\n','One of two methods can be used to calculate the static COP: by using the Barrowman equations, which is', ...
    'restricted in its rocket design options (no canards!), or by calculating the centre of area, which is a', ...
    'less accurate approximation of the static COP.')
disp(COP_options)

COP_selected = input('Enter a number, corresponding to the above key, for the method of calculating the static COP: ');
     

if      COP_method == 1  
        centreOfPressure = (input('Enter a value for the distance of the COP from the nose tip: '))*10^-metricLength;
    
elseif  COP_method == 0
    
    if      COP_selected == COP_options.Barrowman_Equations
        
        if  freeformFins == true
            warning('%s\n','You''ve chosen to use the Barrowman method but with a freeform fin design. ',...
                           'To estimate the fin COP, its shape will be approximated to an ellipse.');
        end
        if  finSetNUM > 1
            warning('%s\n','You''ve chosen to use the Barrowman method but with multiple sets of fins. ',...
                           'Be warned: this will likely produce erroneous COP results.');
        end
            COP_numerator       = nose_normalForce*nose_COP + sum(trans_normalForce.*trans_COP) + sum(fin_normalForce.*fin_COP);
            COP_denominator     = nose_normalForce + sum(trans_normalForce) + sum(fin_normalForce);
            
            centreOfPressure    = COP_numerator/COP_denominator;
        
    elseif  COP_selected == COP_options.Centre_of_Area
            centreOfPressure    = centreOfPlanformArea;
                                % equivalent to centre of projected/planform area
            warning('%s\n','The movement of the COP at large angles of attack cannot be predicted when using the ',...
                           'planform area method. The Centre of Pressure will be assumed constant for all angles of attack.');
    end
end

%%%
end % {static stability analysis}


%% COM Calculations [com08]
%%
%
% currently accounts for longitudinal axis only
%

if staticStability == true
%%%

fprintf('%s\n','Now we will calculate the COM (Centre of Mass) of the rocket.')%,...
       %'Alternatively, if you know the COM location you can enter it directly.');

%COM_method = max(strcmp((input('Enter COM location directly? [Y/N]: ','s')),yes));
COM_method = 0;

if      COM_method == 1
    
        centreOfMass = (input('Enter a value for the distance of the COM from the nose tip AT LAUNCH: '))*10^-metricLength;
        centreOfMass_noProp = (input('Enter a value for the distance of the COM from the nose tip WITHOUT PROPELLANT: '))*10^-metricLength;
    
        
elseif  COM_method == 0
    
        rocketCompsFromFile =  max(strcmp((input('Do you want to import files containing data for rocket design components? [Y/N]: ','s')),yes));
	if  rocketCompsFromFile == false

        rocketCompsName = [rocketName,' (components)'];

        disp([  'You will enter the masses of individual components of the rocket,', ...
                ' and their distances from the nose tip.'   ]);
       
        keep_adding = true;
        
        namesFileID = fopen([rocketName,' (component names).txt'],'wt');

        while keep_adding == true
    
            componentCount = componentCount+1;
            
            component_name{componentCount}      =  input('Type a name for the component: ','s');
            component_mass(componentCount)      = (input('Enter the mass of the component: '))*10^-metricMass;
            component_position(componentCount)  = (input('Enter the distance (mm) of the component from the nose tip: '))*10^-metricMass;
            component_isProp(componentCount) 	= max(strcmp((input('Is this component a propellant {solid (fuel) or liquid (oxidiser)}? [Y/N]: ','s')),yes));
            
            dlmwrite([rocketCompsName,'.txt'],[component_mass(componentCount),...
                component_position(componentCount),component_isProp(componentCount)],...
                '-append','delimiter',','); % 'roffset',0
            
            fprintf(namesFileID,'%s\n',component_name{componentCount});
            
    
            keep_adding = max(strcmp((input('Do you want to add more components? [Y/N]: ','s')),yes));
            
            if  keep_adding == true
                component_mass      = [component_mass,0];
                component_position  = [component_position,0];
                component_isProp 	= [component_isProp,0];
            end
    
        end
        
        fclose(namesFileID);
        
    elseif rocketCompsFromFile == true
        
        if      exist([rocketName,' (components).txt']) == 2
                rocketCompsName = [rocketName,' (components)'];
        elseif  exist([rocketName,' (components).txt']) == 0
                error(['There is no components data file associated with the rocket named ',rocketName,'.'])
        end
        rocketCompsData = dlmread([rocketCompsName,'.txt']);
        rocketCompsDataSize = size(rocketCompsData);
        
        namesFileData = textread([rocketName,' (component names).txt'],'%s','delimiter','\n');
        
        component_name{1}       = namesFileData(1);
        component_mass(1)       = rocketCompsData(1,1);
        component_position(1)   = rocketCompsData(1,2);
        component_isProp(1)     = rocketCompsData(1,3);

        for rocketCompsCount = 2:rocketCompsDataSize(1)
            if  length(component_name)              < rocketCompsCount
                component_mass                      = [component_mass,0];
                component_position                  = [component_position,0];
                component_isProp                    = [component_isProp,0];
            end
            component_name(rocketCompsCount)      = namesFileData(rocketCompsCount);
            component_mass(rocketCompsCount)      = rocketCompsData(rocketCompsCount,1);
            component_position(rocketCompsCount)  = rocketCompsData(rocketCompsCount,2);
            component_isProp(rocketCompsCount)    = rocketCompsData(rocketCompsCount,3);
        end
	end
    
        componentCount = 0;

        centreOfMass = sum(component_mass.*component_position)/sum(component_mass);
        
        if (simRunNum == 1) || ( (simRunNum > 1) && (sameRocket == false) )
        disp('Your rocket components will be listed below, with their moment of inertia about the COM at launch:')
        for component_count = 1:length(component_name)
            component_MOI{component_count} = component_mass(component_count)*(component_position(component_count)-centreOfMass)^2;
            fprintf('%s\n',...
                ['Name: ',              char(component_name{component_count})],...
                ['Mass = ',             num2str(component_mass(component_count)),       ' kilograms'],...
                ['Position = ',         num2str(component_position(component_count)),   ' metres'],...
                ['Moment of inertia = ',num2str(component_MOI{component_count}),        ' kg.m^2'],...
                [])
            tic; while toc < 1.5; end
        end
        end
        
        disp(['The COM at launch is located ',num2str(centreOfMass*10^3),' mm from the nose tip.']);

	component_isProp = ~component_isProp;
	centreOfMass_noProp = sum(component_isProp.*component_mass.*component_position)/sum(component_isProp.*component_mass);

	disp(['The COM at propellant depletion is located ',num2str(centreOfMass_noProp*10^3),' mm from the nose tip.']);
    
    
    disp('The rocket''s mass moment of inertia (MOI) can be calculated from the provided component data.')%, or entered manually.')
            %MOI_method = max(strcmp((input('Enter a value for the rocket''s MOI directly? [Y/N]: ','s')),yes));
            MOI_method = 0;
    if      MOI_method == 1
            MOInertia           = (input('Enter a value for the rocket''s MOI when AT LAUNCH: '))*(10^-metricMass)*(10^-metricLength^2);
            MOInertia_noProp    = (input('Enter a value for the rocket''s MOI when it''s WITHOUT PROPELLANT: '))*(10^-metricMass)*(10^-metricLength^2);
    
    elseif  MOI_method == 0
            MOInertia           = sum(cell2mat(component_MOI));
            MOInertia_noProp    = sum(cell2mat(component_MOI).*component_isProp);
            fprintf('%s\n',['The rocket''s MOI at launch is ',num2str(MOInertia),' kg.m^2'],['The rocket''s MOI at propellant depletion is ',num2str(MOInertia_noProp),' kg.m^2'])
    end

end

%%%
end % {static stability analysis}


%% STABILITY [sta09]
%%
%

if staticStability == true
%%%

%% Rocket Outline Drawing
%       
xEnd = [ x_iArray(end), x_iArray(end), x_iArray(end) ];
yEnd = [ y_i(end), 0, -y_i_projected(end) ];

for rocketDraw_count = 2:length(xAxisRange)
    if  xAxisRange(rocketDraw_count) ~= xAxisRange(rocketDraw_count-1)
        vertLineCounter = vertLineCounter + 1;
        rocketDraw_vertLineX{vertLineCounter} = [ x_iArray(rocketDraw_count),x_iArray(rocketDraw_count),x_iArray(rocketDraw_count) ];
        rocketDraw_vertLineY{vertLineCounter} = [ y_i(rocketDraw_count),0,-y_i_projected(rocketDraw_count) ];
    end
end

figure
plot(x_iArray,y_i,'-k',x_iArray,-y_i_projected,'-k',xEnd,yEnd,'-k');

axis([0,R_length,(-R_length/2),(R_length/2)])
axis square
xlabel('(metres)'); ylabel('(metres)')
set(gcf,'Units','normalized','Position',[0.01,0.15,0.975,0.75]);


for rocketDraw_count = 1:vertLineCounter
    hold on
    plot(rocketDraw_vertLineX{rocketDraw_count},rocketDraw_vertLineY{rocketDraw_count},'-k')
end

    finEndIndex = 0;
for finHorzLineCount = 1:finSetNUM
    for rocketDraw_count = finEndIndex+1:length(xAxisRange)
        if xAxisRange(rocketDraw_count) == 4;
            finStartIndex = rocketDraw_count;
            break
        end
    end
    for rocketDraw_count = finStartIndex:length(xAxisRange)
        if xAxisRange(rocketDraw_count) ~= 4;
            finEndIndex = rocketDraw_count-1;
            break
        end
    end
    hold on
    plot([x_iArray(finStartIndex),x_iArray(finEndIndex)], ...
        [bodyTubeDiameter(finTubeAssign(finHorzLineCount))/2, ...
        bodyTubeDiameter(finTubeAssign(finHorzLineCount))/2],'-k')
    hold on
    if finAngle(finHorzLineCount) ~= pi/2
        plot([x_iArray(finStartIndex),x_iArray(finEndIndex)], ...
        [-finAngle(finHorzLineCount)*(bodyTubeDiameter(finTubeAssign(finHorzLineCount))/2), ...
        -finAngle(finHorzLineCount)*(bodyTubeDiameter(finTubeAssign(finHorzLineCount))/2)],'-k')
    else
        plot([x_iArray(finStartIndex),x_iArray(finEndIndex)], ...
        [-(bodyTubeDiameter(finTubeAssign(finHorzLineCount))/2), ...
        -(bodyTubeDiameter(finTubeAssign(finHorzLineCount))/2)],'-k')
    end
end
  
hold on
plot(y_i_projOverlay(:,1),-y_i_projOverlay(:,2),'-k')

%

hold on
COPpoint1 = plot(centreOfPressure,0,'Marker','o','Color','r');
hold on
COPpoint2 = plot(centreOfPressure,0,'Marker','.','Color','r');
text(centreOfPressure*1.01,0,'COP')

hold on
COMpoint1 = plot(centreOfMass,0,'Marker','o','Color','b','DisplayName','COM');
hold on
COMpoint2 = plot(centreOfMass,0,'Marker','+','Color','b','DisplayName','COM');
hold on
COMpoint3 = plot(centreOfMass,0,'Marker','x','Color','b','DisplayName','COM');
COMpointText = text(centreOfMass*1.01,0,'COM');

saveas(gcf,['Rocket_Design_(',num2str(simRunNum),').jpg']);
tic; while toc < 2; end
commandwindow
pauseFix = input('The current graph represents the rocket''s stability at launch. Hit enter to view it at propellant depletion.');
%drawnow;pause;drawnow % fix for hanging on pause issue with matlab r2016a
figure(gcf)

set(COMpoint1,'Visible','off')
set(COMpoint2,'Visible','off')
set(COMpoint3,'Visible','off')
set(COMpointText,'Visible','off')

hold on
plot(centreOfMass_noProp,0,'Marker','o','Color','b','DisplayName','COM');
hold on
plot(centreOfMass_noProp,0,'Marker','+','Color','b','DisplayName','COM');
hold on
plot(centreOfMass_noProp,0,'Marker','x','Color','b','DisplayName','COM');
text(centreOfMass_noProp*1.01,0,'COM')


%% Summary

tic; while toc < 2; end
commandwindow
bodyAvgDiameter     = sum(bodyTubeDiameter)/bodyTubeNUM;
staticMargin        = centreOfPressure - centreOfMass;
stabCaliber         = staticMargin / bodyAvgDiameter;

fprintf('%s\n',      'At launch, the static stability is as such: ', ...
                    ['Static Margin = ',num2str(staticMargin*10^3),' mm.'], ...
                    ['Stability = ',num2str(stabCaliber),' calibers.'])

staticMargin_noProp     = centreOfPressure - centreOfMass_noProp;
stabCaliber_noProp      = staticMargin_noProp / bodyAvgDiameter;
                
fprintf('%s\n',      'At complete propellant depletion, the static stability becomes: ', ...
                    ['Static Margin = ',num2str(staticMargin_noProp*10^3),' mm.'], ...
                    ['Stability = ',num2str(stabCaliber_noProp),' calibers.'])

close all
%%%
end % {static stability analysis}


%% Control Systems [con10]
%%

moveableFinSet = 1; 
[controlSystem,finDeflection,finRotationRate,...
    maxDeflection,gimbalingRate] = deal(0);

sim_thrustAngle = zeros(1,(expFlightDur/timeStep));

control.Moveable_Fins = 1; control.Gimbaled_Thrust = 2;
control.Vernier_Rocket = 3; control.Thrust_Vane = 4;


if staticStability == true
    
fprintf('%s\n',['As of alpha ',num2str(alphaVer),', this program can simulate certain rocket control mechanisms.'],...
    'These exist to manipulate the rocket''s flight path angle to remain vertical, thus optimising apogee height.',...
    'The options are as such:','1. Moveable fins/canards','2. Gimbaled thrust',...
    '3. Vernier rockets - NOT IMPLEMENTED','4. Thrust vanes - NOT IMPLEMENTED',...
    'Since thrust vanes are functionally similar to gimbaled thrust, you may select that option to simulate the effects.',...
    'Sensitivity of and response times for sensors that would be used for control are not accounted for.')

    control_system =  max(strcmp((input('Does your rocket employ a control system from the list of implemented options? [Y/N]: ','s')),yes));
if  control_system == true
    [controlSystem] = changeInput(simRunNum,'type of control system',...
                        'Enter a value from the above key for the control system to use [1,2]: ',controlSystem,[1,2]);
	switch controlSystem
        case 1
            %disp('The first set of fins/canards being moveable, and a rate of fin rotation of 10 deg./second, are assumed.')
            moveableFinSet = 0;
            while isempty(find((1:finSetNUM) == moveableFinSet)) == true
            [moveableFinSet] = changeInput(simRunNum,'set number for moveable fins/canards',...
                        ['Enter a value (from numbers ',num2str(1:finSetNUM),') for the fin/canard set that is moveable: '],moveableFinSet,(1:finSetNUM));
            end
            [maxDeflection] = changeInput(simRunNum,'maximum fin/canard deflection angle',...
                        'Enter a value for the maximum fin/canard deflection angle (+/- degrees): ',maxDeflection);
            [finRotationRate] = changeInput(simRunNum,'angular rate of fin/canard rotation',...
                        'Enter a number value for the angular rate of fin/canard rotation (degrees per second): ',finRotationRate);
                    finRotationRate = finRotationRate*(pi/180);
        case 2
            %disp('A maximum thrust deflection of +/- 10 degrees, and a rate of gimbal rotation of 10 deg./second, are assumed.')
            [maxDeflection] = changeInput(simRunNum,'maximum thrust deflection angle',...
                        'Enter a value for the maximum thrust deflection angle (+/- degrees): ',maxDeflection);
            [gimbalingRate] = changeInput(simRunNum,'angular rate of gimbaling',...
                        'Enter a value for the angular rate of gimbaling (degrees per second): ',gimbalingRate);
                    gimbalingRate = gimbalingRate*(pi/180);
        otherwise
            error('You entered a number outside the range of the key provided!')
	end
                    maxDeflection = maxDeflection*(pi/180);
end

end


%% FLIGHT SIM [fli11]
%% 
%
if  flight_analysis == true
%%%
%
disp('Please wait while simulation runs...')
tic
%
%% Initial Conditions

dt = timeStep;
initial_mass = grossLiftoffMass;

sim_flightPhase = [1,sim_flightPhase];


simCounter = 0;
firingCounter = 1;

[       sim_mass1(1),   sim_range1(1),  sim_altitude1(1),   sim_velocity1(1),   sim_angle1DEG(1)] = ...
deal(   initial_mass,   int_range,      int_altitude,       int_velocity,       pitchAngle);
%[       sim_mass1(1),   sim_range1(1),  sim_altitude1(1),   sim_velocity1(1),      sim_angle1DEG(1),   sim_angle1DEG_perp(1)] = ...
%deal(   initial_mass,   int_range,      int_altitude,       int_velocity,          pitchAngle,         yawAngle);

sim_componentMass   = component_mass;
component_isProp    = ~component_isProp;
sim_componentMOI    = cell2mat(component_MOI);


motorFailTime = burnout(1)*rand(1);


%% Begin simulation...

descentPhase = false;

while simCounter < numOfSteps
    
    if  (simCounter > 1) && ...
        (sim_time(simCounter) > burnout(end)) && ...
        (sim_altitude1(simCounter) == 0)
        break
    end
    
    
    if  ((failureMode == 7) || (failureMode == 8)) && (simCounter > 1) ...
            && (sim_altitude2(simCounter) - sim_altitude2(simCounter-1) < 0)
        disp(['Failure mode no.',num2str(failureMode),' experienced at T = ',num2str(sim_time(simCounter)),'s.'])
        break
    end

    
    simCounter = simCounter+1;
    
    sim_time(simCounter) = simCounter*timeStep;
    
    
    %% Phase Info
    
    sim_phaseEnd(simCounter) = phase_endTime(phase(sim_flightPhase(simCounter)));
    
    if  sim_time(simCounter) > sim_phaseEnd(simCounter)
        sim_flightPhase(simCounter+1) = sim_flightPhase(simCounter) + 1;
        % special rule for sim_flightPhase, since it needs to start at 1
    else
        sim_flightPhase(simCounter+1) = sim_flightPhase(simCounter);
    end
    
    sim_phaseMass(simCounter) = phase_mass(phase(sim_flightPhase(simCounter+1)));
    
    %% Pre-Iterations

    if      (sim_time(simCounter) >= ignition(firingCounter)) && ...
            (sim_time(simCounter) <= burnout(firingCounter))
        
        if  numOfEngines > 0
            sim_thrust(simCounter) = interp1( thrustData{firingCounter}(:,1),...
                                              thrustData{firingCounter}(:,2),...
                                              sim_time(simCounter),'linear' );

            if  sim_time(simCounter) >= separation(firingCounter)
                firingCounter = firingCounter + 1;
            end
            % Next stage engine ignition occurs after previous stage separation
        
        else
            sim_thrust(simCounter) = phase_thrust(phase(sim_flightPhase(simCounter+1)));
        end
    end
    
        if (failureMode == 4 || failureMode == 6) ...
                && sim_time(simCounter) >= motorFailTime

            if  sim_time(simCounter-1) < motorFailTime
                disp(['Failure mode no.',num2str(failureMode),' experienced at T = ',...
                    num2str(sim_time(simCounter)),'s.'])
            end
            
            sim_thrust(simCounter) = 0;
        end
    
    
    sim_massFlow(simCounter) = phase_massFlow(phase(sim_flightPhase(simCounter+1)));
    
    sim_area(simCounter) = phase_area(phase(sim_flightPhase(simCounter+1)));
    
    if  simCounter > 1
        if  sim_phaseMass(simCounter) < 0
            sim_mass1(simCounter) = sim_mass2(simCounter-1);
        else
            sim_mass1(simCounter) = sim_phaseMass(simCounter);
        end
    end
    
    sim_mass2(simCounter) = sim_mass1(simCounter) - sim_massFlow(simCounter)*dt;

    if  simCounter > 1
        sim_range1(simCounter) = sim_range2(simCounter-1);
        %sim_range1_perp(simCounter) = sim_range2_perp(simCounter-1);
    end
    
    if  simCounter > 1
        if  sim_altitude2(simCounter-1) < 0
            sim_altitude1(simCounter) = 0;
        else
            sim_altitude1(simCounter) = sim_altitude2(simCounter-1);
        end
    end
    
        if  sim_altitude1(simCounter) <= 150 
            % edge of boundary layer for calculating ground winds (R. Newlands)
            sim_windVelocity(simCounter) = wind_velocity*(sim_altitude1(simCounter)/10)^(1/7);
        else
            sim_windVelocity(simCounter) = wind_velocity*(150/10)^(1/7);
        end
        % assumes launch over open, flat terrain (land, not sea)
        % valid up to edge of planetary boundary layer, whose height varies
        
            %{
            sim_COM(simCounter) = interp1([initial_mass,initial_mass-sum(stage_propMass)],...
                                          [centreOfMass,centreOfMass_noProp],   sim_mass2(simCounter));
            sim_MOI(simCounter) = interp1([initial_mass,initial_mass-sum(stage_propMass)],...
                                          [MOInertia,MOInertia_noProp],   sim_mass2(simCounter));
            %}
        
            for massCount = 1:length(sim_componentMass)
                if  component_isProp(massCount) == true
                    sim_componentMass(massCount) = sim_componentMass(massCount) - sim_massFlow(simCounter)*dt * ...
                                        ( component_mass(massCount)/sum(component_mass.*component_isProp) );
                end
            end
            sim_COM(simCounter) = sum(sim_componentMass.*component_position)/sum(sim_componentMass);
            
            for massCount = 1:length(sim_componentMass)
                sim_componentMOI(massCount) = sim_componentMass(massCount)*(component_position(massCount)-sim_COM(simCounter))^2;
            end
            sim_MOI(simCounter) = sum(sim_componentMOI);
            
            
    
    %% Atmospheric Properties

    sim_atmosDensity(simCounter)    = atmosDensity(sim_altitude1(simCounter));
    sim_atmosTemp(simCounter)       = atmosTemp(sim_altitude1(simCounter));
    sim_soundSpeed(simCounter)      = 20.04*sqrt(sim_atmosTemp(simCounter));
    
    %% Pre-Iterations cont.
    
    sim_gravity(simCounter) = 398600441800000/(6375450+sim_altitude1(simCounter))^2;
    
    if  simCounter > 1
        if  sim_altitude2(simCounter-1) < 0
            sim_velocity1(simCounter) = 0;
        else
            sim_velocity1(simCounter) = sim_velocity2(simCounter-1);
        end
    end

    if  simCounter > 1
        if  (sim_altitude1(simCounter) < (towerHeight*sin(pitchAngle*(pi/180)))) && (sim_time(simCounter) < burnout(end))
            sim_angle1DEG(simCounter) = pitchAngle;
            %sim_angle1DEG_perp(simCounter) = yawAngle;
        else
            sim_angle1DEG(simCounter) = sim_angle2DEG(simCounter-1);
            %sim_angle1DEG_perp(simCounter) = sim_angle2DEG_perp(simCounter-1);
        end
    end
    
    sim_angle1RAD(simCounter) = sim_angle1DEG(simCounter)*(pi()/180);
    %sim_angle1RAD_perp(simCounter) = sim_angle1DEG_perp(simCounter)*(pi()/180);


    %% Recovery

    if (numOfChutes > 0) && (simCounter > 3)
        
        if  (  chuteDeployMethod(chuteCounter(simCounter-1)) == parachute.Apogee ...
               && ...
             ( sim_altitude2(simCounter-1) - sim_altitude2(simCounter-2) < 0 )  ) || ...
             ...
             ...
          (  ( chuteDeployMethod(chuteCounter(simCounter-1)) == parachute.Time_delay_post_launch || ...
               chuteDeployMethod(chuteCounter(simCounter-1)) == parachute.Time_delay_post_previous_deployment ) ...
               && ...
            ( (sim_time(simCounter-1) > chuteTimeDelay(chuteCounter(simCounter-1))) && ...
              (sim_time(simCounter-3) < chuteTimeDelay(chuteCounter(simCounter-1))) )  ) || ...
            ...
            ...
          (    chuteDeployMethod(chuteCounter(simCounter-1)) == parachute.Set_altitude_during_descent ...
               && ...
             ( sim_altitude2(simCounter-1) < chuteSetAltitude(chuteCounter(simCounter-1)) ) && ...
             ( sim_altitude2(simCounter-3) > chuteSetAltitude(chuteCounter(simCounter-1)) )  );

         
            if  failureMode == 2 || failureMode == 5
                disp(['Failure mode no.',num2str(failureMode),' experienced at T = ',num2str(sim_time(simCounter)),'s.'])
                deployChute = false;
            else
                deployChute = true;
                descentPhase = true;
            end
        else
                deployChute = false;
        end

        
            if  deployChute == true
                chuteDrag{chuteCounter(simCounter-1)} =	((pi*chuteDiameter(chuteCounter(simCounter-1))^2)/4) * chuteDragCoef(chuteCounter(simCounter-1));
                % new parachute drag (without dynamic pressure)
                        
                deployTime(chuteCounter(simCounter-1)) = sim_time(simCounter);
                
                if  chuteCounter(simCounter-1) < numOfChutes
                    chuteCounter((simCounter):end) = chuteCounter(simCounter)+1;
                    % signals next parachute for remainder of flight
                end
            end
    end

    
    %% Iterations 
    %  to calculate acceleration and velocity at end of time step:
    
        %% FIRST BLOCK
    
    sim_machNum(simCounter) = abs(sim_velocity1(simCounter) / sim_soundSpeed(simCounter));
    
    
    if (simCounter > 1) && (descentPhase == false) && (sim_altitude1(simCounter) > (towerHeight*sin(pitchAngle*(pi/180))))
        switch  controlSystem 
            case control.Moveable_Fins
                if      (sim_angle1DEG(simCounter) > 90) ||...
                        ((sim_angle1DEG(simCounter) == 90) && (finDeflection > 0))
                
                        finDeflection = finDeflection - finRotationRate*dt;
                        
                elseif  (sim_angle1DEG(simCounter) < 90) ||...
                        ((sim_angle1DEG(simCounter) == 90) && (finDeflection < 0))
                    
                        finDeflection = finDeflection + finRotationRate*dt;
                end
                
                    if  abs(finDeflection) > maxDeflection
                        finDeflection = finDeflection*abs(maxDeflection/finDeflection);
                    end
                    
                    
            case control.Gimbaled_Thrust
            if  (sim_time(simCounter) < burnout(end))
                
                if      (sim_angle1DEG(simCounter) < 90) ||...
                        ((sim_angle1DEG(simCounter) == 90) && (sim_thrustAngle(simCounter-1) > 0))
                
                        sim_thrustAngle(simCounter) = sim_thrustAngle(simCounter-1) - gimbalingRate*dt;
                        
                elseif  (sim_angle1DEG(simCounter) > 90) ||...
                        ((sim_angle1DEG(simCounter) == 90) && (sim_thrustAngle(simCounter-1) < 0))
                    
                        sim_thrustAngle(simCounter) = sim_thrustAngle(simCounter-1) + gimbalingRate*dt;
                end
                
                    if  abs(sim_thrustAngle(simCounter)) > maxDeflection
                        sim_thrustAngle(simCounter) = sim_thrustAngle(simCounter)*abs(maxDeflection/sim_thrustAngle(simCounter));
                    end
            end
        end
    end
    
    
            %% Parasite (Zero-Lift) Drag Coefficient
            
if  staticStability == true

    if  ~(((nose_length/R_length) >= 0.2) && ((nose_length/R_length) <= 0.6))
    
            [~,warnID] = lastwarn;
            if  (strcmp(warnID,CDcalcWarnID) == 0) && (sim_dragCoefO(1) == 0)
                warning(CDcalcWarnID,CDcalcWarnMSG)
                setDragCoef = max(strcmp((input('Set the parasite drag coefficient to a constant value (e.g. 0.75)? [Y/N]: ','s')),yes));
                if  setDragCoef == true
                    sim_dragCoefO(1) = input('Enter a value to set the parasite drag to: ');
                end
            end
                    sim_dragCoefO(simCounter+1) = sim_dragCoefO(simCounter);
    end

        if  (sim_machNum(simCounter) > 0) && (setDragCoef == false)
            [sim_dragCoefO(simCounter),critRatioError] = calcParasiteDrag(  sim_altitude1(simCounter),...
                                                            sim_soundSpeed(simCounter),sim_machNum(simCounter),...
                                                            R_length,skinRoughness,bodyMaxDiameter,...
                                                            surfaceArea,fin_rootChord,fin_tipChord,...
                                                            finThickness,fin_span,numOfFins,...
                                                            bodybaseDiameter,nose_length,...
                                                            R_effLength,R_aftLength,finSetNUM,...
                                            numOfLaunchLugs,protub_length,protub_dist,protub_crossArea,protub_surfArea);
                                        
                                        if  isreal(sim_dragCoefO(simCounter)) == false % if parasite drag function fails
                                            if  exist('defaultDragCoefO','var') == false
                                                warning('The parasite drag calculation function has failed to correctly estimate the zero-lift drag coefficient.')
                                                defaultDragCoefO = input('Enter a value for the default drag coefficient to use when the function fails: ');
                                            end
                                            sim_dragCoefO(simCounter) = defaultDragCoefO;
                                        end
                                                        
            if      critRatioError == true
                
                    [~,warnID] = lastwarn;
                    if  strcmp(warnID,critWarnID) == 0
                        warning(critWarnID,critWarnMSG)
                    end
            end
        end

    
elseif  staticStability == false
            
        [~,warnID] = lastwarn;
        if  strcmp(warnID,rocketDataWarnID) == 0
            warning(rocketDataWarnID,rocketDataWarnMSG)
            sim_dragCoefO(1) = input('Enter a value to set the parasite drag to: ');
        end
        
            sim_dragCoefO(simCounter+1) = sim_dragCoefO(simCounter);
end

            %% Lift Coefficient & Force
            
if  (numOfChutes == 0) || ( (numOfChutes > 0) && (descentPhase == false) )
    
    if (sim_altitude1(simCounter) > (towerHeight*sin(pitchAngle*(pi/180)))) && (sim_time(simCounter) < burnout(end))
        resultantFlow = sqrt( sim_velocity1(simCounter)^2 + (sim_windVelocity(simCounter)*sin(sim_angle1RAD(simCounter)))^2 );
        sim_angleOfAttack(simCounter) = acos( min(sim_velocity1(simCounter),resultantFlow) / max(sim_velocity1(simCounter),resultantFlow) );
        
            moveableFinAOA = finDeflection + sim_angleOfAttack(simCounter);
        
        K = 1; % constant between 1.1 & 1.5, but 1.0 tends to be suitable
        COPvsAOA = sim_angleOfAttack(simCounter)*(  ( 4 * K * (sum(y_iArray*x_increment)-fin_area) ) / (pi * (2*nose_baseRadius)^2)  );
        if  COP_selected == COP_options.Barrowman_Equations
            sim_COP(simCounter) = ( COP_numerator   + COPvsAOA*centreOfPlanformArea )...
                                / ( COP_denominator + COPvsAOA );
        else
            sim_COP(simCounter) = centreOfPressure;
        end
        
            moveableFinCOP = fin_COP(moveableFinSet);
        
        
        if  abs(sim_angleOfAttack(simCounter)*(180/pi)) > 10
            sim_liftCoef(simCounter) = liftCoef(sim_angleOfAttack(simCounter));
        else
            sim_liftCoef(simCounter) = (2*pi)*sim_angleOfAttack(simCounter);
        end
        
                if  abs(moveableFinAOA*(180/pi)) > 10
                    moveableFinLiftCoef = liftCoef(moveableFinAOA);
                else
                    moveableFinLiftCoef = (2*pi)*moveableFinAOA;
                end
        
        sim_windLift(simCounter) = 0.5 * sim_atmosDensity(simCounter) * resultantFlow^2 * sum(finArea) * sim_liftCoef(simCounter);
        if controlSystem == control.Moveable_Fins
            moveableFinLift = fin_normalForce(moveableFinSet) - ...
                0.5 * sim_atmosDensity(simCounter) * resultantFlow^2 * finArea(moveableFinSet) * moveableFinLiftCoef;
        else
            moveableFinLift = 0;
        end
        %{
        if numOfChutes > 0
            for liftCount = 1:numOfChutes+1
                if      liftCount == 1
                        angularAccel = ( sim_windLift(simCounter)*(sim_COP(simCounter)-sim_COM(simCounter)) ) / MOInertia;
                elseif  (deployTime(liftCount-1) ~= 0) && (sim_time(simCounter) >= deployTime(liftCount-1))
                        chuteNormalForce = ( chuteDrag{liftCount-1} )* 0.5*sim_atmosDensity(simCounter)*...
                                            ( sim_windVelocity(simCounter)*sin(sim_angle1RAD(simCounter)) )^2;
                        angularAccel = sum([(chuteNormalForce*(sim_COM(simCounter)-nose_length)),...
                                            -sim_windLift(simCounter)*(sim_COP(simCounter)-sim_COM(simCounter))]) / MOInertia;
                end
            end
        else
        %}
                angularAccel = ( (sim_windLift(simCounter) * (sim_COP(simCounter)-sim_COM(simCounter)))...
                               + (sim_thrust(simCounter)*sin(sim_thrustAngle(simCounter))) * (R_length-sim_COM(simCounter))...
                               + (moveableFinLift * (moveableFinCOP-sim_COM(simCounter))) )... 
                                / MOInertia;  % if staticMargin is positive, should rotate anticlockwise
        %end
        
        angleDiff = angularVelo*dt + 0.5*angularAccel*dt^2;
        angularVelo = angularVelo + angularAccel*dt;
    else
        sim_COP(simCounter) = centreOfPressure;
        angleDiff = 0;
    end
    
    sim_stability(simCounter) = ( sim_COP(simCounter) - sim_COM(simCounter) ) / bodyAvgDiameter;
    
    
            %% Lift-Induced Drag Coefficient
            
	for findFins = 1:length(part_identity)
        if  part_identity{1,findFins} == rocket.Fins
            if      part_identity{2,findFins} == finShape.Trapezoidal
                    sweepAngle = fin_leadAngle;
                    % acceptable for Raymer method; Shevell would require
                    % sweep to be taken at quarter-chord
            elseif  part_identity{2,findFins} == finShape.Elliptical
                    sweepAngle = 0;
            else
                    sweepAngle = 0; %%%DEV:NB%%%
            end
            break
        end
	end
    
    aspectRatio = (finSpanTotal(1)^2)/finArea(1);
    ozFactor = oswaldFactor_mod(aspectRatio,sweepAngle,'Average',...
                sim_dragCoefO(simCounter),(bodyAvgDiameter/finSpanTotal(1)),1);
            
    sim_dragCoefI(simCounter) = (sim_liftCoef(simCounter)^2) / ...
                                (pi*aspectRatio*ozFactor);
                            
else
    sim_COP(simCounter) = centreOfPressure;
    angleDiff = 0;
end
    
            %% Remaining calculations...
            
    if numOfChutes > 0
        saveVelocity = sim_velocity1(simCounter);
        
        for dragCount = 1:numOfChutes+1
                sim_dynPressure1(simCounter) = 0.5*sim_atmosDensity(simCounter)*( sim_velocity1(simCounter) + sim_windVelocity(simCounter)*cos(sim_angle1RAD(simCounter)) )^2;
            
                if      dragCount == 1
                        savedynPressure = sim_dynPressure1(simCounter);
                        sim_dragForce(simCounter) = ( sim_area(simCounter)*(sim_dragCoefO(simCounter)+sim_dragCoefI(simCounter)) )* sim_dynPressure1(simCounter);
                        sim_accel1(simCounter) = (sim_thrust(simCounter)*cos(sim_thrustAngle(simCounter))-sim_dragForce(simCounter)-(sim_mass1(simCounter)*sim_gravity(simCounter)*sin(sim_angle1RAD(simCounter))))/sim_mass1(simCounter);
                        sim_velocity1(simCounter) = sim_velocity1(simCounter) + sim_accel1(simCounter)*dt;
                elseif  (deployTime(dragCount-1) ~= 0) && (sim_time(simCounter) >= deployTime(dragCount-1))
                    switch dragCalcMethod
                        case 1
                            sim_dragForce(simCounter) = chuteDrag{dragCount-1} * sim_dynPressure1(simCounter);
                        case 2
                            sim_dragForce(simCounter) = (chuteDrag{dragCount-1} * sim_dynPressure1(simCounter)) + sim_dragForce(simCounter);
                            sim_velocity1(simCounter) = sim_velocity1(simCounter) - ( (chuteDrag{dragCount-1} * sim_dynPressure1(simCounter))/sim_mass1(simCounter) )*dt;
                    end
                end
        end
        
        sim_dynPressure1(simCounter) = savedynPressure;
        sim_velocity1(simCounter)= saveVelocity;
        sim_accel1(simCounter) = (sim_thrust(simCounter)*cos(sim_thrustAngle(simCounter))-sim_dragForce(simCounter)-(sim_mass1(simCounter)*sim_gravity(simCounter)*sin(sim_angle1RAD(simCounter))))/sim_mass1(simCounter);
    else
        sim_dynPressure1(simCounter) = 0.5*sim_atmosDensity(simCounter)*(sim_velocity1(simCounter) - sim_windVelocity(simCounter)*cos(sim_angle1RAD(simCounter)))^2;
        
        sim_dragForce(simCounter) = sim_area(simCounter) * (sim_dragCoefO(simCounter)+sim_dragCoefI(simCounter)) * sim_dynPressure1(simCounter);
                            
        sim_accel1(simCounter) = (sim_thrust(simCounter)*cos(sim_thrustAngle(simCounter))-sim_dragForce(simCounter)-(sim_mass1(simCounter)*sim_gravity(simCounter)*sin(sim_angle1RAD(simCounter))))/sim_mass1(simCounter);
    end
    

        sim_velocity2(simCounter) = sim_velocity1(simCounter) + sim_accel1(simCounter)*dt;
    if ( (numOfChutes > 0) && (descentPhase == true) )
        angleDiff = atan( (sim_windVelocity(simCounter)*sin(sim_angle1RAD(simCounter))) / sim_velocity2(simCounter) );
        sim_velocity2(simCounter) = sqrt( sim_velocity2(simCounter)^2 + (sim_windVelocity(simCounter)*sin(sim_angle1RAD(simCounter)))^2 );
    end
    
    
    sim_deltaAngle(simCounter) = angleDiff + (sim_gravity(simCounter)*dt*cos(sim_angle1RAD(simCounter)))/sim_velocity2(simCounter);
    %sim_deltaAngle_perp(simCounter) = (sim_gravity(simCounter)*dt*cos(sim_angle1RAD(simCounter)))/sim_velocity2(simCounter);
    
    
	if  simCounter > 1
        if  (sim_altitude1(simCounter) < (towerHeight*sin(pitchAngle*(pi/180)))) && (sim_time(simCounter) < burnout(end))
            sim_angle2RAD(simCounter) = sim_angle1RAD(simCounter);
            %sim_angle2RAD_perp(simCounter) = sim_angle1RAD_perp(simCounter);
        else
            sim_angle2RAD(simCounter) = sim_angle1RAD(simCounter) - sim_deltaAngle(simCounter);
            %sim_angle2RAD_perp(simCounter) = sim_angle1RAD_perp(simCounter) - sim_deltaAngle_perp(simCounter);
        end
	end
    
    
        %% SECOND|THIRD|FOURTH|FIFTH BLOCKS
    
    for iterationCount = 2:5
        
        sim_machNum(simCounter) = abs(sim_velocity2(simCounter)/sim_soundSpeed(simCounter));
        
        
            %% Parasite (Zero-Lift) Drag Coefficient
            
if  staticStability == true

    if  ~(((nose_length/R_length) >= 0.2) && ((nose_length/R_length) <= 0.6))
            
            [~,warnID] = lastwarn;
            if  (strcmp(warnID,CDcalcWarnID) == 0) && (sim_dragCoefO(1) == 0)
                warning(CDcalcWarnID,CDcalcWarnMSG)
                setDragCoef = max(strcmp((input('Set the parasite drag coefficient to a constant value (e.g. 0.75)? [Y/N]: ','s')),yes));
                if  setDragCoef == true
                    sim_dragCoefO(1) = input('Enter a value to set the parasite drag to: ');
                end
            end
                    sim_dragCoefO(simCounter+1) = sim_dragCoefO(simCounter);
    end

        if  (sim_machNum(simCounter) > 0) && (setDragCoef == false)
            [sim_dragCoefO(simCounter),critRatioError] = calcParasiteDrag(  sim_altitude1(simCounter),...
                                                            sim_soundSpeed(simCounter),sim_machNum(simCounter),...
                                                            R_length,skinRoughness,bodyMaxDiameter,...
                                                            surfaceArea,fin_rootChord,fin_tipChord,...
                                                            finThickness,fin_span,numOfFins,...
                                                            bodybaseDiameter,nose_length,...
                                                            R_effLength,R_aftLength,finSetNUM,...
                                            numOfLaunchLugs,protub_length,protub_dist,protub_crossArea,protub_surfArea);
                                        
                                        if  isreal(sim_dragCoefO(simCounter)) == false % if parasite drag function fails
                                            if  exist('defaultDragCoefO','var') == false
                                                warning('The parasite drag calculation function has failed to correctly estimate the zero-lift drag coefficient.')
                                                defaultDragCoefO = input('Enter a value for the default drag coefficient to use when the function fails: ');
                                            end
                                            sim_dragCoefO(simCounter) = defaultDragCoefO;
                                        end
                                                        
            if      critRatioError == true
                
                    [~,warnID] = lastwarn;
                    if  strcmp(warnID,critWarnID) == 0
                        warning(critWarnID,critWarnMSG)
                    end
            end
        end

    
elseif  staticStability == false
    
        [~,warnID] = lastwarn;
        if  strcmp(warnID,rocketDataWarnID) == 0
            warning(rocketDataWarnID,rocketDataWarnMSG)
            sim_dragCoefO(1) = input('Enter a value to set the parasite drag to: ');
        end
        
            sim_dragCoefO(simCounter+1) = sim_dragCoefO(simCounter);
end

            %% Lift Coefficient & Force
            
if  (numOfChutes == 0) || ( (numOfChutes > 0) && (descentPhase == false) )
            
        if sim_altitude1(simCounter) > (towerHeight*sin(pitchAngle*(pi/180)))
            resultantFlow = sqrt( sim_velocity2(simCounter)^2 + (sim_windVelocity(simCounter)*sin(sim_angle2RAD(simCounter)))^2 );
            sim_angleOfAttack(simCounter) = acos( min(sim_velocity2(simCounter),resultantFlow) / max(sim_velocity2(simCounter),resultantFlow) );

                moveableFinAOA = finDeflection + sim_angleOfAttack(simCounter);
            
            COPvsAOA = sim_angleOfAttack(simCounter)*(  ( 4 * K * (sum(y_iArray*x_increment)-fin_area) ) / (pi * (2*nose_baseRadius)^2)  );
            if  COP_selected == COP_options.Barrowman_Equations
                sim_COP(simCounter) = ( COP_numerator   + COPvsAOA*centreOfPlanformArea )...
                                    / ( COP_denominator + COPvsAOA );
            else
                sim_COP(simCounter) = centreOfPressure;
            end
            
                moveableFinCOP = fin_COP(moveableFinSet);
        
            
            if  abs(sim_angleOfAttack(simCounter)*(180/pi)) > 10
                sim_liftCoef(simCounter) = liftCoef(sim_angleOfAttack(simCounter));
            else
                sim_liftCoef(simCounter) = (2*pi)*sim_angleOfAttack(simCounter);
            end
            
                if  abs(moveableFinAOA*(180/pi)) > 10
                    moveableFinLiftCoef = liftCoef(moveableFinAOA);
                else
                    moveableFinLiftCoef = (2*pi)*moveableFinAOA;
                end
            
            sim_windLift(simCounter) = 0.5 * sim_atmosDensity(simCounter) * resultantFlow^2 * sum(finArea) * sim_liftCoef(simCounter);
            
            if controlSystem == control.Moveable_Fins
                moveableFinLift = fin_normalForce(moveableFinSet) - ...
                    0.5 * sim_atmosDensity(simCounter) * resultantFlow^2 * finArea(moveableFinSet) * moveableFinLiftCoef;
            else
                moveableFinLift = 0;
            end
            %{
            if numOfChutes > 0
                for liftCount = 1:numOfChutes+1
                    if      liftCount == 1
                            angularAccel = ( sim_windLift(simCounter)*(sim_COP(simCounter)-sim_COM(simCounter)) ) / MOInertia;
                    elseif  (deployTime(liftCount-1) ~= 0) && (sim_time(simCounter) >= deployTime(liftCount-1))
                            chuteNormalForce = ( chuteDrag{liftCount-1} )* 0.5*sim_atmosDensity(simCounter)*...
                                                ( sim_windVelocity(simCounter)*sin(sim_angle2RAD(simCounter)) )^2;
                            angularAccel = sum([(chuteNormalForce*(sim_COM(simCounter)-nose_length)),...
                                                -sim_windLift(simCounter)*(sim_COP(simCounter)-sim_COM(simCounter))]) / MOInertia;
                    end
                end
            else
            %}
            angularAccel = ( (sim_windLift(simCounter) * (sim_COP(simCounter)-sim_COM(simCounter)))...
                            + (sim_thrust(simCounter)*sin(sim_thrustAngle(simCounter))) * (R_length-sim_COM(simCounter))...
                            + (moveableFinLift * (moveableFinCOP-sim_COM(simCounter))) )...
                            / MOInertia;  % if staticMargin is positive, should rotate anticlockwise
            %end
            
            
            angleDiff = angularVelo*dt + 0.5*angularAccel*dt^2;
            angularVelo = angularVelo + angularAccel*dt;
        else
            sim_COP(simCounter) = centreOfPressure;
            angleDiff = 0;
        end
        
        sim_stability(simCounter) = ( sim_COP(simCounter) - sim_COM(simCounter) ) / bodyAvgDiameter;
        
        
            %% Lift-Induced Drag Coefficient
            
        ozFactor = oswaldFactor_mod(aspectRatio,sweepAngle,'Average',...
                    sim_dragCoefO(simCounter),(bodyAvgDiameter/finSpanTotal(1)),1);
            
        sim_dragCoefI(simCounter) = (sim_liftCoef(simCounter)^2) / ...
                                    (pi*aspectRatio*ozFactor);
                                
else
    sim_COP(simCounter) = centreOfPressure;
    angleDiff = 0;
end

            %% Remaining calculations...
            
        if numOfChutes > 0
            saveVelocity = sim_velocity2(simCounter);
            
            for dragCount = 1:numOfChutes+1
                sim_dynPressure2(simCounter) = 0.5*sim_atmosDensity(simCounter)*(sim_velocity2(simCounter) + sim_windVelocity(simCounter)*cos(sim_angle2RAD(simCounter)))^2;
            
                if      dragCount == 1
                        savedynPressure = sim_dynPressure2(simCounter);
                        sim_dragForce(simCounter) = ( sim_area(simCounter)*(sim_dragCoefO(simCounter)+sim_dragCoefI(simCounter)) )* sim_dynPressure2(simCounter);
                        sim_accel2(simCounter) = (sim_thrust(simCounter)*cos(sim_thrustAngle(simCounter))-sim_dragForce(simCounter)-(sim_mass2(simCounter)*sim_gravity(simCounter)*sin(sim_angle2RAD(simCounter))))/sim_mass2(simCounter);
                        sim_velocity2(simCounter) = sim_velocity2(simCounter) + sim_accel2(simCounter)*dt;
                elseif  (deployTime(dragCount-1) ~= 0) && (sim_time(simCounter) >= deployTime(dragCount-1))
                    switch dragCalcMethod
                        case 1
                            sim_dragForce(simCounter) = chuteDrag{dragCount-1} * sim_dynPressure2(simCounter);
                        case 2
                            sim_dragForce(simCounter) = (chuteDrag{dragCount-1} * sim_dynPressure2(simCounter)) + sim_dragForce(simCounter);
                            sim_velocity2(simCounter) = sim_velocity2(simCounter) - ( (chuteDrag{dragCount-1} * sim_dynPressure2(simCounter))/sim_mass2(simCounter) )*dt;
                    end
                end
            end
            
            sim_dynPressure2(simCounter) = savedynPressure;
            sim_velocity2(simCounter)= saveVelocity;
            sim_accel2(simCounter) = (sim_thrust(simCounter)*cos(sim_thrustAngle(simCounter))-sim_dragForce(simCounter)-(sim_mass2(simCounter)*sim_gravity(simCounter)*sin(sim_angle2RAD(simCounter))))/sim_mass2(simCounter);
        else
            sim_dynPressure2(simCounter) = 0.5*sim_atmosDensity(simCounter)*(sim_velocity2(simCounter) - sim_windVelocity(simCounter)*cos(sim_angle2RAD(simCounter)))^2;
        
            sim_dragForce(simCounter) = sim_area(simCounter) * (sim_dragCoefO(simCounter)+sim_dragCoefI(simCounter)) * sim_dynPressure2(simCounter);
                            
            sim_accel2(simCounter) = (sim_thrust(simCounter)*cos(sim_thrustAngle(simCounter))-sim_dragForce(simCounter)-(sim_mass2(simCounter)*sim_gravity(simCounter)*sin(sim_angle2RAD(simCounter))))/sim_mass2(simCounter);
        end
        

            sim_velocity2(simCounter) = sim_velocity1(simCounter) + (0.5*(sim_accel1(simCounter)+sim_accel2(simCounter))*dt);
        if ( (numOfChutes > 0) && (descentPhase == true) )
            angleDiff = atan( (sim_windVelocity(simCounter)*sin(sim_angle1RAD(simCounter))) / sim_velocity2(simCounter) );
            sim_velocity2(simCounter) = sqrt( sim_velocity2(simCounter)^2 + (sim_windVelocity(simCounter)*sin(sim_angle1RAD(simCounter)))^2 );
        end

        
        if  iterationCount == 4
            convergeTestV2 = sim_velocity2(simCounter);
        end
        
        if iterationCount == 5
            break            
        else
            sim_deltaAngle(simCounter) = angleDiff + (sim_gravity(simCounter)*dt*cos(sim_angle1RAD(simCounter)))/sim_velocity2(simCounter);
            %sim_deltaAngle_perp(simCounter) = (sim_gravity(simCounter)*dt*cos(sim_angle2RAD(simCounter)))/sim_velocity2(simCounter);
            
            if  (sim_altitude1(simCounter) < (towerHeight*sin(pitchAngle*(pi/180)))) && (sim_time(simCounter) < burnout(end))
                sim_angle2RAD(simCounter) = sim_angle1RAD(simCounter);
                %sim_angle2RAD_perp(simCounter) = sim_angle1RAD_perp(simCounter);
            else
                sim_angle2RAD(simCounter) = sim_angle1RAD(simCounter) - sim_deltaAngle(simCounter);
                %sim_angle2RAD_perp(simCounter) = sim_angle1RAD_perp(simCounter) - sim_deltaAngle_perp(simCounter);
            end
        end       
    end

    sim_convergence(simCounter) = convergeTestV2/sim_velocity2(simCounter);

    %% Post-Iterations

    if  sim_altitude2(simCounter) < 0
        sim_accel2(simCounter) = 0;
    end
    
    sim_deltaAngle(simCounter) = ( sim_gravity(simCounter) * dt * cos(sim_angle2RAD(simCounter)) ) / sim_velocity2(simCounter);
    %sim_deltaAngle_perp(simCounter) = ( sim_gravity(simCounter) * dt * cos(sim_angle2RAD_perp(simCounter)) ) / sim_velocity2(simCounter);

    
	if  (sim_altitude1(simCounter) < (towerHeight*sin(pitchAngle*(pi/180)))) && (sim_time(simCounter) < burnout(end))
        sim_angle2RAD(simCounter) = sim_angle1RAD(simCounter);
        %sim_angle2RAD_perp(simCounter) = sim_angle1RAD_perp(simCounter);
    else
        sim_angle2RAD(simCounter) = sim_angle1RAD(simCounter) - sim_deltaAngle(simCounter);
        %sim_angle2RAD_perp(simCounter) = sim_angle1RAD_perp(simCounter) - sim_deltaAngle_perp(simCounter);
	end
    
    
    sim_angle2DEG(simCounter) = sim_angle2RAD(simCounter)*(180/pi);
    %sim_angle2DEG_perp(simCounter) = sim_angle2RAD_perp(simCounter)*(180/pi);
    
    sim_angleAvg(simCounter) = ( sim_angle1RAD(simCounter) + sim_angle2RAD(simCounter) )/2;
    %sim_angleAvg_perp(simCounter) = ( sim_angle1RAD_perp(simCounter) + sim_angle2RAD_perp(simCounter) )/2;
    
    
    sim_deltaDisplace(simCounter) = (dt*sim_velocity1(simCounter))+(0.25*(dt^2)*(sim_accel1(simCounter)+sim_accel2(simCounter)));
    
    sim_range2(simCounter) = sim_range1(simCounter)+(sim_deltaDisplace(simCounter)*cos(sim_angleAvg(simCounter)));
    
    sim_altitude2(simCounter) = sim_altitude1(simCounter)+(sim_deltaDisplace(simCounter)*sin(sim_angleAvg(simCounter)));
    
    %sim_range2_perp(simCounter) = sim_range1_perp(simCounter) + ((sim_range2(simCounter)-sim_range1(simCounter))*tan(sim_angleAvg_perp(simCounter)));
    
    if  staticStability == true
        sim_dynViscosity(simCounter) = interp1(dynViscosity(:,1),dynViscosity(:,2),sim_atmosTemp(simCounter),'linear')*10^-5;
        sim_Reynolds(simCounter) = ( sim_atmosDensity(simCounter)*sim_velocity2(simCounter)*R_length )/sim_dynViscosity(simCounter);
    end

    
    %% Energy    
    
    sim_potEnergy(simCounter) = (sim_mass1(simCounter)*sim_gravity(simCounter)*sim_altitude1(simCounter))/1000000;
    
    sim_kinEnergy(simCounter) = (0.5*sim_mass1(simCounter)*(sim_velocity1(simCounter)^2))/1000000;
    
    sim_totalEnergy(simCounter) = sim_potEnergy(simCounter) + sim_kinEnergy(simCounter);


    %% Debugging %%%DEV:NB%%%
    %{
    if  ~isreal(sim_...(simCounter)) || abs() == Inf || isnan() == true
        error('Something went wrong...')
    end
	%}
    
end

sim_flightPhase = sim_flightPhase(2:end);
clear defaultDragCoefO

elapsedTime = toc;
disp(['Time elapsed whilst running flight simulation: ',num2str(elapsedTime),' seconds.'])

if  (simCounter == numOfSteps) && (sim_altitude2(end) ~= 0)
    warning('%s\n','The simulation ended at the predicted flight time, but the rocket was still moving.',...
                   'You may want to re-run the simulation with an increased flight time prediction.')
end

%%%
end % {flight analysis}


%% RESULTS [res12]
%%
%
if  flight_analysis == true
%%%

for simZeroCount = 1:length(sim_altitude1)
    if  sim_altitude2(simZeroCount) < 0
        sim_altitude2(simZeroCount) = 0;
    end
    if  sim_velocity2(simZeroCount) < 0
        sim_velocity2(simZeroCount) = 0;
    end
end

for simEndCount = length(sim_altitude1):-1:2
    if  (sim_altitude1(simEndCount) <= 0) && (sim_altitude1(simEndCount-1) > 0)
        simEndIndex = simEndCount;
        break
    end
end

    if  exist('simEndIndex','var') == 0
        simEndIndex = length(sim_altitude1);
    end
    

    sim_time        = sim_time(1:simEndIndex);
    sim_thrust      = sim_thrust(1:simEndIndex);
    sim_massFlow    = sim_massFlow(1:simEndIndex);
    sim_area        = sim_area(1:simEndIndex);
    sim_mass1       = sim_mass1(1:simEndIndex);
    sim_mass2       = sim_mass2(1:simEndIndex);
    sim_gravity     = sim_gravity(1:simEndIndex);
    sim_velocity1   = sim_velocity1(1:simEndIndex);
    sim_angle1RAD   = sim_angle1RAD(1:simEndIndex);
    sim_angle1DEG   = sim_angle1DEG(1:simEndIndex);
    sim_range1      = sim_range1(1:simEndIndex);
    sim_altitude1   = sim_altitude1(1:simEndIndex);
    sim_accel2      = sim_accel2(1:simEndIndex);
    sim_velocity2   = sim_velocity2(1:simEndIndex);
    sim_machNum     = sim_machNum(1:simEndIndex);
    sim_deltaAngle  = sim_deltaAngle(1:simEndIndex);
    sim_angle2RAD   = sim_angle2RAD(1:simEndIndex);
    sim_angle2DEG   = sim_angle2DEG(1:simEndIndex);
    sim_angleAvg    = sim_angleAvg(1:simEndIndex);
    sim_deltaDisplace = sim_deltaDisplace(1:simEndIndex);
    sim_range2      = sim_range2(1:simEndIndex);
    sim_altitude2   = sim_altitude2(1:simEndIndex);
    sim_atmosDensity= sim_atmosDensity(1:simEndIndex);
    sim_atmosTemp   = sim_atmosTemp(1:simEndIndex);
    sim_soundSpeed  = sim_soundSpeed(1:simEndIndex);
    sim_dragCoefO   = sim_dragCoefO(1:simEndIndex);
    sim_dragCoefI   = sim_dragCoefI(1:simEndIndex);
    sim_angleOfAttack = sim_angleOfAttack(1:simEndIndex);
    sim_liftCoef    = sim_liftCoef(1:simEndIndex);
    sim_dynPressure1= sim_dynPressure1(1:simEndIndex);
    sim_dynPressure2= sim_dynPressure2(1:simEndIndex);
    sim_dragForce   = sim_dragForce(1:simEndIndex);
    sim_accel1      = sim_accel1(1:simEndIndex);
    sim_convergence = sim_convergence(1:simEndIndex);
    sim_phaseEnd    = sim_phaseEnd(1:simEndIndex);
    sim_flightPhase = sim_flightPhase(1:simEndIndex);
    sim_phaseMass   = sim_phaseMass(1:simEndIndex);
    sim_potEnergy   = sim_potEnergy(1:simEndIndex);
    sim_kinEnergy   = sim_kinEnergy(1:simEndIndex);
    sim_totalEnergy = sim_totalEnergy(1:simEndIndex);
    sim_dynViscosity= sim_dynViscosity(1:simEndIndex);
    sim_Reynolds    = sim_Reynolds(1:simEndIndex);
    sim_windLift    = sim_windLift(1:simEndIndex);
    sim_windVelocity= sim_windVelocity(1:simEndIndex);
    sim_COP         = sim_COP(1:simEndIndex);
    sim_COM         = sim_COM(1:simEndIndex);
    sim_MOI         = sim_MOI(1:simEndIndex);
    sim_stability   = sim_stability(1:simEndIndex);

    
[ApoValue,ApoIndex] = max(sim_altitude2);
[RanValue,RanIndex] = max(abs(sim_range2));
[SpeValue,SpeIndex] = max(sim_velocity2);
[MacValue,MacIndex] = max(sim_machNum);
[AccValue,AccIndex] = max(sim_accel2);
[DynValue,DynIndex] = max(sim_dynPressure2);
[DraValue,DraIndex] = max(sim_dragForce);
[PotValue,PotIndex] = max(sim_potEnergy);
[kinValue,kinIndex] = max(sim_kinEnergy);
[totValue,totIndex] = max(sim_totalEnergy);


Results.Apogee               = [num2str(ApoValue/10^3), ' kilometres @ ', ...
                                num2str(sim_time(ApoIndex)),    ' seconds'];
Results.max_Range            = [num2str(RanValue), ' metres @ ', ...
                                num2str(sim_time(RanIndex)),    ' seconds'];
Results.max_Speed            = [num2str(SpeValue), ' m/s @ ', ...
                                num2str(sim_time(SpeIndex)),    ' seconds'];
Results.max_Mach_Number      = [num2str(MacValue), ' @ ', ...
                                num2str(sim_time(MacIndex)),    ' seconds'];
Results.max_Acceleration     = [num2str(AccValue), ' m/s^2 @ ', ...
                                num2str(sim_time(AccIndex)),    ' seconds'];
Results.max_Dynamic_Pressure = [num2str(DynValue/10^6), ' MPascals @ ', ...
                                num2str(sim_time(DynIndex)),    ' seconds'];
Results.max_Drag_Force       = [num2str(DraValue/10^3), ' KNewtons @ ', ...
                                num2str(sim_time(DraIndex)),    ' seconds'];
Results.max_Potential_Energy = [num2str(PotValue), ' MJoules @ ', ...
                                num2str(sim_time(PotIndex)),    ' seconds'];
Results.max_Kinetic_Energy   = [num2str(kinValue), ' MJoules @ ', ...
                                num2str(sim_time(kinIndex)),    ' seconds'];
Results.max_Total_Energy     = [num2str(totValue), ' MJoules @ ', ...
                                num2str(sim_time(totIndex)),    ' seconds'];
Results.Time_of_Flight       = [num2str(sim_time(end)),        ' seconds'];

Results.Ground_Impact_Speed  = [num2str(sim_velocity2(end-1)),    ' m/s'];


pauseFix = input('Hit enter to view the results of the simulation.');
%drawnow;pause;drawnow % fix for hanging on pause issue with matlab r2016a
disp(Results)
fprintf('%s\n',     'Maximum output values - with associated times during flight - are summarised above,', ...
                    'and have been saved to a text file named ''Results.txt''.' )

maxResult_Names     = {  'Apogee = ';'Max. Range = ';'Max. Speed = ';'Max. Mach Number = ';'Max. Acceleration = '; ...
                         'Max. Dynamic Pressure = ';'Max. Drag Force = ';'Max. Potential Energy = '; ...
                         'Max. Kinetic Energy = ';'Max. Total Energy = ';'Time of Flight = ';'Ground Impact Speed = '};
                     
maxResult_Cell(:,1)                         = maxResult_Names;
maxResult_Cell(:,2)                         = struct2cell(Results);
maxResult_Table                             = cell2table(maxResult_Cell);
maxResult_Table.Properties.VariableNames{1} = 'Variable';
maxResult_Table.Properties.VariableNames{2} = 'Value_at_time';
writetable(maxResult_Table,['Results_(',num2str(simRunNum),')'])


totalResults = [sim_time.',sim_thrust.',sim_massFlow.',sim_mass2.',sim_gravity.', ...
                sim_altitude2.',sim_range2.',sim_angle2DEG.',sim_velocity2.',sim_accel2.', ...
                sim_atmosDensity.',sim_atmosTemp.',sim_soundSpeed.',sim_machNum.', ...
                sim_dragCoefO.',sim_dynPressure2.',sim_dragForce.',sim_convergence.', ...
                sim_potEnergy.',sim_kinEnergy.',sim_totalEnergy.',sim_Reynolds.',sim_stability.'];
            
resultNames = { 'Time_(s)','Thrust_(N)','Mass_Flow_(kg/s)','Mass_(kg)','Gravity_(m/s/s)', ...
                'Altitude_(m)','Lateral_Range_(m)','Flight_Path_Angle_(degrees)','Velocity_(m/s)','Acceleration_(m/s/s)', ...
                'Atmospheric_Density_(kg/m3)','Atmospheric_Temperature_(K)','Speed_of_Sound_(m/s)','Mach_Number', ...
                'Drag_Coefficient','Dynamic_Pressure_(Pa)','Drag_(N)','Evidence_of_Convergence', ...
                'Potential_Energy_(MJ)','Kinetic_Energy_(MJ)','Total_Energy_(MJ)','Reynolds_Number','Stability Calibers'};

excelCellRange = {  'A2:A','B2:B','C2:C','D2:D','E2:E','F2:F','G2:G','H2:H', ...
                    'I2:I','J2:J','K2:K','L2:L','M2:M','N2:N','O2:O','P2:P', ...
                    'Q2:Q','R2:R','S2:S','T2:T','U2:U','V2:V','W2:W' };

disp('Saving output data, please wait...')
tic
dlmwrite(['All_Results_(',num2str(simRunNum),').csv'],totalResults,'delimiter',',','precision',4);
%
warning('off','MATLAB:xlswrite:NoCOMServer')
xlswrite(['All_Results_(',num2str(simRunNum),')'],resultNames);
    for resultsCount = 1:(size(totalResults,2))
        xlswrite( ['All_Results_(',num2str(simRunNum),')'], totalResults((1:end),resultsCount), [cell2mat(excelCellRange(resultsCount)),num2str(length(sim_time))] );
    end
[warnMSG,warnID] = lastwarn;
    if  strcmp(warnID,'MATLAB:xlswrite:NoCOMServer')==1
        warning('Could not save results data to spreadsheet file - you need MS Excel installed! Results saved to a comma-separated values (.csv) instead.')
    else
        disp('All results data have been saved to a spreadsheet named ''All_Results.xls''.')
    end
warning('on','MATLAB:xlswrite:NoCOMServer')
%
elapsedTime = toc;
disp(['Time elapsed whilst saving results data: ',num2str(elapsedTime),' seconds.'])

%%%
end % {flight analysis}


%% MICROGRAVITY [mic13]
%%
%
if  flight_analysis == true
    microG_analysis =  max(strcmp((input('Do you want an analysis of the microgravity experienced by the rocket? [Y/N]: ','s')),yes));
if  microG_analysis == true
%
microG_array = 0*sim_time;
microG = (10^-6)*9.807;
% high quality microgravity value, using sea-level gravitational constant

%% Simulation

for microG_count = 1:length(sim_time)

    if      abs(sim_accel2(microG_count)) <= microG*10^0
            microG_array(microG_count) = microG*10^0;
            microG_6_time = microG_6_time + timeStep;
            
    elseif  abs(sim_accel2(microG_count)) <= microG*10^1
            microG_array(microG_count) = microG*10^1;
            microG_5_time = microG_5_time + timeStep;
            
    elseif  abs(sim_accel2(microG_count)) <= microG*10^2
            microG_array(microG_count) = microG*10^2;
            microG_4_time = microG_4_time + timeStep;
            
    elseif  abs(sim_accel2(microG_count)) <= microG*10^3
            microG_array(microG_count) = microG*10^3;
            microG_3_time = microG_3_time + timeStep;
            
    elseif  abs(sim_accel2(microG_count)) <= microG*10^4
            microG_array(microG_count) = microG*10^4;
            microG_2_time = microG_2_time + timeStep;
            
    elseif  abs(sim_accel2(microG_count)) <= microG*10^5
            microG_array(microG_count) = microG*10^5;
            microG_1_time = microG_1_time + timeStep;
    end
    
end

Quality = ['1x10^-6g';'1x10^-5g';'1x10^-4g';'1x10^-3g';'1x10^-2g';'1x10^-1g'];
Duration = [microG_6_time;microG_5_time;microG_4_time; ...
            microG_3_time;microG_2_time;microG_1_time];
microG_table = table(Quality,Duration);


%% Results

close all
commandwindow
pauseFix = input('Hit enter to view microgravity results.');
%drawnow;pause;drawnow % fix for hanging on pause issue with matlab r2016a

fprintf('%s\n',[]); disp(microG_table)
writetable(microG_table,['Microgravity_results(',num2str(simRunNum),')'])
fprintf('%s\n',['The length of time (s) and associated quality of microgravity ', ...
    'the rocket will experience is summarised above: '],['(each lower quality ',...
    'is calculated exclusively from the higher quality preceding it)'])

pauseFix = input('Hit enter to view microgravity graphs.');
%drawnow;pause;drawnow % fix for hanging on pause issue with matlab r2016a

figure(1)

microGlim6  =     microG*10^0 * ones(1,length(sim_accel2));
microGlim5  =     microG*10^1 * ones(1,length(sim_accel2));
microGlim4  =     microG*10^2 * ones(1,length(sim_accel2));
microGlim3  =     microG*10^3 * ones(1,length(sim_accel2));
microGlim2  =     microG*10^4 * ones(1,length(sim_accel2));
microGlim1  =     microG*10^5 * ones(1,length(sim_accel2));
microGlim0  =     microG*10^6 * ones(1,length(sim_accel2));

Microgravity = semilogy(sim_time,abs(sim_accel2),'-k',...
    sim_time,microGlim6,':k',sim_time,microGlim5,':k',...
    sim_time,microGlim4,':k',sim_time,microGlim3,':k',...
    sim_time,microGlim2,':k',sim_time,microGlim1,':k',...
    sim_time,microGlim0,':k'    );

title('Microgravity (I)')
axis([ 0, max(sim_time), 0, 1.1*microG*10^6 ])
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
xlabel('Time /s')
ylabel('Net Acceleration (log scale) /m/s/s')

text(sim_time(simEndIndex),microG*10^5,' 1x10^{-1}g')
text(sim_time(simEndIndex),microG*10^4,' 1x10^{-2}g')
text(sim_time(simEndIndex),microG*10^3,' 1x10^{-3}g')
text(sim_time(simEndIndex),microG*10^2,' 1x10^{-4}g')
text(sim_time(simEndIndex),microG*10^1,' 1x10^{-5}g')
text(sim_time(simEndIndex),microG*10^0,' 1x10^{-6}g')

microgFig = gcf; % gcf = get current figure - needed for figures with multiple graphs
saveas(microgFig,['Microgravity_I_(',num2str(simRunNum),').jpg'])
fprintf('%s\n','Figure 1: Net acceleration values between the horizontal lines indicate microgravity',...
    ' at varying levels of quality.')

tic; while toc < 2; end

figure(2)
Microgravity_II = semilogy(sim_time,microG_array,'-k');
Microgravity_II.LineWidth=1.25;
title('Microgravity (II)')
axis([0,max(sim_time),0,microG*10^6])
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [],'OuterPosition',[0,0,0.8,1]);
xlabel('Time /s')
ylabel('Microgravity Quality (log scale) ')

text(sim_time(simEndIndex),microG*10^5,' Net Accel. <= 1x10^{-1}g')
text(sim_time(simEndIndex),microG*10^4,' Net Accel. <= 1x10^{-2}g')
text(sim_time(simEndIndex),microG*10^3,' Net Accel. <= 1x10^{-3}g')
text(sim_time(simEndIndex),microG*10^2,' Net Accel. <= 1x10^{-4}g')
text(sim_time(simEndIndex),microG*10^1,' Net Accel. <= 1x10^{-5}g')
text(sim_time(simEndIndex),microG*10^0,' Net Accel. <= 1x10^{-6}g')

saveas(Microgravity_II,['Microgravity_II_(',num2str(simRunNum),').jpg'])
fprintf('%s\n','Figure 2: Microgravity quality over the trajectory (acceleration values are',...
    'normalized to a quality band limit','Both graphs use absolute values for acceleration.')

tic; while toc < 2; end

close all

%%%
end % {microgravity analysis...}
end % {...only if flight analysis chosen}


%% GRAPHS [gra14]
%%
%
if  flight_analysis == true
%%% %%%

commandwindow
pauseFix = input('All graphs will be saved as JPEG image files. Hit enter to view graphs of all other results.');
%drawnow;pause;drawnow % fix for hanging on pause issue with matlab r2016a

makeGraph('Line',simRunNum,'Trajectory Profile','Range (Lateral Distance) (m)','Altitude (m)',abs(sim_range2),sim_altitude2)
makeGraph('Line',simRunNum,'Altitude','Time (s)','Altitude (m)',sim_time,sim_altitude2)
makeGraph('Line',simRunNum,'Speed','Time (s)','Speed (m/s)',sim_time,sim_velocity2)
makeGraph('Line',simRunNum,'Acceleration','Time (s)','Acceleration (m/s/s)',sim_time,sim_accel2)
makeGraph('Line',simRunNum,'Mass','Time (s)','Mass (kg)',sim_time,sim_mass2)
makeGraph('Line',simRunNum,'Mach Number','Time (s)','Mach number',sim_time,sim_machNum)
makeGraph('Line',simRunNum,'Zero-Lift Drag Coefficient','Time (s)','Drag coefficient',sim_time,sim_dragCoefO)
makeGraph('Line',simRunNum,'Drag Force','Time (s)','Drag force (N)',sim_time,sim_dragForce)
makeGraph('Line',simRunNum,'Dynamic Pressure','Time (s)','Dynamic Pressure (Pa)',sim_time,sim_dynPressure2)
makeGraph('Line',simRunNum,'Atmospheric Density','Time (s)','Density (kg/m3)',sim_time,sim_atmosDensity)
makeGraph('Line',simRunNum,'Velocity Vector Angle','Time (s)','Pitch angle (deg)',sim_time,sim_angle2DEG)
makeGraph('Line',simRunNum,'Energy','Time (s)','Energy (MJ)','Potential','Kinetic','Total', ...
          sim_time,sim_potEnergy,sim_time,sim_kinEnergy,sim_time,sim_totalEnergy)
makeGraph('Line',simRunNum,'Static Margin','Time (s)','Margin (m)',sim_time,sim_stability*bodyAvgDiameter)
makeGraph('Line',simRunNum,'Stability Calibers','Time (s)','Calibers',sim_time,sim_stability)

    %%
    disp('Lift coefficient calculated for a flat plate (symmetric airfoil) based on NACA-0015 results.')
    makeGraph('Scatter',simRunNum,'Lift Coefficient vs. Angle of Attack','Angle of Attack (deg)','Lift Coefficient',sim_angleOfAttack*(180/pi),sim_liftCoef)
    makeGraph('Scatter',simRunNum,'Lift-to-Drag Ratio (Drag Polar)','Drag (coefficient, parasite+induced)','Lift (coefficient)',(sim_dragCoefO+sim_dragCoefI),sim_liftCoef)
    makeGraph('Scatter',simRunNum,'Drag Curve (Parasite vs. Lift-Induced)','Velocity (m/s)','Drag Coefficient','Parasitic Drag','Lift-Induced Drag','Total Drag Coefficient', ...
              sim_velocity2,sim_dragCoefO,sim_velocity2,sim_dragCoefI,sim_velocity2,(sim_dragCoefO+sim_dragCoefI))

	%%
if  (wind_velocity ~= 0) && (windEffect_method == true)
    
    figure
    plot(abs(sim_windVelocity(1:ApoIndex)),sim_altitude1(1:ApoIndex),abs(windVelocity_preAlt(1:ApoIndex)),sim_altitude1(1:ApoIndex))
    axis([0,max(abs(sim_windVelocity(1:ApoIndex))),0,max(sim_altitude1(1:ApoIndex))])
    title('Wind Speed Irrespective of Altitude vs. WRT Altitude');
	xlabel('Wind speed (m/s)'); ylabel('Altitude (m)')
	legend('With Respect to Altitude','Irrespective of Altitude')
    saveas(gcf,['Wind Speed Irrespective of Altitude vs. WRT Altitude','_(',num2str(simRunNum),').jpg'])
    tic; while toc < 1.5; end; close
    
    
    fprintf('%s\n','Wind graphs use velocity data taken irrespective of altitude ',...
        '(the simulation takes the wind velocity gradient with altitude into account)')
	makeGraph('Line',simRunNum,'Randomized Wind Speed Data','Simulation time (s)','Wind speed (m/s)',...
        'Randomized Rayleigh distribution wind speeds',...
        ['Average wind speed based on random data = ',num2str(sum(abs(windVelocity_preAlt))/length(abs(windVelocity_preAlt))),' m/s'],...
        ['Actual average wind speed from user input = ',num2str(abs(wind_velocity)),' m/s'],...
          sim_time(1:simEndIndex),abs(windVelocity_preAlt(1:simEndIndex)),...
          sim_time(1:simEndIndex),ones(1,simEndIndex).*(sum(abs(windVelocity_preAlt))/length(abs(windVelocity_preAlt))),...
          sim_time(1:simEndIndex),ones(1,simEndIndex).*(abs(wind_velocity)))
      
      windSpeedAxis = linspace(0,max(abs(windVelocity_preAlt)));
      %probDens_Fnct = pdf((makedist('Weibull','a',scaleParameter,'b',shapeParameter)),windSpeedAxis); 
      %%%DEV:NB%%% built-in matlab function fails under certain circumstances
      probDens_Fnct = (shapeParameter/scaleParameter).*((windSpeedAxis./scaleParameter).^(shapeParameter-1)).*exp(-windSpeedAxis/scaleParameter).^shapeParameter;
      
      clear probDens_Data
      for   countA = 1:length(windSpeedAxis)
            counter = 0;
            for countB = 1:length(windVelocity_preAlt)
                if  ( (countA == 1)&& (abs(windVelocity_preAlt(countB)) <= windSpeedAxis(countA)) ) ...
                 || ( (countA > 1) && (abs(windVelocity_preAlt(countB)) > windSpeedAxis(countA-1))  ...
                                   && (abs(windVelocity_preAlt(countB)) <= windSpeedAxis(countA)) )
                    counter = counter + 1;
                end
            end
            probDens_Data{countA} = counter;
      end
            probDens_Data =     cell2mat(probDens_Data)/sum(cell2mat(probDens_Data));
      
            figure(1)
            bar((windSpeedAxis-(windSpeedAxis(2)-windSpeedAxis(1))/2),probDens_Data,1); hold on
            plot((windSpeedAxis-(windSpeedAxis(2)-windSpeedAxis(1))/2),probDens_Fnct/2,'Color','c','LineWidth',2); hold on
            plot([(sum(abs(windVelocity_preAlt))/length(windVelocity_preAlt)),(sum(abs(windVelocity_preAlt))/length(windVelocity_preAlt))],...
                [0,max(probDens_Fnct)],'LineWidth',2,'Color','r'); hold on
            plot([abs(wind_velocity),abs(wind_velocity)],[0,max(probDens_Fnct)],'LineWidth',2,'Color','y')
            axis([0,max(windSpeedAxis),0,max(probDens_Fnct)/2])
            set(gca,'Color',[0.8 0.8 0.8]);
            title('Wind Speed Probability Density Comparison');
            xlabel('Wind speed (m/s)'); ylabel('Probability Density')
            legend('From Randomized Data','From Weibull Function',...
                ['Average speed (calculated) = ',num2str(sum(abs(windVelocity_preAlt))/length(windVelocity_preAlt)),' m/s'],...
                ['Average speed (user input) = ',num2str(abs(wind_velocity)),' m/s'])
            saveas(gcf,['Wind Speed Probability Density Comparison','_(',num2str(simRunNum),').jpg'])
            tic; while toc < 1.5; end; close

            
elseif (wind_velocity ~= 0) && (windEffect_method == false)
    
    windSpeedAxis_ext = linspace(0,gust_velocity,1000);
    cumuProbDens_Fnct = 1-exp(-(windSpeedAxis_ext./scaleParameter).^shapeParameter);
    sigma = ( (gust_velocity-wind_velocity)/3 )/wind_velocity;
    windStandardDev = wind_velocity*sigma; % assumes gust is 3 S.D.s from mean
    
    %{
    descentMass = sum(stage_mass-stage_propMass);
    descentGravity = sum(sim_gravity(ApoIndex:end))/length(sim_gravity(ApoIndex:end));
    descentAtmosDens = sum(sim_atmosDensity(ApoIndex:end))/length(sim_atmosDensity(ApoIndex:end));
    terminalVelocity = sqrt((2*descentMass*descentGravity)/(chuteDrag{end}*descentAtmosDens));
    actualVelocity = sqrt(terminalVelocity^2 + wind_velocity^2);
    ballisticLanding = 2*sim_range2(ApoIndex);
    angleSubtended = atan(terminalVelocity/wind_velocity);
    landingPosition = ( sim_altitude2(ApoIndex)/tan(angleSubtended) ) + sim_range2(ApoIndex);
    %}
    
    landingPosition = max(abs(sim_range2));
    
    circleAngles = linspace(0,2*pi);
    circlePoints = [cos(circleAngles);sin(circleAngles)];
    plot(landingPosition,0,'Marker','*'); hold on
    for stdDevCount = 1:3
        SDplusIndex{stdDevCount} = max(find(windSpeedAxis_ext <= (wind_velocity+stdDevCount*windStandardDev)));
        SDminusIndex{stdDevCount} = max(find(windSpeedAxis_ext <= (wind_velocity-stdDevCount*windStandardDev)));
        if ~isempty(SDminusIndex{stdDevCount})
            SDprob{stdDevCount} = 100*( cumuProbDens_Fnct(SDplusIndex{stdDevCount})-cumuProbDens_Fnct(SDminusIndex{stdDevCount}) );
        else
            SDprob{stdDevCount} = 100*( cumuProbDens_Fnct(SDplusIndex{stdDevCount}) );
        end
        
        %landingPos{stdDevCount} = ( sim_altitude2(ApoIndex) / ...
        %                            (terminalVelocity/(wind_velocity+windStandardDev*stdDevCount)) )...
        %                           + sim_range2(ApoIndex);
        %landingRadius{stdDevCount} = landingPos{stdDevCount} - landingPosition;
        landingRadius{stdDevCount} = landingPosition*((stdDevCount*windStandardDev)/wind_velocity);
        plot(landingPosition+circlePoints(1,:)*landingRadius{stdDevCount},...
             circlePoints(2,:)*landingRadius{stdDevCount})
        hold on
    end
    plot([0,landingPosition],[0,0],'--')
    %axis([0,landingPosition*2,-landingPosition,landingPosition])
    axis square; xlabel('Lateral Range (m)'); ylabel('Longitudinal Range (m)')
    legend('Expected Landing Position',...
          ['1 S.D. from expectation, ',num2str(SDprob{1}),'% of data'],...
          ['2 S.D.s from expectation, ',num2str(SDprob{2}),'% of data'],...
          ['3 S.D.s from expectation, ',num2str(SDprob{3}),'% of data'],...
           'Ground Track'); title('Probable Landing Zones')
    saveas(gcf,['Probable Landing Zones','_(',num2str(simRunNum),').jpg'])
    tic; while toc < 1.5; end; close
end
      
    %%
	maxMachIndex = find(sim_machNum == max(sim_machNum));
makeGraph('Line',simRunNum,'Mach Number vs. Parasite Drag Coefficient, up to max. Mach No.','Mach number',...
    'Zero-Lift Drag Coefficient',sim_machNum(1:maxMachIndex),sim_dragCoefO(1:maxMachIndex))

    sim_ballisticCoef = sim_mass2./(sim_dragCoefO.*sim_area);
if  numOfChutes == 0
    makeGraph('Line',simRunNum,'Ballistic Coefficient','Time (s)','Ballistic Coefficient (kg/m^2)',sim_time,sim_ballisticCoef)
else
    makeGraph('Line',simRunNum,'Ballistic Coefficient, for Ascent','Time (s)','Ballistic Coefficient (kg/m^2)',sim_time(1:ApoIndex),sim_ballisticCoef(1:ApoIndex))
end

    if  staticStability == true
        makeGraph('Line',simRunNum,'Reynolds Number','Time (s)','Reynolds Number',sim_time,sim_Reynolds)
        
        maxReynoldsIndex = find(sim_Reynolds == max(sim_Reynolds));
        semilogx(sim_Reynolds(1:maxReynoldsIndex),sim_dragCoefO(1:maxReynoldsIndex))
        title('Reynold''s Number vs. Parasite Drag Coefficient, up to max. Reynolds No.')
        xlabel('Reynold''s Number (log scale)')
        ylabel('Zero-Lift Drag Coefficient')
        saveas(gcf,['Reynold''s Number vs. Parasite Drag Coefficient, for Ascent_(',num2str(simRunNum),').jpg'])
        tic; while toc < 1.5; end; close
    end % {static stability}
    
%%% %%%
end % {flight analysis}
    
    
%% Animations [ani15]
%%
%
if  flight_analysis == true
%%%

fprintf('%s\n','Graph production will now pause. Next is an animated trajectory plot',...
    'with points plotted at a rate to match the speed of the rocket,','with a current margin of error of 1%.')
commandwindow
close all

if (max(strcmp((input('Play animated trajectory plot? Yes to play, no to skip. [Y/N]: ','s')),yes))) == true
%%% %%%
%
disp('The flight time at regular points during the simulation will be displayed in the command window.')

    if  numOfChutes > 0
        %skipDescentAnim = max(strcmp((input('Skip animating the descent? [Y/N]: ','s')),yes));
        fprintf('%s\n','When using parachutes, descent time will be long and so the animation',...
        'would slow down considerably. Consequently, descent will not be animated.')
        skipDescentAnim = true;
    end

trajFigure = figure('Name','Trajectory','NumberTitle','off');
set(trajFigure, 'Units', 'normalized', 'Position', [0.01,0.15,0.975,0.75])

if  wind_velocity ~= 0
    axis( [ -1*(max(max([sim_altitude2(:), sim_range2(:)])))/2,...
                max(max([sim_altitude2(:), sim_range2(:)]))/2,...
        0, 1.05*max(max([sim_altitude2(:), sim_range2(:)])) ] )
else
    axis( [ 0,...
                max(max([sim_altitude2(:), sim_range2(:)])),...
        0, 1.05*max(max([sim_altitude2(:), sim_range2(:)])) ] )
end
    axis square

anim_line   = animatedline('LineStyle','none','Marker','.');
burnLine    = animatedline('Color','r','Marker','*');
recoverLine = animatedline('Color','b','Marker','o');

legend('Engine OFF','Engine ON','Chute(s) deployed')
xlabel('Range (Lateral Distance) (m)')
ylabel('Altitude (m)')

clock1 = clock;
tic
while	animCounter <= length(sim_time)
    
            clock2      =  clock;
        if  clock2(6)   <  clock1(6)
            time_past   = (clock2(6)+60) - clock1(6);
        else
            time_past   =  clock2(6)     - clock1(6);
        end
        
        
        tCounter(2) = animCounter;
        if      (rem(sim_time(animCounter),1) == 0) && (tCounter(1) ~= tCounter(2))
                disp(['Simulation time = ',num2str(sim_time(animCounter)),' seconds. Real time elapsed = ',num2str(realTime),' seconds.'])
                tCounter(1) = animCounter;
        end
        
        
        if      animCounter == 1 
                addpoints(anim_line,abs(sim_range2(animCounter)),sim_altitude2(animCounter));
                drawnow
                animCounter = animCounter + 1 ;
                
        elseif  animCounter > 1
            
            if	( time_past <= timeStep ) && ( time_past > (1-errorMargin)*timeStep )
            
                if          ( (failureMode == 4 || failureMode == 6) && (sim_time(animCounter) < motorFailTime) ) ...
                        ||  ( sim_time(animCounter) < burnout(length(burnout)) )
                    
                            addpoints(burnLine,abs(sim_range2(animCounter)),sim_altitude2(animCounter))
                    
                elseif  (numOfChutes > 0) && (sim_time(animCounter) > deployTime(1))
                    
                    if      skipDescentAnim == false
                            addpoints(recoverLine,abs(sim_range2(animCounter)),sim_altitude2(animCounter))
                    elseif  skipDescentAnim == true
                            break
                    end
                    
                else
                        addpoints(anim_line,abs(sim_range2(animCounter)),sim_altitude2(animCounter));
                end
                
                drawnow
                animCounter = animCounter + 1 ;
                realTime    =  realTime + time_past;
                clock1 = clock;
                
            elseif  time_past > 1.00*timeStep
                    disp([  'The system cannot currently handle a ',num2str(errorMargin*100), ...
                            '% error margin on the clock timer. Increasing by 1% and continuing.'  ])
                    errorMargin = errorMargin+0.01;
                    clock1 = clock; 
            else
                continue
            end
        end
end
elapsedTime = toc;

if  skipDescentAnim == true
    hold on
    
    deployArray1 = find(sim_time <= deployTime(1));
    deployArray2 = find(sim_time >= deployTime(1));
    
    for deployCount1 = 1:length(deployArray1)
        for deployCount2 = 1:length(deployArray2)
            if deployArray2(deployCount2) == deployArray1(deployCount1)
                deployIndex = deployArray2(deployCount2);
            end
        end
    end
    
	plot( sim_range2( deployIndex : length(sim_time) ) , ...
          sim_altitude2( deployIndex : length(sim_time) ) , ...
          'Color','b','Marker','o' )
end


saveas(gcf,['Animated_Trajectory_(',num2str(simRunNum),').jpg'])
fprintf('%s\n', ['Simulation time at end of animation = ',num2str(sim_time(animCounter-1)),' seconds'],...
                ['Actual time elapsed whilst displaying trajectory animation: ',num2str(elapsedTime),' seconds.'],...
                ['Therefore, processing limitations resulted in the animation taking about ',num2str(round(elapsedTime/sim_time(animCounter-1),1)),'x longer than expected to play.'])

commandwindow
close all
%
%%% %%%
end % {animated trajectory}

%%
if staticStability == true

realTime    = 0;
errorMargin = 0.01;
tCounter    = zeros(1,2);
fireCounter = 1;

%%
fprintf('%s\n','Next is an animated visualization of the rocket as its flight path angle changes,',...
    'once more with points plotted at a rate to match the speed of the rocket,','and with a current margin of error of 1%.',...
    'This animation does not support the separation of multiple stages.')
if (max(strcmp((input('Play rocket stability/control visualisation? Yes to play, no to skip. [Y/N]: ','s')),yes))) == true
%%% %%%
%
fprintf('%s\n','The flight time at regular points during the simulation will be displayed in the command window,',...
        'in addition to the altitude at that time.')

%fprintf('%s\n','If using parachutes, descent time may be long and so the animation may slow down.',...
%        'You can skip the rest of the visualization after a short delay following parachute deployment.')

    if  numOfChutes > 0
        %skipDescentAnim = max(strcmp((input('Skip animating most of the descent? [Y/N]: ','s')),yes));
        skipDescentAnim = true;
        fprintf('%s\n','The animation is not representative of the reality when descending with parachutes,',...
        'so will end near the time of chute deployment.')
    end
    
clock1 = clock;
tic
while	fpaCounter <= length(sim_time)

            clock2      =  clock;
        if  clock2(6)   <  clock1(6)
            time_past   = (clock2(6)+60) - clock1(6);
        else
            time_past   =  clock2(6)     - clock1(6);
        end
        
        
        tCounter(2) = fpaCounter;
        if      (rem(sim_time(fpaCounter),1) == 0) && (tCounter(1) ~= tCounter(2))
                disp(['Simulation time = ',num2str(sim_time(fpaCounter)),' seconds.',' Altitude (of nose tip) at this time = ',num2str(sim_altitude2(fpaCounter)),' metres.']) % Real time elapsed = ',num2str(realTime),' seconds.'])
                if  ~(numOfChutes > 0 && sim_time(fpaCounter) > deployTime(1)) == true
                    disp(['Stability at this time = ',num2str(sim_stability(fpaCounter)),' calibers.'])
                end
                tCounter(1) = fpaCounter;
        end
        

            if	( time_past <= timeStep ) && ( time_past > (1-errorMargin)*timeStep )
            
                %
                if  (numOfChutes > 0) && (sim_time(fpaCounter) > deployTime(1))
                    
                    TransformMatrix = [ cos(sim_angle2RAD(fpaCounter)+pi),-sin(sim_angle2RAD(fpaCounter)+pi); ...
                                        sin(sim_angle2RAD(fpaCounter)+pi),cos(sim_angle2RAD(fpaCounter)+pi) ];
                                    
                    xChute = linspace(-chuteDiameter(1)/2,chuteDiameter(1)/2,100);
                    yChute = sqrt((chuteDiameter(1)/2)^2 - (xChute).^2);
                    TxyChute{1} = TransformMatrix * [(xChute(end/2:end)+sim_COM(fpaCounter));yChute(end/2:end)];
                    TxyChute{2} = TransformMatrix * [(xChute(end/2:end)+sim_COM(fpaCounter));-yChute(end/2:end)];
                    TxyCords{1} = TransformMatrix * [(xChute(1:end/2)+sim_COM(fpaCounter));fliplr(xChute(1:end/2))];
                    TxyCords{2} = TransformMatrix * [(xChute(1:end/2)+sim_COM(fpaCounter));-fliplr(xChute(1:end/2))];
                    
                    plot(TxyChute{1}(1,:),TxyChute{1}(2,:),'-b',TxyChute{2}(1,:),TxyChute{2}(2,:),'-b',...
                         TxyCords{1}(1,:),TxyCords{1}(2,:),'-.b',TxyCords{2}(1,:),TxyCords{2}(2,:),'-.b'); hold on
                    
                else
                    TransformMatrix = [ cos(sim_angle2RAD(fpaCounter)),-sin(sim_angle2RAD(fpaCounter)); ...
                                        sin(sim_angle2RAD(fpaCounter)),cos(sim_angle2RAD(fpaCounter)) ];
                end
                %
                %{
                    TransformMatrix = [ cos(sim_angle2RAD(fpaCounter)),-sin(sim_angle2RAD(fpaCounter)); ...
                                        sin(sim_angle2RAD(fpaCounter)),cos(sim_angle2RAD(fpaCounter)) ];
                                    
                if  (numOfChutes > 0) && (sim_time(fpaCounter) > deployTime(1))
                                    
                    xChute = linspace(-chuteDiameter(1)/2,chuteDiameter(1)/2,100);
                    yChute = sqrt((chuteDiameter(1)/2)^2 - (xChute).^2);
                    TxyChute{1} = TransformMatrix * [(xChute(end/2:end)+sim_COM(fpaCounter));yChute(end/2:end)];
                    TxyChute{2} = TransformMatrix * [(xChute(end/2:end)+sim_COM(fpaCounter));-yChute(end/2:end)];
                    TxyCords{1} = TransformMatrix * [(xChute(1:end/2)+sim_COM(fpaCounter));fliplr(xChute(1:end/2))];
                    TxyCords{2} = TransformMatrix * [(xChute(1:end/2)+sim_COM(fpaCounter));-fliplr(xChute(1:end/2))];
                    
                    plot(TxyChute{1}(1,:),TxyChute{1}(2,:),'-b',TxyChute{2}(1,:),TxyChute{2}(2,:),'-b',...
                         TxyCords{1}(1,:),TxyCords{1}(2,:),'-.b',TxyCords{2}(1,:),TxyCords{2}(2,:),'-.b'); hold on
                end
                %}                
                    
                
                Ty_i            = TransformMatrix * [fliplr(x_iArray)-(R_length-sim_COM(fpaCounter));y_i];
                
                Ty_i_projected	= TransformMatrix * [fliplr(x_iArray)-(R_length-sim_COM(fpaCounter));-y_i_projected];
                
                TxyEnd          = TransformMatrix * [0*xEnd-(R_length-sim_COM(fpaCounter));yEnd];
                
                plot(Ty_i(1,:),Ty_i(2,:),'-k',...
                    Ty_i_projected(1,:),Ty_i_projected(2,:),'-k',...
                    TxyEnd(1,:),TxyEnd(2,:),'-k')
                
                
                if  sim_time(fpaCounter) < burnout(length(burnout))
                   
                    xEnd_firing = [-R_length,0,0,-R_length,-R_length];
                    plumeSpreadFctr = sqrt(R_length/bodyAvgDiameter);
                    fireCounter=-1*fireCounter;
                    
                    if      fireCounter < 0
                            yEnd_firing = [plumeSpreadFctr*y_i(end),y_i(end),-y_i_projected(end),plumeSpreadFctr*-y_i_projected(end),plumeSpreadFctr*y_i(end)];
                    elseif  fireCounter > 0
                            yEnd_firing = [   y_i(end),y_i(end),-y_i_projected(end),   -y_i_projected(end),   y_i(end)];
                    end
                    

                    TxyEnd_firing = TransformMatrix * [ xEnd_firing-(R_length-sim_COM(fpaCounter)) ; yEnd_firing ];
                    distToOriginX = sum(TxyEnd_firing(1,2:3))/2;
                    distToOriginY = sum(TxyEnd_firing(2,2:3))/2;
                    for gimbalCount = 1:5
                        TxyEnd_firing(:,gimbalCount) = TxyEnd_firing(:,gimbalCount) - [distToOriginX;distToOriginY]; %translate to origin
                        TxyEnd_firing(:,gimbalCount) = [ cos(sim_thrustAngle(fpaCounter)),-sin(sim_thrustAngle(fpaCounter));...
                                                          sin(sim_thrustAngle(fpaCounter)),cos(sim_thrustAngle(fpaCounter)) ]...
                                                     * TxyEnd_firing(:,gimbalCount); % rotate about origin
                        TxyEnd_firing(:,gimbalCount) = TxyEnd_firing(:,gimbalCount) + [distToOriginX;distToOriginY]; % translate back
                    end
                    
                    
                    hold on; fill(TxyEnd_firing(1,:),TxyEnd_firing(2,:),[1,0.3067,0]);
                end
                
                
                hold on;            plot(0,0,'Marker','+','Color','b')
                hold on;            plot(0,0,'Marker','x','Color','b')
                hold on; COMpoint = plot(0,0,'Marker','o','Color','b',...
                                       'DisplayName','Centre of Mass');
                                         
                Transformed_COP = TransformMatrix*[(sim_COM(fpaCounter)-sim_COP(fpaCounter));0];
                                
                hold on;            plot(Transformed_COP(1),Transformed_COP(2),'Marker','.','Color','r');
                hold on; COPpoint = plot(Transformed_COP(1),Transformed_COP(2),'Marker','o','Color','r',...
                                           'DisplayName','Centre of Pressure');
                
                hold on; plot([R_length,-R_length],[0,0],'--k')
                hold on; plot([0,0],[R_length,-R_length],'--k')
                
 
                hold off
                set(gcf,'Units','normalized','Position',[0.01,0.175,0.975,0.725]);

                xticklabel = linspace(  sim_range1(fpaCounter)-(sim_COM(fpaCounter)+centreOfPressure)/2,...
                                        sim_range1(fpaCounter)+(sim_COM(fpaCounter)+centreOfPressure)/2, 5);
                                
                yticklabel = linspace(  sim_altitude1(fpaCounter)-(sim_COM(fpaCounter)+centreOfPressure)/2,...
                                        sim_altitude1(fpaCounter)+(sim_COM(fpaCounter)+centreOfPressure)/2, 5);

                xyTick = [-(sim_COM(fpaCounter)+centreOfPressure)/2,-(sim_COM(fpaCounter)+centreOfPressure)/4, 0,...
                           (sim_COM(fpaCounter)+centreOfPressure)/4, (sim_COM(fpaCounter)+centreOfPressure)/2];
                       
                set(gca,'XTick',xyTick)
                set(gca,'XTickLabel',xticklabel);
                set(gca,'XTickLabelRotation',45)
                xlabel('Range (m)')
                    
                set(gca,'YTick',xyTick)
                set(gca,'YTickLabel',yticklabel)
                ylabel('Altitude (m)')

                axis( [ -(centreOfMass+centreOfPressure)/2, (centreOfMass+centreOfPressure)/2,...
                        -(centreOfMass+centreOfPressure)/2, (centreOfMass+centreOfPressure)/2 ] )
                axis square
                legend([COMpoint,COPpoint])
                
                
                drawnow
                fpaCounter = fpaCounter + 1 ;
                realTime    =  realTime + time_past;
                clock1 = clock;
                
                frame = getframe(1);
                frameImage{fpaCounter} = frame2im(frame);
                
                
                if (numOfChutes > 0) && (skipDescentAnim == true) && ...
                   ( sim_time(fpaCounter-1) >= (sim_time(ApoIndex) + 0.015*sim_time(end)) )
                    break
                end % ends visualization after short time delay after parachute deployment
                
                
            elseif  time_past > 1.00*timeStep
                    disp([  'The system cannot currently handle a ',num2str(errorMargin*100), ...
                            '% error margin on the clock timer. Increasing by 1% and continuing.'  ])
                    errorMargin = errorMargin+0.01;
                    clock1 = clock; 
            else
                continue
            end
end
elapsedTime = toc;
fprintf('%s\n', ['Simulation time at end of visualization = ',num2str(sim_time(fpaCounter-1)),' seconds'],...
                ['Actual time elapsed whilst displaying stability visualization: ',num2str(elapsedTime),' seconds.'],...
                ['Therefore, processing limitations resulted in the animation taking about ',num2str(round(elapsedTime/sim_time(fpaCounter-1),1)),'x longer than expected to play.'])

            
disp('Please wait... (saving stability visualization as GIF sequential image file)')
filename = [rocketName,'_animation_(',num2str(simRunNum),').gif'];
    for gifCount = 2:fpaCounter
        [indexedImage,colourMap] = rgb2ind(frameImage{gifCount},256);
        if      gifCount == 2
                imwrite(indexedImage,colourMap,filename,'gif','LoopCount',Inf,...
                    'DelayTime',1);
        elseif  gifCount == fpaCounter
                imwrite(indexedImage,colourMap,filename,'gif','WriteMode','append',...
                    'DelayTime',1);
        else
                imwrite(indexedImage,colourMap,filename,'gif','WriteMode','append',...
                    'DelayTime',timeStep);
        end
    end % saves stability visualization as a GIF file
clear frameImage
            
%%% %%%
end % {rocket visualization}

%%% %%%
end % {static stability}

%%%
end % {flight analysis}

%%
%%%
commandwindow
reRunSim = max(strcmp((input('Re-run the simulation with different inputs? [Y/N]: ','s')),yes));
close all
diary off
end % un-comment this line when integrating scripts %%%DEV:NB%%%
warning('on','MATLAB:hg:AutoSoftwareOpenGL')
%%%
%
%%