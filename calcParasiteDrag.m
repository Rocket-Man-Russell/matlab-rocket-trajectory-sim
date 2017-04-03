function [ParasiteDrag,critRatioError] = calcParasiteDrag(sim_altitude1,sim_soundSpeed,sim_machNum,...
                                           R_length,skinRoughness,bodyMaxDiameter,...
                                           surfaceArea,fin_rootChord,fin_tipChord,...
                                           finThickness,fin_span,numOfFins,...
                                           bodybaseDiameter,nose_length,...
                                           R_effLength,R_aftLength,finSetNUM,...
          numOfLaunchLugs,protub_length,protub_dist,protub_crossArea,protub_surfArea)
                                       % Assumes constant fin thickness
                                       % (square fin cross-section)
%% 1 Friction Drag
%
%%   1.1 Body Friction Drag
%
%       Speed of Sound
%{
if      sim_altitude1   <= (37000/3.28084) % 1 metre = 3.28084 feet
        sim_soundSpeed  = ( -0.004*(sim_altitude1(simCounter)*3.28084) + 1116.45 ) / 3.28084;
        
elseif  sim_altitude1   <= (64000/3.28084)
        sim_soundSpeed  = ( 968.08 ) / 3.28084;
    
elseif  sim_altitude1   > (64000/3.28084)
        sim_soundSpeed  = ( 0.0007*(sim_altitude1(simCounter)*3.28084) + 924.99 ) / 3.28084;
end
%}

%       Kinematic Viscosity
if      sim_altitude1       <= (15000/3.28084)
        sim_kineViscosity   = 0.000157*exp(0.00002503*(sim_altitude1*3.28084));
        
elseif  sim_altitude1       <= (30000/3.28084)
        sim_kineViscosity   = 0.000157*exp(0.00002760*(sim_altitude1*3.28084) - 0.03417);
        
elseif  sim_altitude1       > (30000/3.28084)
        sim_kineViscosity   = 0.000157*exp(0.00004664*(sim_altitude1*3.28084) - 0.6882);
end

%       Compressible Reynolds Number
        body_Reynolds       = ( (sim_soundSpeed*sim_machNum*R_length)/(12*sim_kineViscosity) ) * ...
                              (1 + 0.0283*sim_machNum - 0.043*sim_machNum^2 + 0.2107*sim_machNum^3 ...
                                 - 0.03829*sim_machNum^4 + 0.002709*sim_machNum^5);
                             
%       Incompressible skin friction coefficient
    %if body_Reynolds ~= 0
        body_skinFrictionI   = 0.037036*body_Reynolds^-0.155079;
    %else
        %body_skinFrictionI = 0;
    %end
        
%       Compressible skin friction coefficient
        body_skinFrictionII    = ( body_skinFrictionI ) * ...
                              (1 + 0.00798*sim_machNum - 0.1813*sim_machNum^2 + 0.0632*sim_machNum^3 ...
                                 - 0.00933*sim_machNum^4 + 0.000549*sim_machNum^5);
                             
%       Incompressible skin friction coefficient with roughness
        body_skinFrictionIII   = 1/( 1.89 * 1.62*log10(R_length/skinRoughness) )^2.5;

%       Compressible skin friction coefficient with roughness
        body_skinFrictionIV  = body_skinFrictionIII / (1 + 0.2044*sim_machNum^2);
        
%       Final skin friction coefficient
        if      body_skinFrictionII >= body_skinFrictionIV
                body_skinFrictionV   = body_skinFrictionII;
                
        elseif  body_skinFrictionII  < body_skinFrictionIV
                body_skinFrictionV   = body_skinFrictionIV;
        end
        
%       Body coefficient of drag due to friction
        body_skinFrictionTOTAL = body_skinFrictionV * ( 1 + 60/((R_length/bodyMaxDiameter)^3) + 0.0025*(R_length/bodyMaxDiameter) ) * ( (4*surfaceArea)/(pi*bodyMaxDiameter^2) );


%%  1.2 Fin Friction Drag

for finsetCount = 1:finSetNUM
%       Compressible Reynolds Number
        fin_ReynoldsI       = ( (sim_soundSpeed*sim_machNum*fin_rootChord(finsetCount))/(12*sim_kineViscosity) ) * ...
                              (1 + 0.0283*sim_machNum - 0.043*sim_machNum^2 + 0.2107*sim_machNum^3 ...
                                 - 0.03829*sim_machNum^4 + 0.002709*sim_machNum^5);
                             
%       Incompressible skin friction coefficient
    %if fin_ReynoldsI ~= 0
        fin_skinFrictionI   = 0.037036*fin_ReynoldsI^-0.155079;
    %else
        %fin_skinFrictionI = 0;
    %end
        
%       Compressible skin friction coefficient
        fin_skinFrictionII  = ( fin_skinFrictionI ) * ...
                              (1 + 0.00798*sim_machNum - 0.1813*sim_machNum^2 + 0.0632*sim_machNum^3 ...
                                 - 0.00933*sim_machNum^4 + 0.000549*sim_machNum^5);
                             
%       Incompressible skin friction coefficient with roughness
        fin_skinFrictionIII = 1/( 1.89 * 1.62*log10(fin_rootChord(finsetCount)/skinRoughness) )^2.5;

%       Compressible skin friction coefficient with roughness
        fin_skinFrictionIV  = fin_skinFrictionIII / (1 + 0.2044*sim_machNum^2);
        
%       Final skin friction coefficient
        if      fin_skinFrictionII >= fin_skinFrictionIV
                fin_skinFrictionV   = fin_skinFrictionII;
                
        elseif  fin_skinFrictionII  < fin_skinFrictionIV
                fin_skinFrictionV   = fin_skinFrictionIV;
        end
        
%       Incompressible Reynolds Number
        fin_ReynoldsII      = (sim_soundSpeed*sim_machNum*fin_rootChord(finsetCount))/(12*sim_kineViscosity);
        %if fin_ReynoldsII < 0, fin_ReynoldsII = 1; end %%%DEV:NB%%% bug fix
        
%       Ratio of fin tip chord to root chord
        fin_taperRatio = fin_tipChord(finsetCount)/fin_rootChord(finsetCount);
        %if fin_taperRatio < 0, fin_taperRatio = 1; end %%%DEV:NB%%% bug fix
        
%       Average flat plate skin friction coefficient for each fin panel
        if  fin_taperRatio == 0
            fin_skinFrictionVI = fin_skinFrictionV * ( 1 + (0.5646)/log10(fin_ReynoldsII) );
        else
            fin_skinFrictionVI = fin_skinFrictionV * ( ((log10(fin_ReynoldsII))^2.6)/((fin_taperRatio^2)-1) ) * ...
                               (   ( (fin_taperRatio^2)/(log10(fin_ReynoldsII*fin_taperRatio))^2.6 ) ...
                                   - 1/(log10(fin_ReynoldsII))^2.6 ...
                                   + 0.5646 * (  ( (fin_taperRatio^2)/(log10(fin_ReynoldsII*fin_taperRatio))^3.6 ) ...
                                                 - 1/(log10(fin_ReynoldsII))^3.6  )   );

        end
        
%       Coefficient of friction drag for all fins
        X = 0; % Distance from fin leading edge to maximum thickness
        fin_surfaceArea = (fin_span(finsetCount)/2)*(fin_rootChord(finsetCount)+fin_tipChord(finsetCount)); % Total wetted area of each fin
        fin_skinFrictionTOTAL{finsetCount} = fin_skinFrictionVI * ( 1 + 60*(finThickness(finsetCount)/fin_rootChord(finsetCount))^4 + 0.8*(1+5*X^2)*(finThickness(finsetCount)/fin_rootChord(finsetCount)) ) ...
                                                   * ((4*numOfFins(finsetCount)*fin_surfaceArea)/(pi*bodyMaxDiameter^2));
end
        

%%  1.3 Protuberance Friction Drag

    ProtuberanceDrag = 0;
%
        %{
        Need to define:
            Protuberance length,
                protub_length;
            Distance from rocket nose to front edge of 
                protuberance, protub_dist;
            Maximum cross-section area of protuberance,
                protub_crossArea;
            Wetted surface area of protuberance,
                protub_surfArea;
            number of protuberances (launch lugs)
                numOfLaunchLugs
        %}
        
    if numOfLaunchLugs ~= 0
    %%%
    for protubCount = 1:numOfLaunchLugs
    %%% %%%
    

        %       Compressible Reynolds Number
        protub_ReynoldsI       = ( (sim_soundSpeed*sim_machNum*protub_length{protubCount})/(12*sim_kineViscosity) ) * ...
                              (1 + 0.0283*sim_machNum - 0.043*sim_machNum^2 + 0.2107*sim_machNum^3 ...
                                 - 0.03829*sim_machNum^4 + 0.002709*sim_machNum^5);
                             
%       Incompressible skin friction coefficient
    %if protub_ReynoldsI ~= 0
        protub_skinFrictionI   = 0.037036*protub_ReynoldsI^-0.155079;
    %else
        %protub_skinFrictionI = 0;
    %end
        
%       Compressible skin friction coefficient
        protub_skinFrictionII  = ( protub_skinFrictionI ) * ...
                              (1 + 0.00798*sim_machNum - 0.1813*sim_machNum^2 + 0.0632*sim_machNum^3 ...
                                 - 0.00933*sim_machNum^4 + 0.000549*sim_machNum^5);
                             
%       Incompressible skin friction coefficient with roughness
        protub_skinFrictionIII = 1/( 1.89 * 1.62*log10(protub_length{protubCount}/skinRoughness) )^2.5;

%       Compressible skin friction coefficient with roughness
        protub_skinFrictionIV  = protub_skinFrictionIII / (1 + 0.2044*sim_machNum^2);
        
%       Final skin friction coefficient
        if      protub_skinFrictionII >= protub_skinFrictionIV
                protub_skinFrictionV   = protub_skinFrictionII;
                
        elseif  protub_skinFrictionII  < protub_skinFrictionIV
                protub_skinFrictionV   = protub_skinFrictionIV;
        end
        
%       Friction coefficient of protuberance
    % if (protub_dist{protubCount}/protub_length{protubCount}) ~= 0
        protub_FrictionCoef = 0.8151*protub_skinFrictionV*(protub_dist{protubCount}/protub_length{protubCount})^-0.1243;
    % else
    %   protub_FrictionCoef = 0;
    % end
        
%       Drag coefficient of protuberance due to friction
        protub_DragCoef = protub_FrictionCoef * ...
            (1+1.798*(sqrt(protub_crossArea{protubCount})/protub_length{protubCount})^(3/2)) * ...
            ((4*protub_surfArea{protubCount})/(pi*bodyMaxDiameter^2));

        ProtuberanceDrag = ProtuberanceDrag + protub_DragCoef;
    
    %%% %%%
    end
    %%%
    end
%

        
%%  1.4 Drag due to Excrescencies
        if      sim_machNum < 0.78
                Ke = 0.00038;
        elseif  sim_machNum < 1.04
                Ke = -0.4501*sim_machNum^4 + 1.5954*sim_machNum^3 - 2.1062*sim_machNum^2 + 1.2288*sim_machNum - 0.26717;
        else
                Ke = 0.0002*sim_machNum^2 - 0.0012*sim_machNum + 0.0018;
        end
        
        Excrescencies_Drag = Ke * (4*surfaceArea)/(pi*bodyMaxDiameter^2); % scratches, holes etc.

        
%%  1.5 Total Friction & Interference Drag Coefficient

        Kf = 1.04; % Mutual interference factor of fins and launch lug with body
        FrictionDragTOTAL = body_skinFrictionTOTAL  + Kf*sum(cell2mat(fin_skinFrictionTOTAL)) ...
                          + Kf*ProtuberanceDrag     + Excrescencies_Drag;


%% 2 Base Drag Coefficient
%
%%  2.1 Base Drag Coefficient for M < 0.6

        Kb = 0.0274*atan((R_aftLength/bodyMaxDiameter)+0.0116);
        if  R_aftLength ~= 0
            n = 3.6542*(R_aftLength/bodyMaxDiameter)^-0.2733;
        else
            n = 0;
        end
        baseDrag = Kb * ( ((bodybaseDiameter/bodyMaxDiameter)^n)/sqrt(FrictionDragTOTAL) );


%%  2.1 Base Drag Coefficient for M > 0.6
    if  sim_machNum >= 0.6
        if      sim_machNum <= 1
                fb = 1 + 215.8*(sim_machNum-0.6)^6;
        elseif  sim_machNum <= 2
                fb = 2.0881*(sim_machNum-1)^3 - 3.7938*(sim_machNum-1)^2 + 1.4618*(sim_machNum-1) + 1.883917;
        elseif  sim_machNum > 2
                fb = 0.297*(sim_machNum-2)^3 - 0.7937*(sim_machNum-2)^2 - 0.1115*(sim_machNum-2) + 1.64006;
        end
        baseDrag = baseDrag*fb;
    end


%% 3 Transonic Wave Drag Coefficient

%       ratio of nose length to effective rocket length:
        criticalRatio = nose_length/R_effLength;

	if  criticalRatio  > 0.6
        critRatioError = true;
    else
        critRatioError = false;
	end

%       Transonic drag divergence Mach Number
        MD = -0.0156*(nose_length/bodyMaxDiameter)^2 + 0.136*(nose_length/bodyMaxDiameter) + 0.6817;

%       Final Mach Number of Transonic Region
        if  criticalRatio < 0.2
            a = 2.4;
            b = -1.05;
        elseif criticalRatio < 0.6
            a = -321.94*criticalRatio^2 + 264.07*criticalRatio - 36.348;
            b = 19.634*criticalRatio^2 - 18.369*criticalRatio + 1.7434;
        else
            [a,b] = deal(0);
        end
            
        MF = a*(nose_length/bodyMaxDiameter)^b + 1.0275;

%       Maximum drag rise over transonic region
        c = 50.676*(nose_length/R_length)^2 - 51.734*(nose_length/R_length) + 15.642;
        g = -2.2538*(nose_length/R_length) + 1.3108*(nose_length/R_length) - 1.7344;
        
        if      (R_effLength/bodyMaxDiameter) >= 6
                deltaCDtransMAX = c*(R_effLength/bodyMaxDiameter)^g;
        elseif  (R_effLength/bodyMaxDiameter) < 6
                deltaCDtransMAX = c*6^g;
        end

%       Transonic drag rise for given Mach Number ‘M’
        x = (sim_machNum - MD)/(MF - MD);
        F = -8.3474*x^5 + 24.543*x^4 - 24.946*x^3 + 8.6321*x^2 + 1.1195*x;
        
        if      (MD <= sim_machNum) && (sim_machNum <= MF)
                deltaCDtrans = deltaCDtransMAX * F;
        else
                deltaCDtrans = 0;
        end


%% 4 Supersonic Wave Drag Coefficient

%       Supersonic drag rise for given Mach Number ‘M’
        if      sim_machNum >= MF
                deltaCDsuper = deltaCDtransMAX;
        elseif  sim_machNum < MF
                deltaCDsuper = 0;
        end

        
%% 5 TOTAL DRAG COEFFICIENT
    
    ParasiteDrag = FrictionDragTOTAL + baseDrag + deltaCDtrans + deltaCDsuper; 
    % This is the zero-lift drag coefficient for the rocket.

%%
end
%{
source: web.aeromech.usyd.edu.au/AERO2705/Resources/Research/Drag_Coefficient_Prediction.pdf
%}