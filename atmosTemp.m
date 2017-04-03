function [sim_atmosTemp] = atmosTemp(sim_altitude1)
    if      sim_altitude1 < 11000
            sim_atmosTemp = (-0.0064895*sim_altitude1) + 288.19;
    elseif  sim_altitude1 < 20000
            sim_atmosTemp = 217;
    elseif  sim_altitude1 < 50000
            sim_atmosTemp = (-0.0000000000025*(sim_altitude1^3))+(0.0000002942*(sim_altitude1^2))-(0.0089406*sim_altitude1) + 299.9;
    elseif  sim_altitude1 < 85000
            sim_atmosTemp = (-0.0024429*sim_altitude1) + 391.64;
    else
            sim_atmosTemp = (0.0008667*sim_altitude1) + 113.33;
    end
end
% sim_atmosTemp(simCounter) = atmosTemp(sim_altitude1(simCounter))