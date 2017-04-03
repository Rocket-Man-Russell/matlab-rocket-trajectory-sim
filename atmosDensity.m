function [sim_atmosDensity] = atmosDensity(sim_altitude1)
    if  sim_altitude1 < 30000
        sim_atmosDensity = (-0.0000000000000395*(sim_altitude1^3))+(0.0000000037*(sim_altitude1^2))-(0.0001151*sim_altitude1)+1.225;
    else
        sim_atmosDensity = 1.0776*(2.71828^(-0.0001387*sim_altitude1));
    end
end
% sim_atmosDensity(simCounter) = atmosDensity(sim_altitude1(simCounter))