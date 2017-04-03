function [sim_liftCoef] = liftCoef(sim_angleOfAttack)
	if      (abs(sim_angleOfAttack*(180/pi)) <= 15)
            sim_liftCoef = 1.1;
            
	elseif  (abs(sim_angleOfAttack*(180/pi)) <= 20)
            sim_liftCoef = (-1.824*pi*sim_angleOfAttack + 2.583);
            
    elseif  (abs(sim_angleOfAttack*(180/pi)) <= 36)
            sim_liftCoef = (0.472*pi*sim_angleOfAttack + 0.069);
            
	elseif  (abs(sim_angleOfAttack*(180/pi)) <= 55)
            sim_liftCoef = 1.0;
            
	elseif  (abs(sim_angleOfAttack*(180/pi)) <= 125) % <= 90)
            sim_liftCoef = (-0.474*pi*sim_angleOfAttack + 2.445);
    else
            sim_liftCoef = 0;
            warning('Extremely high angle of attack encountered.')
	end
end
% sim_liftCoef(simCounter) = liftCoef(sim_angleOfAttack(simCounter))