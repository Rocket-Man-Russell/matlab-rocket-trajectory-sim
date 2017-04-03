function e = oswaldFactor_mod(AR,sweep,method,CD0,df_over_b,u)
% This is a modified version of Sky Sartorius' original function
% for estimating the Oswald efficiency factor, simply shortened by 
% removing some comments (version history etc.) - F. Russell


% Set defaults.
if nargin < 6
    u = 1;
end
if nargin < 5 || isempty(df_over_b)
    df_over_b = 0;
end
if nargin >= 4 && isempty(method)
    method = 'shevell';
end
if nargin < 3 || isempty(method)
    method = 'raymer';
end
if nargin < 2 || isempty(sweep)
    sweep = 0;
end


if isnumeric(method) || strcmpi(method,'average')
    if strcmpi(method,'average')
        avgweight = 1/2;
    elseif (all(method(:) >= 0)) && (all(method(:) <= 1))
        avgweight = method;
        method = 'average';
    else
        error('All elements of a numeric METHOD must be between 0 and 1.')
    end
end

if ~strcmpi('shevell',method) %need to calculate raymer at least
    if any(AR(:)>10)
        warning('Raymer method not recommended for high aspect ratios.')
    end
    % with sweep
    e_ray_sweep = 4.61*(1-0.045*AR.^.68).*(cos(sweep)).^0.15-3.1;
    % no sweep
    e_ray_no_sweep = 1.78*(1-0.045*AR.^.68)-0.64;
    
    weight = abs(sweep)/(30*pi/180);
    e_ray = weight.*e_ray_sweep+(1-weight).*e_ray_no_sweep;
    
    hisweepind = sweep > 30*pi/180;
    e_ray(hisweepind) = e_ray_sweep(hisweepind);
    if strcmpi('raymer',method)
        e = e_ray;
        return
    end
end

%calculate shevell
s = 1 - 1.9316 * df_over_b.^2;
k = 0.375 * CD0;
sweep_corr = 1  -  0.02 * AR.^0.7 .* sweep.^2.2;

e_shev = sweep_corr./(pi*AR.*k + 1./(u.*s));

switch lower(method)
    case 'shevell'
        e = e_shev;
    case 'pessimistic'
        e = min(e_shev,e_ray);
    case 'optimistic'
        e = max(e_shev,e_ray);
    case 'average'
        e = avgweight.*e_shev+(1-avgweight).*e_ray;
    otherwise
        error('Bad method string.')
end

end

%{
   e = OSWALDFACTOR(AR,sweep,method,C_D0,df_b,u)

   AR     - Wing effective aspect ratio (simply aspect ratio for single
            wing w/o winglets).
   sweep  - Wing sweep in radians. Default sweep = 0.
   method - Choose between different methods. Options listed below.
            Default method = 'Shevell' if C_D0 is provided, otherwise
            default method = 'Raymer'.
       'Raymer'      - Raymer's "Aicraft Design: A Conceptual Approach"
                       equations 12.48/49 in fifth edition.
       'Shevell'     - Adapted from Shevell's "Fundamentals of Flight."
       'Pessimistic' - The worst of the two.
       'Optimistic'  - The best of the two.
       numeric 0-1   - method×Shevell + (1-method)×Raymer.
       'Average'     - Same as numeric method = 0.5.
   C_D0   - Parasite drag coefficient (drag independent of lift).
   df_b   - Ratio of fuselage diameter (or width) to wing span. Default 
            df_b = 0.
   u      - Planform efficiency — usually 0.98 < u < 1. Default u = 1.
%}

%{
    Copyright 2012-2013 Sky Sartorius
    Author contact: mathworks.com/matlabcentral/fileexchange/authors/101715

References:
    Shevell: "Fundamentals of Flight" 2ed. 1989
             See also http://adg.stanford.edu/aa241/drag/induceddrag.html
             for some info on the Shevell theory and assumptions leading to
             coefficients.
    Raymer:  "Aircraft Design: A Conceptual Approach" 5ed. 2012
             See Raymers leading-edge-suction method (section 12.6.2 in 5th
             ed.) for an alternate method for high-AR cases.
%}
