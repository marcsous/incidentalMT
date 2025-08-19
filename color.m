function c = color(k)
%
% Return the default MATLAB color codes
validateattributes(k,{'numeric'},{'scalar','integer'});

switch mod(k-1,7)
    case 0; c = [0.0000 0.4470 0.7410];
    case 1; c = [0.8500 0.3250 0.0980]; 
    case 2; c = [0.9290 0.6940 0.1250]; 
    case 3; c = [0.4940 0.1840 0.5560];
    case 4; c = [0.4660 0.6740 0.1880];
    case 5; c = [0.3010 0.7450 0.9330];
    case 6; c = [0.6350 0.0780 0.1840];
end
