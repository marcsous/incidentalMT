% code to optimize TR for a given T1 to null and T1 of interest
clear

T1null = 4000; % T1 to null
T1int = 790;   % T1 of interest
TElast = 100;  % time of last echo (fast spin echo)

% avoid singularity at zero
initial_guess = 3.57*T1int;

% find the minimum of -efficiency
myfun = @(TR)-efficiency(TR,T1null,T1int,TElast);
TRopt = fminsearch(myfun,initial_guess,optimset('display','off'));

% display results
[E TI M C] = efficiency(TRopt,T1null,T1int,TElast);
fprintf('TRopt=%f TI=%f M=%+f E=%f C=%f\n',TRopt,TI,M,E,C);

%% function to calculate efficiency/contrast
function [E TI M C] = efficiency(TR,T1null,T1int,TElast)

% TI needed for T1null
TI = T1null * (log(2) - log(1+exp(-(TR-TElast)/T1null)));

% magnetization at T1int
M = 1 - 2*exp(-TI/T1int) + exp(-(TR-TElast)/T1int);

% efficiency normalized to spin echo
E = abs(M) / TR^(1/2);
E = E / (0.64 / T1int^(1/2));

% contrast at T1int normalized to spin echo
C = abs(exp(-(TR-TElast)/T1int)*(TR-TElast) - 2*exp(-TI/T1int)*TI) / T1int^2;
C = C / (2*exp(-2) / T1int);

end
