function g = SuperLorentz_LineShape(freq,T2b)
% Returns the SuperLorentzian LineShape (spherical averaging)
%  freq in Hz: frequency offset with respect to the center of the distribution
%  T2b in s: T2 of the underlying distribution

g = zeros(size(freq));

for k = 1:numel(freq)
    omega = 2*pi*freq(k); % (rad/s)
    g(k) = (1/sqrt(2*pi))*quadgk(@(theta)SphericalLineShape_function(omega,T2b,theta),0,pi/2);
end

%% Function for Spherical lineshape integration
function dg = SphericalLineShape_function(omega,T2b,theta)

% Include neighbors' contribution to remove singularity at the magic angle.
% Ref:  Pampel et al. NeuroImage 114 (2015) 136–146

R2res = 31.4; % (s^-1)
R2b = abs(3*cos(theta).^2-1)./(2*T2b);
R2tot = hypot(R2res,R2b);
dg = (1./R2tot).*exp(-(1/2)*omega^2./R2tot.^2).*sin(theta);
