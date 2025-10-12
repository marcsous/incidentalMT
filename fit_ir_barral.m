function [T1 x phi ci95] = fit_ir_barral(TI,data,T1,lambda)
%[T1 x phi ci95] = fit_ir_barral(TI,data,T1,lambda)
%
% Fits inversion recovery data to estimate T1.
% Model: [x(1) + x(2)*exp(-TI/T1)] * exp(i*phi)
%
% Arguments:
%  TI is a vector of inversion times (Nx1)
%  data is an array of complex data points (NxM)
%  T1 (optional) is an initial T1 estimate (scalar)
%  lambda (optional) is regularization term (scalar)
%
% Notes:
% Signed or complex data only. Not magnitude.
% Multiple coupled RHS okay (for multi-coil).
% Lambda penalizes x(1)+x(2)/2 (should be 0).
%
% Ref. Barral (2010) doi.org/10.1002/mrm.22497
%
%% check inputs
if isreal(data) && all(data>=0)
    warning('data should be signed or complex');
end
if length(data)~=length(TI)
    error('data and TI size mis-match')
end
if ~exist('T1','var') || isempty(T1)
    R1 = 1/mean(TI);
else
    R1 = 1/T1;
end
if ~exist('lambda','var') || isempty(lambda)
    lambda = 0;
end

% shape and type consistency
R1 = reshape(double(R1),1,1);
TI = reshape(double(TI),[],1);
lambda = reshape(double(lambda),1,1);
data = reshape(double(data),numel(TI),[]);

% no. RHS vectors
nrhs = size(data,2); 

%% least squares fitting

% optimization
opts = optimset('display','off','TolFun',1e-16);
[R1,~,~,~,~,H] = fminunc(@(R)myfun(R,TI,data,lambda),R1,opts);

% 95% confidence interval
[resnorm f x phi] = myfun(R1,TI,data,lambda);
if isreal(data)
    v = numel(data)-1-2*nrhs; % [T1 A1 A2]
else
    v = 2*numel(data)-1-3*nrhs; % [T1 A1 A2 phi]
end
err = resnorm/v;
cov = inv(H)*err;
stderr = sqrt(diag(cov));
ci95 = 1.96*stderr;

% calculate r2
r2 = 1 - resnorm/sum(sum(abs(data-mean(data)).^2));

% return arguments (T1 = 1/R1)
T1 = 1/R1;
ci95 = ci95/R1^2; 

%% display
plot(TI,real(exp(-i*phi).*data),'o','color',[0 0.447 0.741]);
hold on
if ~isreal(data)
    plot(TI,imag(exp(-i*phi).*data),'o','color',[0.85 0.325 0.098]);    
    plot(TI,imag(exp(-i*phi).*f),'color',[0.85 0.325 0.098]);
end
plot(TI,real(exp(-i*phi).*f),'color',[0 0.447 0.741]);
hold off
title('Inversion recovery'); xlabel('TI'); ylabel('Signal'); grid on

% text box
str{1} = sprintf('T_1 = %.3f ± %.3f',T1,ci95);
pow10 = ceil(log10(norm(x,'fro')));
str{2} = sprintf('x_1 = [%s].10^{%.0f}',num2str(x(1,:)/10^pow10,' %+.3f'),pow10);
str{3} = sprintf('x_2 = [%s].10^{%.0f}',num2str(x(2,:)/10^pow10,' %+.3f'),pow10);
if ~isreal(data)
    str{4} = sprintf('\\phi_0 = [%s] rad',num2str(phi,' %+.3f'));
end
str{end+1} = sprintf('r^2 = %.9f',r2);
text(0.03,0.83,str,'Units','Normalized','FontName','FixedWidth')
fprintf('Fitted T1: %.3f ± %.3f\n',T1,ci95);

%% function: [x(1) + x(2)*exp(-R*t)] * exp(i*phi)
function [resnorm f x phi] = myfun(R,ti,b,lambda)

% multiple coils (nc)
[nti nc] = size(b);

% doi.org/10.1002/mrm.22497
A = [ones(nti,1) exp(-R*ti)];

% penalize x(1)+x(2)/2
L = [1 1/2]*lambda;

% doi.org/10.1016/j.laa.2010.07.011
Ab = A'*b; M = A'*A+L'*L;
if isreal(Ab)
    % linear
    x = pinv(M)*Ab;
    phi = 0; % unused
else
    % phase constrained
    invM = pinv(real(M)); 
    phi = angle(sum(Ab.*(invM*Ab),1))/2;
    x = invM*real(Ab.*exp(-i*phi));
end

% function and residual norm
f = A*(x.*exp(i*phi));
resnorm = norm(f-b,'fro')^2;

% cosmetic change of sign to make x(2) negative
if nargout>2
    for c = 1:nc
        if x(2,c)>0
            x(:,c) = -x(:,c);
            if phi(c)>0
                phi(c) = phi(c) - pi;
            else
                phi(c) = phi(c) + pi;
            end
        end
    end
end
