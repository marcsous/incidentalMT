% Bloch simulation of continuous wave saturation of Mb
%
% Mf is free pool magnetization
% Mb is bound pool magnetization
%
%  R1f     kf    R1b+RFb
% <--- Mf <==> Mb --->
%          kb
%
clear

set(gcf,'position',[56 1006 1300 342]);

%% constants

T1 = 0.5/log(2); % baseline T1 value (s)
T2 = 10e-6;      % T2 of the bound pool (s)

R1f = 1/T1;      % free spin-lattice relaxation (s^-1)
R1b = 1/T1;      % bound spin-lattice relaxation (s^-1)

kf = 1;          % forward rate (s^-1)
kb = 10;         % backward rate (s^-1)

TR = 5;          % repetition time (s)
dt = 1e-6;       % simulation time step (s)
t = 0:dt:TR;     % time after inversion (s)

Mf_0 = 1;        % free pool size
Mb_0 = kf/kb;    % bound pool size

%% saturating b1 field 

b1 = logspace(-7,-4,51); % (T)
delta = 1000; % frequency (Hz)

%% lineshape (doi.org/10.1002/mrm.29071)

shape = 'Gaussian';

switch shape
    case 'Lorentzian';
        g = (T2/pi)/(1+(2*pi*delta*T2)^2);
    case 'Gaussian';
        g = (T2/sqrt(2*pi))*exp(-(2*pi*delta*T2)^2/2);
    case 'superLorentzian';
        f = @(z)abs(3*z.^2-1);
        g = (T2*sqrt(2/pi))*integral(@(z)exp(-2*(2*pi*delta*T2./f(z)).^2)./f(z),0,1);
end

%% Bloch simulations

T1pre = zeros(numel(b1),1); % fitted T1 pre-zero crossing (s)
T1post = zeros(numel(b1),1); % fitted T1 post-zero crossing (s)
Zcross = zeros(numel(b1),1); % location of zero crossing (s)

for j = 1:numel(b1)

    % loop 1: initialize
    % loop 2: steady state
    for loop = 1:2

        % longitudinal magnetization
        N = numel(t);
        M = zeros(N,2); 

        % inversion pulse (invert Mf, saturate Mb)
        if loop==1
            Mf =-Mf_0;
            Mb = 0;
        else
            Mf =-Mf;
            Mb = 0;
        end

        % time evolution
        for n = 1:N

            M(n,1) = Mf;
            M(n,2) = Mb;

            % RF saturation rate (doi:10.1016/j.neuroimage.2015.03.068)
            gamma = 2.675e8; % rad/s/T
            w1 = gamma*b1(j);
            RFb = pi*w1^2*g;

            % equilibration and T1 relaxation (doi:10.1002/mrm.10386)
            dMfdt = -R1f*(Mf - Mf_0) - kf*Mf + kb*Mb;
            dMbdt = -R1b*(Mb - Mb_0) - kb*Mb + kf*Mf - RFb*Mb;

            % update magnetization
            Mf = Mf + dMfdt*dt;
            Mb = Mb + dMbdt*dt;

        end

    end

    %% calculate parameters

    % zero crossing
    fprintf('B1 = %.1eT\n',b1(j)); 
    [~,k] = min(abs(M(:,1)));
    Zcross(j) = interp1(M(k-1:k+1,1),t(k-1:k+1),0);
    fprintf('Mz zero crossing: %.3fs\n',Zcross(j));

    % analytical values
    T1S = 2/(R1f+R1b+kf+kb+sqrt((R1f-R1b+kf-kb).^2+4*kf*kb)); % (s)
    T1L = 2/(R1f+R1b+kf+kb-sqrt((R1f-R1b+kf-kb).^2+4*kf*kb)); % (s)
    T1CW = 1/(R1f+kf); % (s)

    % fitted T1s
    range = t<Zcross(j); % pre zero-crossing
    T1pre(j) = fit_ir_barral(t(range),M(range,1),T1L);
    range = t>Zcross(j); % post zero-crossing
    T1post(j) = fit_ir_barral(t(range),M(range,1),T1L);

    %% plot signal
    subplot(1,3,1);    
    range = t<2.5; % display range
    plot(1e3*t(range),M(range,1),'color',color(2));
    hold on
    plot(1e3*t(range),5*M(range,2),'color',color(5));
    hold off
    xlabel('Time (ms)'); ylabel('Magnetization');
    grid on; axis square
    ylim([-1 1]);
    xticks((0:5)*500);
    h1 = gca;
    h1.YMinorGrid='on';
    h1.YAxis.MinorTickValues=-1:0.25:1;
    h1.XMinorGrid='off';
    hold on
    M0 = 1-2*exp(-t/T1)+exp(-TR/T1); % monoexponential
    plot(1e3*t(range),M0(range),':','color',color(1));
    hold off
    text(800,-0.375,'CW saturation');
    legend({'M_f','M_b × 5'},'location','northwest');

    % plot observed T1
    subplot(1,3,2);
    h2 = semilogx(1e6*b1,1e3*[T1post T1pre]);
    h2(2).LineStyle='--';
    grid on; xlim([1e-1 1e2]);
    xticks([1e-1 1e0 1e1 1e2]);
    xticklabels({'0.1','1','10','100'});
    axis square
    ylabel('Fitted T_1 (ms)');
    xlabel('B_1 (µT)');
    text(0.13,970*T1L,sprintf('T_{1L}'));    
    text(33,1045*T1CW,sprintf('T_{1CW}'));
    legend({'Post zero crossing','Pre  zero crossing'});
    %set(get(gca(),'XAxis'),'MinorTickValues',10.^(-1:2));
    drawnow;

end

%% generate dsir curve

TI = [0.350; 0.500]; % (s)
T1 = 0.400:dt:0.900; % (s)
T1CW = 1./(1./T1+kf-0.7768893696437); % chosen s.t. ΔT1~100ms (s)

% signals
S1 = 1 - 2*exp(-TI./T1);
S2 = 1 - 2*exp(-TI./T1CW);

dSIR1 = -diff(abs(S1))./sum(abs(S1));
dSIR2 = -diff(abs(S2))./sum(abs(S2));

% apparent shift of T1
[~,j] = min(dSIR1);
[~,k] = max(dSIR1);
deltaT1 = T1(k)-interp1(dSIR1(j:k),T1(j:k),interp1(T1,dSIR2,T1(k)));
fprintf('dSIR1 at T1=%.3fs: %+.5f\n',T1(k),interp1(T1,dSIR1,T1(k)));
fprintf('dSIR2 at T1=%.3fs: %+.5f\n',T1(k),interp1(T1,dSIR2,T1(k)));
fprintf('T1 at dSIR1=%+.5f: %.3fs\n',interp1(T1,dSIR2,T1(k)),T1(k)-deltaT1);
fprintf('Apparent ΔT1: %.3fs\n',deltaT1);
deltaT1-0.1

% plots
subplot(1,3,3);
plot(1e3*T1,dSIR1,'color',color(1));
hold on
plot(1e3*T1,dSIR2,'--','color',color(2));
hold off
grid on
xlabel('T_1 (ms)');
ylabel('dSIR signal');
line([1e3 1e3]*T1(k),[-1 1],'linestyle',':','color','black','linewidth',1.25);
line([1e3 1e3]*interp1(dSIR1(j:k),T1(j:k),interp1(T1,dSIR2,T1(k))),[-1 1],'linestyle',':','color','black','linewidth',1.25);
legend({'MT off','MT on'},'location','northwest');
axis square
text(1e3*interp1(dSIR1(j:k),T1(j:k),interp1(T1,dSIR2,T1(k)))+1e3*deltaT1/4.5,0.25,'\DeltaT_1');

% arrow overlays
subplot(1,3,1); annotation(gcf,'arrow',[0.220 0.212],[0.420 0.518]);
subplot(1,3,3); annotation(gcf,'doublearrow',[0.789 0.826],[0.585 0.585]);
