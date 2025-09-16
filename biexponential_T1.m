% Bloch simulation of biexponential T1 with MT
%
% Mf is free pool magnetization
% Mb is bound pool magnetization
%
%  R1f     kf     R1b
% <--- Mf <==> Mb --->
%          kb
%
clear

set(gcf,'position',[56 1006 1300 342]);

%% constants

T1 = 0.5/log(2); % baseline T1 value (s)

R1f = 1/T1;      % free spin-lattice relaxation (s^-1)
R1b = 1/T1;      % bound spin-lattice relaxation (s^-1)

kf = 1;          % forward rate (s^-1)
kb = 10;         % backward rate (s^-1)

TR = 5;          % repetition time (s)
dt = 1e-6;       % simulation time step (s)
t = 0:dt:TR;     % time after inversion (s)

Mf_0 = 1;        % free pool size
Mb_0 = kf/kb;    % bound pool size

%% "soft" inversion (inverts Mf, saturates Mb)

saturate = 0; % 0=fully saturate 1=do nothing
fprintf('Saturation factor = %.1f\n',saturate);

%% Bloch simulations

% loop 1: exchange off
% loop 2: exchange on
% loop 3: incidental RF on
% loop 4: steady state on
for loop = 1:4

    % longitudinal magnetization
    N = numel(t);
    M = zeros(N,2);

    % inversion pulse (invert Mf, saturate Mb)
    if loop==1 || loop==2 || loop==3
        Mf =-Mf_0;
        Mb = Mb_0 * saturate;
    elseif loop==4
        Mf =-Mf;
        Mb = Mb * saturate;
    end

    % time evolution
    for n = 1:N

        M(n,1) = Mf;
        M(n,2) = Mb;

        % incidental RF
        if loop==3 || loop==4
            if mod(n*dt,0.125)<eps
                Mb = Mb * saturate;
            end
        end

        % equilibration and T1 relaxation (doi:10.1002/mrm.10386)
        dMfdt = -R1f*(Mf - Mf_0) - (kf*Mf - kb*Mb) * (loop>1);
        dMbdt = -R1b*(Mb - Mb_0) - (kb*Mb - kf*Mf) * (loop>1);

        % update magnetization
        Mf = Mf + dMfdt*dt;
        Mb = Mb + dMbdt*dt;

    end

    % location of zero crossing (s)
    [~,k] = min(abs(M(:,1)));
    Zcross(loop) = interp1(M(k-1:k+1,1),t(k-1:k+1),0);
    fprintf('loop %i: zero crossing at %.3fs\n',loop,Zcross(loop));

    % display
    if loop==1

        M0 = M; % save for overplotting on next loops
    
    else

        range = t<2.5; % display range

        subplot(1,3,loop-1);
        plot(1e3*t(range),M(range,1),'color',color(2));
        hold on
        plot(1e3*t(range),5*M(range,2),'color',color(5));
        hold off
        xlabel('Time (ms)');
        if loop==2; ylabel('Magnetization'); end
        grid on;
        ylim([-1 1]);
        xticks((0:5)*500);

        h(loop-1) = gca;
        h(loop-1).YMinorGrid='on';
        h(loop-1).YAxis.MinorTickValues=-1:0.25:1;
        h(loop-1).XMinorGrid='off';
        
        hold on
        plot(1e3*t(range),M0(range,1),':','color',color(1));
        hold off

        legend({'M_f','M_b × 5'},'location','northwest');

        switch loop
            case 2; text(400,-0.625,'Equilibration');
            case 3; text(400,-0.625,'Incidental RF');
            case 4; text(400,-0.625,'Steady state');
        end
        axis square; drawnow

    end

end

fprintf('Long-term magnetization: %f\n',Mf/Mf_0);

% arrow overlays
subplot(1,3,1); annotation(gcf,'arrow',[0.162 0.139],[0.295 0.427]);
subplot(1,3,2); annotation(gcf,'arrow',[0.441 0.423],[0.295 0.427]);
subplot(1,3,3); annotation(gcf,'arrow',[0.724 0.695],[0.271 0.230]);