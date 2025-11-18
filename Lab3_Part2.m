%% Lab 3 Part 2
clc;
close all;
clear;

% Part 2 Task 1 - Prandtl Lifting Line Code

%wing characteristics
b= 100; %span [feet]
a0_t = 6.3; %cross-sectional lift slope at the tips [radians]
a0_r = 6.5; % cross-section lift slope at the root [radians]
c_t = 8; % chord at tips [ft]
c_r = 10; % chord at root [ft]
aero_t = 0; % zero-lift AoA at tips [radians]
aero_r = -2*pi/180; % zero-lift AoA at root [radians]
geo_t = 5*pi/180; % geometric AoA (geometric twist + alpha) at tips [radians]
geo_r = 7*pi/180; % geometric AoA (geometric twist + alpha) at root [radians]

N = 5; % # of odd terms for circulation


[e,c_L,c_Di,A_matrix] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);


function [e,c_L,c_Di,A_matrix] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
    theta = (1:N)' * pi / (2*N);
    c = c_r + (c_t-c_r) .* cos(theta);
    a0 = a0_r + (a0_t - a0_r) .* cos(theta);
    alphaL0 = aero_r + (aero_t - aero_r) .* cos(theta);
    alpha_geo = geo_r + (geo_t - geo_r) .* cos(theta);
    b_i = (alpha_geo - alphaL0);
    Mij = zeros(N);
    for i = 1:N
        ac_Theta = theta(i);
        mu = 4 * b / (a0(i) * c(i));
        for j = 1:N
            n = 2 * j - 1;
            Mij(i,j) = (mu + (n / sin(ac_Theta))) * sin(n * ac_Theta);
        end
    end
    A_matrix = Mij\b_i;
    
    S = b*((c_r +c_t)/2);
    AR = b^2/S;
    %coefficent of lift
    c_L = A_matrix(1)*pi*AR;

    C_Di_total = 0;
    for j=1:N
        n = 2 * j - 1;
        C_Di_total = C_Di_total + n * (A_matrix(j).^2);
    end

    c_Di = pi * AR * C_Di_total;

    %e
    e = A_matrix(1).^2 ./C_Di_total;
    %coefficent of drag
    c_Di = (c_L^2)/(pi*AR*e);

    %taper = linspace(0,(c_t/c_r),length(delta));
    % figure()
    % hold on;
    % plot(taper, delta,'k');
    % xlabel('Taper ratio c_t/c_r');
    % ylabel('delta');

end

    % 
    % for i=1:length(y_0)
    %     sec_lift = rho_inf*v_inf *circ(y_0);
    % end
