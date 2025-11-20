%% Lab 3 Part 2
clc;
close all;
clear;

%% Part 2 Task 1 - Prandtl Lifting Line Code

%wing characteristics
b= 100; %span [feet]
a0_t = 2*pi; %cross-sectional lift slope at the tips [radians]
a0_r = 2*pi; % cross-section lift slope at the root [radians]
c_t = 8; % chord at tips [ft]
c_r = 10; % chord at root [ft]
aero_t = 0; % zero-lift AoA at tips [radians]
aero_r = 0; % zero-lift AoA at root [radians]
geo_t = 5*pi/180; % geometric AoA (geometric twist + alpha) at tips [radians]
geo_r = 5*pi/180; % geometric AoA (geometric twist + alpha) at root [radians]

N = 5; % # of odd terms for circulation


[e,c_L,c_Di,A_matrix] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N); % Task 1 Function Call

%% Part 2 Task 2 - Effect of Taper 

N = 50; % Redefine terms to match textbook figure

AR = [4; 6; 8; 10]; % fix AR to match textbook figure
count = 50; % arb count for vec length
scale = ones(1, count);
deltaVec = zeros(count, length(AR));
c_t_fix = linspace(0, c_r, count); % linear space for varying tip chord vals
c_r_fix = c_r * scale; % all factors except c_t and AR are fixed to match conditions in above test case
ratioVec = c_t_fix ./ c_r_fix;
b_fix = AR .* ((c_r_fix + c_t_fix)./2 .* scale); 
a0_t_fix = a0_t .* scale; 
a0_r_fix = a0_r .* scale; 
aero_t_fix = aero_t .* scale;
aero_r_fix = aero_r .* scale;
geo_t_fix = geo_t .* scale; 
geo_r_fix = geo_r .* scale;

legend_entry = cell(1, length(AR)); % legend entries setup

figure; % figure setup
hold on;
for j = 1:length(AR) % iterate for aspect ratios
    for i = 1:count % iterate delta for different chord ratio
        [e_2,~,~,~] = PLLT(b_fix(j, i), a0_t_fix(i), a0_r_fix(i), c_t_fix(i), c_r_fix(i), aero_t_fix(i), aero_r_fix(i), geo_t_fix(i), geo_r_fix(i), N);
        deltaVec(i, j) = (1-e_2)/e_2;
    end
    legend_entry{j} = ['AR = ', num2str(AR(j))];
end
plot(ratioVec, deltaVec, 'LineWidth', 2); % plot
xlabel('Chord Ratio (c_t/c_r)');
ylabel('Induced Drag Factor (\delta)');
title('Effect of Chord Ratio on Induced Drag Efficiency');
legend(legend_entry, 'Location','northeast');
grid on;
box on;
filename = 'task2delta';  % save figures
ax = gca; 
set(ax, 'LooseInset', get(ax, 'TightInset'));
print(filename,'-r500','-dpng')



function [e,c_L,c_Di,A_matrix] = PLLT(b, a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
    theta = (1:N)' * pi / (2*N); % setup and treating data
    c = c_r + (c_t-c_r) .* cos(theta);
    a0 = a0_r + (a0_t - a0_r) .* cos(theta);
    alphaL0 = aero_r + (aero_t - aero_r) .* cos(theta);
    alpha_geo = geo_r + (geo_t - geo_r) .* cos(theta);
    b_i = (alpha_geo - alphaL0);
    Mij = zeros(N);
    for i = 1:N % iterate
        ac_Theta = theta(i);
        mu = 4 * b / (a0(i) * c(i));
        for j = 1:N % iterate
            n = 2 * j - 1;
            Mij(i,j) = (mu + (n / sin(ac_Theta))) * sin(n * ac_Theta);
        end
    end
    A_matrix = Mij \ b_i; % matrix out
    
    S = b*((c_r +c_t)/2); % geometric setup
    AR = b^2/S;
    c_L = A_matrix(1)*pi*AR; %lift coeff

    sum_n_An2 = 0; % iterative solve
    for j=1:N
        n = 2 * j - 1;
        sum_n_An2 = sum_n_An2 + n * (A_matrix(j).^2);
    end

   c_Di = pi * AR * sum_n_An2; % induced drag

   if sum_n_An2 == 0
        e = 0;
    else
        e = A_matrix(1)^2 / sum_n_An2;
    end
    
    
end


