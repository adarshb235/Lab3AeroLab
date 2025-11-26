clc;
clear;
close all;


alpha = (-10:15) .* (pi/180);

tipCode = 0012;
panelCount = 64;
tipVec = digi(tipCode);
[x_b_tip, y_b_tip] = genNACA(1, tipVec(1), tipVec(2), str2double(strcat(num2str(tipVec(3)), num2str(tipVec(4)))), (panelCount/2) + 2, 0);

rootCode = 2412;
panelCount_2 = 68;
rootVec = digi(rootCode);
[x_b_root, y_b_root] = genNACA(1, rootVec(1), rootVec(2), str2double(strcat(num2str(rootVec(3)), num2str(rootVec(4)))), (panelCount_2/2) + 2, 0);

alpha_2D_deg = -4:2:4; % Degrees for 2D analysis
alpha_2D_rad = alpha_2D_deg .* (pi/180); 
Cl_tip_2D = zeros(1, length(alpha_2D_deg));
Cl_root_2D = zeros(1, length(alpha_2D_deg));

for i = 1:length(alpha_2D_deg)
    Cl_tip_2D(i)  = Vortex_Panel(x_b_tip, y_b_tip, alpha_2D_deg(i));
    Cl_root_2D(i) = Vortex_Panel(x_b_root, y_b_root, alpha_2D_deg(i));
end

p_tip = polyfit(alpha_2D_rad, Cl_tip_2D, 1);
a0_t = p_tip(1);           % Lift slope [per rad]
aero_t = -p_tip(2)/p_tip(1); % Zero-lift AoA [rad] (x-intercept)

p_root = polyfit(alpha_2D_rad, Cl_root_2D, 1);
a0_r = p_root(1);          % Lift slope [per rad]
aero_r = -p_root(2)/p_root(1); % Zero-lift AoA [rad] (x-intercept)

%wing characteristics
b= 36; %span [feet]
c_t = 5 + 4/12; % chord at tips [ft]
c_r = 3 + 7/12; % chord at root [ft]


N = 5; % # of odd terms for circulation

A_mat = zeros(length(alpha), N);
liftingLine = struct();
Cl_ref = struct();


for i = 1:length(alpha)
    alpha_1 = alpha(i);
    geo_t = alpha(i) - (2 * pi / 180); % geometric AoA (geometric twist + alpha) at tips [radians]
    geo_r = alpha(i); % geometric AoA (geometric twist + alpha) at root [radians]
    Cl_ref.tip(i) = Vortex_Panel(x_b_tip, y_b_tip, alpha(i));
    Cl_ref.root(i) = Vortex_Panel(x_b_root, y_b_root, alpha(i));
    [liftingLine.e(i),liftingLine.c_L(i),liftingLine.c_Di(i), A_mat(i, :)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r, N);

end

figure(); % Task 1 Plot
hold on;
plot(alpha * 180 / pi, liftingLine.c_L);
xlabel('Angle of Attack (deg)');
ylabel('Lift Coefficient (C_L)');
title('Lift Coefficient vs Angle of Attack');
grid on;
box on;
legend('Cessna 180 (NACA 2412 Root, NACA 0012 Tip)', 'Location', 'best')






% Part 1 NACA Code Generator

function [x_b, y_b] = genNACA(c, m, p, t, panel, graph)

m = m/100; 
p = p/10;
t = t/100;

theta = linspace(0, pi, panel);

x = c./2 * (1-cos(theta));
x_c = x/c;

y_t = (t/0.2) * c * ((0.2969 * sqrt(x_c)) - (0.126 * x_c) - (0.3516 * (x_c).^2) + (0.2843 * (x_c).^3) - (0.1036 * (x_c).^4));

y_c = zeros(1, length(x));

if (m == 0 || p == 0)
    dyc_dx = 0;
else
    for i = 1:length(x)
        if (x(i) >= 0) && (x(i) < p*c)
            y_c(i) = m * (x(i) / p^2) * (2*p - x_c(i));
        elseif (x(i) >= p*c) && (x(i) <= c)
            y_c(i) = m * ((c - x(i)) / (1 - p)^2) * (1 + (x_c(i)) - 2*p);
        end
    end
    dyc_dx = gradient(y_c, x);
end


xi = atan(dyc_dx);

x_u = x - y_t .* sin(xi);
x_l = x + y_t .* sin(xi);
y_u = y_c + y_t .* cos(xi);
y_l = y_c - y_t .* cos(xi);

x_b = [fliplr(x_l), x_u(2:end)];
y_b = [fliplr(y_l), y_u(2:end)];

if (graph == 1)
    figure;
    hold on;
    axis equal;
    grid on;
    box on;
    plot(x_b, y_b, 'k');
    plot(x, y_c)
    title("NACA" + m*100 + p*10 + t*100)
    legend('Foil Outline', 'Mean Camber Line', 'Location', 'best');
    filename = append('airfoilNACA', num2str(m*100), num2str(p*10), num2str(t*100)); 
    ax = gca; 
    set(ax, 'LooseInset', get(ax, 'TightInset'));
    print(filename,'-r500','-dpng')
end

end

% Part 2 PLLT Code

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

% Aux Function To Supplement Use of Part 1 Work

function [vec] = digi(input)
    vec = num2str(input);
    if length(vec) < 4
        % Pad with leading zeros if it's shorter than 4
        vec = pad(vec, 4, 'left', '0'); 
    elseif length(vec) > 4
        % Truncate to 4 characters if it's longer
        vec = vec(1:4);
    end
    vec = vec - '0';
end
