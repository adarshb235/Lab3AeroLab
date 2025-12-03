clc;
clear;
close all;


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
b = 36; %span [feet]
c_r = 5 + 4/12; % chord at root [ft]
c_t = 3 + 7/12; % chord at tips [ft]


N = 5; % # of odd terms for circulation

alpha_exp = [-18.5, -18.25, -18.0, -17.75, -17.5, -17.25, -17.0, -16.75, -16.5, -16.25, -16.0, -15.75, -15.5, -15.25, -15.0, -14.75, -14.5, -14.25, -14.0, -13.75, -13.5, -13.25, -13.0, -12.75, -12.5, -12.25, -12.0, -11.75, -11.5, -11.25, -11.0, -10.75, -10.5, -10.25, -10.0, -9.75, -9.5, -9.25, -9.0, -8.75, -8.5, -8.25, -8.0, -7.75, -7.5, -7.25, -7.0, -6.75, -6.5, -6.25, -6.0, -5.75, -5.5, -5.25, -5.0, -4.75, -4.5, -4.25, -4.0, -3.75, -3.5, -3.25, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75, 18.0] * pi/180;
p_drag = [0.09929, 0.09187, 0.08451, 0.07744, 0.07074, 0.06446, 0.05874, 0.05351, 0.04856, 0.04283, 0.038, 0.0342, 0.03122, 0.02885, 0.0269, 0.02518, 0.02363, 0.02223, 0.02099, 0.01991, 0.01896, 0.01816, 0.0167, 0.01576, 0.01504, 0.01438, 0.01376, 0.01321, 0.01274, 0.01175, 0.01118, 0.01069, 0.01023, 0.00982, 0.00912, 0.00868, 0.00833, 0.00801, 0.00747, 0.00712, 0.00682, 0.00642, 0.00612, 0.00587, 0.00552, 0.00527, 0.005, 0.00475, 0.00446, 0.00419, 0.00395, 0.00371, 0.00348, 0.00304, 0.00284, 0.00266, 0.00248, 0.00232, 0.00217, 0.00204, 0.00191, 0.0018, 0.0017, 0.0016, 0.00151, 0.00144, 0.00137, 0.00131, 0.00126, 0.00122, 0.00118, 0.00116, 0.00115, 0.00114, 0.00115, 0.00116, 0.00118, 0.00122, 0.00126, 0.00131, 0.00137, 0.00144, 0.00151, 0.0016, 0.0017, 0.0018, 0.00191, 0.00204, 0.00218, 0.00232, 0.00248, 0.00265, 0.00284, 0.00304, 0.00348, 0.00371, 0.00395, 0.00419, 0.00446, 0.00475, 0.00499, 0.00527, 0.00552, 0.00587, 0.00612, 0.00642, 0.00682, 0.00712, 0.00747, 0.00801, 0.00832, 0.00868, 0.00912, 0.00982, 0.01023, 0.01069, 0.01118, 0.01175, 0.01274, 0.01321, 0.01376, 0.01438, 0.01504, 0.01576, 0.0167, 0.01815, 0.01895, 0.01991, 0.02099, 0.02222, 0.02361, 0.02516, 0.02687, 0.02885, 0.03121, 0.03418, 0.038, 0.0428, 0.04854, 0.0534, 0.05856, 0.06434, 0.0706, 0.07731, 0.08442, 0.09174, 0.09922];

% experimental data attained from airfoiltools.com for NACA 0012 ^

A_mat = zeros(length(alpha_exp), N);
liftingLine = struct();
Cl_ref = struct();

for i = 1:length(alpha_exp)
    alpha_1 = alpha_exp(i);
    geo_r = alpha_1 + (2 * pi / 180); % geometric AoA (geometric twist + alpha) at tips [radians]
    geo_t = alpha_1; % geometric AoA (geometric twist + alpha) at root [radians]
    Cl_ref.tip(i) = Vortex_Panel(x_b_tip, y_b_tip, alpha_1);
    Cl_ref.root(i) = Vortex_Panel(x_b_root, y_b_root, alpha_1);
    [liftingLine.e(i), liftingLine.CL(i), liftingLine.CD_i(i), liftingLine.CD_prof(i), ~] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N, p_drag);    
    liftingLine.CD_tot(i) = liftingLine.CD_prof(i) + liftingLine.CD_i(i);
end
S = (c_r + c_t)/2 * b;


figure(); % Task 1 Plot
hold on;
grid on;
box on;
plot(alpha_exp * 180/pi, liftingLine.CL, 'LineWidth', 1.5);
xlabel('Angle of Attack (deg)');
ylabel('Lift Coefficient (C_L)');
xlim([-10, 15])
title('Lift Coefficient vs Angle of Attack');
legend('Cessna 180 (NACA 2412 Root, NACA 0012 Tip)', 'Location', 'best')

figure; % task 2
hold on;
grid on;
box on;
xlim([-16 16])
plot(alpha_exp * 180/pi, liftingLine.CD_prof, 'LineWidth', 2);
plot(alpha_exp * 180/pi, p_drag, 'LineWidth', 1.5)
xlabel('Angle of Attack (deg)');
ylabel('Drag Coefficient (C_D)');
title('Drag Coefficient vs Angle of Attack');
legend('Calc Profile Drag', 'Experimental Data', 'Location', 'best')

figure; % task 3 plot
hold on;
grid on;
box on;
xlim([-15 15]);
plot(alpha_exp * 180/pi, liftingLine.CD_tot, 'k', 'LineWidth', 2);
plot(alpha_exp * 180/pi, liftingLine.CD_i, 'b--', 'LineWidth', 1.5);
plot(alpha_exp * 180/pi, liftingLine.CD_prof, 'r-.', 'LineWidth', 1.5);
xlabel('Angle of Attack (deg)');
ylabel('Drag Coefficient');
legend('Total Drag (C_D)', 'Induced Drag (C_{Di})', 'Profile Drag (C_{D0})', 'Location', 'best');
title('Drag Components vs Angle of Attack');

rho = 0.001756; % density [slugs/ft^3] at 10,000' standard day
W = 2500; % aircraft weight [lbs]

V = zeros(1, length(alpha_exp));
T_req = zeros(1, length(alpha_exp));

for i = 1:length(alpha_exp)
    if liftingLine.CL(i) > 0
        V_1 = sqrt((2*W)/(rho * S * liftingLine.CL(i)));
        V(i) = 0.592484 * V_1;
        T_req(i) = 0.5 * rho * V_1^2 * S * liftingLine.CD_tot(i);
    else
        V(i) = NaN;
        T_req(i) = NaN;
    end
end

figure(); % task 4 plot
plot(V, T_req, 'b-', 'LineWidth', 2);
xlabel('Airspeed (knots)');
ylabel('Thrust Required (lbs)');
title('Thrust Required for Steady Level Flight (10,000 ft)');
legend('Cessna 180 (NACA 2412 Root, NACA 0012 Tip)', 'Location', 'northwest')
grid on;

filenamesVec = ["task1liftpolar", "task2dragvsangle", "task3dragpolar", "task4velvsthrust"];

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

function [e, c_L, c_Di, C_D0, A_matrix] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N, p_drag)
    theta = (1:N)' * pi / (2*N); 
    
    c = c_r + (c_t - c_r) .* cos(theta);              
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
    
    A_matrix = Mij \ b_i; 
    
    S = (c_r + c_t)/2 * b; 
    AR = b^2/S;
    c_L = A_matrix(1) * pi * AR; 
    
    sum_n_An2 = 0; 
    for j = 1:N
        n = 2 * j - 1;
        sum_n_An2 = sum_n_An2 + n * (A_matrix(j).^2);
    end
    c_Di = pi * AR * sum_n_An2; 
    
    if sum_n_An2 == 0
        e = 0;
    else
        e = A_matrix(1)^2 / sum_n_An2;
    end

    alpha_i_vec = zeros(size(theta));
    for i = 1:length(theta)
        sum_val = 0;
        for j = 1:N
            n = 2*j - 1;
            sum_val = sum_val + n * A_matrix(j) * sin(n * theta(i));
        end
        alpha_i_vec(i) = sum_val / sin(theta(i));
    end

    alpha_eff = alpha_geo - alpha_i_vec;
    
    alpha_exp = [-18.5, -18.25, -18.0, -17.75, -17.5, -17.25, -17.0, -16.75, -16.5, -16.25, -16.0, -15.75, -15.5, -15.25, -15.0, -14.75, -14.5, -14.25, -14.0, -13.75, -13.5, -13.25, -13.0, -12.75, -12.5, -12.25, -12.0, -11.75, -11.5, -11.25, -11.0, -10.75, -10.5, -10.25, -10.0, -9.75, -9.5, -9.25, -9.0, -8.75, -8.5, -8.25, -8.0, -7.75, -7.5, -7.25, -7.0, -6.75, -6.5, -6.25, -6.0, -5.75, -5.5, -5.25, -5.0, -4.75, -4.5, -4.25, -4.0, -3.75, -3.5, -3.25, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0, 14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75, 17.0, 17.25, 17.5, 17.75, 18.0] * pi/180;
    cd_local = interp1(alpha_exp, p_drag, alpha_eff, 'pchip', 'extrap');

    y_stations = (b/2) * cos(theta);
    integrand = cd_local .* c;
    trapz_int = abs(trapz(y_stations, integrand)); 
    
    C_D0 = (2/S) * trapz_int;
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
