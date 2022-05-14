%% Simulation of a cold room filled with apples
% Implicit method is used with constant values of lambda.
% Proposta mallat:
% 5 - reforços
% 100 - poma i poliuretà

%% Fisical
l_air = 2.5;
l_poli = 0.2;
l_reforc = 0.01;
l_palot = 0.5;
T_ext = 15;
T_air = 0;
T_aux = 10;
alpha_ext = 10;
alpha_air = 10;
alpha_aux = 10;
h = 2.4;    % Height

%% Numerical
t_inic = 0;
n_hours = 24;
t_max = n_hours*3600;
inc_t = 60;
n_el_poma = 100;
n_el_poli = 100;
n_el_reforc = 10;

n_nod_poma_n = n_el_poma + 1;
n_nod_poli_n = n_el_poli;
n_nod_reforc_n = n_el_reforc+1; % One node is shared with the poliutheran core.


%% Some geometrical values.

ax_reforc = l_reforc/n_el_reforc;
ax_poli = l_poli/n_el_poli;
ax_poma = l_palot/n_el_poma;

%% Numerating nodes.
total_el = n_el_poli + n_el_reforc*2 + n_el_poma;
total_nod = n_nod_poma_n + n_nod_poli_n + n_nod_reforc_n*2-1;

nod_reforc_1 = 1:1:n_nod_reforc_n;
nod_poli = (n_nod_reforc_n+1):1:(n_nod_reforc_n+1+n_nod_poli_n-1);
nod_reforc_2 = (n_nod_reforc_n+1+n_nod_poli_n):1:(n_nod_reforc_n+1+n_nod_poli_n+n_nod_reforc_n-2);
nod_poma = (n_nod_reforc_n+1+n_nod_poli_n+n_nod_reforc_n-1):1:...
    (n_nod_reforc_n+1+n_nod_poli_n+n_nod_reforc_n+n_nod_poma_n-2);

%% Setting node coordinates
[coord_total] = ...
    node_coord(nod_reforc_1,nod_poli, nod_reforc_2, nod_poma, ax_reforc, ax_poli, ax_poma, l_air, total_nod);

%% Calculating coefficients ([W/K]) and solving

T = zeros(1, total_nod);

[ap,ae, aw, bp] = coefficient_calc(coord_total, total_nod, alpha_ext,...
    alpha_air, nod_reforc_1, nod_poli, nod_reforc_2, nod_poma, T_ext, T_air, T_aux, alpha_aux);
[P,R] = matrix_elements(ap,ae, aw, bp, total_nod, nod_poma);
T(1,:) = temp_field_calc(P, R, total_nod);


%% POSTPROCESSING

y = [0 2.4];

% T_air must be included in the 
T_plot = zeros(size(T,1), size(T,2)+1);
T_plot(1,1:nod_reforc_2(end)) = T(1, 1:nod_reforc_2(end));
T_plot(1, nod_reforc_2(end)+1) = T_air;
T_plot(1, nod_reforc_2(end)+2:end) = T(1, nod_poma(1):end);

[X,Y] = meshgrid(coord_total,y);

figure
surf(X,Y, zeros(size(X)),T_plot(1,2:end), 'edgecolor','none');
xlabel('x [m]');
ylabel('y [m]');
title('Mapa de temperatures');
xlim([0 3.22]);
ylim([0 2.4]);
view(0,90);
colorbar

% shading interp

% figure
% imagesc(coord_total, [0 2.4], T(20,:));
