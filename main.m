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
alpha_ext = 10;
alpha_air = 10;

%% Numerical
t_inic = 0;
t_max = 1800;
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

T = zeros(t_max/inc_t+1, total_nod);

T_inic = T_ext;
T(1,:) = T_inic;


%% Calculating coefficients ([W/K])
% 
i = inc_t;
j = 1;

while i<=t_max
    [ap,ae, aw, bp] = coefficient_calc(coord_total, total_nod, alpha_ext, alpha_air ,T(j,:), nod_reforc_1, nod_poli, nod_reforc_2, nod_poma, inc_t);
    [P,R] = matrix_elements(ap,ae, aw, bp, total_el);
    T(j+1,:) = temp_field_calc(P, R, n);
    
end

