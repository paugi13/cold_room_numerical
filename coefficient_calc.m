function [ap,ae, aw, bp] = coefficient_calc(coord_total, total_nod, alpha_ext,...
    alpha_air, T, nod_reforc_1, nod_poli, nod_reforc_2, nod_poma, inc_t, T_ext, T_air, j)
% Function to calculate all the coefficients through the fin.
% They are returned in vector format.

ap = zeros(total_nod, 1);
aw = zeros(total_nod, 1);
ae = zeros(total_nod, 1);
bp = zeros(total_nod, 1);

% Lambda values
lambda_poma = 0.5115;
lambda_poli = 0.028;
lambda_alu = 204;

% Aluminium, poliutheran and apple's densities.
p_alu = 2698.4;
p_poli = 50;
p_poma = 947;

% Values for cp.
cp_alu = 896;
cp_poli = 1674;
cp_poma = 3.7628e+03;

% First conductivity constants must be calculated
% Critical nodes: Boundaries 
%   nod_reforc_1(1), nod_reforc_1(end), nod_poli(end), nod_reforc_2(end),
%   nod_poma(1), nod_poma(end)

% Lambda calculus

for i = 2:(size(coord_total,2)-1)
lambda_e = 0;
lambda_w = 0;

    if i < nod_reforc_1(end)
        lambda_w = lambda_alu;
        lambda_e = lambda_poli;
        p = p_alu;
        cp = cp_alu;
    elseif i == nod_reforc_1(end)
        lambda_w = lambda_alu;
        lambda_e = lambda_poli;
        p = (p_alu + p_poli)/2;
    elseif i < nod_poli(end)
        lambda_w = lambda_poli;
        lambda_e = lambda_poli;
        p = p_poli;
        cp = cp_poli;
    elseif i == nod_poli(end)
        lambda_w = lambda_poli;
        lambda_e = lambda_alu;
        p = (p_poli + p_alu)/2;
    elseif i < nod_reforc_2(end)
        lambda_w = lambda_alu;
        lambda_e = lambda_alu;
        p = p_alu;
        cp = cp_alu;
    elseif i == nod_reforc_2(end)
        lambda_w = lambda_alu;
        p = p_alu;
    elseif i == nod_poma(1)
        lambda_e = lambda_poma;
        p = p_poma;
    elseif i < nod_poma(end)
        lambda_w = lambda_poma;
        lambda_e = lambda_poma;
        p = p_poma;
        cp = cp_poma;
    end
    d_PW = coord_total(1,i)-coord_total(1,i-1);
    d_PE = coord_total(1,i+1)-coord_total(1,i);
    aw(i) = lambda_w/d_PW;
    ae(i) = lambda_e/d_PE;
    ap(i) = ae(i) + aw(i) + p*cp*(d_PW/2+d_PE/2)/inc_t;
    bp(i) = T(j-1,i)*p*cp*(d_PW/2+d_PE/2)/inc_t;
end


%% Nodes with convection flow
% First node
d_PE = coord_total(1,2)-coord_total(1,1);
aw(1) = 0;
ae(1) = lambda_alu/(d_PE);
ap(1) = ae(1) + alpha_ext + p_alu*cp_alu*(d_PE/2)/inc_t;
bp(1) = alpha_ext*T_ext + p_alu*cp_alu*(d_PE/2)*T(1)/inc_t;

%Last reforc node
d_PW = coord_total(1,nod_reforc_2(end))-coord_total(1,nod_reforc_2(end)-1);
aw(nod_reforc_2(end)) = lambda_alu/(coord_total(1,2)-coord_total(1,1));
ae(nod_reforc_2(end)) = 0;
ap(nod_reforc_2(end)) = aw(nod_reforc_2(end)) + alpha_air + p_alu*cp_alu*(d_PW/2)/inc_t;
bp(nod_reforc_2(end)) = alpha_air*T_air + p_alu*cp_alu*(d_PE/2)*T(nod_reforc_2(end))/inc_t;

% First apple node
d_PE = coord_total(1,nod_poma(1)+1)-coord_total(1,nod_poma(1));
aw(nod_poma(1)) = 0;
ae(nod_poma(1)) = lambda_poma/(d_PE);
ap(nod_poma(1)) = ae(1) + alpha_air + p_poma*cp_poma*(d_PE/2)/inc_t;
bp(nod_poma(1)) = alpha_air*T_air + p_poma*cp_poma*(d_PE/2)*T(nod_poma(1))/inc_t;

% Last apple node: Adiabatic end is assumed because of symmetric conditions
aw(nod_poma(end)) = 0;
ae(nod_poma(end)) = 1;
ap(nod_poma(end)) = 1;
bp(nod_poma(end)) = 0;
