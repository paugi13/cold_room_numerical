function [ap,ae, aw, bp] = coefficient_calc(coord_total, total_nod, alpha_ext, alpha_air, T, nod_reforc_1, nod_poli, nod_reforc_2, nod_poma, inc_t)
% Function to calculate all the coefficients through the fin.
% They are returned in vector format.

ap = zeros(total_nod, 1);
aw = zeros(total_nod, 1);
ae = zeros(total_nod, 1);
bp = zeros(total_nod, 1);


% First conductivity constants must be calculated
% Critical nodes: Boundaries 
%   nod_reforc_1(1), nod_reforc_1(end), nod_poli(end), nod_reforc_2(end),
%   nod_poma(1), nod_poma(end)

% Lambda calculus

for i = 2:(size(coord_total,2)-1)
lambda_e = 0;
lambda_w = 0;

    if i < nod_reforc_1(end)
        lambda_w = 204;
        lambda_e = 0.028;
    elseif i == nod_reforc_1(end)
        lambda_w = 204;
        lambda_e = 0.028;
    elseif i < nod_poli(end)
        lambda_w = 0.028;
        lambda_e = 0.028;
    elseif i == nod_poli(end)
        lambda_w = 0.028;
        lambda_e = 204;
    elseif i < nod_reforc_2(end)
        lambda_w = 204;
        lambda_e = 204;
    elseif i == nod_reforc_2(end)
        lambda_w = 204;
    elseif i == nod_poma(1)
        lambda_e = 0.5115;
    elseif i < nod_poma(end)
        lambda_w = 0.5115;
        lambda_e = 0.5115;
    end
    d_PW = coord_total(1,i)-coord_total(1,i-1);
    d_PE = coord_total(1,i+1)-coord_total(1,i);
    aw(i) = lambda_w/d_PW;
    ae(i) = lambda_e/d_PE;
    ap(i) = ae + aw + p(i)*cp*(d_PW/2+d_PE/2)/inc_t;
    bp(i) = T(i)*p(i)*cp*(d_PW/2+d_PE/2)/inc_t;
    
end


%% Nodes with convection flow
% First node
d_PE = coord_total(1,2)-coord_total(1,1);
aw(1) = 0;
ae(1) = 204/(d_PE);
ap(1) = ae(1) + alpha_ext + p(1)*cp*(d_PE/2)/inc_t;
bp(1) = alpha_ext*Text + p(1)*cp*(d_PE/2)*T(1)/inc_t;

%Last reforc node
d_PW = coord_total(1,nod_reforc_2(end))-coord_total(1,nod_reforc_2(end)-1);
aw(nod_reforc_2(end)) = 204/(coord_total(1,2)-coord_total(1,1));
ae(nod_reforc_2(end)) = 0;
ap(nod_reforc_2(end)) = aw(nod_reforc_2(end)) + alpha_air + p(1)*cp*(d_PW/2)/inc_t;
bp(nod_reforc_2(end)) = alpha_air*T_air + p(nod_reforc_2(end))*cp*(d_PE/2)*T(nod_reforc_2(end))/inc_t;

% First apple node
d_PE = coord_total(1,nod_poma(1)+1)-coord_total(1,nod_poma(1));
aw(nod_poma(1)) = 0;
ae(nod_poma(1)) = 0.5115/(d_PE);
ap(nod_poma(1)) = ae(1) + alpha_air + p(nod_poma(1))*cp*(d_PE/2)/inc_t;
bp(nod_poma(1)) = alpha_air*Text + p(nod_poma(1))*cp*(d_PE/2)*T(nod_poma(1))/inc_t;





% aigua = 0.84*(0.57109+0.0017625*T-0.0000067036*T^2);
% carbo = 0.136*(0.20141+0.0013874*T-0.0000043312*T^2);
% fibra = 0.024*(0.18331+0.0012497*T-0.0000031683*T^2);





