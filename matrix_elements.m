function [P,R] = matrix_elements(ap,ae, aw, bp, total_nod, nod_poma)
% Critical nodes: Boundaries 
%   nod_reforc_1(1), nod_reforc_1(end), nod_poli(end), nod_reforc_2(end),
%   nod_poma(1), nod_poma(end)

P = zeros(total_nod,1);
R = zeros(total_nod,1);

P(1) = ae(1)/ap(1);
R(1) = bp(1)/ap(1);

for i = 2:total_nod
    P(i) = ae(i)/(ap(i) - aw(i)*P(i-1));
    R(i) = (bp(i) + aw(i)*R(i-1))/(ap(i)-aw(i)*P(i-1));
    if i == nod_poma(1)
    % nod_poma(1) has no node at its left
        P(i) = ae(i)/ap(i);
        R(i) = bp(i)/ap(i);
    end
end



