function [T] = temp_field_calc(P, R, total_nod)

T = zeros(total_nod,1);

T(total_nod) = R(total_nod);

for i = total_nod:-1:1
    T(i) = P(i)*T(i+1) + R(i);
end
