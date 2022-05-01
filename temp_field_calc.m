function [T] = temp_field_calc(P, R, n)

T = zeros(n+1,1);

T(n+1) = R(n+1);

for i = n:-1:1
    T(i) = P(i)*T(i+1) + R(i);
end
