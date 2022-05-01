function [P,R] = matrix_elements(ap,ae, aw, bp, n)

P = zeros(n+1,1);
R = zeros(n+1,1);

P(1) = ae(1)/ap(1);
R(1) = bp(1)/ap(1);

for i = 2:(n+1)
    P(i) = ae(i)/(ap(i) - aw(i)*P(i-1));
    R(i) = (bp(i) + aw(i)*R(i-1))/(ap(i)-aw(i)*P(i-1));
end



