aigua = 0.84*(0.57109+0.0017625*T-0.0000067036*T^2);
carbo = 0.136*(0.20141+0.0013874*T-0.0000043312*T^2);
fibra = 0.024*(0.18331+0.0012497*T-0.0000031683*T^2);

% Proposta mallat:
% 5 - reforços
% 100 - poma i poliuretà


%% Calculating coefficients ([W/K])
[ap,ae, aw, bp, node] = coefficient_calc(Rext,Rint,lambda,n, ef, alpha_ext, Text, alpha_end, Twall);
[P,R] = matrix_elements(ap,ae, aw, bp, n);


%% Initiation
[T] = temp_field_calc(P, R, n);


%Postprocessing
text_line = 0:0.5:Rext;
Text_vec = zeros(size(text_line,2), 1);
for i = 1:size(Text_vec,1)
   Text_vec(i) = Text;
end
figure
plot(node, T, 'r');
xlabel('r [m]');
ylabel('T [K]');
title('Temperature along the circular fin');
grid on

hold on
plot(text_line, Text_vec, 'black');
hold off
