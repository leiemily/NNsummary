function tumfrac = cellmovemodel_3(Yp)
% 应用于cellmodel_3函数
a1 = 0.0174001;
b1 = 1.32887;
a2 = 12.0809;
b2 = -5.83762;
c1 = 2892.12;

nf = 8;
omega = 6;

lambda_f = a1 * exp (b1 * Yp);
mu_f = a2 * exp(-(b2 - Yp).^4 / c1);


P_run = zeros(1,length(Yp));
for i = omega : nf
    P_ccwi = nchoosek(nf, i) * (mu_f ./(mu_f + lambda_f)).^i .* ( lambda_f ./(mu_f + lambda_f)).^(nf - i);
    P_run = P_run + P_ccwi;
end
P_tumble = 1 - P_run;
tumfrac = P_tumble;
% lambda = omega * lambda_f .* P_ccww ./ P_run;
% mu = (nf - omega + 1) * mu_f .* P_ccww_1 ./P_tumble;
% tumfrac = lambda ./ (lambda + mu);

end





