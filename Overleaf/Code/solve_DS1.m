function err = solve_DS1(Ds1, p, z, d,exp1)
x = (d.*z.^(-1/p.eta_h-1)./Ds1).^(exp1);
x = min(max(x, p.xmin), p.xmax);
err = Ds1 - (  (p.Mu( x )   .* z.^(-1/p.eta_h-1) .*d).^(exp1) * ones(p.N,1)).^(1/exp1) ;

