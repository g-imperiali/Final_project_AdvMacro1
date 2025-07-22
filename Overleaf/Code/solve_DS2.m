function err = solve_DS2(Ds2, p, z, mu,exp2)
g = (mu.*z.^(1/p.eta-1)./Ds2).^(exp2);
g = min(max(g, p.gmin), p.gmax);
err = Ds2 - (  (p.D( g )   .* z.^(1/p.eta-1) .*mu ).^(exp2) * ones(p.N,1)).^(1/exp2) ;

