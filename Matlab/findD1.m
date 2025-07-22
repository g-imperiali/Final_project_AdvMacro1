function [err_omega, mu, omega, err_mu]  = findD1(D1, p, A1, exp1, delta)

    f    = (delta .* A1 ./ D1).^exp1; 
    mu   = p.Mu(log(f));
    err_omega  = sum(mu.^exp1 .* f, 2) - 1;
    
    if nargout > 2
    
        omega = mu.^exp1 .* f;
        err_mu = 1 - 1./mu - 1/p.eta - (1/p.theta - 1/p.eta) .* omega;

    end
    
end
