function [err_alpha, delta, alpha, err_delta] = findD2(D2, p, A2, exp2, mu)

    g     = (mu .* A2 ./ D2).^exp2;
    delta = p.Delta(log(g));

    err_alpha = sum(delta.^exp2 .* g, 2) - 1;

    if nargout > 2

        alpha = delta.^exp2 .* g;

        err_delta = delta - 1 - 1/p.etaL - (1/p.thetaL - 1/p.etaL) .* alpha;

    end

end