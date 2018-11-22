function phi_inv = quantile(p)
    phi_inv = sqrt(2) * 1/erf(2 * p - 1);
end
