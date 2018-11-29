function phi_inv = quantile(p)
    phi_inv = sqrt(2) * erfinv(2 * p - 1);
end
