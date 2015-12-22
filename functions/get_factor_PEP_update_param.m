function E = get_factor_PEP_update_param(dist_sqr, beta_sr, beta_rd, g, sigma_sqr_d, sigma_sqr_r, x)

tmp_1 = x * sigma_sqr_r + beta_sr * dist_sqr;
tmp_2 = x * sigma_sqr_d ./ (g ^ 2 * beta_rd * tmp_1);
E = x * sigma_sqr_r ./ tmp_1 + beta_sr * dist_sqr .* tmp_2 ./tmp_1 .* exp(tmp_2) .* expint(tmp_2); % Note the difference between the exponential integral function definition in the reference and Matlab

