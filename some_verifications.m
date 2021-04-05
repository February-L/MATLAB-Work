clear;

M = 1.0; timelike = true;
r_h1_func = @(a) (M - sqrt(M.^2.0 - a.^2.0));
r_h2_func = @(a) (M + sqrt(M.^2.0 - a.^2.0));
r_ss2 = 2.0 * M;
n1 = 0; n_join1 = 0; n2 = 0; n_join2 = 0; n3 = 0; n_join3 = 0;
a_vals = 0.0: 0.01: 1.0; L_vals = 0.2: 0.2: 50.0;
container1 = zeros(length(a_vals), length(L_vals));
container2 = zeros(length(a_vals), length(L_vals));
container3 = zeros(length(a_vals), length(L_vals));
for i_a = 1: 1: length(a_vals)
    a = a_vals(i_a);
    r_h1 = r_h1_func(a); r_h2 = r_h2_func(a);
    for i_L = 1: 1: length(L_vals)
        L = L_vals(i_L);
        [r_ex_a, ex_p_a, ad_step_a] = effp_extreme_s([r_h2, 100.0], 10.0.^(-2), (L >= 0.0), timelike, L, M, a);
        the_i = find(ex_p_a == (-1.0 + 2.0 * double(L >= 0.0)));
        if (~isempty(r_ex_a) && ~isempty(the_i))
            n_join1 = n_join1 + 1;
            r0 = r_ss2 * (1 - a*eff_potential(r_ex_a(the_i), (L >= 0.0), timelike, L, M, a, false)/L);
            container1(i_a, i_L) = r_ex_a(the_i) - r0;
            if (r_ex_a(the_i) < r0)
                n1 = n1 + 1;
            end
        end
        [r_ex_b, ex_p_b, ad_step_b] = effp_extreme_s([r_h2, 100.0], 10.0.^(-2), (L <= 0.0), timelike, L, M, a);
        the_i = find(ex_p_b == (-1.0 + 2.0 * double(L <= 0.0)));
        if (~isempty(r_ex_b) && ~isempty(the_i))
            n_join2 = n_join2 + 1;
            r0 = r_ss2 * (1 - a*eff_potential(r_ex_b(the_i), (L <= 0.0), timelike, L, M, a, false)/L);
            container2(i_a, i_L) = r_ex_b(the_i) - r0;
            if (r_ex_b(the_i) < r0)
                n2 = n2 + 1;
            end
        end
        [r_ex_c, ex_p_c, ad_step_c] = effp_extreme_s([0.0, r_h1], 10.0.^(-4), (L <= 0.0), timelike, L, M, a);
        the_i = 1;
        if (~isempty(r_ex_c))
            n_join3 = n_join3 + 1;
            r0 = r_ss2 * (1 - a*eff_potential(r_ex_c(the_i), (L <= 0.0), timelike, L, M, a, false)/L);
            container3(i_a, i_L) = r_ex_c(the_i) - r0;
            if (r_ex_c(the_i) > r0)
                n3 = n3 + 1;
            end
        end
    end
end
fprintf("\njoin1: %g, n1 = %g\njoin2: %g, n2 = %g\njoin3: %g, n3 = %g\n", n_join1, n1, n_join2, n2, n_join3, n3);
