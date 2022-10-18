function [outputArg1] = computeError(u, coords, d)

    eps = [-1/sqrt(3) 1/sqrt(3)];
    ws = [1 1];
    Fe = 0;
    for g=1:2
        N = 1/2 * [1-eps(g) 1+eps(g)];
        ue = subs(u, N * coords');
        uh = N * d;
        q = ue -uh;
        Fe = Fe + ws(g) * (q ^ 2);
    end
    outputArg1 = Fe;
end