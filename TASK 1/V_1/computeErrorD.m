function [outputArg1] = computeErrorD(up, coords, d, he)

    eps = [-1/sqrt(3) 1/sqrt(3)];
    ws = [1 1];
    B = 1/he *[-1 1];
    Fe = 0;
    for g=1:2
        N = 1/2 * [1-eps(g) 1+eps(g)]
        ue = subs(up, N * coords');
        uh = B * d;
        q = ue -uh;
        Fe = Fe + ws(g) * (q ^ 2);
    end
    outputArg1 = Fe;
end