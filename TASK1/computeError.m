function [outputArg1] = computeError(u, coords, d)
    syms x f(x)
    eps = [-1/sqrt(3) 1/sqrt(3)];
    ws = [1 1];
    l = coords(2) - coords(1);
    m = (l)/(d(2)-d(1));
    f(x) = m*(x-coords(1)) + d(1);
    eps = (eps / 2 + 0.5 )*l + coords(1);
    Fe = 0;
    for g=1:2
        uh = f(eps(g));
        ue = u(eps(g));
        q = (ue -uh)^2;
        Fe = Fe + ws(g) * q;
    end
    outputArg1 = Fe;
end