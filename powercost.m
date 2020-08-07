function p = powercost(h, s)
    p = 1.28 * (10^(50 * s/200) - 1)/h/h;
end
