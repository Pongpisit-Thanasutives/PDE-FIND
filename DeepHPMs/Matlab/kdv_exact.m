function sol = kdv_exact(x,c,a)
    sol = (c/2)*sech(sqrt(c)/2*(x-a))^2;
end
