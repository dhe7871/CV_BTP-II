%Calculating numerical values of all the derivatives at some position 'p'
function df = calculateDerivatives(fm, p)
    a = p(1);
    b = p(2);
    
    f = fm{1};
    fx = fm{2};
    fy = fm{3};
    fxx = fm{4};
    fyy = fm{5};
    fxy = fm{6};
    
    df = [f(a, b), fx(a, b), fy(a, b), fxx(a, b), fyy(a, b), fxy(a, b)];
end
