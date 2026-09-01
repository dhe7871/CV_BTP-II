clear; clc; close all;
syms x y;


p0 = [15, 20]';
kappa0 = 10;
s = -1;

o = [0, -20]';
kappa0Obs = 0.005;

% alpha = 20;
% 
% fn = (x./alpha).^2 + (y./alpha).^2 ...
%    + 0.6*exp(-((x./alpha - 0.5).^2 + (y./alpha + 0.3).^2)) ...
%    + 0.4*exp(-((x./alpha + 0.7).^2 + (y./alpha - 0.6).^2)) ...
%    + 0.35*exp(-((x./alpha + 0.2).^2 + (y./alpha + 0.8).^2)) ...
%    + 0.25*exp(-((x./alpha - 0.9).^2 + (y./alpha - 0.1).^2)) ...
%    + 0.15*sin(3*x./alpha).*cos(2*y./alpha) ...
%    + 0.12*sin(5*y./alpha + 0.7) ...
%    + 0.10*cos(4*x./alpha - 1.1) ...
%    + 0.08*(x./alpha).*(y./alpha) - 1.8;
fn = x^2 - 4 * y^2 * (1 - y^2);
fm = matlabFuncTillDoubleDiff(fn);

r1 = 25;
r2 = 10;
[a, b] = meshgrid(linspace(-50, 50, 100), linspace(-50, 50, 100));
u = zeros(size(a)); v = zeros(size(a));


gn = (x - o(1))^2 + (y - o(2))^2 - r2^2;
gm = matlabFuncTillDoubleDiff(gn);
dg = calculateDerivatives(gm, p0);
rObs0 = dg(1);

kappaObs = kappa0Obs;
if(rObs0 < 0.001)
    kappaObs = kappaObs/abs(rObs0);
end


df = calculateDerivatives(fm, p0);
r0 = df(1);
kappa = kappa0;
if(r0 > 0.001)
    kappa = kappa/abs(r0);
end

for i = 1:size(a, 1)
    for j = 1:size(a, 2)
        p = [a(i, j), b(i, j)]';

        df = calculateDerivatives(fm, p);
        r = df(1);

        vc = [df(2);df(3);];
        vc = vc./sqrt(vc'*vc);
        vs = [df(3);-df(2);];
        vs = vs./sqrt(vs'*vs);
        vd = -tanh(kappa*r)*vc + s*sech(kappa*r)*vs;

        dg = calculateDerivatives(gm, p);
        rObs = dg(1);

        voc = [dg(2);dg(3);];
        voc = voc./sqrt(voc'*voc);
        vos = [dg(3);-dg(2);];
        vos = vos./sqrt(vos'*vos);
        vod = tanh(kappaObs*rObs)*voc + s*sech(kappaObs*rObs)*vos;

        if(norm(p - o) - 1.5*r2 > 0)
            wtf = sech(norm(p - o) - 1.5*r2);
        else
            wtf = 1;
        end
        u(i, j) = (1 - wtf) * vd(1)/norm(vd) + wtf * vod(1)/norm(vod);
        v(i, j) = (1 - wtf) * vd(2)/norm(vd) + wtf * vod(2)/norm(vod);
    end
end


figure;
hold on;
grid on; grid minor; axis equal;
fimplicit(fn, 2*[-50, 50, -50, 50], "LineStyle", "--", "LineWidth", 1.5, "Color", "#ff0055");
fimplicit(gn, 2*[-50, 50, -50, 50], "LineStyle", "--", "LineWidth", 1.5);

quiver(a, b, u, v);

scatter(u(abs(u) < 0.1 & abs(v) < 0.1), v(abs(u) < 0.1 & abs(v) < 0.1), "Marker", ".", "SizeData", 500);

xlabel("x[m]");
ylabel("y[m]");
hold off;



%Helper functions
%Calculating the MatlabFunctions till double derivatives
function fnm = matlabFuncTillDoubleDiff(fn)
syms x y;
fx = diff(fn, x);
fy = diff(fn, y);
fxx = diff(fx, x);
fyy = diff(fy, y);
fxy = diff(fx, y);

fm = matlabFunction(fn, "Vars", [x, y]);
fxm = matlabFunction(fx, "Vars", [x, y]);
fym = matlabFunction(fy, "Vars", [x, y]);
fxxm = matlabFunction(fxx, "Vars", [x, y]);
fyym = matlabFunction(fyy, "Vars", [x, y]);
fxym = matlabFunction(fxy, "Vars", [x, y]);

fnm = {fm, fxm, fym, fxxm, fyym, fxym};
end

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
