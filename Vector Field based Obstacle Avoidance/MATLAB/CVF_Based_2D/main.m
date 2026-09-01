%% CVF Path Following and CVF Obstacle
clc; clear; close all;
projectPath = "C:\Users\Dheeraj\Academics\BTP\1-MATLAB\Vector Field based Obstacle Avoidance\CVF_Based\2D\v1";
addpath(genpath(projectPath + "\lib"));

simTime = 300; params.simTime = simTime; %s
dt = 0.01; params.dt = dt;



choice = questdlg("Do you want to generate the curve points?", "Confirm", "Yes", "No", "No");
usePrevCurves = (choice == "No");

params.usePrevCurves = usePrevCurves; %Generate the new curve points or just use the previously saved points

syms x y;
fn = (x - 10)^2 + (y - 15)^2 - 400^2; kappa0 = 20;
% fn = (y - 100*sin(x/100))/100; kappa0 = 10;
% fn = (x/200)^2  + (y/150)^2 - 1; kappa0 = 8;
% fn = y - 100*cos(x/100); kappa0 = 5;
% k = 400;
% fn = (x./k).^2 + (y./k).^2 ...
%     + 0.6*exp(-((x./k - 0.5).^2 + (y./k + 0.3).^2)) ...
%     + 0.4*exp(-((x./k + 0.7).^2 + (y./k - 0.6).^2)) ...
%     + 0.35*exp(-((x./k + 0.2).^2 + (y./k + 0.8).^2)) ...
%     + 0.25*exp(-((x./k - 0.9).^2 + (y./k - 0.1).^2)) ...
%     + 0.15*sin(3*x./k).*cos(2*y./k) ...
%     + 0.12*sin(5*y./k + 0.7) ...
%     + 0.10*cos(4*x./k - 1.1) ...
%     + 0.08*(x./k).*(y./k) - 1.8; kappa0 = 20;

params.fn = fn;

%Differentiated matlab functions
fm = matlabFuncTillDoubleDiff(fn); params.fm = fm;

%UAV Specifications
s = 1; params.uav.s = s;                    %direction for rotational VF
va = 15; params.uav.va = va;                %m/s
w = [6, 8]'; params.uav.w = w;              %m/s

p0 = [-600, 600]'; params.uav.p0 = p0;      %m
chi0 = -pi/2; params.uav.chi0 = chi0;     %radians

g = 9.81;                                   %m/s^2
phiMax = pi/4; params.uav.phiMax = phiMax;  %Maximum Bank Angle (radians)
chiDotMax = g*tan(phiMax)/va; params.uav.chiDotMax = chiDotMax;



%Path Following Gains/Parameters
%"kappa0" defined earlier with the curve
kChi = 1.5; params.uav.kChi = kChi;

df = calculateDerivatives(fm, p0);
r0 = df(1);
kappa = kappa0;
if(abs(r0) > 0.001)
    kappa = kappa/abs(r0);
end
params.uav.kappa = kappa;

%Obstacle Parameters
o = [-800, 0];
so = -s; params.obs.s = so;

%Obstacle to Avoid
% gn = (x - o(1))^2 + (y - o(2))^2 - 100^2; kappa0Obs = 5;
gn = ((x - o(1))/100)^2 + ((y - o(2))/50)^2 - 1; kappa0Obs = 1;
% gn = y + x - 100; kappa0Obs = 0.01;
% gn = ((x - o(1))/100)^2 + ((y - o(2))/50)^2 + 0.2*sin(8*atan(y, x)) + 0.1*cos(12*atan(y, x)) - 1; kappa0Obs = 1;
% k = 100;
% gn = ((x - o(1))./k).^2 + ((y - o(2))./k).^2 ...
%     + 0.6*exp(-(((x - o(1))./k - 0.5).^2 + ((y - o(2))./k + 0.3).^2)) ...
%     + 0.4*exp(-(((x - o(1))./k + 0.7).^2 + ((y - o(2))./k - 0.6).^2)) ...
%     + 0.35*exp(-(((x - o(1))./k + 0.2).^2 + ((y - o(2))./k + 0.8).^2)) ...
%     + 0.25*exp(-(((x - o(1))./k - 0.9).^2 + ((y - o(2))./k - 0.1).^2)) ...
%     + 0.15*sin(3*(x - o(1))./k).*cos(2*(y - o(2))./k) ...
%     + 0.12*sin(5*(y - o(2))./k + 0.7) ...
%     + 0.10*cos(4*(x - o(1))./k - 1.1) ...
%     + 0.08*((x - o(1))./k).*((y - o(2))./k) - 1.8; kappa0Obs = 0.8;

params.gn = gn;

% [x1, y1, x2, y2] = normLenBwCurves(fn, gn);
% x1, y1
% x2, y2

gm = matlabFuncTillDoubleDiff(gn); params.gm = gm;

dg = calculateDerivatives(gm, p0);
rObs0 = dg(1);

kappaObs = kappa0Obs;
if(abs(rObs0) > 0.001)
    kappaObs = kappaObs/abs(rObs0);
end
params.obs.kappaObs = kappaObs;

%Obstacle Distance Margin (Safety Margin)
delta = 2*va/chiDotMax; params.obs.delta = delta;

%Smooth Field Switch factor (can't be zero)
%(on increasing switching distance increases)
alpha = 25; params.obs.alpha = alpha;


%Generating curve points
if(params.usePrevCurves && exist(fullfile(projectPath, "temp\curves.mat"), "file"))
    load(fullfile(projectPath, "temp\curves.mat"));
    params.curveF = curveF; params.curveG = curveG;
else
    curveF = curvePoints2D(fm{1}, p0, 40000, "Direction", "forward"); params.curveF = curveF;
    curveG = curvePoints2D(gm{1}, p0, 50000, "Direction", "both"); params.curveG = curveG;
    save(fullfile(projectPath, "temp\curves.mat"), "curveF", "curveG");
end

%find CurveParams i.e., Center of Mass, xLim, yLim etc.
params = find2DCurveParams(params);
params.obsPos = obsPosWRTPath(params);

xComO = params.obs.com(1); yComO = params.obs.com(2);
xClose = params.uav.minObsComToPath.pathPt(1); yClose = params.uav.minObsComToPath.pathPt(2);
obsPos = params.obsPos;

[a, b] = meshgrid(linspace(params.xLim(1), params.xLim(2), 100), linspace(params.yLim(1), params.yLim(2), 100));
[u, v] = calculateVF(params, a, b);



%State Initializations
%state = [x, y, dx, dy, vg, chi, chid, dchid, dchic, r = f(x, y)]
state = zeros(10, simTime/dt + 1);
state(:, 1) = [p0', zeros(1, 2), 0, chi0, 0, zeros(1, 3)]';

writerObj = VideoWriter(fullfile(projectPath,"\media\animation.avi"));
open(writerObj);

figure;
hold on;
grid on; grid minor;
axis equal;

plot(curveF(1, :), curveF(2, :), "LineStyle", "--", "LineWidth", 1.5, "Color", "#ff0055");
plot(curveG(1, :), curveG(2, :), "LineStyle", "-", "LineWidth", 1.5)
% fimplicit(fn, [xLim, yLim], "LineStyle", "--", "LineWidth", 1.5, "Color", "#ff0055");
% fimplicit(gn, [xLim, yLim], "LineStyle", "-", "LineWidth", 1.5);
quiver(a, b, u, v);

scatter(p0(1), p0(2), "LineWidth", 1.5, "Marker", "o", "SizeData", 100);

h1 = plot(NaN, NaN, "LineWidth", 1.5, "Color", "#0099aa");
h2 = plot(NaN, NaN, "LineWidth", 1.5, "Marker", ">", "MarkerSize", 8, "Color", "#0099aa");

legend("Curve to follow","Obstacle to Avoid", "Vector Fields", "Initial Position",...
    "Actual trajectory", "Current Position");
title("2D CVF Path following Guidance and Obstacle Avoidance");



i = 2;
progressbar("Initialize", [i, size(state, 2)]);
for t = dt:dt:simTime
    p = state(1:2, i - 1);
    df = calculateDerivatives(fm, p);

    chi = state(6, i - 1);
    vg = (w(1)*cos(chi) + w(2)*sin(chi)) + sqrt(va^2 - (w(1)*sin(chi) - w(2)*cos(chi))^2);
    r = df(1);
    state(10, i - 1) = r;
    
    guessPt = p;
    [xn, yn] = normLen(gm{1}, gm{2}, gm{3}, guessPt);
    nlen = norm(p - [xn, yn]');
    if(nlen - delta > 0)
        wtf = sech((nlen - delta)/alpha);
    else
        wtf = 1;
    end

    vc = [df(2); df(3)];
    vc = vc ./ sqrt(vc' * vc);
    vs = [df(3); -df(2)];
    vs = vs ./ sqrt(vs' * vs);
    vp = -tanh(kappa*r)*vc + s*sech(kappa*r)*vs;
    normvp = vp./sqrt(vp'*vp);

    %obstacle Avoidance
    dg = calculateDerivatives(gm, p);
    rObs = dg(1);

    % %Experimental Code, 'Comment' or 'Delete' below code after experiment
    % % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % vco = [dg(2);dg(3);];
    % normvco = vco ./ norm(vco);
    % 
    % vso = [dg(3); -dg(2)];
    % normvso = vso ./ norm(vso);
    % 
    % vo = sech(kappaObs * rObs) .* normvco;
    % normvo = vo ./ norm(vo);
    % 
    % vp90o = normvp - (normvp' * normvo) .* normvo;
    % 
    % vo = wtf * normvo + (1 - wtf) * vp90o;
    % 
    % normvo = vo ./ norm(vo);
    % %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    %Uncomment Below Code when experiment is over.......
    % ---------------------------------------------------------------
    vco = [dg(2); dg(3)];
    vco = vco ./ sqrt(vco' * vco);
    vso = [dg(3); -dg(2)];
    vso = vso ./ sqrt(vso' * vso);
    vo = sech(kappaObs*rObs)*vco + so*tanh(kappaObs*rObs)*vso;
    normvo = vo ./ norm(vo);
    % -----------------------------------------------------------------

    chig = atan2(normvp(2), normvp(1));
    chio = atan2(normvo(2), normvo(1));

    if(obsPos == "out" || obsPos == "in")
        zeta = wrapTo2Pi(atan2(yComO - yClose, xComO - xClose));
        chig = zeta + wrapTo2Pi(chig - zeta);
        chio = zeta + wrapTo2Pi(chio - zeta);
    elseif(obsPos == "intersect")
        if(r >= 0)
            if(so == 1 && chig < chio)
                chig = chig + 2*pi;
            elseif (so == -1 && chio < chig)
                chio = chio + 2*pi;
            end
        else
            if(so == -1 && chig < chio)
                chig = chig + 2*pi;
            elseif (so == 1 && chio < chig)
                chio = chio + 2*pi;
            end
        end
        
    end

    chid = wrapToPi((1 - wtf) * chig ...
        + wtf * chio);


    %previous algorithm to combine
    chie = wrapToPi(chi - chid);
    
    % Saturated PD Controller
    omega1 = -kChi*chie;
    chidDot = omega1;
    chicDot = max(min(omega1, chiDotMax), -chiDotMax);
    
    state(5, i - 1) = vg;
    state(7:9, i - 1) = [chid, chidDot, chicDot]';

    chi = wrapToPi(chi + chicDot*dt);

    xDot = vg*cos(chi);
    yDot = vg*sin(chi);
    
    state(1:2, i) = p + dt*[xDot, yDot]';
    state(3:4, i) = [xDot, yDot]';
    state(6, i) = chi;

    if(~mod(i, 100))
        set(h1, "XData", state(1, 1:i)', "YData",state(2, 1:i)');
        set(h2, "XData", state(1, i), "YData", state(2, i));

        drawnow limitrate;
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end

    progressbar(i);
    i = i + 1;
end

hold off; %End of first figure
close(writerObj);

plotFigures(params, state);



%Helper functions
function [x1n, y1n, x2n, y2n] = normLenBwCurves(fn, gn)
    syms x y x1 y1 x2 y2;

    f1 = subs(fn, [x y], [x1, y1]);
    f2 = subs(gn, [x y], [x2, y2]);

    df1dx1 = diff(f1, x1); df1dy1 = diff(f1, y1);
    df2dx2 = diff(f2, x2); df2dy2 = diff(f2, y2);

    eq1 = f1 == 0;
    eq2 = f2 == 0;
    eq3 = df1dx1 * df2dy2 - df1dy1 * df2dx2 == 0;
    eq4 = (y2 - y1) * df1dx1 - (x2 - x1) * df1dy1 == 0;

    sol = solve([eq1, eq2, eq3, eq4], [x1, y1, x2, y2]);

    x1 = double(sol.x1); 
    y1 = double(sol.y1); 
    x1n = x1((imag(x1) == 0) & (imag(y1) == 0));
    y1n = y1((imag(x1) == 0) & (imag(y1) == 0));

    x2 = double(sol.x2);
    y2 = double(sol.y2);
    x2n = x2((imag(x2) == 0) & (imag(y2) == 0));
    y2n = y2((imag(x2) == 0) & (imag(y2) == 0));
end