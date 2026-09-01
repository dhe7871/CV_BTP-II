clc; clear; close all;

dt = 0.001;
simTime = 75;


wpts = [
    100, -100;
    300, 600;
    200, 800;
    ];
nwpts = size(wpts, 1);

alpha1 = 0.6;
oCenter1 = wpts(1, :) + alpha1*(wpts(2, :) - wpts(1, :));
ro1 = 50;

alpha2 = 0.6;
oCenter2 = wpts(2, :) + alpha2*(wpts(3, :) - wpts(2, :));
ro2 = 20;


va = 20; %m/s
g = 9.81; %m/s^2
phiMax = pi/4;
psiDotMax = g*tan(phiMax)/va;

k1 = 5; 
k2 = 0.3;
k3 = 1.5;
k4 = 10/ro1;

p0 = [-200, 100]';
psi0 = pi/3;

figure;
hold on; grid on; grid minor; axis equal;
plot(wpts(:, 1), wpts(:, 2), "LineWidth", 1.5, "LineStyle", ":");

%obstacle
thetao = 0:0.01:2*pi;
xo = oCenter1(1) + ro1*cos(thetao);
yo = oCenter1(2) + ro1*sin(thetao);
plot(xo, yo, "LineWidth", 1.5, "LineStyle", "--");

% xo = oCenter2(1) + ro2*cos(thetao);
% yo = oCenter2(2) + ro2*sin(thetao);
% plot(xo, yo, "LineWidth", 1.5, "LineStyle", "--");


scatter(p0(1), p0(2), "LineWidth", 1.5);

state = zeros([6, simTime/dt + 1]);
state(:, 1) = [p0', psi0, 0, 0, 0]';

%line no. to follow
l = 1;
if(nwpts > 2)
    unitL1 = (wpts(2, :) - wpts(1, :))./norm(wpts(2, :) - wpts(1, :));
    unitL2 = (wpts(3, :) - wpts(2, :))./norm(wpts(3, :) - wpts(2, :));
    thetal = pi - cos(sum(unitL1.*unitL2));
    switchRadius = va/(psiDotMax*tan(thetal/2));
end

i = 2;
for t = dt:dt:simTime
    x = state(1, i-1); y = state(2, i - 1);
    psi = wrapToPi(state(3, i-1));

    if(nwpts > l + 1)
        if((x - wpts(l + 1, 1))^2 + (y - wpts(l + 1, 2))^2 - (switchRadius)^2 <= 0)
            l = l + 1;
            disp(l)
            if(nwpts > l + 1)
                unitL1 = (wpts(l + 1, :) - wpts(l, :))./norm(wpts(l + 1, :) - wpts(l, :));
                unitL2 = (wpts(l + 2, :) - wpts(l + 1, :))./norm(wpts(l + 2, :) - wpts(l + 1, :));
                thetal = pi - cos(sum(unitL1.*unitL2));
                switchRadius = va/(psiDotMax*tan(thetal/2));
            end
        end
    end
    x1 = wpts(l, 1); y1 = wpts(l, 2);
    x2 = wpts(l + 1, 1); y2 = wpts(l + 1, 2);

    unitLOS = [x2-x1, y2-y1]'./sqrt((x2 - x1)^2 + (y2 - y1)^2);
    thetaLOS = atan2(unitLOS(2), unitLOS(1));
    % disp(fprintf("thetaLOS: %f", thetaLOS))

    thetaP = atan2(y2 - y, x2 - x);
    
    d = cross([unitLOS', 0], [x - x1, y - y1, 0]);
    d = d(3);
    
    % disp(k1*(thetaP - psi));
    % disp(k2*(-d)*(thetaLOS - sign(d)*pi/2 - psi))
    psiDot = k1*wrapToPi(thetaP - psi) + k2*abs(d)*wrapToPi(thetaLOS - sign(d)*pi/2 - psi);
    % disp(psiDot)

    psiDot = max(min(psiDot, psiDotMax), -psiDotMax);
    
    state(4:6, i - 1) = [psiDot, thetaP, d]';
    
    psi = wrapToPi(psi + psiDot*dt);

    %Obstacle Avoidance Algorithm
    vo = oCenter1' - state(1:2, i - 1);
    vw = [x2, y2]' - state(1:2, i - 1);
    phivo = wrapToPi(atan2(vo(2), vo(1)) - atan2(vw(2), vw(1)));
    lvo = norm(vo);

    wtf = sech(k3*phivo)*sech(k4*(lvo - ro1));
    psi = psi*(1 - wtf) + (atan2(vw(2), vw(1)) + pi/2*W(phivo))*wtf;
    

    state(1:2, i) = state(1:2, i - 1) + [cos(psi), sin(psi)]'*va*dt;
    state(3, i) = psi;

    i = i + 1;
end

plot(state(1, :), state(2, :), "LineWidth", 1.5);
hold off;

figure;
hold on; grid on; grid minor;
plot(0:dt:simTime, state(3, :)*180/pi, "LineWidth", 1.5);
xlabel("Time[s]");
ylabel("\psi[deg]");
hold off;

figure;
hold on; grid on; grid minor;
plot(0:dt:simTime, state(4, :), "LineWidth", 1.5);
xlabel("Time[s]");
ylabel("\psi_D_o_t[rad/s]");
hold off;

figure;
hold on; grid on; grid minor;
plot(0:dt:simTime, state(5, :)*180/pi, "LineWidth", 1.5);
xlabel("Time[s]");
ylabel("\theta_P[deg]");
hold off;


function w = W(x)
    if(x == 0)
        w = 1;
    else
        w = -sign(x);
    end
end