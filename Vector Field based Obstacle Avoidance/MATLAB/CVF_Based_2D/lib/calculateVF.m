%Calculate 2D VF
function [u, v] = calculateVF(params, a, b)
    if(isprop(params, "fm"))
        fm = params.fm;
    else
        fm = matlabFuncTillDoubleDiff(params.fn);
        params.fm = fm;
    end

    if(isprop(params, "gm"))
        gm = params.gm;
    else
        gm = matlabFuncTillDoubleDiff(params.gn);
        params.gm = gm;
    end
    kappa = params.uav.kappa; kappaObs = params.obs.kappaObs;
    delta = params.obs.delta; 
    alpha = params.obs.alpha;
    % alpha = 100;
    s = params.uav.s; so = params.obs.s;

    xComO = params.obs.com(1); yComO = params.obs.com(2);
    fprintf("\n\nObstactle COM: [%0.2f, %0.2f] m\n", xComO, yComO);

    xClose = params.uav.minObsComToPath.pathPt(1); yClose = params.uav.minObsComToPath.pathPt(2);
    obsPos = params.obsPos;

    zeta = wrapTo2Pi(atan2(yComO - yClose, xComO - xClose));
    fprintf("Bifurcating Angle: %0.2f°\n\n", zeta*180/pi);

    u = zeros(numel(a), 1); v = zeros(numel(a), 1);
    parfor k = 1:(size(a, 1) * size(a, 2))
        i = ceil(k / size(a, 2));
        j = mod(k - 1, size(a, 2)) + 1;
    
        p = [a(i, j), b(i, j)]';
    
        df = calculateDerivatives(fm, p);
        r = df(1);
        
        dg = calculateDerivatives(gm, p);
        rObs = dg(1);
        if(rObs < 0)
            continue;
        end

        guessPt = p;
        [xn, yn] = normLen(gm{1}, gm{2}, gm{3}, guessPt);
        nlen = norm(p - [xn, yn]');
        if(nlen - delta > 0)
            wtf = sech((nlen - delta)/alpha);
        else
            wtf = 1;
        end
    
        % %Experimental Code, 'Comment' or 'Delete' below code after experiment
        % % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % vc = [df(2);df(3);];
        % vc = vc./sqrt(vc'*vc);
        % vs = [df(3);-df(2);];
        % vs = vs./sqrt(vs'*vs);
        % vp = -tanh(kappa*r)*vc + s*sech(kappa*r)*vs;
        % normvp = vp ./ norm(vp);
        % 
        % vco = [dg(2);dg(3)];
        % normvco = vco ./ norm(vco);
        % 
        % %projection vector of 'normvp' on to perpendicular basis of 'normvo'
        % gamma = normvp' * normvco;
        % vp0vco = gamma .* normvco;
        % vp90vco = normvp - vp0vco;
        % 
        % wtf1 = wtf * 1 + (1 - wtf) * gamma; 
        % 
        % vo = wtf1 * normvco + (1 - wtf) * vp90vco;
        % 
        % normvo = vo ./ norm(vo);
        % 
        % %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




        %UNCOMMENT, after experiment
        %-----------------------------------------------------------------
        %Global Path Following Vector Field
        vc = [df(2);df(3);];
        vc = vc./sqrt(vc'*vc);
        vs = [df(3);-df(2);];
        vs = vs./sqrt(vs'*vs);
        vp = -tanh(kappa*r)*vc + s*sech(kappa*r)*vs;
        normvp = vp ./ norm(vp);

        %Obstacle Vector Field
        voc = [dg(2);dg(3);];
        voc = voc./sqrt(voc'*voc);
        vos = [dg(3);-dg(2);];
        vos = vos./sqrt(vos'*vos);

        vo = sech(kappaObs*rObs)*voc + so*tanh(kappaObs*rObs)*vos;
        normvo = vo ./ norm(vo);
        %-----------------------------------------------------------------

        chig = wrapTo2Pi(atan2(normvp(2), normvp(1)));
        chio = wrapTo2Pi(atan2(normvo(2), normvo(1)));

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
        % elseif(obsPos == "in")
        %     if(so == -1 && chig < chio)
        %         chig = chig + 2*pi;
        %     elseif (so == 1 && chio < chig)
        %         chio = chio + 2*pi;
        %     end
        end



        chid = wrapToPi((1 - wtf) * chig ...
            + wtf * chio);

        u(k) = cos(chid);
        v(k) = sin(chid);

        % u(k) = normvo(1);
        % v(k) = normvo(2);
    end
    
    u = reshape(u, flip(size(a)))';
    v = reshape(v, flip(size(a)))';
end