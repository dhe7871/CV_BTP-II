function verdict = obsPosWRTPath(params) %output: "in", "out", "intersect"
    fMatFunc = params.fm; gMatFunc = params.gm;
    fm = fMatFunc{1};

    delta = params.obs.delta;

    %coordinates of point "delta" distance away from the Obstacle
    % only using 1% of the total points for determination
    curveGPt = params.curveG(:, 1:100:end);
    posMap = zeros(size(curveGPt, 2), 1);

    parfor k = 1:size(curveGPt, 2)
        dg = calculateDerivatives(gMatFunc, curveGPt(:, k));

        gradG = [dg(2), dg(3)]';
        normGradG = sqrt(gradG(1)^2 + gradG(2)^2);
        if(normGradG < 1e-6)
            continue;
        end
        
        extCurveGPt = curveGPt(:, k) + (gradG ./ normGradG) .* delta;

        r = fm(extCurveGPt(1), extCurveGPt(2));
        posMap(k) = (r > 0) - (r < 0);
    end
    if(~all(posMap) || (~all(posMap + 1) && ~all(posMap - 1)))
        verdict = "intersect";
        return;
    elseif(~any(posMap - 1))
        verdict = "out";
        return;
    elseif(~any(posMap + 1))
        verdict = "in";
        return;
    end
    verdict = "";
end