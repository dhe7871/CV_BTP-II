function params = find2DCurveParams(params)
    curveF = params.curveF; curveG = params.curveG;
    p0 = params.uav.p0;

    %Finding Curve Params
    params.uav.xMin = min(curveF(1, :)); params.uav.xMax = max(curveF(1, :));
    params.uav.yMin = min(curveF(2, :)); params.uav.yMax = max(curveF(2, :));
    
    params.obs.xMin = min(curveG(1, :)); params.obs.xMax = max(curveG(1, :));
    params.obs.yMin = min(curveG(2, :)); params.obs.yMax = max(curveG(2, :));
    
    params.xMin = min([params.uav.xMin, params.obs.xMin, p0(1)]);
    params.xMax = max([params.uav.xMax, params.obs.xMax, p0(1)]);
    
    params.yMin = min([params.uav.yMin, params.obs.yMin, p0(2)]);
    params.yMax = max([params.uav.yMax, params.obs.yMax, p0(2)]);

    params.uav.xLim = [
        params.uav.xMin - 0.2 * (params.uav.xMax - params.uav.xMin), ...
        params.uav.xMax + 0.2 * (params.uav.xMax - params.uav.xMin)
        ];
    params.uav.yLim = [
        params.uav.yMin - 0.2 * (params.uav.yMax - params.uav.yMin), ...
        params.uav.yMax + 0.2 * (params.uav.yMax - params.uav.yMin)
        ];

    params.obs.xLim = [
        params.obs.xMin - 0.2 * (params.obs.xMax - params.obs.xMin), ...
        params.obs.xMax + 0.2 * (params.obs.xMax - params.obs.xMin)
        ];
    params.obs.yLim = [
        params.obs.yMin - 0.2 * (params.obs.yMax - params.obs.yMin), ...
        params.obs.yMax + 0.2 * (params.obs.yMax - params.obs.yMin)
        ];

    params.xLim = [
        params.xMin - 0.2 * (params.xMax - params.xMin), ...
        params.xMax + 0.2 * (params.xMax - params.xMin)
        ];
    params.yLim = [
        params.yMin - 0.2 * (params.yMax - params.yMin), ...
        params.yMax + 0.2 * (params.yMax - params.yMin)
        ];

    params.uav.com = sum(curveF, 2) ./ size(curveF, 2);
    params.obs.com = sum(curveG, 2) ./ size(curveG, 2);
    
    guessPt = params.obs.com;
    [x, y] = normLen(params.fm{1}, params.fm{2}, params.fm{3}, guessPt); 
    params.uav.minObsComToPath.pathPt = [x, y]';
    params.uav.minObsComToPath.nlen = norm(params.obs.com - [x, y]');

    guessPt = params.uav.com;
    [x, y] = normLen(params.gm{1}, params.gm{2}, params.gm{3}, guessPt); 
    params.obs.minUavComToObs.obsPt = [x, y]';
    params.obs.minUavComToObs.nlen = norm(params.uav.com - [x, y]');
end