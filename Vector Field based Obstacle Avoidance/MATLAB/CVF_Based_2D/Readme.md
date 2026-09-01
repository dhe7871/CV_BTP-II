Previous Combining algorithm

    vd = (1 - wtf) .* vd./norm(vd) + wtf .* vdo./norm(vdo);
    chid = atan2(vd(2), vd(1));

which results in a combined vector field which may have **Singularities**.

But this new way of combining the fields, in which we are essentially **combining the only the angles**
makes the resulting Combined Vector Field, **Singularity free**.

    chig = atan2(normvd(2), normvd(1)); chio = atan2(normvdo(2), normvdo(1));
    chid = wrapToPi((1 - wtf) * chig + wtf * chio);
