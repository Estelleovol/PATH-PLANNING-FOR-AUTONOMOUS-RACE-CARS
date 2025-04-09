function f = evaluate2(xx)
    load('data.mat')
    xx = xx';
    xx = [(trackData(:,5)-trackData(:,3)).*xx+trackData(:,3),(trackData(:,6)-trackData(:,4)).*xx+trackData(:,4)];
    xx = [xx;xx(1,:)];
    dist = zeros(nseg,1);
    f = 0;
    for i = 1:nseg
        dist(i) = sqrt(sum((xx(i,:)-xx(i+1,:)).^2));
    end
    xresSP = xx(:,1);
    yresSP = xx(:,2);
    dx = gradient(xresSP);%速度
    d2x = gradient(dx);%加速度
    dy = gradient(yresSP);
    d2y = gradient(dy);
    k = sum(abs(dx.*d2y-dy.*d2x)./((dx.^2+dy.^2).^(3/2)));%曲率

    v_origin = sqrt(dx.^2+dy.^2);
    a_origin = gradient(v_origin);
    a_origin_min = min(a_origin);
    a_origin_max = max(a_origin);

    scale_v = vmax/max(v_origin);
    scale_a = min(amax/a_origin_max,amin/a_origin_min);
    scale = min(scale_a,scale_v);
    v = v_origin*scale;
    f = sum(dist./v(2:end));
end