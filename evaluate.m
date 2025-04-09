function f = evaluate(xx)
    load('data.mat')
    xx = xx';
    xx = [(trackData(:,5)-trackData(:,3)).*xx+trackData(:,3),(trackData(:,6)-trackData(:,4)).*xx+trackData(:,4)];
    xx = [xx;xx(1,:)];
    f = 0;
    for i = 1:nseg
        f = f+sqrt(sum((xx(i,:)-xx(i+1,:)).^2));
    end
end