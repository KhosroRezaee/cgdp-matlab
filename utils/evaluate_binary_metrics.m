function m = evaluate_binary_metrics(y, p)
%EVALUATE_BINARY_METRICS AUROC + AUPRC + ECE-ish quick stats
y = y(:);
p = p(:);

[m.auroc] = auroc(y,p);
[m.auprc] = auprc(y,p);

% simple calibration error (10 bins)
[m.ece] = ece10(y,p);

end

function a = auroc(y,p)
[~,ord] = sort(p,'descend');
y = y(ord);
tp = cumsum(y==1);
fp = cumsum(y==0);
tpr = tp / max(1,sum(y==1));
fpr = fp / max(1,sum(y==0));
a = trapz(fpr, tpr);
end

function a = auprc(y,p)
[~,ord] = sort(p,'descend');
y = y(ord);
tp = cumsum(y==1);
fp = cumsum(y==0);
prec = tp ./ max(1,(tp+fp));
rec = tp / max(1,sum(y==1));
a = trapz(rec, prec);
end

function e = ece10(y,p)
bins = linspace(0,1,11);
e = 0;
n = numel(y);
for i=1:10
    mask = (p>=bins(i) & p<bins(i+1));
    if i==10
        mask = (p>=bins(i) & p<=bins(i+1));
    end
    if ~any(mask), continue; end
    acc = mean(y(mask));
    conf = mean(p(mask));
    e = e + (sum(mask)/n)*abs(acc-conf);
end
end
