function m = train_nnpu_logreg(X, y, pi_prior, cfg)
%TRAIN_NNPU_LOGREG nnPU logistic regression with minibatch SGD.
% y: 1 for positives, 0 for unlabeled/background.

Xp = X(y==1,:);
Xu = X(y==0,:);

[d1, D] = size(Xp); %#ok<ASGLU>
[d2, ~] = size(Xu); %#ok<ASGLU>

% Initialize
rng(cfg.model.seed);
w = 0.01 * randn(D,1);
b = 0.0;

lr = cfg.model.lr;
bs = cfg.model.batch_size;
epochs = cfg.model.epochs;
wd = cfg.model.weight_decay;
beta = cfg.model.nnpu_beta;

for ep=1:epochs
    % minibatch sampling
    ip = randi(size(Xp,1), [bs,1]);
    iu = randi(size(Xu,1), [bs,1]);
    Xpb = Xp(ip,:);
    Xub = Xu(iu,:);

    fp = Xpb*w + b;
    fu = Xub*w + b;

    % logistic losses
    % positive loss: log(1+exp(-f))
    lp_pos = mean(log1pexp(-fp));
    % negative loss: log(1+exp(f))
    lp_neg = mean(log1pexp(fp));
    lu_neg = mean(log1pexp(fu));

    % nnPU risk
    risk_neg = lu_neg - pi_prior * lp_neg;
    if risk_neg < -beta
        risk = pi_prior*lp_pos - beta;
        % gradient from positive part only (stop gradient on negative)
        [gw, gb] = grad_pos_part(Xpb, fp, pi_prior);
    else
        risk = pi_prior*lp_pos + risk_neg;
        [gw_pos, gb_pos] = grad_pos_part(Xpb, fp, pi_prior);
        [gw_neg, gb_neg] = grad_neg_part(Xpb, fp, Xub, fu, pi_prior);
        gw = gw_pos + gw_neg;
        gb = gb_pos + gb_neg;
    end

    % weight decay
    gw = gw + wd*w;

    % update
    w = w - lr*gw;
    b = b - lr*gb;

    % mild LR decay
    if mod(ep, 500)==0
        lr = lr * 0.7;
    end

    if mod(ep, 500)==0
        fprintf('[nnPU] ep=%d risk=%.4f\n', ep, risk);
    end
end

m = struct('w', w, 'b', b, 'pi', pi_prior);
end

function y = log1pexp(x)
% stable log(1+exp(x))
y = max(x,0) + log(1 + exp(-abs(x)));
end

function [gw, gb] = grad_pos_part(Xp, fp, pi_prior)
% gradient of pi*E_p[ log(1+exp(-f)) ]
sp = sigmoid(-fp);  % d/d f of log(1+exp(-f)) = -sigmoid(-f)
% Actually derivative wrt f: -1/(1+exp(f)) = -sigmoid(-f)
df = -sp;
gb = pi_prior * mean(df);
gw = pi_prior * (Xp' * df) / size(Xp,1);
end

function [gw, gb] = grad_neg_part(Xp, fp, Xu, fu, pi_prior)
% gradient of (E_u[log(1+exp(f))] - pi*E_p[log(1+exp(f))])
su = sigmoid(fu); % derivative of log(1+exp(f)) wrt f = sigmoid(f)
sp = sigmoid(fp);

df_u = su;
df_p = sp;

gb = mean(df_u) - pi_prior*mean(df_p);
gw = (Xu' * df_u)/size(Xu,1) - pi_prior*(Xp' * df_p)/size(Xp,1);
end
