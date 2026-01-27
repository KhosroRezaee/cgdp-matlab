function z = predict_logit(m, X)
%PREDICT_LOGIT raw logit
z = X*m.w + m.b;
end
