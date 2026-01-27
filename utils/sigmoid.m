function y = sigmoid(x)
%SIGMOID Numerically stable sigmoid
y = 1 ./ (1 + exp(-max(min(x, 50), -50)));
end
