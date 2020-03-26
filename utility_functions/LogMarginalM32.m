
function [logpdf, gradlogpdf] = LogMarginalM32(logl,X,Y,g1,g2,sig_e)
sig_f = 1;
l = exp(logl);
X = scaleInput(X);

% the likelihood model
N = length(Y);
d = sqrt((X-X.').^2);
K = sig_f^2 * exp(-sqrt(3)*d/l) .* (1 + sqrt(3)*d/l).* ((g2 - g1) .*(g2.' - g1.'));
Ky = K +  sig_e^2 * eye(N);

[Sy, flag] = chol(Ky);

while flag
    Ky = Ky + 0.00001*eye(N);
   [Sy, flag] = chol(Ky); 
    warning('Chol is ill conditioned') 
end

v = Sy.'\Y;
logGP = -1/2*(v.'*v) - sum(log(diag(Sy))) - N/2*log(2*pi);

% compute gradients
alpha = Sy\v;
dKdl = sig_f^2 * 3 * d.^2 / l^3 .*exp(-sqrt(3)*d/l).* ((g2 - g1) .*(g2.' - g1.'));
dGPdl = 0.5 * alpha.' * dKdl * alpha - 0.5 * trace(Sy\(Sy.'\dKdl));

% add it all up
logpdf = -logGP;
gradlogpdf = -l*dGPdl;
end
