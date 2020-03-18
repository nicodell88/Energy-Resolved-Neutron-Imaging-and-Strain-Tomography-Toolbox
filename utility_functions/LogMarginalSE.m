function [logpdf, gradlogpdf] = LogMarginalSE(logl,X,Y,g1,g2,sig_e)
% log marginal likelihood for the squared exponenetial covariance function

sig_f = 1;
l = exp(logl);

% the likelihood model
N = length(Y);
d = (X-X.').^2;
K = sig_f^2 * exp(-0.5*d / l^2).* ((g2 - g1) .*(g2.' - g1.'));
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
dKdl = (d.*K)/l^3;
dGPdl = 0.5 * alpha.' * dKdl * alpha - 0.5 * trace(Sy\(Sy.'\dKdl));

% add it all up
logpdf = -logGP;
gradlogpdf = -l*dGPdl;
end

