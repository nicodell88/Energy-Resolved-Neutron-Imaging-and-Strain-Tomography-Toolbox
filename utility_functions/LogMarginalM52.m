function [logpdf] = LogMarginalM52(logl,X,Y,g1,g2,sig_e)
sig_f = 1;
l = exp(logl);

% the likelihood model
N = length(Y);
d = sqrt((X-X.').^2);
d2 = (X-X.').^2;
K = sig_f^2 * (exp(-sqrt(5)*d/l) .* (1 + sqrt(5)*d/l + 5*d2/l^2/3)).* ((g2 - g1) .*(g2.' - g1.'));
Ky = K +  sig_e^2 * eye(N);

[Sy, flag] = chol(Ky);

while flag
    Ky = Ky - 2*min(eig(Ky))*eye(N);
   [Sy, flag] = chol(Ky); 
    warning('Chol is ill conditioned') 
end

v = Sy.'\Y;
logGP = -1/2*(v.'*v) - sum(log(diag(Sy))) - N/2*log(2*pi);

% compute gradients (there is an error in here somewhere
% alpha = Sy\v;
% grad1 = (sqrt(5)*d/l^2).*exp(-sqrt(5)*d/l);   % gradient of part 1
% grad2 = -sqrt(5)*d/l^2.*exp(-sqrt(5)*d/l) + (5*d2/l^3).*exp(-sqrt(5)*d/l); % gradient of part 2
% grad3 = -10*d2/3/l^3.*exp(-sqrt(5)*d/l) + 5*sqrt(5)*d2*d/3/l^4*exp(-sqrt(5)*d/l);
% dKdl = (grad1+grad2+grad3).* ((g2 - g1) .*(g2.' - g1.'));
% dGPdl = 0.5 * alpha.' * dKdl * alpha - 0.5 * trace(Sy\(Sy.'\dKdl));

% add it all up
logpdf = -logGP;
% gradlogpdf = -l*dGPdl;
end


