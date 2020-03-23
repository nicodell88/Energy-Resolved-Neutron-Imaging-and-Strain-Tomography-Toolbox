function [Phi,Phi_T,SLambda,lambdas,dPhi_T] = hilbert_approxM32(l,sig_f,m,L,x_test,x_obs)
%   l: length scales
%   sig_f: prior std
%   m: number of basis functions
%   x_test: test points to estimate at
%   x_obs: observation points
lambdas = pi * [0:m-1]/2/L;
Phi =  1/sqrt(L)*sin(bsxfun(@times,x_obs+L,lambdas));
Phi_T = 1/sqrt(L)*sin(bsxfun(@times,x_test+L,lambdas));
dPhi_T = lambdas .* cos(bsxfun(@times,x_test+L,lambdas))/sqrt(L);
SLambda = sig_f^2*4 * sqrt(3^3) /l^3 ./(3/l^2 + lambdas.^2).^2;
end