function D = iss2dd(L,A,C,K)

% Calculate dynamical dependence of projection L for
% innovations-form state-space model with parameters A,C,K.
%
% NOTE 1: assumes uncorrelated residuals
% NOTE 2: projection L MUST be orthonormal!!!

% Calculate residuals covariance matrix V of projected model (solve DARE)

[~,V,rep] = mdare(A,L'*C,K*K',[],K*L); % mdare from MVGC2

if rep < 0  || rep > 1e-08 % DARE failed
	D = NaN;
	return
end

% D = log-determinant of residuals covariance matrix V

[R,p] = chol(V);
if p == 0
	D = 2*sum(log(diag(R)));
else
	D = NaN; % fail: V not positive-definite
end
