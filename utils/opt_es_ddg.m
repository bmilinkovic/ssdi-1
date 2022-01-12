function [dopt,Lopt,converged,sig,iters,dhist] = opt_es_ddg(H,Lopt,maxiters,sig,ifac,nfac,tol,hist)

% Assumptions
%
% 1 - Lopt is orthonormal
% 2 - Residuals covariance matrix is identity

if isscalar(tol)
	stol = tol;
	dtol = tol;
	gtol = tol;
else
	stol = tol(1);
	dtol = tol(2);
	gtol = tol(3);
end

[n,m] = size(Lopt);

% Calculate dynamical dependence of initial projection

dopt  = trfun2dd(Lopt,H);
grad  = trfun2ddgrad(Lopt,H);               % dynamical dependence gradient
mgrad = sqrt(sum(grad(:).^2));              % magnitude of gradient vector

if hist
	dhist = zeros(maxiters,3);
	dhist(1,:) = [dopt sig mgrad];
else
	dhist = [];
end

% Optimise

converged = 0;
for iters = 2:maxiters

	% Move (hopefully) down gradient and orthonormalise

	grad  = trfun2ddgrad(Lopt,H);               % dynamical dependence gradient
	mgrad = sqrt(sum(grad(:).^2));              % magnitude of gradient vector
	Ltry = orthonormalise(Lopt-sig*grad/mgrad); % gradient descent

	% Calculate dynamical dependence of mutated projection

	dtry = trfun2dd(Ltry,H);

	% If dynamical dependence smaller, accept move

	if dtry < dopt
		Lopt = Ltry;
		dopt = dtry;
		sig  = ifac*sig;
	else
		sig  = nfac*sig;
	end

	if hist
		dhist(iters,:) = [dopt sig mgrad];
	end

	% Test convergence

	if sig < stol
		converged = 1;
		break
	end

	if dopt < dtol
		converged = 2;
		break
	end

	if mgrad < gtol
		converged = 3;
		break
	end

end

if hist
	dhist = dhist(1:iters,:);
end
