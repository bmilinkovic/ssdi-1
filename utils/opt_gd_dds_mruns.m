function [dopt,Lopt,conv,iopt,sopt,cput,ohist] = opt_gd_dds_mruns(H,L0,niters,sig0,gdls,gdtol,hist,pp)

if nargin < 9 || isempty(pp), pp = false; end

nruns = size(L0,3);

dopt = zeros(1,nruns);
Lopt = zeros(size(L0));
conv = false(1,nruns);
iopt = zeros(1,nruns);
sopt = zeros(1,nruns);
cput = zeros(1,nruns);
if hist
	ohist = cell(nruns,1);
else
	ohist = [];
end

% DD optimisation (gradient descent)

if pp
	parfor k = 1:nruns
		fprintf('GD/ES optimisation parallel run %4d of %4d : ',k,nruns);
		tcpu = cputime;
		[dopt(k),Lopt(:,:,k),conv(k),sopt(k),iopt(k),ohist{k}] = opt_gd_dds(H,L0(:,:,k),niters,sig0,gdls,gdtol,hist);
		cput(k) = cputime-tcpu;
		fprintf('dd = %.4e : sig = %.4e : ',dopt(k),sopt(k));
		if conv(k) > 0, fprintf('converged(%d)',conv(k)); else, fprintf('unconverged '); end
		fprintf(' in %4d iterations : CPU secs = %6.2f\n',iopt(k),cput(k));
	end
else
	for k = 1:nruns
		fprintf('GD/ES optimisation serial run %4d of %4d : ',k,nruns);
		tcpu = cputime;
		[dopt(k),Lopt(:,:,k),conv(k),sopt(k),iopt(k),ohist{k}] = opt_gd_dds(H,L0(:,:,k),niters,sig0,gdls,gdtol,hist);
		cput(k) = cputime-tcpu;
		fprintf('dd = %.4e : sig = %.4e : ',dopt(k),sopt(k));
		if conv(k) > 0, fprintf('converged(%d)',conv(k)); else, fprintf('unconverged '); end
		fprintf(' in %4d iterations : CPU secs = %6.2f\n',iopt(k),cput(k));
	end
end

% Sort everything by dynamical dependence

[dopt,sidx] = sort(dopt);
Lopt = Lopt(:,:,sidx);
iopt = iopt(sidx);
conv = conv(sidx);
sopt = sopt(sidx);
cput = cput(sidx);
if hist
	ohist = ohist(sidx);
end
