% Startup stuff

fprintf('[ssdi startup] Initialising ssdi\n');

% Add ssdi root dir + appropriate subdirs to path

global ssdi_root;
ssdi_root = fileparts(mfilename('fullpath')); % directory containing this file
addpath(ssdi_root);
fprintf('[ssdi startup] Added path %s\n',ssdi_root);

sdir = 'sims';
addpath(fullfile(ssdi_root,sdir));
fprintf('[ssdi startup] Added sub-path %s\n',sdir);

sdir = 'utils';
addpath(fullfile(ssdi_root,sdir));
fprintf('[ssdi startup] Added sub-path %s\n',sdir);

sdir = 'tests';
addpath(fullfile(ssdi_root,sdir));
fprintf('[ssdi startup] Added sub-path %s\n',sdir);

clear sdir

% Initialize mvgc library with gpmat and gvmat APIs

mvgc_path  = getenv('MVGC2_PATH');
gpmat_path = getenv('GPMAT_PATH');
gvmat_path = getenv('GVMAT_PATH');
assert(exist(mvgc_path,'dir') == 7,'bad MVGC path: ''%s'' does not exist or is not a directory',mvgc_path);
cd(mvgc_path);
startup;
cd(ssdi_root);
clear mvgc_path

fprintf('[ssdi startup] Initialised (you may re-run `startup'' at any time)\n');