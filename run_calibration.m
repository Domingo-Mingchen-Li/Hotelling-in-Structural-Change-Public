function report = run_calibration(resume_file)
root=setup_solver();
addpath(fullfile(root,'calibration'));
c=calibration_config();
if nargin<1,resume_file='';end
if strcmp(resume_file,'resume')
 files=dir(fullfile(c.output_root,'*','progress.mat'));
 assert(~isempty(files),'cal:Resume','No progress.mat found under configured output_root');
 [~,k]=max([files.datenum]);
 resume_file=fullfile(files(k).folder,files(k).name);
end
report=calibrate_transition(c,resume_file);
end
