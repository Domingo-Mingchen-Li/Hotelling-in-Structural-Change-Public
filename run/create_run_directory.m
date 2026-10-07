function [out,started] = create_run_directory(cfg)
if ~isfield(cfg,'run_kind'),cfg.run_kind='experiments';end
assert(any(strcmp(cfg.run_kind,{'experiments','validation','robustness'})), ...
    'v2:Output','run_kind must be experiments, validation or robustness');
if ~isfield(cfg,'output_tag'),cfg.output_tag='';end
assert(ischar(cfg.output_tag)&&(isrow(cfg.output_tag)||isempty(cfg.output_tag)), ...
    'v2:Output','output_tag must be a character row');
if ~isfolder(cfg.output_root),mkdir(cfg.output_root);end
clock_value=now;started=datestr(clock_value,'yyyy-mm-dd HH:MM:SS.FFF');
stamp=datestr(clock_value,'yyyymmdd_HHMMSS_FFF');
horizons=sprintf('%d_',cfg.horizons);horizons(end)=[];
name=sprintf('%s_%s_T%s_%s',cfg.run_kind,cfg.model,horizons,stamp);
tag=regexprep(cfg.output_tag,'[^A-Za-z0-9_-]','_');
if ~isempty(tag),name=[name '_' tag];end
out=fullfile(cfg.output_root,name);count=1;
while exist(out,'file')~=0
    count=count+1;out=fullfile(cfg.output_root,sprintf('%s_%03d',name,count));
end
[ok,message]=mkdir(out);assert(ok,'v2:Output','Cannot create output directory: %s',message);
end
