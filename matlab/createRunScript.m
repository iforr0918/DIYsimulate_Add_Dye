%%%
%%% createRunScript.m
%%%
%%% Convenience function to create an execution script for model runs.
%%%
function createRunScript (  local_home_dir, ...
                            run_name, ...
                            model_code_dir, ...
                            exec_name, ...
                            use_intel, ...
                            use_pbs, ... 
                            use_cluster, ...
                            uname, ...
                            cluster_addr, ...
                            cluster_home_dir)

  %%% File/directory names    
  local_run_dir = fullfile(local_home_dir,run_name);    
  sfname = fullfile(local_run_dir,'Run.sh');
  cleanfname = fullfile(local_run_dir,'Clean.sh');
  upfname = fullfile(local_run_dir,'Upload.sh');
  downfname = fullfile(local_run_dir,'Download.sh');  
  
  %%% Create a script file to run the code
  sfid = fopen(sfname,'w');
  if (sfid == -1)
    error(['Could not open ',sfname]);
  end

  %%% Copy code files over
  if (~copyfile(['../',strtrim(model_code_dir),'/*'],local_run_dir))
    error(['Could not copy code files to ',local_run_dir]);
  end

  %%% To use intel compilers, need to add executables and libraries to system path
  if (use_intel)
    fprintf(sfid,'source /opt/intel/bin/compilervars.sh intel64\n');
  end

  %%% This line executes the C code
  run_str = ['./',exec_name,' ',run_name,'_in ','. \n'];    
  
  %%% PBS submission requires a PBS script to be prepared
  if (use_pbs)
    
    pbsfid = fopen(fullfile(local_run_dir,'run_pbs'),'w');
    
    pbsstr = ...
    ['#!/bin/sh \n' ...
    '# \n' ...
    '# Your job name \n' ...
    '#$ -N onelayer \n' ...
    '# \n' ...
    '# Use current working directory \n' ...
    '#$ -cwd \n' ...
    '# \n' ...
    '# Join stdout and stderr \n' ...
    '#$ -j y \n' ...
    '# \n' ...
    '# pe request for MPI. Set your number of processors here. \n' ...
    '# Make sure you use the "impi" parallel environment. \n' ...
    '#$ -pe impi 20 \n' ...
    '# \n' ...
    '# Run job through bash shell \n' ...
    '#$ -S /bin/bash \n' ...
    '# \n' ...
    '## Output file \n' ...
    '#$ -o ./output.txt \n' ...
    '# \n' ...
    '\n'];
  
    fprintf(pbsfid,'%s',pbsstr);
    fprintf(pbsfid,run_str);
    fclose(pbsfid);
    
    fprintf(sfid,'qsub run_pbs');   

  else
    
    fprintf(sfid,run_str);  
    
  end

  fclose(sfid);  
  
  %%% Create a file to clean the run folder of output files
  cleanstr = ['rm ./*n=*.dat'];
  cleanfid = fopen(cleanfname,'w');
  if (cleanfid == -1)
    error(['Could not open ',cleanfname]);
  end
  fprintf(cleanfid,cleanstr);
  fclose(cleanfid);
  
  if (use_cluster)
    
    %%% Create a file to upload the run to the cluster    
    upstr = ['rsync -av --update ',fullfile('..',run_name),' ',uname,'@',cluster_addr,':',cluster_home_dir];
    upfid = fopen(upfname,'w');
    if (upfid == -1)
      error(['Could not open ',upfname]);
    end
    fprintf(upfid,upstr);
    fclose(upfid);
    
    %%% Create a file to download the run from the cluster
    downstr = ['rsync -av --update ',uname,'@',cluster_addr,':',cluster_home_dir,'/',run_name,' ../'];
    downfid = fopen(downfname,'w');
    if (downfid == -1)
      error(['Could not open ',downfname]);
    end
    fprintf(downfid,downstr);
    fclose(downfid);
    
  end

end


