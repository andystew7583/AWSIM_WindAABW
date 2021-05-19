%%%
%%% findLastOutput.m
%%%
%%% Finds the last output file in an experiment directory.
%%%
function last_idx = findLastOutput (local_home_dir, run_name)

  last_idx = -1;
  listing = dir(fullfile(local_home_dir,run_name,'H0_n=*.dat'));
  for m = 1:length(listing)
    filename = listing(m).name;
    filenum = filename(6:end-4);
    idx = str2num(filenum);
    if (idx > last_idx) 
      last_idx = idx;
    end
  end

end

