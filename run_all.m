%% RUN_ALL  Reproduce the results and figures of the malaria-hookworm fuzzy manuscript.
%   Press Run. The analysis prints the key numbers in the Command Window and the seven
%   figures of the manuscript then open on screen (they are also saved as PDF in ./figures).
clear; close all; clc;

RUN_ANALYSIS = true;    % true: run the full analysis (about 30-60 minutes)
                        % false: only draw the figures from the saved results (seconds)
WRITE_CSV    = false;   % true: also save CSV tables of all results in ./results
RUN_CHECKS   = false;   % true: also run the robustness checks (several hours)

here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here); addpath(fullfile(here, 'src'));
setappdata(0, 'fmh_write_csv', WRITE_CSV);
setappdata(0, 'fmh_run_checks', RUN_CHECKS); setappdata(0, 'fmh_root', here);

if RUN_ANALYSIS
    run(fullfile(here, 'main_fuzzy_malaria_hookworm.m'));   % note: this script clears the workspace
end

here = getappdata(0, 'fmh_root');
make_figures();                                             % reads results/figure_data.mat

if getappdata(0, 'fmh_run_checks')
    run(fullfile(here, 'checks', 'conv_static.m'));
    for N = [200 500 1000]
        alv = [0 0.25 0.5 0.75];
        run(fullfile(here, 'checks', 'conv_scen.m'));
    end
    run(fullfile(here, 'checks', 'alt_shapes.m'));
    run(fullfile(here, 'checks', 'alt_shapes_fuzzy.m'));
end
