% Copyright (c) 2022-2023, Intel Corporation.  All rights reserved.
% TWINS method for SLM phase response measurement
% Please refer to the paper below for details:
% "TWINS: improved spatial and angular phase calibration for holography,"
% Appl. Opt. 62, 7860-7867 (2023)
% Pre-print: https://doi.org/10.1364/opticaopen.24029844
close all;
ssp = 12; sw = 4; % slit spacing, slit width
slits_cfg = [3 20 ssp sw]; % half-period, full periods count, spacinf, width
slits_coord = {[931 531],[931-300 531-150]}; % measurement, reference
fnum = 1;
% this call is used inside to get SLM parameters
hs = mk_hs(); % change the function accordingly
dPall = cell(2,3,2); % scan full set (v,h;ph;l/u)
dts = convertStringsToChars(string(datetime('now','TimeZone','local','Format','yyyyMMdd_hhmmss')));
mkdir('mdata');
spath = ['mdata\ph_r', num2str(ssp), '_gr3_h_mfvib_' dts ];

% Measurement mode selection: 
% sim=0 - real setup with SLM; sim=1 - SLM and hologram simulation
sim=1;

if(sim==1) 
	%simulation
	pwnd_lu = {[1278 10 390 390+300-1 4]',[1416 10 390 390+300-1 4 ]'}; %{[upper measurement ROI], [upper reference ROI]}   
	pwnd_ll = {[642 10 390 390+300-1 4]',[458 10 390 390+300-1 4 ]'};   %{[lower measurement ROI], [lower reference ROI]}
    % use next 2 lines instead for visual ROI calibration
    % just add 6th element to the measurment ROI to enable ROI calibration
    % it is very useful for real camera setup as well
	%pwnd_lu = {[1278 10 390 390+300-1 4 8]',[1416 10 390 390+300-1 4 ]'}; %{[upper measurement ROI], [upper reference ROI]}   
	%pwnd_ll = {[642 10 390 390+300-1 4 8]',[458 10 390 390+300-1 4 ]'};   %{[lower measurement ROI], [lower reference ROI]}
    [dP, dPf, dPnf, dPn, Ic] = scan_slit_phase_full_avg(slits_cfg, slits_coord, {pwnd_lu,pwnd_ll}, [16666 18 60], fnum, [spath, '_smp.mat'],zeros(1080,1920));
else
	%%real setup with SLM 
	pwnd_lu = {[3370-3 10 977 977+1120-1 4 ]',[3567-3 10 977 977+1120-1 4 ]'};
	pwnd_ll = {[2502-18 10 977 977+1120-1 4 ]',[2263-21 10 977 977+1120-1 4 ]'};
	[dP, dPf, dPnf, dPn, Ic] = scan_slit_phase_full_avg(slits_cfg, slits_coord, {pwnd_lu,pwnd_ll}, [16666 18 60], fnum, [spath, '_smp.mat']);
end
figure; hold on;
Gg = 0:(hs.Glmax(2)/(hs.Glmax(1)-1)):hs.Glmax(2); % gray levels
plot(Gg ,dP{1}*2,"red"); plot(Gg,dP{2}*2,"blue");
hold off;
dPall{2,3,1} = dP{2}; dPall{2,3,2} = dP{1}; % dPu, dPl scan full set
save([spath, '.mat'], 'dP', 'dPf', 'dPn', 'dPnf','-v7.3');
save(['mdata\ph_r' num2str(ssp) '_mfvib_' dts '.mat'], 'dPall', '-v7.3');
