Description

This code is used for measuring phase modulation using programmable double slit method called TWINS. Please refer to the paper for the method details.

1. Main measurement function:
scan_slit_phase_full_avg(gparam, pos0, row_nr_df0, exp_gain, fnum, outpath, Psim):     
    % gparam = [ph=grating 1/2-period, ng=# of periods, sp=slit spacing, w=slit width]
    % pos0 - slit position, if cell array then first element is measurement
    % slit position, and the second is reference slit position. If
    % row_nf_df elements are transposed, then pos0 are coordinates
    % in transposed SLM.
    % row_nr_df = (1,2) cell array of parameter arrays 
    % [row = approximate row location, nr = #of rows to avg,
    % x0, x1, df = min shift from center frequency]. 
    % Or (1,2) cell array of 2 such cell arrays.
    % If parameter arrays are transposed (size=[5,1]) then slits
    % are horizontal and all coordinates are in trasposed camera 
    % image.
    % See readme.md for more details.
    % exp_gain - [exposure gain]
    % fnum - # of measurements to average
    % outpath - path wherevto save results
    % Psim - simulation or pre-captured images (see readme.md)

gparam - double-slit definition, all values are in SLM pixels
pos0 - double-slit position ([y x]), top left corner of the first slit in SLM pixels. If pos0 is cell array then the first element is position of measurement slit and the 2nd element is position of reference slit. If row_nr_df0 elements are transposed, then coordinates in transposed SLM ([x y]). 
row_nr_df0 - cell array specifying camera image ROI and peak/carrier frequency search ROI. (1,2) cell array containing parameters of top and bottom (-1 and +1 order) fringe patterns on the camera image: {uparam, lparam}. uparam and lparam can be (1,2) cell arrays or arrays of parameters.
	If uparam and lparam are cell arrays then element {1,1} is ROI for measurement slit fringe pattern and element {1,2} is ROI for pattern of measurement slit. That way we can use -1, +1 orders of reference slit.
	For measurements each parameter array has structure [row, nr, x0, x1, df]. Camera image ROI defines rectangle in the camera image where fringe pattern is located (defined by row,nr,x0,x1).
	row - top row of the ROI (img pixels)
	nr - height of the ROI (img pixels), the signal will be averaged over this height, resulting in 1 row
	x0 - left column of the ROI, x1 - right column of the ROI (img pixels)
	df - peak frequency search ROI start (fft frequency bins). The frequncies between central and df (between 0 and df-1) will be excluded from the search for carrier peak. Typically required for real experiments to filter out non-uniform background noise.
	Calibration mode: if first array of uparam and lparam contains 6 elements [row, nr, x0, x1, df, gray_level], then the function will display captured image with fringe pattern for measurement (and reference) slit(s) and rectangle(s) corresponding to the specified ROI(s). gray_level will be used for 'on' pixels of measurement slit.
	Horizontal slits: if individual arrays [row, nr, x0, x1, df (, gray_level)] are transposed (size=(5,1) or (6,1)) then slits will be horizontal and coordinates are in transposed camera image. Calibraton mode will display tranposed camera images for easy ROI adjustment.
exp_gain - imaging camera paramters, array of up to 4 values: [exposure gain framerate hw trigger delay]. Please refer to cam_init description.
outpath - path to store intermediate measurement results for each frame
Psim - simulation mode with initial SLM phase profile (use zeroes(SLM dimensions) to simulate ideally flat SLM). If (fnum, group_num, measurement_per_group) cell array with pre-captured camera images, then they will be used insted of live capture.

Attention: to use horizontal slits instead of vertical pls use transposed row_nr_df0 array, the code will internally transpose SLM image and camera image. The pos0 and row_nr_df0 need to be specified in transposed coordinates (correspondent to transposed SLM and camera image).

Attention: make sure FScreen(<display index>, ...) call in the code refers to the SLM display.   

The function mk_hs() is used to define parameters for the system. Some important parameters are: 
slm_pix - defines SLM image resolution
wvl     - defines wavelength of the light source
f       - defines propagation distance

2.scan_all.m is an example script using the measurement function with some default values

Requirements

1. Compile FScreen.cpp to make a mex file using the following command in MATLAB (refer to MATLAB instructions on builifding mex files if there are any propblems):
mex FScreen.cpp -R2018a
This file uses Windows APIs to display an input in a full-screen window on one of the displays.
We provide only Windows version of the mex function.

2. The code uses Image Acquisition Toolbox to get images from the camera. Only limited set of input parameter is supported, the rest are hardcoded (e.g. we always request first 'gentl' camera supporting 12-bit monochrome mode).
All parameters are located in one file cam_init.m. Please refer to the source if modifications are required.
cam_init function takes an array of up to 4 parameters: exposure, gain, framerate and hw trigger delay. If hw trigger delay is specified the camera will be instructed to use hw trigger on line 2.