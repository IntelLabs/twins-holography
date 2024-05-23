function [ hs ] = mk_hs( device_name, im_shift_pix )
    if(nargin < 1)
        device_name = 'himax';
    end
    hs.device_name = device_name; %'holoeye';%
    % target image shift in shortest wvl pixels
    if(nargin <2)
        hs.s0 = [0 0];
    else
        hs.s0 = im_shift_pix; 
    end
    hs.f = 1.30; % image distance in m

    hs.wvl = 638e-9;% red 520e-9; % green 450;% blue fisba doc

    hs.wvl_rgb = [638e-9 520e-9 450e-9];% [r g b] fisba doc

    %spectrometer
    hs.wvl_rgb = [637.2875e-9 521.76e-09 443.61e-09];
    %hs.wvl_rgb = [637.2875e-9 540e-09 460e-09];

    hs.wvl = hs.wvl_rgb(1);



    hs.slm_pix = [1080 1920];
    hs.lut_ti = [];
    hs.clut_ti = [];
    hs.Glmax = [256 255]; % [levels_count max_level] supported gray levels that phase maps to
    if( strcmp(hs.device_name, 'holoeye') )
        hs.slm_pp = 6.4e-6;%6.393e-6; %6.4e-6;%4e-6; % 4 micron % holoeye LETOII SLM
        hs.slm_ff = 0.93; % slm fill factor
        hs.Glmax = [256 255]; % supported gray levels that phase maps to
    elseif( strcmp(hs.device_name, 'compound') )
        hs.slm_pp = 3.015e-6; %3e-6; % compound photonics SLM spec
        hs.slm_ff = 0.935; % compound photonics slm fill factor
        hs.Glmax = [256 255]; % supported gray levels that phase maps to
    elseif( strcmp(hs.device_name, 'himax') )
        hs.slm_pp = 4.25e-6; %3e-6; % compound photonics SLM spec
        hs.slm_ff = 0.935; % compound photonics slm fill factor
        hs.Glmax = [128 254]; % supported gray levels that phase maps to
    elseif( strcmp(hs.device_name, 'ti30') )
        hs.slm_pp = 10.8e-6; % ti SLM spec
        hs.slm_ff = 0.935; % ti slm fill factor
        hs.slm_pix = [800 1280];
        hs.Glmax = [16 15]; % supported gray levels that phase maps to
    elseif( strcmp(hs.device_name, 'ti30d') )
        hs.slm_pp = 10.8e-6; % ti SLM spec
        hs.slm_ff = 0.935; % ti slm fill factor
        hs.slm_pix = [800 1280];
        hs.Glmax = [16 15]; % supported gray levels that phase maps to
        [hs.lut_ti, hs.clut_ti] = phase_lut_ti(0.05);
                hs.lut_ti = hs.lut_ti(1,:);
    elseif( strcmp(hs.device_name, 'ti60') )
        hs.slm_pp = 10.8e-6; % ti SLM spec
        hs.slm_ff = 0.935; % ti slm fill factor
        hs.slm_pix = [800 1358];
        hs.Glmax = [16 15]; % supported gray levels that phase maps to
        [hs.lut_ti, hs.clut_ti] = phase_lut_ti();
    elseif( strcmp(hs.device_name, 'ti30q') )
        hs.slm_pp = 10.8e-6; % ti SLM spec
        hs.slm_ff = 0.935; % ti slm fill factor
        hs.slm_pix = [800 1280];
        hs.Glmax = [16 15]; % supported gray levels that phase maps to
        hs.lut_ti = phase_lut_ti();
        hs.lut_ti = hs.lut_ti(1,:);
    else %basler
        hs.slm_pp = 2.2e-06;
        hs.slm_pix = [1944 2592];
        hs.Glmax = [256 255]; % supported gray levels that phase maps to        
    end
    
    hs.Asinc = ones(hs.slm_pix(1), hs.slm_pix(2), 3);
    hs.slm_pa = hs.slm_pp-0.2e-6;
    hs.r = 0.045;%0.05; % light source distance for double step fresnel
    hs.z = 1.30;%0.45; %1;% distance from slm plane to screen plane

    if nargin>=6
        hs.slm_ff = slm_ff;
    end

    if nargin>=5
        hs.f = f;
    end
    if nargin>=4
        hs.wvl = wvl;
    end
    if nargin>=3
        hs.slm_pp = slm_pp;
    end


    hs.slm_hdim = [(hs.slm_pix(1)-1) (hs.slm_pix(2)-1) 0];
    hs.slm_hdim(3) = sqrt(hs.slm_hdim(1)*hs.slm_hdim(1)+hs.slm_hdim(2)*hs.slm_hdim(2));
    hs.slm_hdim = hs.slm_hdim./2;

end

