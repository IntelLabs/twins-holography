% post-process scan results
scale = 1.14;
Gg = [0:2:254];
hs = mk_hs();
names = [
    "..\HoloHUD_data\ph_r12_gr1_v_mfvib_031623_smp.mat"
    "..\HoloHUD_data\ph_r12_gr1_v_mfvib_031723_smp.mat"
    "..\HoloHUD_data\ph_r12_gr1_v_mfvib_031823_smp.mat"
    "..\HoloHUD_data\ph_r12_gr2_v_mfvib_031623_smp.mat"
    "..\HoloHUD_data\ph_r12_gr2_v_mfvib_031723_smp.mat"
    "..\HoloHUD_data\ph_r12_gr2_v_mfvib_031823_smp.mat"
    "..\HoloHUD_data\ph_r12_gr3_v_mfvib_031623_smp.mat"
    "..\HoloHUD_data\ph_r12_gr3_v_mfvib_031723_smp.mat"
    "..\HoloHUD_data\ph_r12_gr3_v_mfvib_031823_smp.mat"];

names0 = [
    "..\HoloHUD_data\ph_r12_gr1_h_mfvib_031623_smp.mat"
    "..\HoloHUD_data\ph_r12_gr1_h_mfvib_031723_smp.mat"
    "..\HoloHUD_data\ph_r12_gr1_h_mfvib_031823_smp.mat"
    "..\HoloHUD_data\ph_r12_gr2_h_mfvib_031623_smp.mat"
    "..\HoloHUD_data\ph_r12_gr2_h_mfvib_031723_smp.mat"
    "..\HoloHUD_data\ph_r12_gr2_h_mfvib_031823_smp.mat"
    "..\HoloHUD_data\ph_r12_gr3_h_mfvib_031623_smp.mat"
    "..\HoloHUD_data\ph_r12_gr3_h_mfvib_031723_smp.mat"
    "..\HoloHUD_data\ph_r12_gr3_h_mfvib_031823_smp.mat"];

names0 = [
    "..\HoloHUD_data\ph_r12_gr1_v_mfvib_051123_smp.mat"
    "..\HoloHUD_data\ph_r12_gr2_v_mfvib_051123_smp.mat"
    "..\HoloHUD_data\ph_r12_gr3_v_mfvib_051123_smp.mat"];

names0 = [
    "..\HoloHUD_data\ph_r12_gr1_h_mfvib_051123_smp.mat"
    "..\HoloHUD_data\ph_r12_gr2_h_mfvib_051123_smp.mat"
    "..\HoloHUD_data\ph_r12_gr3_h_mfvib_051123_smp.mat"];

%     gr1 = load_measurements(names(2:3));
%     gr2 = load_measurements(names(5:6));
%     gr3 = load_measurements(names(8:9));
    gr1 = load_measurements(names(1:3));
    gr2 = load_measurements(names(4:6));
    gr3 = load_measurements(names(7:9));
%     gr1 = load_measurements(names(1:1));
%     gr2 = load_measurements(names(2:2));
%     gr3 = load_measurements(names(3:3));

    gr1.dPl = decode_groups(gr1.dPln, numel(Gg));
    gr1.dPl = unwrap(gr1.dPl)*2*scale;
    gr1.dPu = decode_groups(gr1.dPun, numel(Gg));
    gr1.dPu = unwrap(gr1.dPu)*2*scale;
    
    gr2.dPl = decode_groups(gr2.dPln, numel(Gg));
    gr2.dPl = unwrap(gr2.dPl)*2*scale;
    gr2.dPu = decode_groups(gr2.dPun, numel(Gg));
    gr2.dPu = unwrap(gr2.dPu)*2*scale;
    
    gr3.dPl = decode_groups(gr3.dPln, numel(Gg));
    gr3.dPl = unwrap(gr3.dPl)*2*scale;
    gr3.dPu = decode_groups(gr3.dPun, numel(Gg));
    gr3.dPu = unwrap(gr3.dPu)*2*scale;
    
    figure; hold on; 
    plot(Gg, gr1.dPl, 'r'); plot(Gg, gr1.dPu, 'r');
    plot(Gg, gr2.dPl, 'g'); plot(Gg, gr2.dPu, 'g');
    plot(Gg, gr3.dPl, 'b'); plot(Gg, gr3.dPu, 'b');
    %plot(Gg, (gr1.dPl+gr1.dPu+gr2.dPl+gr2.dPu+gr3.dPl+gr3.dPu)/6, 'm.');
    hold off;
    
    Py = interpolate_y({gr1.dPu gr1.dPl; gr2.dPu gr2.dPl; gr3.dPu gr3.dPl} , hs);
    Ply = lines(Py, Gg/Gg(end)*2*pi); % lines are normalized to take values from 0..2*pi
    
    figure; hold on;
    for y = 1:40:size(Py,1)
        plot(Gg/Gg(end)*2*pi, Py(y,:)); 
        plot(Gg/Gg(end)*2*pi, polyval(Ply(y,:), Gg/Gg(end)*2*pi));
    end
    hold off;
    %save(['..\HoloHUD_data\Ply_r_' num2str(scale*100) '.mat'], 'Ply');

function [] = reproduce_graph(gr, N, Navg, Gg, scale)
    figure; hold on;
    for i=1:N
        dPln = mean_dpn(gr.dPlnc, gr.lvalid, Navg);
        dPl = decode_groups(dPln, numel(Gg));
        dPl = unwrap(dPl)*2*scale;
        plot(Gg, dPl, 'b');
    end
    hold off;
    figure; hold on;
    for i =1:N
        dPun = mean_dpn(gr.dPunc, gr.uvalid, Navg);
        dPu = decode_groups(dPun, numel(Gg));
        dPu = unwrap(dPu)*2*scale;
        plot(Gg, dPu, 'b');
    end
    hold off;
end

function [Pl]= lines(Py, Gx)
    patch = true;
    Pl = zeros(size(Py,1), 2); % P = Ax+b;
    for y = 1:size(Py,1)
        l = polyfit(Gx, Py(y, :), 1);
        Pl(y,:) = l;
        if patch
            %v = polyval(l, [Gx(1) Gx(end)]);
            v = [Py(y,1) Py(y,end)];
            % Gx(1)->0 Gx(end)->l(Gx(end));
            A = (v(2)-v(1))/(Gx(end)-Gx(1)); 
            b = v(1)-A*Gx(1);
            Pl(y,:) = [A b];
        end
    end
end 
function [Py]= interpolate_y(Pgr, hs)
    % create per-y response curves by interpolating from 6 angles
    % Pgr{grperiod/2, (ord+3)/2}
    % negative angles are up, positive-down to make consistent with the
    % target image
    % known angles
    theta0d = [grating_angle(2, -1, 0), grating_angle(4, -1, 0), grating_angle(6, -1, 0), ...
    grating_angle(6, 1, 0), grating_angle(4, 1, 0), grating_angle(2, 1, 0)];
    L = numel(Pgr{1,1}); % phase levels
    % extract response based on angle
    Pin = zeros(6, L);
    Pin(1,:) = Pgr{1, 1}; Pin(2,:) = Pgr{2, 1}; Pin(3,:) = Pgr{3, 1}; % upper angles
    Pin(4,:) = Pgr{3, 2}; Pin(5,:) = Pgr{2, 2}; Pin(6,:) = Pgr{1, 2}; % lower angles
    H = hs.slm_pix(1); % height
    Py = zeros(H,L); % each y has phase response
    % required angles (y coord mapped to angle assuming normal incidence)
    ytan=(-H/2:1:H/2-1)/H*(hs.wvl_rgb(1)/hs.slm_pp); % (-H+1:2:H-1)/2
    thetay = atan(ytan); % asin
    for l=1:L
        Py(:,l) = interp1(theta0d, Pin(:,l)', thetay, 'linear');
    end
    figure; hold on;
    for i=1:20:H
        plot(Py(i,:));
    end
    hold off;
    figure; hold on;
    plot(Py(H/2,:)); plot(Py(H/2+1,:));
    hold off;
    figure; hold on;
    plot(thetay, Py(:, end)');
    plot(theta0d, Pin(:, end)');
    hold off;
end

function [gr] = load_measurements(names)
    sz = 0;
    for i=1:numel(names)
        load(names(i), 'dPlnc');
        sz = sz + size(dPlnc,1);
    end
    gr.dPlnc = cell(sz,1);
    gr.dPunc = cell(sz,1);
    gr.Fmaxlc = cell(sz,1);
    gr.Fmaxuc = cell(sz,1);
    p = 1;
    for i=1:numel(names)
        load(names(i), 'dPlnc', 'dPunc', 'Fmaxlc', 'Fmaxuc');
        sz0 = size(dPlnc,1);
        gr.dPlnc(p:p+sz0-1, 1 ) = dPlnc(:,1);
        gr.dPunc(p:p+sz0-1, 1 ) = dPunc(:,1);
        gr.Fmaxlc(p:p+sz0-1, 1 ) = Fmaxlc(:,1);
        gr.Fmaxuc(p:p+sz0-1, 1 ) = Fmaxuc(:,1);
        p = p+sz0;
    end
    lvalid = valid_mask(gr.Fmaxlc);
    uvalid = valid_mask(gr.Fmaxuc);
    [gr.dPlnc, luvalid] = unwrap_mask(gr.dPlnc);
    [gr.dPunc, uuvalid] = unwrap_mask(gr.dPunc);
    lvalid = lvalid & luvalid;
    uvalid = uvalid & uuvalid;
    gr.lvalid = lvalid;
    gr.uvalid = uvalid;
    gr.dPln = mean_dpn(gr.dPlnc, lvalid);
    gr.dPun = mean_dpn(gr.dPunc, uvalid);
end

function [mf] = meanf(Fmax)
    N = 0;
    mf = 0;
    for i = 1:size(Fmax,1)
        mf = mf+sum(Fmax{i,1}{1}(:,:,1), 'all');
        N = N+numel(Fmax{i,1}{1}(:,:,1));
    end
    mf = mf/N;
end

function [msk]=valid_mask(Fmax)
    fmax = meanf(Fmax);
    disp(['Max f = ' num2str(fmax) '; fi = ' num2str(round(fmax))]);
    fmax = round(fmax);
    
    sz = size(Fmax,1);
    msk = false(sz, 1);
    for i = 1:sz
        msk(i) = all(Fmax{i,1}{1}(:,:,1)==fmax, 'all');
        if(~msk(i))
            disp(['ATTN: pos=' num2str(i) ' peak frequency error!' ]);
        end
    end
end

function [dPout, msk]=unwrap_mask(dPin)   
    sz = size(dPin,1);
    msk = false(sz, 1);
    dPout = dPin;
    for i = 1:sz
        dPout{i,1}{1}=unwrap(dPin{i,1}{1}, [], 2);
        msk(i) = all(dPout{i,1}{1}==dPin{i,1}{1}, 'all'); % did unwrap change anything?
        if(~msk(i))
            disp(['ATTN: pos=' num2str(i) ' delta phase discontinuity error!' ]);
        end
    end
end

function [dPn]=mean_dpn(dPnc, msk, N)   
    sz = size(dPnc,1);
    i_all = 1:1:sz;
    i_msk = i_all(msk);
    if nargin>2
        N = min(N, numel(i_msk));
        i_N = randperm(numel(i_msk));
        i_N = i_N(1:N);
        i_msk = i_msk(i_N);
    end
    sz = numel(i_msk);
    dPn = zeros(size(dPnc{1,1}{1}));
    for ii = 1:sz
        i = i_msk(ii);
        dPn=dPn + dPnc{i,1}{1};
    end
    dPn = dPn/sz; % mean
end

function [] = plot_dpn(dPnc, i_msk)
    sz = numel(i_msk);
    figure; hold on;
    for ii = 1:sz
        i = i_msk(ii);
        plot(dPnc{i,1}{1}');
    end
    hold off;
end

function [P]=decode_groups(dPn, glevels, Pcorr)
    if nargin < 3
        Pcorr = zeros(size(dPn));
    end
    P = zeros(1, glevels);
    ngroups = size(dPn, 1);
    gn = size(dPn,2)-1;
    % remove tail of the last group by patching the last group size
    lastgs = gn-(ngroups*gn-glevels);
    gsizes = ones(1, ngroups)*gn; gsizes(ngroups) = lastgs;
    % make 0-based pieces of response curve
    for g=1:ngroups
        dPn(g, :)= dPn(g, :)-dPn(g, 1); % subtract base level
    end
    dPn = dPn-Pcorr; % should be safe to subtract here because Pcorr(:,1)=0;
    % first group starts the curve
    gs = gsizes(1);
    P(1:gs)=dPn(1, 2:(2+gs-1));
    % all other groups continue the curve
    for g=2:ngroups
        gs = gsizes(g);
        gbase = (g-1)*gn;
        P((gbase+1):(gbase+gs)) = dPn(g, 2:(2+gs-1))+P(gbase);
    end
end
