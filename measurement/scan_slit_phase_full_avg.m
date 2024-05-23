function [dP, dPf, dPnf, dPn, Ic] = scan_slit_phase_full_avg(gparam, pos0, row_nr_df0, exp_gain, fnum, outpath, Psim)
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
    % exp_gain - [exposure gain fps] - camera parameters
    % fnum - # of measurements to average
    % outpath - path to save detailed measurement information
    % Psim - simulation or image re-use mode. If Psim is an array of SLM 
    % dimensions, then simulation mode is on, and Psim is a static phase
    % curvature on simulated SLM. If Psim is cell array then image re-use
    % mode is on and Psim cell should contain (fnum, group_num, measurements_per_group)
    % pre-captured camera images.

    % Query SLM data
    hs = mk_hs();
    Gg = 0:(hs.Glmax(2)/(hs.Glmax(1)-1)):hs.Glmax(2);

    
    sim=0;
    Imem = [];
    if(nargin>6)
        if(iscell(Psim))
            Imem=Psim;
        else 
            sim = 1;
        end
    end

    if(nargin < 4 || isempty(exp_gain))
        exp_gain = [16666 18 60];
    end
    exp_gain_cam = exp_gain;
    if(sim==0 && isempty(Imem))
        exp_gain_cam = cell(1,3);
        [exp_gain_cam{2}, exp_gain_cam{3}] = cam_init(exp_gain);
        exp_gain_cam{1} = [];
    end
    dPlnc = cell(fnum,1);
    dPlnfc = cell(fnum,1);
    Fmaxlc = cell(fnum,1);
    Fmaxlfc = cell(fnum,1);
    dPunc = cell(fnum,1);
    dPunfc = cell(fnum,1);
    Fmaxuc = cell(fnum,1);
    Fmaxufc = cell(fnum,1);

    calib_roi = (numel(row_nr_df0{1}{1})>=6); % scan_slit_... returns early, Ic is empty
    for fr =1:fnum
        if(sim==0)
            if(isempty(Imem)) % live capture mode
                [dPlnc{fr}, Fmaxlc{fr}, dPlnfc{fr}, Fmaxlfc{fr}, Ic] = scan_slit_phase_full_one(gparam, pos0, row_nr_df0{1}, exp_gain_cam, fnum, []);
                if ~calib_roi
                    [dPunc{fr}, Fmaxuc{fr}, dPunfc{fr}, Fmaxufc{fr}, ~] = scan_slit_phase_full_one(gparam, pos0, row_nr_df0{2}, exp_gain_cam, fnum, [], Ic);
                else % Ic is empty when we use calibration mode
                    [dPunc{fr}, Fmaxuc{fr}, dPunfc{fr}, Fmaxufc{fr}, ~] = scan_slit_phase_full_one(gparam, pos0, row_nr_df0{2}, exp_gain_cam, fnum, []);
                end                    
            else % image re-use mode
                [dPlnc{fr}, Fmaxlc{fr}, dPlnfc{fr}, Fmaxlfc{fr}, ~] = scan_slit_phase_full_one(gparam, pos0, row_nr_df0{1}, exp_gain_cam, fnum, [], Imem(fr));
                [dPunc{fr}, Fmaxuc{fr}, dPunfc{fr}, Fmaxufc{fr}, ~] = scan_slit_phase_full_one(gparam, pos0, row_nr_df0{2}, exp_gain_cam, fnum, [], Imem(fr));
            end
        else % simulation mode
            [dPlnc{fr}, Fmaxlc{fr}, dPlnfc{fr}, Fmaxlfc{fr}, Ic] = scan_slit_phase_full_one(gparam, pos0, row_nr_df0{1}, exp_gain_cam, fnum, [], Psim);
            if ~calib_roi
                [dPunc{fr}, Fmaxuc{fr}, dPunfc{fr}, Fmaxufc{fr}, ~] = scan_slit_phase_full_one(gparam, pos0, row_nr_df0{2}, exp_gain_cam, fnum, [], Ic);
            else % Ic is empty when we use calibration mode
                [dPunc{fr}, Fmaxuc{fr}, dPunfc{fr}, Fmaxufc{fr}, ~] = scan_slit_phase_full_one(gparam, pos0, row_nr_df0{2}, exp_gain_cam, fnum, [], Psim);
            end

        end
    end
    if(~isempty(outpath))
        save(outpath, 'dPlnc','Fmaxlc','dPlnfc','Fmaxlfc', 'dPunc','Fmaxuc','dPunfc','Fmaxufc', '-v7.3');
    end
    nwins = size(dPlnc{1}, 2);
    % find carrier frequency
    % currently uses first frame
    % but can be easily modified to avg over i

    Fmaxl = Fmaxlc{1}; Fmaxu = Fmaxuc{1};
    Fmaxlf = Fmaxlfc{1}; Fmaxuf = Fmaxufc{1};
    for i=2:fnum
        for w=1:nwins
            Fmaxl{1,w} = Fmaxl{1,w} + Fmaxlc{1}{1,w};
            Fmaxu{1,w} = Fmaxu{1,w} + Fmaxuc{1}{1,w};
            Fmaxlf{1,w} = Fmaxlf{1,w} + Fmaxlfc{1}{1,w};
            Fmaxuf{1,w} = Fmaxuf{1,w} + Fmaxufc{1}{1,w};
        end
    end
    for w=1:nwins
        Fmaxl{1,w} = round(Fmaxl{1,w}/fnum);
        Fmaxu{1,w} = round(Fmaxu{1,w}/fnum);
        Fmaxlf{1,w} = round(Fmaxlf{1,w}/fnum);
        Fmaxuf{1,w} = round(Fmaxuf{1,w}/fnum);
    end
    dPln = avgg(dPlnc, Fmaxlc, Fmaxl);
    dPlnf = avgg(dPlnfc, Fmaxlfc, Fmaxlf);
    dPun = avgg(dPunc, Fmaxuc, Fmaxu);
    dPunf = avgg(dPunfc, Fmaxufc, Fmaxuf);

    [dPl, Pcorrl] = decode_corr(dPln, Fmaxl, numel(Gg));
    [dPlf, Pcorrlf] = decode_corr(dPlnf, Fmaxlf, numel(Gg));
    [dPu, Pcorru] = decode_corr(dPun, Fmaxu, numel(Gg));
    [dPuf, Pcorruf] = decode_corr(dPunf, Fmaxuf, numel(Gg));
    dP = {dPl, dPu};
    dPf = {dPlf, dPuf};
    dPn = {dPln, dPun};
    dPnf = {dPlnf, dPunf};

    if(iscell(exp_gain_cam))
        vid = exp_gain_cam{2};
        exp_gain_cam = {};
        delete(vid);
    end
end

function dPn = avgg(dPnc, Fmaxc, Fmax)
    fnum = size(dPnc,1);
    dPn = dPnc{1};
    nwins = size(dPn,2);
    for fr = 2:fnum
    	for w=1:nwins
            dPn{w} = dPn{w}+dPnc{fr,1}{w};
            if(~all(Fmaxc{fr,1}{w}==Fmax{w}))
                disp('Error: wrong peaks found');
            end
        end
    end
    for w=1:nwins
        dPn{w} = dPn{w}/fnum;
    end
end

function [dP, Pcorr] = decode_corr(dPn, Fmax, numGg)
    %decode dP
    Pcorr = zeros(size(dPn{1}));
    nwins = size(dPn,2);
    if nwins>1 % 1 or 2 windows for now
        Pcorr = reference_groups(dPn{nwins});
        Pcorr = phase_to_phase(Pcorr, Fmax{nwins}, Fmax{1});
        dPcorr = flatten_groups(Pcorr);
    end
        
    dP = decode_groups(dPn{1}, numGg, Pcorr);
    dP = unwrap(dP);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dPn, Fmax, dPnf, Fmaxf, Ic] = scan_slit_phase_full_one(gparam, pos0, row_nr_df0, exp_gain_cam, fnum, outpath, Psim)
    ng0=60;
    sim=0;
    Imem = [];
    if(nargin>6)
        if(iscell(Psim))
            Imem=Psim;
        else 
            sim = 1;
            % DEBUG: enable to check mean phase difference across slits
            % in the presence of static curvature.
            % Download 2D phase unwrapping code unwrap_phase.m from 
            % https://github.com/mfkasim1/unwrap_phase before enabling.
            % Pu = unwrap_phase(Psim);
        end
    end
    % gparam = [ph=grating 1/2-period, ng=# of periods, sp=slit spacing,
    % w=slit width]
    % pos - slit position
    % row_nr_df = [row = approximate row location, nr = #of rows to avg, x0, x1, df = min shift from center frequency]
    % exp_gain - [exposure gain]
    % fnum - # of frames to average
    % rot_angle - introduce later?
    % if it is a cell then we use reference slit
    if(iscell(row_nr_df0))
        row_nr_df = row_nr_df0; 
    else % legacy/no reference
        row_nr_df = {row_nr_df0};
    end
    if(iscell(pos0))
        pos = pos0;  % get position for measurement slit
    else % legacy/no refence
        pos = {pos0};
    end
    
    wnd = 20;
    bkg_on = 0;
    if(isempty(gparam))
        ph = 1; ng = ng0; sp = 12; w = 8; % sp+w=20 - optimal
    else
        ph = gparam(1); ng = gparam(2); sp = gparam(3); w = gparam(4);
        if numel(gparam)>4
            bkg_on = gparam(5);
        end
    end
    if(numel(row_nr_df{1})<2) % simulation mode only
        row_nr_df{1} = [row_nr_df(1) 1 0 0 7];
    end
    xpose = 0; % transpose? are we doing horizontal slits/vertical pattern
    if(size(row_nr_df{1},1)>1)
        % transposed slit
        %row_nr_df{1} = row_nr_df{1}';
        xpose = 1;
    end
    if(sim==0 && isempty(Imem) && numel(exp_gain_cam)~=3)
        disp(['Error: wierd camera parameters']);
        return;
    end
    
    gdim = [ph*2*ng sp+w];
    hs = mk_hs();
    slm_dim = [hs.slm_pix(1) hs.slm_pix(2)];
    if(xpose)
      slm_dim = [slm_dim(2) slm_dim(1)];
    end
    if(numel(pos{1})<2)
        pos{1} = slm_dim/2-[ph*2*ng sp+2*w]/2 + [0 750]*0; %ATTN
        disp(['Setting slit pos to center: y:', num2str(pos{1}(1)), ' x:', num2str(pos{1}(2))]);
    end
    Gg = 0:(hs.Glmax(2)/(hs.Glmax(1)-1)):hs.Glmax(2);
    Gstep = (hs.Glmax(2)/(hs.Glmax(1)-1));
    Gc = Gg(numel(Gg)/2+1); % contrast phase delta
    Gwrap = hs.Glmax(2)+Gstep; % one level beyond max level
    dP = zeros(1,numel(Gg));
    dPf = zeros(1,numel(Gg));
    Md = zeros(1,numel(Gg));
    PG = zeros(1,numel(Gg));
    gn = 8; %ATTN:!!!! 
    ngroups = ceil(numel(Gg))/gn;
    lastgs = gn-(ngroups*gn-numel(Gg)); % actual size of the last group
    Gsz=zeros(1,ngroups); Gsz(1:ngroups-1)=gn; Gsz(ngroups) = lastgs; % set group sizes
    Gn = zeros(ngroups, gn+1); % per-group grayscales with base as the first element
    % measurements
    nwins = numel(row_nr_df);
    dPn = cell(1, nwins); dPn(1,:) = {zeros(ngroups, gn+1)}; % per-group measurements with base - first element
    dPnf = cell(1, nwins); dPnf(1,:) = {zeros(ngroups, gn+1)}; % per-group measurements with base - first element
    Fmax = cell(1, nwins); Fmax(1,:) = {ones(ngroups, gn+1, 2)}; % max frequency index, width
    Fmaxf = cell(1, nwins); Fmaxf(1,:) = {ones(ngroups, gn+1, 2)}; % max frequency index, width
    Ic=cell(ngroups, gn+1);     % images
    Nc=cell(ngroups, gn+1);     % image file names
    Pcontr = zeros(ngroups, gn+1);
    % initialize per-group gray levels
    Gn(1,1) = Gg(1); Gn(2:ngroups,1)=Gg(gn:gn:gn*(ngroups-1)); % base levels
    for g = 1:ngroups %per-group levels
        gs = Gsz(g);
        Gn(g,2:(2+gs-1))=Gg(((g-1)*gn+1):((g-1)*gn+gs));
    end
    for g = 1:ngroups
        gend = Gsz(g)+1;
        for i=1:gend
            if(numel(row_nr_df{1})>=6) % calibration mode
                Gp=row_nr_df{1}(6); % get specific step
                i = mod(Gp/2, gn)+1+1;
                g = floor(Gp/2/gn)+1;
            end
            Gp = Gn(g, i); % current gray level
            Gbase = Gn(g, 1); % take base (last element from previous group)
            % standard gbase
%             G0 = [Gbase Gp+128]; G0 = mod(G0, 256);
%             G1 = [Gp Gp+128]; G1 = mod(G1, 256);
            % gbase contrast - gbc
            G0 = [Gbase+Gc Gbase ]; G0 = mod(G0, Gwrap);
            G1 = [Gbase+Gc Gp ]; G1 = mod(G1, Gwrap);
            % ATTN: test G0 = [128 0]; G1 = [128 0];
            P = grating_slit2([ng w], pos{1}, pos{1}+[0 w+sp], [ph ph], G0, G1, mod(Gbase+Gc, Gwrap)*bkg_on, xpose);%ATTN! %figure; imshow(P);
            if(numel(pos)>1) % use reference slits?
                ph1 = 1*(ph==2)+2*(ph~=2); % ATTN: move reference pattern next to ph-pattern
                %ph1 = ph; % ATTN !!!!
                % apply extra slits (up to 2)
                for s2=2:numel(pos) % TODO: fix this part to support different set of slits vs windows
                    % reference and measurement slits height must be the same
                    % to have the same brightness
                    ng1 = ceil(ng*ph/ph1);
                    xpose1 = (size(row_nr_df{s2},1)>1);
                    P = grating_slit2_bkg(P, [ng1 w], pos{s2}, sp, ph1, [0 Gc], 0, 120, xpose1);
                end
            end
            %P = mod(P, 256);

            % one frame start
            if(sim == 0)
                if(isempty(Imem))
                    Pim = uint8(P);
                    Id = get_frame(Pim, 1, exp_gain_cam);
                else
                    Id = double(Imem{g,i})/4095;
                end
            else
                % various phase mapping functions for simulation
                PGP = @(x) 2*pi*x/Gg(end);
                %PGP = @(x) (tan(-pi/4+pi/2*x/Gg(end))+1)*pi;
                %PGP = @(x) ((-1+2*x/Gg(end)).^3+1)*pi;
                PG = PGP(Gg);
                Id = far_field(PGP(P)+Psim, 1, 0.7)*1000; 
                %ctr = floor(size(Id)/2)+1;
                %Id(ctr(1), ctr(2)) = 0;
                % vibration simulation
                Svib = 0;
                Id = imtranslate(Id,Svib*[normrnd(0,1), normrnd(0,1)],'FillValues',0); 
                % DEBUG: enable to check mean phase difference across slits
                % in the presence of static curvature
                %Md = mean_slit_x(Pu, [ng w], pos{1}, pos{1}+[0 w+sp], [ph ph]);
                Pcontr(g, i) = PGP(mod(Gp+Gc,Gwrap))-PGP(Gp);
            end
            % ATTN
            %Id = imrotate(Id, atan2(1087-1117,1843-152)/pi*180*0.5, 'bicubic', 'crop');
            % check calibration
            if(numel(row_nr_df{1})>=6) % calibration mode
                %show_calib1(Id, gr_y, avg_y, [gr_x0 gr_x1], wnd);
                show_calib(Id*2, row_nr_df, wnd);
                if(iscell(exp_gain_cam))
                    vid = exp_gain_cam{2};
                    exp_gain_cam = {};
                    delete(vid);
                end
                return;
            end
            
            for s2 = 1:numel(row_nr_df)
                % process window
                gr_y=row_nr_df{s2}(1); avg_y = row_nr_df{s2}(2);
                gr_x0 = max(row_nr_df{s2}(3), 1); % set to 1 if 0
                gr_x1 = row_nr_df{s2}(4); % 0 - don't know the resolution, correct later
                df = row_nr_df{s2}(5);
                % fix end x
                if(gr_x1 <=0)
                    gr_x1 = size(Id,2);
                end

                % transpose
                if(xpose)
                    Id = Id';
                end

                %figure;imshow(Id((gr_y-5):(gr_y+avg_y-1+5), gr_x0:gr_x1)*1);
                %Id = medfilt2(Id, [11 11]);
                Yg = mean(Id(gr_y:(gr_y+avg_y-1), gr_x0:gr_x1), 1); 
                [dPn{s2}(g,i), mf]=get_phase(Yg, df);
                [dPnf{s2}(g,i), mff] =get_phase(filt_data(Yg, wnd), df);
                Fmax{s2}(g,i,:) = [mf numel(Yg)]; % store max frequency and width in pixels for correction
                Fmaxf{s2}(g,i,:) = [mff numel(Yg)];
%             diff_ptrn =@(X, S, c, x0, dp) S*2*(cos(2*pi*c*(X-x0)*(sp+w)/2+dp)+1).*sinc(c*(X-x0)*w/2).^2;
%             SSE = @(B) sum(Yg - diff_ptrn([-numel(Yg)/2:(numel(Yg)/2-1)], B(1), B(2), B(3), B(4))).^2;
%              B00 = [1, numel(Yg)/(numel(Yg)/2+1-mf)/2/pi, 0, dPn(g,i)];
%             B0 = fminunc(SSE, B00);
%             dPnf(g, i) = B0(4);
%            disp(['S: ', num2str(B0(1)), ' c0: ',num2str(B00(2)) , ' c: ', num2str(B0(2)), ' x0:', num2str(B0(3)), ' dP:', num2str(B0(4))]);
            
            
                disp(['G1: ', num2str(G0(1)),', ',num2str(G0(2)),...
                    'G2: ', num2str(G1(1)), ', ', num2str(G1(2)),...
                    ', peak=',num2str(Fmax{s2}(g,i,1)), 'Dp=', num2str(dPn{s2}(g,i))]);
            
                if(xpose)
                    Id = Id';
                end
            end % for over windows
            
            if(nargout >= 5) % save data
                Ic{g,i} = uint16(Id*4095);
            end
            % frame end

        end % for over i - measurements in the gorup
    end % for over g - groups   
    % save data mode
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function [Pref] = reference_groups(Pn)
    Pref = zeros(size(Pn));
    for g = 1:size(Pn,1)
        Pref(g,:) = Pn(g,:)-Pn(g,1);
    end
end

function [P] = flatten_groups(Pg)
    P = Pg(:,2:end)';
    P = P(:)';
end

function [Spix] = phase_to_shift(P, Fmax)
    % convert phase in P to signed pixel shift using
    % Fmax - a table with [freq_idx wnd_size_pix]
    % Fmax was initialized with 1s so we shouldn't
    % have freq_idx = 0
    Size = Fmax(:,:,2);
    F = Fmax(:,:,1); % frequencies
    if any(F(:)==0)
        disp('Zero frequency(s) detected: possible error. Patching with 1.');
        F(F(:,:)==0) = 1;
    end
    Spix = P./(2*pi).*Size./F;
end

function [P] = shift_to_phase(Spix, Fmax)
    % convert signed pixel shift Spix to phase P using
    % Fmax - a table with [freq_idx wnd_size_pix]
    % Fmax was initialized with 1s so we shouldn't
    % have wnd_size_pix = 0
    Size = Fmax(:,:,2);
    F = Fmax(:,:,1); % frequencies
    if any(Size(:)==0)
        disp('Zero window size detected: possible error. Patching with 1.');
        Size(Size(:,:)==0) = 1;
    end
    P = Spix.*F./Size.*(2*pi);
end

function [P1] = phase_to_phase(P0, Fmax0, Fmax1)
    % convert shift corresponding to phase P0 defined by
    % Fmax0 - a table with [freq_idx wnd_size_pix]
    % to phase P1 corresponding to the same pixel shift
    % of the signal defined by Fmax1
    if(all(Fmax0==Fmax1))
        P1 = P0;
        return;
    end
    Spix = phase_to_shift(P0, Fmax0);
    P1 = shift_to_phase(Spix, Fmax1);
end

function [Peq] = equalize_phase(P, Fmax, glevels)
    % check if frequencies are different in different measurements
    % and convert everything to the first frequency
    % Fmax - a table with [freq_idx wnd_size_pix]
    % P a raw phase measured
    Peq = P;
    % patch the tail of the last group with Fmax(1,1,:)
    ngroups = size(P,1);
    gn = size(P,2)-1; % group size
    tailst = gn+1 - (ngroups*gn-glevels) + 1; % start of undefined tail
    for i = tailst:gn
        Fmax(ngroups, i, :) = Fmax(1,1,:);
    end
    for g = 1:ngroups
        for i = 1:size(P,2)
            if any(Fmax(g,i,:)~=Fmax(1,1,:))
                disp(['Patching (g,i):(', num2str(g), ', ', num2str(i), ')']);
                Spix = phase_to_shift(P(g,i), Fmax(g,i,:));
                Peq(g,i) = shift_to_phase(Spix, Fmax(1,1,:));
            end
        end
    end
end


function [P, mf]=get_phase(Y, df)
    Ygf = fftshift(fft(fftshift(Y-mean(Y))));
    [~, mf]=max(abs(Ygf(1:(size(Ygf,2)/2+1-df))));
    P = angle(Ygf(mf));
    mf = numel(Ygf)/2+1-mf; % convert to frequency #
%     wf = 3;
%      P0 = angle(Ygf([mf-wf:1:mf+wf]));
%      W0 = abs(Ygf([mf-wf:1:mf+wf]));
%      W0 = W0.^2;
%      P = sum(P0.*W0)/sum(W0);
%      mf = sum([mf-wf:1:mf+wf].*W0)/sum(W0);
end

function [Xf] = filt_data(X, wnd)
    Xf = filter((1/wnd)*ones(1,wnd), 1, X);
end

function []=show_calib(Id, row_nr_df, wnd)
    % show image with rectangles for all windows
    clut = ['r', 'g', 'b', 'm'];
    figure; imshow(Id);
    hold on;
    for s2 = 1:numel(row_nr_df)
        gr_y = row_nr_df{s2}(1); avg_y=row_nr_df{s2}(2);
        xrng = [row_nr_df{s2}(3) row_nr_df{s2}(4)];
        if(size(row_nr_df{s2},1)==1) % no transpose
            rect = [xrng(1), gr_y, xrng(2)-xrng(1)+1, avg_y];
        else % transpose
            rect = [gr_y, xrng(1), avg_y, xrng(2)-xrng(1)+1];
        end
        rectangle('Position', rect, 'LineWidth', 0.5, 'EdgeColor', clut(s2));
    end
    hold off;
    
    for s2 = 1:numel(row_nr_df)
        gr_y = row_nr_df{s2}(1); avg_y=row_nr_df{s2}(2);
        xrng = [row_nr_df{s2}(3) row_nr_df{s2}(4)];
        df = row_nr_df{s2}(5);
        xpose = (size(row_nr_df{s2},1)>1);
        if xpose
            Id = Id';
        end
        Yg = mean(Id(gr_y:(gr_y+avg_y-1), xrng(1):xrng(2)), 1); 
        figure; plot(Yg);
        Ygf = fftshift(fft(fftshift(Yg-mean(Yg))));   
        figure; plot(abs(Ygf));
        Yf = filt_data(Yg, wnd);
        figure; plot(Yf);
        Yff = fftshift(fft(fftshift(Yf-mean(Yf))));   
        figure; plot(abs(Yff));
        [~, mf] = get_phase(Yg, df);
        [~, mff] = get_phase(Yf, df);
        disp(['Freq.max = ', num2str(mf), ', npix = ', num2str(numel(Yg)/mf), 'FILT Freq.max = ', num2str(mff), ', npix = ', num2str(numel(Yf)/mff)]);
        if xpose
            Id = Id';
        end
    end
end

function []=show_calib1(Id, gr_y, avg_y, xrng, wnd)
    figure; imshow(Id);
    hold on;
    rectangle('Position', [xrng(1), gr_y, xrng(2)-xrng(1)+1, avg_y], 'LineWidth', 1, 'EdgeColor', 'r');
    hold off;
    Yg = mean(Id(gr_y:(gr_y+avg_y-1), xrng(1):xrng(2)), 1); 
    figure; plot(Yg);
    Ygf = fftshift(fft(fftshift(Yg-mean(Yg))));   
    figure; plot(abs(Ygf));
    Yf = filt_data(Yg, wnd);
    figure; plot(Yf);
    Yff = fftshift(fft(fftshift(Yf-mean(Yf))));   
    figure; plot(abs(Yff));
    [~, mf] = get_phase(Yg, 5);
    [~, mff] = get_phase(Yf, 5);
    mf = numel(Yg)/2+1-mf;
    mff = numel(Yf)/2+1-mff;
    disp(['Freq.max = ', num2str(mf), ', npix = ', num2str(numel(Yg)/mf), 'FILT Freq.max = ', num2str(mff), ', npix = ', num2str(numel(Yf)/mff)]);
end

function [Md]= mean_slit_x(P, dim, p0, p1, ph)
    % [ng w], [10 50], [10 50 + w+sp], [ph ph]
    dim0 = [dim(1)*abs(ph(1))*2 dim(2)];
    dim1 = [dim(1)*abs(ph(2))*2 dim(2)];
    roi0 = P(p0(1):(p0(1)+dim0(1)-1), p0(2):(p0(2)+dim0(2)-1));
    roi1 = P(p1(1):(p1(1)+dim1(1)-1), p1(2):(p1(2)+dim1(2)-1));
    M0 = mean(roi0(:));
    M1 = mean(roi1(:));
    roid = roi1-roi0;
    Md = mean(roid(:));
end

function [I] = far_field(P, c, f, lpos)
    if(nargin<2)
        Ig = fftshift(fft2(fftshift(exp(1i*P))));
    else
        hs = mk_hs();
        if(nargin>3)
            P = P+phase_light(c, lpos);
        end
        [ Ig, Sout ] = free_prop( exp(1i*P), f, hs.wvl_rgb(c), hs.slm_pp, +1 );
        Ig = Ig/Sout;
    end
    I = abs(Ig);
    I = I/max(I(:));
    I = I.^2;
end

function [Id] = get_frame(Pim, fnum, exp_gain)
    Pim = repmat(Pim,1,1,3);
    FScreen(Pim, 1);
    pause(0.07);
    % time average
    I = getimage(exp_gain);
    Id = double(I)/4095;
    for f = 1:1:(fnum-1)
        I = getimage(exp_gain);
        Id = Id+double(I)/4095;
    end
    Id = Id/fnum;
end