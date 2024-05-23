function [P]=grating_slit2_bkg(P, Syx, O1yx, sp, Ph, G1, Pbg, Sbg, transp)
    % Sxy - size [ng w], O1yx - origin, sp-spacing, Ph - half-period
    % Sbg - size of background rectangle, Pbg - gray level for background
    if(nargin<=8)
        transp = false;
    end
    if(transp) 
        P = P';
    end
    if(nargin<=6)
        Pbg = 0;
    end
    if(size(Sbg,2)==1)
        Sbg = [Sbg(1) Sbg(1)];
    end
    ng = Syx(1); w = Syx(2);
    sz = [2*Ph*ng, 2*w+sp]; % size of both slits
    Sbg = max(Sbg, sz);
    Obg = O1yx-floor((Sbg-sz)/2); % slit is centered in bg rect
    
    P(Obg(1):(Obg(1)+Sbg(1)-1), Obg(2):(Obg(2)+Sbg(2)-1)) = Pbg; % prepare bg rectangle
    
    sz1 = [2*Ph*ng, w];
    G = bin_grating_y(Syx, Ph);

    P(O1yx(1):(O1yx(1)+sz1(1)-1),O1yx(2):(O1yx(2)+sz1(2)-1)) = G*G1(2)+(1-G)*G1(1); 
    
    
    O2yx = O1yx+[0 w+sp];
    P(O2yx(1):(O2yx(1)+sz1(1)-1),O2yx(2):(O2yx(2)+sz1(2)-1)) = G*G1(2)+(1-G)*G1(1);
    
    if(transp)
        P = P';
    end
end


function [G] = bin_grating_y(Syx, Ph)
    Gy = cat(1, zeros(Ph, 1), ones(Ph, 1));
    Gy = repmat(Gy, [Syx(1), 1]);
    G = repmat(Gy, [1, Syx(2)]);
end
