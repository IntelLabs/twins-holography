function [P]=grating_slit2(Syx, O1yx, O2yx, Ph, G1, G2, Pbkg, transp)
    % Sxy - size [ng w], O1yx, O2yx - origin, Ph - half-period
    hs = mk_hs();
    if(nargin<=6)
        Pbkg = 0;
    end
    
    if(nargin>=8 && transp>0)
        slm_dim = [hs.slm_pix(2) hs.slm_pix(1)];
    else
        slm_dim = hs.slm_pix;
    end
    
    P = zeros(slm_dim)+Pbkg;
    
    sz1 = sz_gr_y(Syx, Ph(1));
    G = bin_grating_y(Syx, Ph(1));

    P(O1yx(1):(O1yx(1)+sz1(1)-1),O1yx(2):(O1yx(2)+sz1(2)-1)) = G*G1(2)+(1-G)*G1(1); 
    
    
    sz2 = sz_gr_y(Syx,Ph(2));
    G = bin_grating_y(Syx, Ph(2));
    
    P(O2yx(1):(O2yx(1)+sz2(1)-1),O2yx(2):(O2yx(2)+sz2(2)-1)) = G*G2(2)+(1-G)*G2(1);
    
    if(transp)
        P = P';
    end
end

function [G] = bin_grating_y(Syx, Ph)
    % Ph <0 special case: no grating
    no_grating  = 0;
    if Ph < 0
        Ph = abs(Ph);
        no_grating = 1;
    end
    Gy = cat(1, zeros(Ph, 1)+no_grating, ones(Ph, 1));
    Gy = repmat(Gy, [Syx(1), 1]);
    G = repmat(Gy, [1, Syx(2)]);
end

function [sz] = sz_gr_y(Syx, Ph)
    Ph = abs(Ph);
    sz = [Syx(1)*2*Ph Syx(2)];
end