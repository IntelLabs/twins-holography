function G = ft2(g, delta1, delta2)
    % function G = ft(g, delta)
    if(nargin<2)
        delta1 = 1;
    end
    if(nargin<3)
        delta2 = delta1;
    end
    G = fftshift(fft2(fftshift(g))) * delta1*delta2;
    %dbg_print_min_max(fftshift(fft2(fftshift(g))));
end

function []=dbg_print_min_max(Z)
    disp(['Re =', num2str(min(min(real(Z)))),', ' num2str(max(max(real(Z)))),...
        'Im =', num2str(min(min(imag(Z)))),', ' num2str(max(max(imag(Z))))]);
end