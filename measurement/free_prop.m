function [ Uout, Sout ] = free_prop( Uin, f, wvl, d1, dir )
% free space propagation via Fresnel integral approximation
% applies virtual lens with focal len f at distance d from the slm
%   Detailed explanation goes here

    pix_loc = 0;
    
    N = size(Uin, 1); % assume rect grid
    M = size(Uin, 2);
	k = 2.0*pi/wvl; % optical wavevector
	fY = (-N/2 : 1 : N/2 - 1) / (N * d1);
	fX = (-M/2 : 1 : M/2 - 1) / (M * d1);
	% observation plane coordinates
	[x2, y2] = meshgrid(wvl * f * fX, wvl * f * fY);
    [xp, yp] = meshgrid(((-M/2:1:M/2-1)+pix_loc)*d1, ((-N/2:1:N/2-1)+pix_loc)*d1);
	clear('fX');

    if (nargin < 5)||(dir>=0)
        % evaluate the Fresnel-Kirchhoff integral with
        % the quadratic phase factor inside 
        Uout = 1/(1i*wvl*f)*ft2(Uin.* exp(1i*k/(2*f)*(xp.^2 + yp.^2)), d1).*exp(1i*k*(f+(x2.^2+y2.^2)/(2*f)));
        Sout = abs(1/(1i*wvl*f))*d1*d1;
    else
        Uout = (1i*wvl*f)*ift2(Uin./exp(1i*k*(f+(x2.^2+y2.^2)/(2*f))), 1/(N * d1), 1/(M * d1))./exp(1i*k/(2*f)*(xp.^2 + yp.^2));
        Sout = abs(1i*wvl*f)*(1/d1)^2;
    end
end

