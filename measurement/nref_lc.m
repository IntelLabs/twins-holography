function [nref] = nref_lc(theta0, alpha, ne, no)
    if (nargin <=2)
        ne = 1.9135; no = 1.535;
    end
    if(numel(theta0)>1 && numel(alpha)>1)
        [theta0, alpha] = meshgrid(theta0, alpha);
    end
    % ne0 = ne(theta0, alpha0, ne, no);
    nref = neff(theta0, alpha, ne, no);
end

function [neff] = neff(theta0, alpha, ne, no)
    neff = 1.0./sqrt((cos(theta0+alpha)./ne).^2 + (sin(theta0+alpha)./no).^2);
end