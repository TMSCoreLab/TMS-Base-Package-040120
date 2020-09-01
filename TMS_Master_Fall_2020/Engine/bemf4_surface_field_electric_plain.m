function [P, E] = bemf4_surface_field_electric_plain(c, Center, Area)
%   Computes potential/continuous electric field on a surface facet due to
%   charges on ALL OTHER facets using plain FMM
%   Self-terms causing discontinuity may not be included
%   To obtain the true field, use E = E/eps0;
%
%   Copyright SNM 2017-2020    
   
    %  FMM 2019   
    %----------------------------------------------------------------------
    %   Fields plus potentials of surface charges (potential not used)
    %   Only FMM (without correction)
    prec             = 1e-1;
    pg              = 2;
    srcinfo.sources = Center';
    srcinfo.charges = (c.*Area)';
    U               = lfmm3d(prec, srcinfo, pg);
    P               = +U.pot'/(4*pi);
    E               = -U.grad'/(4*pi);            
end
