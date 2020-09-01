function Einc = bemf3_inc_field_electric(strcoil, Points, dIdt, mu0)    
%   Computes electric field from the coil via the FMM
%   in terms of the pseudo electric potential evaluated for segment centers
%
%   Copyright SNM 2017-2020

    %   Compute pseudo potentials       
    N  = size(Points, 1);
    segvector   =     (strcoil.Pwire(strcoil.Ewire(:, 2), :) - strcoil.Pwire(strcoil.Ewire(:, 1), :)).*repmat(strcoil.Swire, 1, 3);
    segpoints   = 0.5*(strcoil.Pwire(strcoil.Ewire(:, 1), :) + strcoil.Pwire(strcoil.Ewire(:, 2), :));
    PseudoQx    = segvector(:, 1);
    PseudoQy    = segvector(:, 2);
    PseudoQz    = segvector(:, 3);    
    Einc        = zeros(N, 3);
    const       = mu0*dIdt/(4*pi);  % normalization    
    
    %   This is -dAdt*j 
    %   FMM 2019
    srcinfo.nd      = 3;                            %   three charge vectors    
    srcinfo.sources = segpoints';                   %   source points
    targ            = Points';                      %   target points
    prec            = 1e-1;                         %   precision    
    pg      = 0;                                    %   nothing is evaluated at sources
    pgt     = 1;                                    %   potential is evaluated at target points
    srcinfo.charges(1, :)    = PseudoQx.';     %   charges
    srcinfo.charges(2, :)    = PseudoQy.';     %   charges
    srcinfo.charges(3, :)    = PseudoQz.';     %   charges
    U                        = lfmm3d(prec, srcinfo, pg, targ, pgt);
    Einc(:, 1)               = const*U.pottarg(1, :);
    Einc(:, 2)               = const*U.pottarg(2, :);
    Einc(:, 3)               = const*U.pottarg(3, :); 
end

