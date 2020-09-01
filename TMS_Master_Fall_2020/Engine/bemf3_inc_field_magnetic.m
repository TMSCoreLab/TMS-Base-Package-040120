function Binc = bemf3_inc_field_magnetic(strcoil, Points, I0, mu0)  
%   Computes magnetic field from the coil via the FMM
%   in terms of the pseudo electric potential
%
%   Copyright SNM 2017-2020

    %  Compute dipole positions, directions, and moments
    M  = size(strcoil.Ewire, 1);
    N  = size(Points, 1);
    segvector   = (strcoil.Pwire(strcoil.Ewire(:, 2), :) - strcoil.Pwire(strcoil.Ewire(:, 1), :));
    moments     = I0*segvector.*repmat(strcoil.Swire, 1, 3);
    segpoints   = 0.5*(strcoil.Pwire(strcoil.Ewire(:, 1), :) + strcoil.Pwire(strcoil.Ewire(:, 2), :));   
   
     %  Compute and normalize pseudo normal vectors
    nx(:, 1) = +0*moments(:, 1);
    nx(:, 2) = -1*moments(:, 3);
    nx(:, 3) = +1*moments(:, 2);
    ny(:, 1) = +1*moments(:, 3);
    ny(:, 2) = +0*moments(:, 2);
    ny(:, 3) = -1*moments(:, 1);
    nz(:, 1) = -1*moments(:, 2);
    nz(:, 2) = +1*moments(:, 1);
    nz(:, 3) = +0*moments(:, 3);

    %   FMM 2019
    srcinfo.nd      = 3;                            %   three charge vectors    
    srcinfo.sources = segpoints';                   %   source points
    targ            = Points';                      %   target points
    prec            = 1e-1;                         %   precision    
    pg      = 0;                                    %   nothing is evaluated at sources
    pgt     = 1;                                    %   potential is evaluated at target points
    srcinfo.dipoles(1, :, :)     = nx.';                             
    srcinfo.dipoles(2, :, :)     = ny.';                             
    srcinfo.dipoles(3, :, :)     = nz.';                       
    U                        = lfmm3d(prec, srcinfo, pg, targ, pgt);
    Binc(:, 1)               = mu0*U.pottarg(1, :)/(4*pi);
    Binc(:, 2)               = mu0*U.pottarg(2, :)/(4*pi);
    Binc(:, 3)               = mu0*U.pottarg(3, :)/(4*pi); 
end

