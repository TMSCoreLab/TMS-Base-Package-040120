function p = bemf5_volume_field_potential(Points, c, P, t, Center, Area, normals, Size, eps0, R, flag)
%   Computes electric potential for an array Points anywhere in space (line,
%   surface, volume). This field is due to surface charges at triangular
%   facets only. Includes accurate neighbor triangle integrals for
%   points located close to a charged surface.   
%   R is the dimensionles radius of the precise-integration sphere
%
%   Copyright SNM 2017-2020
%   R = is the local radius of precise integration in terms of average triangle size
    
    %   FMM 2019
    srcinfo.sources = Center';                      %   source points
    targ            = Points';                      %   target points
    prec            = 1e-2;                         %   precision    
    pg      = 0;                                    %   nothing is evaluated at sources
    pgt     = 1;                                    %   potential is evaluated at target points
    srcinfo.charges = c.'.*Area'/eps0;              %   charges
    U               = lfmm3d(prec, srcinfo, pg, targ, pgt);
    p               = +U.pottarg'/(4*pi);  
    
    if flag == 0    % Only the center-point approximation is used
        return;
    end
    %   Contribution of the charge of triangle m to the field at all points is sought                 
    %   Undo the effect of the m-th triangle charge on neighbor obs. points and
    %   add precise integration instead 
    M = size(Center, 1); 
    const = 4*pi*eps0;     
    ineighborlocal   = rangesearch(Points, Center, R*Size, 'NSMethod', 'kdtree'); % over triangles: M by X  
    for m =1:M
        index       = ineighborlocal{m};  % index into points that are close to triangle m   
        if ~isempty(index)
            temp        = repmat(Center(m, :), length(index), 1) - Points(index, :);   %   these are distances to the observation points
            DIST        = sqrt(dot(temp, temp, 2));                                    %   single column                
            I           = Area(m)./DIST;                                               %   center-point integral, standard format    
            p(index)    = p(index) - c(m)*I/const;         
            r1      = P(t(m, 1), :);    %   row
            r2      = P(t(m, 2), :);    %   row
            r3      = P(t(m, 3), :);    %   row           
            I       = potint(r1, r2, r3, normals(m, :), Points(index, :));     %   analytical precise integration MATLAB            
            p(index)= p(index) + c(m)*I/const;
        end    
    end                        
end

function p = PFieldFMMST(points, centers, areas, charges, eps)
%   Computes electric potential at target poiints resulting from charge densities
%   uniformy distributed over individual triangles using FMM (source to
%   source). Ignores self-contributions
%%   Inputs
%   points      - (M, 3) array of target points
%   centers     - (N, 3) array of triangle centers
%   areas       - (N, 1) array of triangle areas
%   charges     - (N, 1) charge densities
%   eps         - dielectric constant
%%  Outputs
%   Potential - (M, 1) array of electric potential values at target points, complex
%%  Operation
%   phi(m) = +1/(4*pi*eps)*...
%   sum(charges(n)*areas(n)*1/norm(centers(n, :)-points(m, :))

%   Copyright SNM 2017-2018
%   The Athinoula A. Martinos Center for Biomedical Imaging, Massachusetts General
%   Hospital & ECE Dept., Worcester Polytechnic Inst.

    %%  Use FMM by Zydrunas Gimbutas and Leslie Greengard
    %   Flags
    iprec = 1;                              %   this precision should be enough
    ifcharge = 1;
    ifdipole = 0;
    ifpot = 0;
    iffld = 0;
    ifpottarg = 1;
    iffldtarg = 0;
    nsource = size(centers, 1);             %   # of source points
    source  = centers';                     %   source points
    dipstr  = zeros(1, nsource);            %   dipole values (must keep even if not used)
    dipvec  = zeros(3, nsource);            %   dipole values (must keep even if not used)
    ntarget = size(points, 1);              %   # of target points
    target  = points';                      %   target points
    [U] = lfmm3dpart(iprec, nsource, source, ifcharge, charges'.*areas'/eps, ...
            ifdipole, dipstr, dipvec, ifpot, iffld, ntarget, target, ifpottarg, iffldtarg);
    p   = U.pottarg/(4*pi);
    p   = real(p.'); 
end


