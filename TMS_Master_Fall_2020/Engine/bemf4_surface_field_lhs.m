function LHS = bemf4_surface_field_lhs(c, Center, Area, contrast, normals, weight, EC)   
%   Computes the left hand side of the charge equation for surface charges
%
%   Copyright SNM 2017-2020
    toc
    %   LHS is the user-defined function of c equal to c - Z_times_c which is
    %   exactly the left-hand side of the matrix equation Zc = b
    [P0, E0]      = bemf4_surface_field_electric_plain(c, Center, Area);    %   Plain FMM result    
    correction  = EC*c;                                                     %   Correction of plain FMM result
    LHS         = +c - 2*correction ...                                     %   This is the dominant (exact) matrix part and the "undo" terms for center-point FMM
                     - 2*(contrast.*sum(normals.*E0, 2)) ...                %   This is the full center-point FMM part
                     + weight*sum(c.*Area)/sum(Area);                       %   This is weight correction (optional)    
end
