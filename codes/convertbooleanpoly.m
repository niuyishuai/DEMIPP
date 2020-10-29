function [g,y] = convertbooleanpoly(f,x,mode)
    %stdbooleanpoly
    %% DESCRIPTION:
    %  Convert a Boolean polynomial function between {0,1} variables and {-1,1} variables
    %% SYNTAX:
    %   [g,y] = stdbooleanpoly(f,x,mode)
    %% INPUTS:
    %  f: Boolean polynomial with 0-1 variables
    %  x: polynomial variables (unused)
    %  mode: 0: from {-1,1} to {0,1}; 1: from {0,1} to {-1,1}. 
    %% OUTPUTS:
    %  g: result Boolean polynomial
    %% COPYRIGHT:
    %  Copyright since 2019, Yi-Shuai NIU. All Rights Reserved.
    %  2020/08/21    Initial Coding, support polylab
    
    if isa(f,'MPOLY')
        y=polylabvar(f.n);
        if mode==0 % from {-1,1} to {0,1}
            g=f.subs(2*y-1); % y is a {0,1} vector
        else % from {0,1} to {-1,1}
            g=f.subs((1+y)/2); % y is a {-1,1} vector
        end
        g=g.simplify;
    elseif isa(f,'sdpvar') || isa(f,'polynomial')
    elseif isa(f,'sym')
    elseif isa(f,'struct')
    end
end