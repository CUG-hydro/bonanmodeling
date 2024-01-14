classdef module_Ipaper
  
  %% module_Radiation
  %
  %% Example:
  % import module_Ipaper.*
  
  properties (Constant)
  end
  
  properties
  end
  
  methods (Static)
    function s1 = merge_struct(s1, s2)
      %% merge two structures
      % if s1 has a field that s2 does not have, then add that field to s2
      names = fieldnames(s2);
      for i = 1:length(names)
        name = names{i};
        if ~isfield(s1, name)
          s1.(name) = s2.(name);
        end
      end
    end
    
    function s = var2struct(varargin)
      %% convert variables to a structure
      % s = var2struct(var1, var2, ...)
      s = struct;
      for var = 1:nargin
        s.(genvarname(inputname(var), fieldnames(s))) = varargin{var};
      end
    end

  end
end
