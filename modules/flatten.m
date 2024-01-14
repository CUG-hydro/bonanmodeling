function varargout = flatten(varargin)
%% flatten cells into variable
% 
%% Examples
% [x, y, z] = flatten(1, 2, 3);
% [x, y, z] = flatten({1, 2, 3});
% r = flatten(1);
if length(varargin) == 1 && iscell(varargin{1}); varargin = varargin{1}; end

varargout = varargin;
end
