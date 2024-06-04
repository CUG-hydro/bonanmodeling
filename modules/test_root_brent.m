% func(2)

% % root_brent (func, xa, xb, tol, varargin)
% [fx, x] = root_brent(@func, -10, 10, 1e-2)
% % [y, x] = func(1, 2)

% function [y, x] = func(x, varargin)
%     y = x.^2 - 4;
% end
% % func = @(x) x.^2 - 4;

[fluxvar, dx] = most(x, physcon, forcvar, surfvar, fluxvar)

open fluxvar
