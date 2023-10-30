function [results] = get_ivim_v4(bvals,data,varargin)
%out = get_ivim(bvals,data,method,bsplit)
%provide this function with:
%bvals - b values vector
%data - 2D (x, b values) 3D (x,y,b values) or 4D (x,y,z,b values) with IVIM data
%method - "1step" Fits all parameters at once
%         "segmented" fit for D and S0 and then for Dstar and f
%         "bayes" 2 step grid search for D and S0 and Dstar

p = inputParser; %branch test

allowedmethods = {'1step','segmented','bayesian'};
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'bvals');
addRequired(p,'data');
addParameter(p,'method','1step',@(x) any(validatestring(x,allowedmethods)));
addParameter(p,'bsplit',250,validScalarPosNum);
addParameter(p,'mask',0);
parse(p,bvals,data,varargin{:});
