function normalizedUnits(x,varargin)
% Parse the arguments
% - Create input parser
p = inputParser;
% - Parameter names must be specified in full
p.PartialMatching = false;
% - These parameters are optional, and can be in any order
addParameter(p,'normalize',false,@islogical);
addParameter(p,'denormalize',false,@islogical);
% - Apply input parser
parse(p,x,varargin{:});
% - Extract parameters
process = p.Results;

if (process.normalize)
    
end

if (process.denormalize)
    
end

end