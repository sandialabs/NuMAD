
% params.Class = 5;
% 
% if ~isequal(params.Class,1) && ~isequal(params.Class,2) && ~isequal(params.Class,3)
%         error('`params.Class` must be equal to 1, 2, or 3');
% end


params.fastsim = 'kelley';

if ~strcmp(params.fastsim,'fast') && ~strcmp(params.fastsim,'fast simulink') && ~strcmp(params.fastsim,'adams')
        error('`params.fastsim` must be equal to "fast", "fast simulink", "adams"');
end   