function Positions = initialization(SearchAgents_no, dim, ub, lb)

% convert bounds to row vectors
ub = ub(:)';  
lb = lb(:)';

% if scalar bounds, expand
if numel(ub) == 1
    ub = repmat(ub, 1, dim);
end
if numel(lb) == 1
    lb = repmat(lb, 1, dim);
end

% safety check
if numel(ub) ~= dim || numel(lb) ~= dim
    error('initialization: ub/lb must be scalar or length = dim');
end

Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;

end