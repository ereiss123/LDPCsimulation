function output = prob_combinations(p)

and_streams = [];
or_streams = [];
%xor_streams = [];

for l=1:numel(p)
    mask = de2bi(l,numel(p));
    idx = find(mask);
    and_streams = [and_streams p(l)^2 prod(p(idx))];
    or_streams = [or_streams orprob([p(l) p(l)]) orprob(p(idx))];
    %xor_streams = [xor_streams xorprob(p(idx))];
end

output = unique([p unique(and_streams) unique(or_streams) unique(1-and_streams) unique(1-or_streams)]'); %unique(xor_streams)]');

end

function output = orprob(p)
  if (numel(p) > 1)
      output = p(1) + (1-p(1))*orprob(p(2:numel(p)));
  else
      output = p;
  end
end

function output = xorprob(p)
  if (numel(p) > 1)
      output = p(1)*(1-xorprob(p(2:numel(p)))) + (1-p(1))*xorprob(p(2:numel(p)));
  else
      output = p;
  end
end