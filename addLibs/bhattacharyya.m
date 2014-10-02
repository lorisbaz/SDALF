function d = bhattacharyya(k,q)
% d = bhattacharyya(k,q)
% Computer the bhattacharyya-coefficient based distance.
%
% Input:
%   k,q : 2 color (or ldg) histogram (same dim.)
%
% Output:
%   d : bhattacharyya metric value

k = normalize(k);
q = normalize(q);

d = real(sqrt(1-sum(sqrt(k.*q)))); % bhattacharyya metric
% d = -log(sum(sqrt(k.*q)));  % bhattacharyya distance
% d = sum(sqrt(k.*q));  % bhattacharyya coefficient