function assert(out, mode)
%ASSERT signal the failure of the condition s
%
% s must be a valid boolean expression
% if mode = 'fail' the program aborts upon failure
% if mode = 'warn' a warnig is raised upon failure

% Author: A. Fusiello, 2005

if nargin == 1
    mode = 'fail';
  end
  
if ~out
  if strcmp(mode,'fail') 
    error('ASSERTION failed')
  elseif strcmp(mode,'warn')
    warning('ASSERTION failed')
  end
  
end
