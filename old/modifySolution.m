% Discretize solution: find rounded solution with minimal squared error with respect to measurements
function [discreteOutput]=modifySolution(fractionalOutput,M,y)

% sort the fractionalOutput from highest to lowest
[junk,index] = sort(fractionalOutput,'descend');
index(find(junk<0.2)) = [];


currXround = zeros(size(fractionalOutput));
errorRound = zeros(1,length(index));

val = ones(1,length(index));
% round the top 1 to i entries in x and for each of them calculate the squared error
for i=1:length(index)
  
  if abs(fractionalOutput(index(i))-1)<abs(fractionalOutput(index(i))-2)
    currXround(index(i)) = 1;
  else
    currXround(index(i)) = 2;
    val(i) = 2;
  end
  errorRound(i) = sqrt(sum((M*currXround-y).^2)); % squared error
end

% find the minimal error
[junk,res] = min(errorRound);
discreteOutput = zeros(size(fractionalOutput));
discreteOutput(index(1:res)) = val(1:res);
