%Adds comma to separate digits for clarity. Shamelessly copied from an
%official reply to
%https://se.mathworks.com/matlabcentral/answers/96131-is-there-a-format-in-
%matlab-to-display-numbers-such-that-commas-are-automatically-inserted-
%into-the
%and slightly altered to avoid adding '.00' to the integers. 

function [str]=num2bank(num)
str = arrayfun(@(x) num2bankScalar(x) , num, 'UniformOutput', false);
end
function [str]=num2bankScalar(num)
%Include commented code to add '.00' to the end of the number. 

%num=floor(num*100)/100;
str = num2str(num);
%k=find(str == '.', 1);
% if(isempty(k))
%     str=[str,'.00'];
% end
%FIN = min(length(str),find(str == '.')-1);
FIN = length(str);
for i = FIN-2:-3:2
    str(i+1:end+1) = str(i:end);
    str(i) = ',';
end
end