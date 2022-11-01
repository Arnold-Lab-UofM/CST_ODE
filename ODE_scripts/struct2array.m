function a = struct2array(s)
%STRUCT2ARRAY Convert structure with doubles to an array.
...dpb snip...
% Convert structure to cell
c = struct2cell(s);
% Construct an array
a = [c{:}];