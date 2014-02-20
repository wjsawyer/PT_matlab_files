%smoothen
function out = smooth(in)
out = bwmorph(in, 'majority', inf);
out = bwmorph(out, 'thin', inf);
%out = bwmorph(out, 'diag');




