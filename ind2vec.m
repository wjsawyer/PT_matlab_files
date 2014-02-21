%smoothen
function [X, Y] = ind2vec(A, width)
width = int16(width);
%defines two vecotrs such that the location of point p inside a is (X(p),
%Y(p))....
%width is number of columns
X = zeros(1,length(A)); %the column number
Y = X; % the row number
for i = 1: length(A)
  Y(i) = rem(A(i),width); %remainder
  X(i) = idivide(A(i), width, 'ceil'); %int divide, rounding up
end
    
    


