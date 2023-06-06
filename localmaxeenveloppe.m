function B = localmaxeenveloppe(A,n)
%B = localmaxeenveloppe(A,n)
%Calculate the enveloppe of a variable by computing for each data point the
%maximum value in a radius of n elements. 
%INPUT:
%- A : data vectorforwhich the envelope must be calculated
%- n : window radius for the maximum search


B = zeros(size(A));
%center
for i = (n+1)/2:numel(A)-(n+1)/2
    B(i) = max(A(i-(n-1)/2:i+(n-1)/2));
end
%edges
for i = 1:(n-1)/2
B(i) = max(A(1:i+i-1));
B(end-i+1) = max(A(end-i+1:end));
end
end
