% Normalize matrix rows 
function [M, z] = normalizeRows(A, dim)

z=sum(A,dim);
s = z + (z==0);
L=size(A,dim);
d=length(size(A));
v=ones(d,1);
v(dim)=L;
c=repmat(s,v');
M=A./c;

