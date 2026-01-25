function [ds,dsmag]=data2source(data,A)

[nchan,n]=size(data);
[nchan,ns,ndum]=size(A);
ds=zeros(ns,ndum,n);

for i=1:ndum
    ds(:,i,:)=(A(:,:,i))'*data;
end

dsmag=sqrt(sum(ds.^2,3));
return;
  