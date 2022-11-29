function [W,zv] = central_weights(al,bet,z0,zv)
num_data = length(zv);
Z = repmat(zv,1,num_data); Z(:,1) = 1; Z = cumprod(Z,2);
EI = @(zk,pow) zk.^(pow+bet+al+1).*(gamma(1+bet)*gamma(pow+al+1))./(gamma(pow+bet+2+al));

[POW,Zk]=meshgrid(0:(num_data-1),z0);
W = EI(Zk,POW)/Z;