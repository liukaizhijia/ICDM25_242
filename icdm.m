clearvars;close all;
dimension=400;
r=10;%W size is D*r
a=rand(dimension,5*dimension);
Sw=a*a';
b=3*rand(dimension,5*dimension);
Sb=b*b';
itr=6;
obj_alg3=zeros(itr,1);
eye_matrix=eye(dimension);
W_ini=eye_matrix(:,1:r);
W_3=W_ini;
for i=1:itr
    obj_alg3(i)=obj_value(Sw,Sb,W_3);
    [V,D]=eig(Sb-Sw*obj_alg3(i));
    eigvals = diag(D);
    [~, idx] = sort(eigvals, 'descend');                   
    W_3 = V(:, idx(1:r)); 
end
figure
semilogy(obj_alg3,"LineWidth",1.5,"Marker","+");
grid on;
function obj = obj_value(Sw,Sb,W)
    obj=trace(W'*Sb*W)/trace(W'*Sw*W);
end