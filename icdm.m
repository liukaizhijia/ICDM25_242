% this code is to implement alg2 and alg3 to minize 
% trace(W'SwW)/trace(W'SbW)with W'W=I constraint
% this is same as Eq. (3) in the submitted paper
% the 2D case follows similarly
clearvars;close all;
dimension=400;
r=10;%W size is D*r
a=rand(dimension,5*dimension);
Sw=a*a';
b=3*rand(dimension,5*dimension);
Sb=b*b';
itr=60;
obj_alg3=zeros(itr,1);
obj_alg2=zeros(itr,1);
eye_matrix=eye(dimension);
W_ini=eye_matrix(:,1:r);
W_3=W_ini;
for i=1:itr
    obj_alg3(i)=obj_value(Sw,Sb,W_3);
    [V,D]=eig(Sw-Sb*obj_alg3(i));
    eigvals = diag(D);
    [~, idx] = sort(eigvals, 'ascend');                   
    W_3 = V(:, idx(1:r)); 
end
W_2=W_ini;
for i=1:itr
    obj_alg2(i)=obj_value(Sw,Sb,W_2);
    grad=(Sw-Sb*obj_alg2(i))*W_2;
    step=1/norm(Sw-Sb*obj_alg2(i));
    [U,~,V]=svd(W_2-step*grad,0);
    W_2=U*V';
end
figure
plot(obj_alg2,"LineWidth",1.5);
hold on;
plot(obj_alg3,"LineWidth",1.5);
legend('alg2','alg3')
function obj = obj_value(Sw,Sb,W)
    obj=trace(W'*Sw*W)/trace(W'*Sb*W);
end