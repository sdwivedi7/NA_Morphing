function [ W, grad, Hess ] = potential(Nux, Nuy, Nvx, Nvy, q, lambda, mu, dofs, force_dofs, Fext, nt, A11, A12, A22, AT, varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if isempty(varargin) == 0
    q(varargin{1})=0;
end
qq=q(dofs);
ux=sum(Nux.*qq,2);
uy=sum(Nuy.*qq,2);
vx=sum(Nvx.*qq,2);
vy=sum(Nvy.*qq,2);
%% Strain & stress
E11=ux+(ux.^2+uy.^2)/2;
E22=vy+(vx.^2+vy.^2)/2;
E12=(uy+vx+ux.*vx+uy.*vy)/2;
S11=(lambda+2*mu)*E11+lambda*E22;
S22=(lambda+2*mu)*E22+lambda*E11;
S12=2*mu*E12;
%% Derivatives
dW_dE11=S11;
dW_dE22=S22;
dW_dE12=S12;
dE11_dq=Nux;
dE22_dq=Nvy;
dE12_dq=(Nuy+Nvx)/2;
%% Potential
W_int=sum((E11.*S11+E12.*S12+E22.*S22))/2;
W_ext=sum(Fext'*q(force_dofs));
W=W_int-W_ext;
%% Gradient
F_int=zeros(length(q),1);
F_ext=F_int;
F_ext(force_dofs)=Fext;
for i=1:nt
    dE11_dq(i,:)=dE11_dq(i,:)+qq(i,:)*reshape(A11(i,:),6,6);
    dE12_dq(i,:)=dE12_dq(i,:)+qq(i,:)*reshape(A12(i,:),6,6);
    dE22_dq(i,:)=dE22_dq(i,:)+qq(i,:)*reshape(A22(i,:),6,6);
end
dW_dq=repmat(dW_dE11,1,6).*dE11_dq+repmat(dW_dE12,1,6).*dE12_dq+repmat(dW_dE22,1,6).*dE22_dq;
for i=1:nt
    F_int(dofs(i,:))=F_int(dofs(i,:))+dW_dq(i,:)';
end
grad=F_int-F_ext;
%%
if isempty(varargin) == 0
    Hess = hessian(Nux, Nuy, Nvx, Nvy, q, lambda, mu, dofs, nt, A11, A12, A22, varargin{2},varargin{3});
end
end
