function [ H ] = hessian( Nux, Nuy, Nvx, Nvy, q, lambda, mu, dofs, nt, A11, A12, A22, iK, jK,varargin)
%UNTITLED5 Summary of this function goes here
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
dE11_dq=Nux;
dE22_dq=Nvy;
dE12_dq=(Nuy+Nvx)/2;
d2E11_dqdq=A11;
d2E22_dqdq=A22;
d2E12_dqdq=A12;
d2W_dqdq11=zeros(nt,36);
d2W_dqdq12=d2W_dqdq11;
d2W_dqdq22=d2W_dqdq11;
for i=1:nt
    dE11_dq(i,:)=dE11_dq(i,:)+qq(i,:)*reshape(A11(i,:),6,6);   %accelerate
    dE12_dq(i,:)=dE12_dq(i,:)+qq(i,:)*reshape(A12(i,:),6,6);
    dE22_dq(i,:)=dE22_dq(i,:)+qq(i,:)*reshape(A22(i,:),6,6);
    k=kron(dE11_dq(i,:)',dE11_dq(i,:));
    d2W_dqdq11(i,:)=k(:);
    k=kron(dE12_dq(i,:)',dE12_dq(i,:));
    d2W_dqdq12(i,:)=k(:);
    k=kron(dE22_dq(i,:)',dE22_dq(i,:));
    d2W_dqdq22(i,:)=k(:);
end
% d2W_dE11dq=(lambda+2*mu)*dE11_dq+lambda*dE22_dq;
% d2W_dE22dq=(lambda+2*mu)*dE22_dq+lambda*dE11_dq;
% d2W_dE12dq=2*mu*dE12_dq;
d2W_dqdq=repmat(S11,1,36).*d2E11_dqdq+repmat(S22,1,36).*d2E22_dqdq+repmat(S12,1,36).*d2E12_dqdq;
sK=d2W_dqdq+d2W_dqdq11+d2W_dqdq12+d2W_dqdq22;
% K=sparse(length(q),length(q));
% for i=1:nt
%     d=dofs(i,:);
%     K(d,d)=K(d,d)+reshape(sK,6,6);
% end
K=sparse(iK,jK,sK);
H=(K+K')/2;
end

