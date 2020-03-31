clear all
close all
%% TO LOAD RECTANGULAR MESH
% load mesh_coarse
%% TO LOAD AIRFOIL MESH
[p,t] = tessellation('N0012_Airfoil.mat');
%% Material properties
E = 7e4; %young's modulus
nu = 0.33; %poissons ration
lambdael = E*nu/(1+nu)/(1-2*nu);
mu = E/2/(1+nu);
Fext = -1000;
%% Domain Properties
t=t(1:3,:)';
figure
simpplot(p',t)
x=p(1,:)';
y=p(2,:)';
x1=x(t(:,1));
x2=x(t(:,2));
x3=x(t(:,3));
y1=y(t(:,1));
y2=y(t(:,2));
y3=y(t(:,3));
AT=.5*(x1.*(y2-y3)+x2.*(y3-y1)+x3.*(y1-y2));
nt=length(AT);
b1=(y2-y3)./AT/2;
b2=(y3-y1)./AT/2;
b3=(y1-y2)./AT/2;
c1=(x3-x2)./AT/2;
c2=(x1-x3)./AT/2;
c3=(x2-x1)./AT/2;
%% Shape Functions
z=zeros(nt,1);
Nux=[b1 z b2 z b3 z];
Nvx=[z b1 z b2 z b3];
Nuy=[c1 z c2 z c3 z];
Nvy=[z c1 z c2 z c3];
dofs=[2*t(:,1)-1 2*t(:,1) 2*t(:,2)-1 2*t(:,2) 2*t(:,3)-1 2*t(:,3)];
n_dofs=max(max(dofs));
% A11b=zeros(6,6,nt);
A11=zeros(nt,36);
A22=zeros(nt,36);
A12=zeros(nt,36);

for i=1:nt
%   A11b(:,:,i)=Nux(i,:)'*Nux(i,:)+Nuy(i,:)'*Nuy(i,:);
    A11(i,:)=reshape(Nux(i,:)'*Nux(i,:)+Nuy(i,:)'*Nuy(i,:),1,36);
    A12(i,:)=0.5*reshape(Nux(i,:)'*Nvx(i,:)+Nuy(i,:)'*Nvy(i,:)+ Nvx(i,:)'*Nux(i,:)+Nvy(i,:)'*Nuy(i,:),1,36);
    A22(i,:)=reshape(Nvx(i,:)'*Nvx(i,:)+Nvy(i,:)'*Nvy(i,:),1,36);
end

%% Boundary Condition
kk=boundary(x,y,0.01);
x_bound=x(kk);
y_bound=y(kk);
figure
plot(x_bound,y_bound,'*'); grid on; hold off;
p_bound=[x_bound y_bound];

% minbound = min(x_bound) + 0.3*(1-min(x_bound)); % Because the max chord is always scaled to 1.
% maxbound = min(x_bound) + 0.55*(1-min(x_bound));
% boundind = [(minbound <= x_bound) & (x_bound <= maxbound)];
% x_bound = x_bound(boundind);
% y_bound = y_bound(boundind);

% force_node=find((abs(x_bound-100)<=1)&(abs(y_bound-30)<=1));
force_node=find(x_bound==1.0);
% imposed_node=find((x_bound >= 0.3) & (x_bound <= 0.55));
imposed_node=find((x_bound >= 0.3) & (x_bound <= 0.65));
% imposed_node=find(abs(x_bound-0)<=1);

force_dofs=2*force_node;
imposed_dofs= unique([2*imposed_node-1; 2*imposed_node]);
n_imposed_dofs=length(imposed_dofs);

% force_node = find((abs(x_bound-80)<=0)&(abs(y_bound-10)<=0));
% imposed_node =find(abs(x_bound-0)<=1);
% 
% forced_dofs = 3*forced_node;
% imposed_dofs = unique([2*imposed_node-1; 2*imposed_node]);
% n_imposed_dofs = length(imposed_dofs);


%Domain properties
nelx=length(x);
nely=length(y);
nodes_number= reshape(1:(1+nely)*(1+nelx),(nely+1),(nelx+1));
% dofVec= reshape(2* nodes_number(1:end-1,1:end-1)+1, nelx*nely,1);
% dofMat= repmat(dofVec,1,8) + repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely,1);

iK= kron(ones(1,6),dofs); % indices of assembling the matrices for the structure
jK= kron(dofs,ones(1,6));

%% Min compliance
q0= zeros(n_dofs,1);
% tic
pot=@(q) potential(Nux, Nuy, Nvx, Nvy, q, lambdael, mu, dofs, force_dofs, Fext, nt, A11, A12, A22, AT);
% toc

% tic
hess=@(q) hessian(Nux, Nuy, Nvx, Nvy, q, lambdael, mu, dofs, nt, A11, A12, A22, iK, jK);
% toc
% [ W0, grad0 ]=potential(Nux, Nuy, Nvx, Nvy, q0, lambda, mu, dofs, force_dofs, Fext, nt, A11, A12, A22, AT);
%% Potential minimization
options = optimoptions('fmincon','GradObj', 'on', 'Display','iter',...
    'Algorithm','Interior-point');%,'CheckGradient',true,'FiniteDifferenceType','central');
Aeq = sparse(imposed_dofs,1,1,n_dofs,1);
Aeq = diag(Aeq);
Aeq = Aeq(imposed_dofs,:);
beq = zeros(length(imposed_dofs),1);
q = fmincon(pot,q0,[],[],Aeq,beq,[],[],[],options); 
u =q(1:2:end-1);
v =q(2:2:end);
pp =p+[u v]';
figure
simpplot(pp',t); hold on; plot(x,y,'*'); hold off;

% options = optimoptions(@fmincon,'Algorithm','interior-point',...
%         'Display','iter','SpecifyObjectiveGradient',true,...
%         'HessianFcn',hess); %on github

% options = optimoptions('fmincon','Algorithm','interior-point','CheckGradients',false,'Display','none',...
%     'SpecifyObjectiveGradient', true,'HessianFcn', hess);

%% Fminunc usage 
q0= zeros(n_dofs,1);
% tic
pot=@(q) potential(Nux, Nuy, Nvx, Nvy, q, lambdael, mu, dofs, force_dofs, Fext, nt, A11, A12, A22, AT,imposed_dofs, iK, jK);
% toc

% tic
hess=@(q) hessian(Nux, Nuy, Nvx, Nvy, q, lambdael, mu, dofs, nt, A11, A12, A22, iK, jK, imposed_dofs);
% toc
% [ W0, grad0 ]=potential(Nux, Nuy, Nvx, Nvy, q0, lambda, mu, dofs, force_dofs, Fext, nt, A11, A12, A22, AT);
%% Fminunc - minimization
options = optimoptions('fminunc', 'Display', 'iter','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective','PlotFcns',@optimplotfval);

options.StepTolerance = 1e-9;
q =fminunc(pot,q0,options);

u1 =q(1:2:end-1);
v1 =q(2:2:end);
pp =p+[u1 v1]';
figure
simpplot(pp',t); hold on; plot(x,y,'*'); hold off;
%% Check by Finite differences(equilibrium point) - 
err=zeros(n_dofs,1);
q00=q;
[W0,grad0] = potential(Nux, Nuy, Nvx, Nvy, q00, lambdael, mu, dofs, force_dofs, Fext, nt, A11, A12, A22, AT);
for i=1:n_dofs
    if ~ismember(i,imposed_dofs)
        q=q00;
        q(i)=1.001*q00(i);
        dq=0.001*q00(i);
        [W,grad] = potential(Nux, Nuy, Nvx, Nvy, q, lambdael, mu, dofs, force_dofs, Fext, nt, A11, A12, A22, AT);
        err(i)=norm((W-W0)-grad0(i)*dq)/norm(W0);
    end
end
figure
plot(err)
%% Hessian function - FD Check
errhes = zeros(n_dofs,1);
q11 = q;
[H0] = hessian(Nux, Nuy, Nvx, Nvy, q11, lambdael, mu, dofs, nt, A11, A12, A22, iK, jK);
for i=1:n_dofs
    if ~ismember(i, imposed_dofs)
        q=q11;
        q(i) = 1.001*q11(i);
        dq = 0.001*q11(i);
        [H] = hessian(Nux, Nuy, Nvx, Nvy, q, lambdael, mu, dofs, nt, A11, A12, A22, iK, jK);
        errhes(i)=norm(H-(H0(i)*dq)/norm(H));
    end
end
figure
plot(errhes)
