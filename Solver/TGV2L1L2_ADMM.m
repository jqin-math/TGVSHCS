function output = TGV2L1L2_ADMM(I,P,B,para)
%% Solve TGV2-L1-L2 model:
% $\frac\beta2\|Ku-b\|_2^2+\lambda\|T(u)\|_1+TGV_\alpha^2(u)$
% Note: T is a general sparsifying transform, which is defined by
%       the input variables: para.kernels and para.trans.
%
% Inputs:
% I - ground truth image
% P - full sampling matrix
% B - full partial k-space data
% para - parameter structure
%     .maxiter
%     .alpha - 2X1 vector
%     .lambda/.beta/
%     .mu - 3X1 vector
%     .gamma - 3X1 vector
%     .correctopt - correction option of results
%     .kernels/trans - sparsifying kernels/transform
%     .initial.* - initial intermediate results
%
% Outputs:
% output.u - reconstructed image
%       .err - relative errors for all iterations
%       .enderr - final relative error
%
% Refer to the paper ``A New Detail-preserving Regularity Scheme".
% For FFST toolbox, please visit http://www.mathematik.uni-kl.de/imagepro/members/haeuser/ffst/ 
% Jing, 4-26-2013

%% Read input parameters
global m n
[m,n]            = size(I);
maxiter          = para.maxiter;
alpha            = para.alpha;
lambda           = para.lambda;
beta             = para.beta;
mu               = para.mu;
gamma            = para.gamma;

a                = alpha(1)*mu(2);
b                = alpha(2)*mu(3);

% Read the sparsifying transform kernels in freq.domain
H                = para.kernels;
N                = size(H,3);

% Read the sparsifying transform
trans            = para.trans;

% Selection of the correction step
if isfield(para,'correctopt')
    correctopt = para.correctopt;
else
    correctopt = 'no';
end

switch correctopt
    case 'real'
        correct = @real;
    case 'abs'
        correct = @abs;
    case 'no'
        correct = @(x)(x);
end

%% Initialize
u  = zeros(m,n);
p  = zeros(m,n,2);
x  = zeros(m,n,N);
du = zeros(m,n,2);
xx = zeros(m,n,N);
yy = zeros(m,n,2);
zz = zeros(m,n,4);

uband = trans(u);

% Initialize the error
err = zeros(maxiter,1);
errpre = 1;

d = zeros(1,n);
d(1) = 2; d(2) = -1; d(n) = -1;
e = ones(m,1);
D1 = abs(kron(fft(d),e));
D2 = D1';
denom1 = D1+D2;
denom2 = D1+D2./2;
denom3 = D1./2+D2;

%% Main loop
h = waitbar(0,'Please wait ...');
for i = 1:maxiter
    u2 = u;
    
    % Step 1: shrinkage of shearlets
    for j = 1:N
        x(:,:,j) = shrink(uband(:,:,j)+xx(:,:,j),1/mu(1));
    end
    
    % Step 2: shrinkage of gradients
    y = shrink2(du-p+yy,1/mu(2));
    
    % Step 3: shrinkage of second order derivatives
    Epp = Ep(p);
    z = shrink2(Epp+zz,1/mu(3));
    
    % Step 4: Solve for u
    RHS1 = 0;
    for j = 1:N
        RHS1 = RHS1 + H(:,:,j).*fft2(x(:,:,j)-xx(:,:,j));
    end
    RHS1 = lambda.*mu(1).*RHS1 + a*(fft2(Dxt(p(:,:,1)+y(:,:,1)-yy(:,:,1)))...
        +fft2(Dyt(p(:,:,2)+y(:,:,2)-yy(:,:,2))))+beta*B;
    LHS1 = lambda*mu(1).*sum(abs(H).^2,3)+a*denom1+beta*P;
    u = ifft2(RHS1./LHS1);
    u = correct(u);
    du(:,:,1) = Dx(u);
    du(:,:,2) = Dy(u);
    
    % Record the error
    err(i) = norm(u-I,'fro')/norm(I,'fro');

    % Step 5: Solve for p1
    RHS2 = a*(du(:,:,1)-y(:,:,1)+yy(:,:,1))+...
        b*Dxt(z(:,:,1)-zz(:,:,1))...
        +b*Dyt(-.5*Dx(p(:,:,2))+z(:,:,3)-zz(:,:,3));
    LHS2 = a+b*denom2;
    p(:,:,1) = ifft2(RHS2./LHS2);
        
    % Step 6: Solve for p2
    RHS3 = a*(du(:,:,2)-y(:,:,2)+yy(:,:,2))+...
        b*Dyt(z(:,:,2)-zz(:,:,2))...
        +b*Dxt(-.5*Dy(p(:,:,1))+z(:,:,3)-zz(:,:,3));
    LHS3 = a+b*denom3;
    p(:,:,2) = ifft2(RHS3./LHS3);

    % Step 7: Bregman update
    uband = trans(u);
    xx = xx + gamma(1)*(uband-x);
    yy = yy + gamma(2)*(du-p-y);
    zz = zz + gamma(3)*(Ep(p)-z);
    
    % Stopping criterion
    output.relerr = norm(u(:)-u2(:))/norm(u(:));
    if output.relerr < 1e-5 
        fprintf('Iter %d, err %4.2e%%\n',i,output.relerr)
        break
    elseif errpre < err(i)
        fprintf('Iter %d, errpre < errcur\n',i)
        break
    else 
        errpre = err(i);
    end
    
    % Display the progress
    waitbar(i/maxiter,h,sprintf('iteration %d',i));
           
end
close(h)

if i<maxiter
    err(i+1:end) = [];
end

%% Output results
output.u        = abs(u);
output.err      = err;
output.enderr   = err(end);

function z = Ep(p)
global m n
z = zeros(m,n,4);
z(:,:,1) = Dx(p(:,:,1));
z(:,:,2) = Dy(p(:,:,2));
z(:,:,3) = (Dy(p(:,:,1))+Dx(p(:,:,2)))./2;
z(:,:,4) = z(:,:,3);

function dxu = Dx(U)
% x-axis forward difference
dxu = U(:,[2:end,1])-U;

function dyu = Dy(U)
% y-axis forward difference
dyu = U([2:end,1],:)-U;           

function dxtu = Dxt(U)
% -(x-axis backward difference)
dxtu = U(:,[end,1:end-1])-U;

function dytu = Dyt(U)
% -(y-axis backward difference)
dytu = U([end,1:end-1],:)-U;
