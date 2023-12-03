function output = TGV2L1L2_ADMM(I,P,B,para)
%% Solve TGV2-(shearlet)L1-L2 model:
% $\frac\beta2\|Ku-b\|_2^2+\lambda\|SH(u)\|_1+TGV_\alpha^2(u)$
% In this version, we solve the (u,p)-subproblem directly.
%
% Inputs:
% I - ground truth image
% P - full sampling matrix
% B - full partial k-space data
% para - parameter structure
%     .maxiter
%     .alpha - 2X1 vector
%     .lambda/.beta/
%     .mu - 3X1 vector storing $\mu_1$, $\mu_2$ and $\mu_3$
%     .gamma - scalar, $\mu$ in the paper
%     .skrinkopt - shrink method
%     .period/.contfact - continuation parameters
%     .correctopt - correction option of results
%     .kernels/trans - sparsifying kernels/transform
%     .initial.* - initial intermediate results
%
% Outputs:
% output.y - reconstructed image
%       .err/.enderr/.x/.y/.z/.xx/.yy/.zz
%
% Refer to the paper ``A New Detail-preserving Regularity Scheme" by 
%     Weihong Guo, Jing Qin and Wotao Yin.
%
% Jing, updated on 07-16-2014

%% Read input parameters
global m n
[m,n]            = size(I);
maxiter          = para.maxiter;
alpha            = para.alpha;
lambda           = para.lambda;
beta             = para.beta;
mu               = para.mu;
gamma            = para.gamma;

a                = alpha(2)*mu(2);
b                = alpha(1)*mu(3);

% Read the sparsifying transform kernels in freq.domain
H                = para.kernels;
N                = size(H,3);

% Read the sparsifying transform
trans            = para.trans;

% Selection of shrinkage method
if isfield(para,'shrinkopt')
    method = para.shrinkopt;
else
    method = 'soft';
end

switch method
    case 'soft'
        shrinkmethod = @shrink;
    case 'nng' % nonnegative garotte shrinkage
        shrinkmethod = @nnshrink;
    case 'hyper'
        shrinkmethod = @hypershrink;
    case 'firm'
        shrinkmethod = @firmshrink;
end

% Selection of correction step
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
if isfield(para,'initial')
    u  = para.initial.u;
    p  = para.initial.p;
    x  = para.initial.x;
    y  = para.initial.y;
    z  = para.initial.z;
    xx = para.initial.xx;
    yy = para.initial.yy;
    zz = para.initial.zz;
else
    u  = zeros(m,n);
    p  = zeros(m,n,2);
    x  = zeros(m,n,N);
    y  = zeros(m,n,2);
    du = zeros(m,n,2);
    z  = zeros(m,n,4);
    xx = zeros(m,n,N);
    yy = zeros(m,n,2);
    zz = zeros(m,n,4);
end

uband = trans(u);

% Initialize the error
err = zeros(maxiter,1);
errpre = 1;

D1 = abs(psf2otf([-1,1],[n,n])).^2;
D2 = abs(psf2otf([-1;1],[n,n])).^2;
d1 = lambda*mu(1)*sum(abs(H).^2,3)+a.*(D1+D2)...
    +beta*P; % note that P'*P = P
d2 = a+b*(D1+.5*D2);
d3 = a+b*(D2+.5*D1);
d4 = psf2otf([1, -1],[n, n]);
d5 = psf2otf([1; -1],[n, n]);
d6 = d5(:,1)*d5(:,1)';

d4 = -a*d4;
d5 = -a*d5;
d6 = .5*d6;

d4t = conj(d4);
d5t = conj(d5);
d6t = conj(d6);

denom = d1.*d2.*d3+d4t.*d6t.*d5+d4.*d5t.*d6...
    -d2.*d5.*d5t-d1.*d6t.*d6-d3.*d4.*d4t;


%% Main loop
h = waitbar(0,'Please wait ...');
for i = 1:maxiter
    u2 = u;
    
    % Step 1: shrinkage of shearlets
    for j = 1:N
        x(:,:,j) = shrinkmethod(uband(:,:,j)+xx(:,:,j),1/mu(1));
    end
    
    % Step 2: shrinkage of gradients
    y = shrink2(du-p+yy,mu(2));
    
    % Step 3: shrinkage of second order derivatives
    Epp = Ep(p);
    z = shrink2(Epp+zz,mu(3));
    
    % Solve the (u,p1,p2)-subproblem
    RHS1 = 0;
    for j = 1:N
        RHS1 = RHS1 + H(:,:,j).*fft2(x(:,:,j)-xx(:,:,j));
    end
    FB1 = lambda.*mu(1).*RHS1 + a*(fft2(Dxt(y(:,:,1)-yy(:,:,1)))...
        +fft2(Dyt(y(:,:,2)-yy(:,:,2))))+beta*B; % note P'*B = B
    FB2 = a*fft2(yy(:,:,1)-y(:,:,1))+b*(fft2(Dxt(z(:,:,1)-zz(:,:,1)))...
        +fft2(Dyt(z(:,:,3)-zz(:,:,3))));
    FB3 = a*fft2(yy(:,:,2)-y(:,:,2))+b*(fft2(Dyt(z(:,:,2)-zz(:,:,2)))...
        +fft2(Dxt(z(:,:,3)-zz(:,:,3))));
    
    % Step 4: Solve for u
    RHS1 = (d2.*d3-d6.*d6t).*FB1-(d3.*d4t-d6.*d5t).*FB2...
        +(d4t.*d6t-d2.*d5t).*FB3;
    u = ifft2(RHS1./denom);
    u = correct(u);
    du(:,:,1) = Dx(u);
    du(:,:,2) = Dy(u);
    
    % Record the error
    err(i) = norm(u-I,'fro')/norm(I,'fro');

    % Step 5: Solve for p1
    RHS2 = (d5.*d6t-d3.*d4).*FB1+(d1.*d3-d5.*d5t).*FB2...
        +(d4.*d5t-d1.*d6t).*FB3;
    p(:,:,1) = ifft2(RHS2./denom);
        
    % Step 6: Solve for p2
    RHS3 = (d4.*d6-d2.*d5).*FB1+(d4t.*d5-d1.*d6).*FB2...
        +(d1.*d2-d4.*d4t).*FB3;
    p(:,:,2) = ifft2(RHS3./denom);


    % Step 7: Bregman update
    uband = trans(u);
    xx = xx + gamma*(uband-x);
    yy = yy + gamma*(du-p-y);
    zz = zz + gamma*(Ep(p)-z);
    
    % Stopping criterion
    output.relerr = norm(u(:)-u2(:))/norm(u(:));
    if output.relerr < 1e-5 
        break
    elseif errpre < err(i)
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
output.x        = x;
output.y        = y;
output.z        = z;
output.xx       = xx;
output.yy       = yy;
output.zz       = zz;
output.shsum    = sum(sum(sum(abs(uband),3)));
output.tv1      = sum(sum(sum(abs(du-p),3)));
output.tv2      = sum(sum(sum(abs(Ep(p)),3)));
output.fid      = norm(P.*fftshift(fft2(u))-B,'fro');

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
