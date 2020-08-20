function x = snr(sig, ref)
% 
% mse = mean((ref(:)-sig(:)).^2);
% if mse == 0; x = inf; return; end
% 
% dv = var(ref(:),1);
% x = 10*log10(dv/mse);

persistent ref_save;
if nargin > 1 
    ref_save = ref; 
    %return; 
end

if nargin <2
    return
end

mse = mean((ref_save(:)-sig(:)).^2);
if mse == 0; x = inf; return; end

dv = var(ref_save(:),1);
x = 10*log10(dv/mse);