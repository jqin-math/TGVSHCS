function Psi = scalesShearsAndSpectra ( l, numOfScales, shearlet_spect, shearlet_arg)
% computes the spectra of the shearlet for given dimensions l and number of
% scales numOfscales with optionally given name of shearlet
% (shearlet_spect) with respective shearlet argument (shearlet_arg).
% The output Psi is a 3-d-matrix with the shearlets in the fourier domain
% ordered with ascending scale and within each scale ordered by the
% direction of the respective shears (see comments below for further details).
%
% INPUT:
%  l				(vector) dimensions of image
%  numOfScales		(int) number of scales
%  shearlet_spect	(string) name of shearlet(spectrum) (optional)
%  shearlet_arg		(var) shearlet argument (optional)
%
% OUTPUT:
%  Psi				(3-d-matrix) spectrum of shearlets
%
%--------------------------------------------------------------------------
% 2012-01-20, v1.0, (c) Sören Häuser

	%% initialization / error trapping
	if(nargin < 4)
		shearlet_arg = [];
	end

	if(nargin < 3)
		shearlet_spect = 'meyerShearletSpect';
	end
	
	%compute number of scales
	if(nargin<2)
		numOfScales = floor(0.5 * log2(max(l)));
		if(numOfScales<1)
			error('image too small!');
		end
	end
	
	%% new construction
% 	new = 0;
    
    %% computation of the shearlets
    
    %for better symmetrie each l should be odd
    l_orig = l;
	lm = logical(1-mod(l,2)); 
	l(lm) =  l(lm) + 1;
	
	% create meshgrid
	%largest value where psi_1 is equal to 1
	X = 2^(2*(numOfScales-1)+1); % = 2^(2*numOfScales - 1)
	%xi_x = linspace(-X,X-1/l(2)*2*X,l(2)); %not exactly symmetric
    xi_x_init = linspace(0,X,(l(2)+1)/2);
	xi_x_init = [-fliplr(xi_x_init(2:end)) xi_x_init];
	[xi_x, xi_y] = meshgrid(xi_x_init,fliplr(xi_x_init));
	
	% number of shears: |-2^j,...,0,...,2^j| = 2 * 2^j + 1 
	% now: inner shears for both cones: 
	% |-(2^j-1),...,0,...,2^j-1| 
	% = 2 * (2^j - 1) + 1
	% = 2^(j+1) - 2 + 1 = 2^(j+1) - 1
	% outer scales: 2 ("one" for each cone)
	% shears for each scale: hor: 2^(j+1) - 1, ver: 2^(j+1) - 1, diag: 2
	%  -> hor + ver + diag = 2*(2^(j+1) - 1) +2 = 2^(j + 2)
	%  + 1 for low-pass
	shearsPerScale = 2.^((0:numOfScales-1)+2);
	numOfAllShears = 1 + sum(shearsPerScale);

	%init
	Psi = zeros([l, numOfAllShears]);
	% frequency domain:
	% k  2^j 0 -2^j
	%
	%     4  3  2  -2^j
	%      \ | /
	%   (5)- x -1  0
	%      / | \
	%              2^j
	%
	%        [0:-1:-2^j][-2^j:1:2^j][2^j:-1:1] (not 0)
	%           hor          ver        hor
	%
	% start with shear -2^j (insert in index 2^j+1 (with transposed
	% added)) then continue with increasing scale. Save to index 2^j+1 +- k,
	% if + k save transposed. If shear 0 is reached save -k starting from
	% the end (thus modulo). For + k just continue.
	% 
	% then in time domain:
	%
	%  2  1  8
	%   \ | /
	%  3- x -7
	%   / | \
	%  4  5  6
	%

	%lowpass
	lp = feval(shearlet_spect, xi_x, xi_y, NaN, NaN, 'scaling');
	Psi(:,:,1) = lp;

	%cone indicator function for appropiate cut-off for the lowest and
	%biggest shear
	H = coneIndicator(l);

	%loop for each scale
	for j=0:numOfScales-1
		%starting index
		idx = 2^j + 1;
		start_index = 1 + sum(shearsPerScale(1:j));
		shift = 1;
		for k = -2^j:1:2^j
			%shearlet spectrum
			P = feval(shearlet_spect, xi_x, xi_y, 2^(-2*j), k*2^(-j), shearlet_arg);
% 			if (new)
% 				P = shearletSpect2( shearlet_spect, 2^(-2*j), k*2^(-j), xi_x, xi_y, shearlet_arg);
% 			else
% 				P = shearletSpect( shearlet_spect, 2^(-2*j), k*2^(-j), xi_x, xi_y, shearlet_arg);
% 			end
			if( k==-2^j )
				Psi(:,:,start_index + idx) = P .* double(H.hor) + P' .* double(H.ver);
			elseif( k==2^j )
				Psi(:,:,start_index + idx+shift) = P .* double(H.hor) + P' .* double(H.ver);
			else
				new_pos = mod(idx-shift,shearsPerScale(j+1));
				if(new_pos == 0)
					new_pos = shearsPerScale(j+1);
				end
				Psi(:,:,start_index + new_pos) = P;
				Psi(:,:,start_index + idx + shift) = P';

				%update shift
				shift = shift + 1;
			end
		end
	end
    
    %generate output with size l
	Psi = Psi(1:l_orig(1), 1:l_orig(2), :);
	
end

