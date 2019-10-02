function [W,w0,Y] = ica(X,W0,w00,max_it,lrate)
% [W,w0,Y] = ica(X,W0,w00,max_it,lrate)
% Function to perform Independent Component Analysis using the method of
% maximization of the entropy of a non-linear function (sigmoid) of W*x,
% where x is a column of X.   
% The method is  described in "An  information maximization approach  to
% blind separation and blind deconvolution," written by Anthony Bell and
% Terrence Sejnowski  (ICASSP-95);  and in  "An information-maximization
% approach to blind separation and blind deconvolution", also written by
% A.  Bell and T.  Seijnowski",(Neural Computation, vol.7, pp.1129-1159,
% 1995).                                                                 
% 
% Inputs:
% 	X  -> Matrix whose rows are vectors, each of them representing
% 	      linear combinations of independent sources.
% 	W0 -> Initial "guess" for the separation matrix W. (Optional)
% 	w00-> Initial "guess" for the bias vector w0. (Optional)
%	max_it-> Maximum number of iterations. (Optional. Default = 100)
%	lrate -> Learning rate. (Optional. Default = 0.01)
% Outputs:
% 	W  -> Separation matrix. 
% 	w0 -> Bias vector.
% 	Y  -> Matrix given by Y = W*X + w0
% 	      W and w0 are such that the components of Y= W*X+w0 are as 
%               independent as possible.

format compact
[N,M] = size(X); 		% M=17408, N=2, for example 
if N > M, X = X'; [N,M] = size(X); end % Signal vectors must be in rows of X
X = X./repmat(std(X')',1,M)*3;  % Set the row of X to have standard dev = 3; 
permute = randperm(M); 		% generate a permutation vector
Xp = X(:,permute);		% time-scrambled inputs for stationarity
% Initialize unmixing matrix and bias vector
idm = eye(N);		% Make identity matrix of dimension N
if nargin < 2,  W = idm;
  elseif isempty(W0) W = idm;
  else  W = W0; end
if nargin < 3,  w0 = 0*ones(N,1)/2;
  elseif isempty(w00) w0 = 0*ones(N,1)/2;
  else w0 = w00; end
if nargin < 4,  max_it = 200; end
if nargin < 5,  lrate = 0.01; end % Learning rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = 300;		% Block size
noblocks = fix(M/B);	% Number_of_blocks = Number_of_vectors/block_size
for i = 1:max_it	% 100 = number of iterations
  for t = 1:B:noblocks*B	% For each block
    xblock = Xp(:,t:t+B-1);	% Extract the block from X
    u = W*xblock+repmat(w0,1,B);% Compute u, a separation trial
    y = 1./(1+exp(-u));		% "Squash" u with sigmoid function
    yy = 1-2*y;			% Compute the derivative
%    dw = inv(W') + (yy*xblock')/B;%Get variation of W in the direction that
				   %makes the components of u more independent
    dw = (idm + yy*u'/B)*W;	% Use the natural gradient instead of the
				% "conventional" gradient used in the 
				% commented line.
    dw0 = mean(yy')';		% Get variation of the bias vector
    W_anterior = W;		% Save previous value of W
    W = W + lrate*dw;		% Get new W and w0 closer to the optimum W and
    w0 = w0 + lrate*dw0;	% w0 that makes the components of u as 
				% independent as possible
  end;
  if ~rem(i,10)			% Check separation performance
    Y = W*X+repmat(w0,1,M);	% Y is the vector of separated components
    Y_cov = cov(Y');		% Covariance matrix of Y
    %
    % If the rows of Y are independent, they are also uncorrelated. In such a
    % case, the product of the elements of the diagonal of Y_cov equals the 
    % determinant of Y_cov. 'perf' gives a distortion measure between the
    % the product of the elements of Y_cov and its determinant. The smaller
    % this distortion, the closer Y_cov is of a diagonal matrix. Since the
    % separation algorithm used here minimizes the dependence among the
    % components of Y, it must also minimize the covariance.
    perf = abs(1 - prod(diag(cov(Y')))/det(cov(Y')));
    % Finally, the variation observed in W in two successive iterations is
    % given by the Frobenius norm of the difference between the current W
    % and the previus W.
    aux = norm(W-W_anterior,'fro')/norm(W,'fro');
    fprintf(1,['Iteration # ',int2str(i),...
	'\n    Deviation from identity matrix: ',num2str(perf,4),...
	'\n    Variation in W from previous iteration: ',num2str(aux,4),'\n']);
    if aux < lrate/1000, break, end % Stop if W practically does not vary
  end
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
