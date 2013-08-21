function phi=capon(y,m,L)
%
% The Capon spectral estimator.
%
% phi=capon(y,m,L);
%
%    y   <- the data vector (length N)
%    m   <- the length of the Capon filter
%    L   <- the number of estimated spectral samples
%    phi -> the estimated spectrum
%

% Copyright 1996 by R. Moses

y=y(:);
N=length(y);       % data length

% form the sample covariance matrix
R=zeros(m+1,m+1);
for i = m+1 : N,
   R=R+y(i:-1:i-m)*y(i:-1:i-m)';
end
R=R/(N-m);

% compute the inverse of R
IR=inv(R);

% compute the spectrum
phi=zeros(L,1);
for k = 1 : L, 
   a=exp(- j*2*pi*(k-1)/L*[0:m].');      % form the a(w) vector
   phi(k)=real(a'*IR*a);
end
phi=(m+1)./phi;
