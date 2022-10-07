function y = rect(N, D)
	% function y = rect(x, D)
	if nargin == 1, D = 1; end
    if D > N, D=N; end
	y = zeros(N);
    y(N/2-D/2:N/2+D/2,N/2-D/2:N/2+D/2)=1;