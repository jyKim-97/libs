function pseudo_x = get_surrogate_data(x)
n = length(x);
if mod(n, 2) == 1
    % odd
    n2 = (n - 1)/2;
else
    % even
    n2 = n/2 - 1;
end
y = fft(x);
pseudo_y = abs(y);
% uniform distrib
% Wikipedia, Surrogate data testing, Algorithm 1
uni_angle = exp(2j*pi*rand(1, n2)); 
pseudo_y(2:n2+1) = pseudo_y(2:n2+1) .* uni_angle;
pseudo_y(end-n2+1:end) = pseudo_y(end-n2+1:end) .* flip(conj(uni_angle));
pseudo_x = ifft(pseudo_y);