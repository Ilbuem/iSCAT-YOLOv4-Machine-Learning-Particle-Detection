
function Out = Bessel(n,x) %Spherical Bessel function of first kind
    % Bessel function of first kind: besselj
    Out = sqrt(pi./(2*x)).*besselj(n+0.5,x);
end