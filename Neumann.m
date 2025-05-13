
function Out = Neumann(n,x) %Spherical Bessel function of second kind
    % Bessel function of second kind: bessely
    Out = sqrt(pi./(2*x)).*bessely(n+0.5,x);
end
