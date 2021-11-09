function [u,X,T,M,N] = HangingChain(g,nu,v0, v1, t, x, lambda, h)
% Rozwi¹zanie analityczne równania 
%     u_tt = g*(x * u_x)_x - nu*u_t
% z warunkami pocz¹tkowymi 
%     u(t0,x) = v0(x)
%     u_t(t0,x) = v1(x)
% i warunkami brzegowymi:
% u(t,L) = v0(x) dla t>0
% u(t,0) < infinity dla t > 0
% Parametry: 
% t = [0, t0] - przedzia³ czasowy
% x = [0, L] - przedzia³ przestrzenny
% h - krok przestrzenny, lambda = h/k, gdzie k - krok czasowy.

k = lambda * h;
X = (x(1) + h/2):h:(x(2) - h/2);    M = length(X);
T = t(1):k:t(2);                    N = length(T);

u = zeros(M+2,N);

% uzupe³nienie warunku pocz¹tkowego dla n = 1:
for m = 2:(M+1)
u(m,1) = feval(v0,X(m-1));
end

%Na lewym brzegu musimy wykorzystaæ warunek brzegowy u(0,n) = 2 u(1,n) -
%u(2,n)
%Teraz robimy dla n = 1.
u(1,1) = 2 * u(2,1) - u(3,1);
%Na prawym brzegu skorzystamy z zerowego warunku brzegowego (przybli¿enia
%wartoœci (u(M,n) + u(M + 1,n))/2 do zera)
%Teraz robimy dla n = 1.
u(M + 2,1) = - u(M + 1,1);

for m = 2:(M + 1)
    u(m,2) = 2 * u(m,1) + 2 * k * feval(v1, X(m-1)) - ((k^2)/(2*h))*g*(u(m+1,1) - u(m-1,1))+...
        lambda^2 * g * X(m-1)*(u(m+1,1)-2*u(m,1) + u(m-1,1)) -...
        nu * feval(v1,X(m-1)); 
end

%Korzystamy z wspomnianych wczeœniej warunków brzegowych dla n = 2.
u(1,2) = 2 * u(2,2) - u(3,2);
u(M + 2,2) = - u(M + 1,2);


% g³ówna pêtla
for n = 2:(N-1)
    for m = 2:(M-1)
        u(m, n + 1) = (2*u(m,n) - u(m, n - 1) +...
            ((k^2)/h)*g*(u(m+1,n)-u(m-1,n)) +...
            g * X(m-1) * lambda^2 * (u(m + 1) -...
            2 * u(m,n) + u(m-1,n)) - ...
            (nu/(2*k))*u(m,n-1))/((1 - (nu * k^2))/(2*k));
    end
    % Warunki brzegowe
    u(1,n + 1) = 2 * u(2,n + 1) - u(3,n + 1); 
    u(M + 2,n + 1) = - u(M + 1,n + 1);   
end