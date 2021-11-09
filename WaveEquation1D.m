function [u, X, T, M, N] = ...
    WaveEquation1D(c, F, u0, v0, fL, fR, boundary, t, x, lambda, h)
% Rownanie u_tt = c^2 u_xx + F(t, x)
% z warunkiem poczatkowym:
% u(t0, x) = u0(x)
% u_t(t0, x) = v0(x)
% i warunkiem brzegowym fL(t) i fR(t) typu Dirichleta lub Neumanna na
% lewym i prawym brzegu.
% Parametry:
% boundary: 'Dirichlet', 'Neumann'
% t = [t0, tk] - przedzia³ czasowy
% x = [x0, xk] - przedzia³ przestrzenny
% h - krok przestrzenny, lambda = k/h, gdzie k - krok czasowy

k = lambda * h;
X = (x(1) + h/2):h:(x(2) - h/2);    M = length(X);
T = t(1):k:t(2);                    N = length(T);

u = zeros(M,N);

% uzupelnienie warunku poczatkowego dla n = 1:
u(:,1) = feval(u0,X);
if strcmp(boundary, 'Dirichlet')
    uL = 2*feval(fL,T(1)) - u(1,1);
    uR = 2*feval(fR,T(1)) - u(M,1);
elseif strcmp(boundary, 'Neumann')
    uL = u(1,1) - h*feval(fL,T(1));
    uR = u(M,1) + h*feval(fR,T(1));
end

% uzupelnienie warunku poczatkowego dla n = 2:
for m = 2: (M - 1)
    u(m,2) = u(m,1) +k*feval(v0,X(m)) + ...
        c^2*lambda^2/2*(u(m+1,1) - 2 * u(m,1) + u(m-1,1)) + k^2*F(T(1),X(m));
end
    % wartosci na brzegu (korzystaja z uL i uR)
    u(1,2) = u(1,1) +k*feval(v0,X(1)) + ...
        c^2*lambda^2/2*(u(1+1,1) - 2 * u(1,1) + uL)+k^2*F(T(1),X(1));
    u(M,2) = u(M,1) +k*feval(v0,X(M)) + ...
        c^2*lambda^2/2*(uR - 2 * u(M,1) + u(M-1,1))+k^2*F(T(1),X(M));
    
%   nowe wartosci uL i uR
if strcmp(boundary, 'Dirichlet')
    uL = 2*feval(fL,T(2)) - u(1,2);
    uR = 2*feval(fR,T(2)) - u(M,2);
elseif strcmp(boundary, 'Neumann')
    uL = u(1,2) - h*feval(fL,T(2));
    uR = u(M,2) + h*feval(fR,T(2));
end

% glowna petla
for n = 2:(N-1)
    for m = 2:(M-1)
        u(m,n+1) = 2*u(m,n) - u(m, n - 1) + ...
            c^2*lambda^2*(u(m+1,n) - 2*u(m,n) + u(m-1,n)) + k^2*F(T(n),X(m));
    end
    % wartosci na brzegu (korzystaja z uL i uR)
    u(1,n+1) = 2*u(1,n) - u(1, n - 1) + ...
            c^2*lambda^2*(u(1+1,n) - 2*u(1,n) + uL) + k^2*F(T(n),X(1));
    u(M,n+1) = 2*u(M,n) - u(M, n - 1) + ...
            c^2*lambda^2*(uR - 2*u(M,n) + u(M-1,n)) + k^2*F(T(n),X(M));   
    % nowe wartosci uL i uR
    if strcmp(boundary, 'Dirichlet')
        uL = 2*feval(fL, T(n+1)) - u(1,n+1);
        uR = 2*feval(fR,T(n+1)) - u(M,n+1);
    elseif strcmp(boundary, 'Neumann')
        uL = u(1, n+1)- h*feval(fL,T(n+1));
        uR = u(M,n+1) + h*feval(fR,T(n+1));
    end
end