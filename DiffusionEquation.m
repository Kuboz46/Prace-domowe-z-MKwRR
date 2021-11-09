function [u, X, T, M, N] = DiffusionEquation(D, u0, fL, fR, boundary, t, x, mu, h)
% Rownanie 
%     u_t = (D(x)*u_x)_x
% z warunkiem poczatkowym
%     u(t0,x) = u0(x)
% W przypadku warunkow brzegowych rozpatrzymy wszystkie mozliwe kombinacje typow warunkow brzegowych 
% 1. War. Dirichleta na obu brzegach.
%     u(t,x0) = fL(t)
%     u(t,x1) = fR(t).
% 2. War. Dirichleta z lewej strony i Neumanna z prawej.
%     u(t,x0) = fL(t)
%     u_x(t,x1 = fR(t).
% 3. War. Neumanna z lewej strony i Dirichleta z prawej.
%     u_x(t,x0) = fL(t)
%     u(t,x1) = fR(t).
% 4. War. Neumanna na obu brzegach.
%     u_x(t,x0) = fL(t)
%     u_x(t,x1) = fR(t).
% Parametry:
% t = [t0, t1] - przedzial czasowy
% x = [x0, x1] - przedzial przestrzenny
% h - krok przestrzeny, mu = k/(h^2), gdzie k - krok czasowy.

k = mu*(h^2);

% punkty siatki (siatka przestrzenna - centralna)
X = (x(1)+h/2):h:(x(2)-h/2);    M = length(X);
T = t(1):k:t(2);                N = length(T);

% zainicjowanie macierzy rozwiazania (uwzgledniamy poza siatka)
u = zeros(M+2,N);

% uzupelnienie warunku poczatkowego dla n = 1

u(2:(M+1),1) = u0(X);

if (strcmp(boundary(1), 'Dirichlet') && strcmp(boundary(2), 'Dirichlet'))
    u(1,1) = 2*feval(fL,T(1)) - u(2,1);
    u(M+2,1) = 2*feval(fR,T(1)) - u(M+1,1);
elseif (strcmp(boundary(1), 'Dirichlet') && strcmp(boundary(2), 'Neumann'))
    u(1,1) = 2*feval(fL,T(1)) - u(2,1);
    u(M+2,1) = u(M+1,1) - h*feval(fR,T(1));
elseif (strcmp(boundary(1), 'Neumann') && strcmp(boundary(2), 'Dirichlet'))
    u(1,1) = u(2,1) - h*feval(fL,T(1));
    u(M+2,1) = 2*feval(fR,T(1)) - u(M+1,1);
else
    u(1,1) = u(2,1) - h*feval(fL,T(1));
    u(M+2,1) = u(M+1,1) - h*feval(fR,T(1));
end

% wspolczynniki ukladu rownan dla punktow wewnetrznych
alpha = -mu/2*feval(D,X-h/2);
beta = 1 + mu*(feval(D,X+h/2)+feval(D,X-h/2))/2;
gamma = -mu/2*feval(D,X+h/2);

if (strcmp(boundary(1), 'Dirichlet') && strcmp(boundary(2), 'Dirichlet'))
    a = [0, alpha, 0.5]; % dolna przekatna macierzy ukladu
    b = [0.5, beta, 0.5]; % glowna przekatna macierzy ukladu
    c = [0.5, gamma, 0]; % gorna przekatna macierzy ukladu
elseif (strcmp(boundary(1), 'Dirichlet') && strcmp(boundary(2), 'Neumann'))
    a = [0, alpha, 1/h]; % dolna przekatna macierzy ukladu
    b = [0.5, beta, -1/h]; % glowna przekatna macierzy ukladu
    c = [1/h, gamma, 0]; % gorna przekatna macierzy ukladu
elseif (strcmp(boundary(1), 'Neumann') && strcmp(boundary(2), 'Dirichlet'))
    a = [0, alpha, 0.5]; % dolna przekatna macierzy ukladu
    b = [-1/h, beta, 0.5]; % glowna przekatna macierzy ukladu
    c = [1/h, gamma, 0]; % gorna przekatna macierzy ukladu
else
    a = [0, alpha, 1/h]; % dolna przekatna macierzy ukladu
    b = [-1/h, beta, -1/h]; % glowna przekatna macierzy ukladu
    c = [1/h, gamma, 0]; % gorna przekatna macierzy ukladu
end

delta = zeros(1,M); 

for n = 1:N-1
    % wspolczynniki wektora "prawej strony" ukladu rownan
    for m = 1:M
        delta(m) = u(m,n)*(mu/2*feval(D,X(m)-h/2)) + ...
            u(m+1,n)*(1-mu/2*feval(D,X(m)+h/2)-mu/2*feval(D,X(m)-h/2)) + ...
            u(m+2,n)*(mu/2*feval(D,X(m)+h/2));
    end
    
    if (strcmp(boundary(1), 'Dirichlet') && strcmp(boundary(2), 'Dirichlet'))
        d = [feval(fL,T(n+1)), delta, feval(fR,T(n+1))];
    elseif (strcmp(boundary(1), 'Dirichlet') && strcmp(boundary(2), 'Neumann'))   
        d = [feval(fL,T(n+1)), delta, h*feval(fR, T(n+1))];
    elseif (strcmp(boundary(1), 'Neumann') && strcmp(boundary(2), 'Dirichlet'))
        d = [-h*feval(fL,T(n+1)), delta, feval(fR,T(n+1))];
    else
        d = [-h*feval(fL,T(n+1)), delta, h*feval(fR, T(n+1))];
    end    
    u(:,n+1) = ThomasSolve(a,b,c,d);
end
