function [u, X, T, M, N] = NumericSolve(a, u0, g, t, x, lambda, h, method)
% Funkcja rozwi¹zuj¹ca numerycznie zagadnienie dla r.
% transportu.

k = lambda*h;
X = x(1):h:x(2); M = length(X);
T = t(1):k:t(2); N = length(T);

u = zeros(M,N); 

u(:,1) = feval(u0, X); %Tu mo¿na daæ feval(u0, X)', ale nie ma to znaczenia.
if a>0
    u(1,:) = feval(g,T);
else
    u(M,:) = feval(g,T);
end

if strcmp(method,'forward-time forward-space')
    for n = 1:(N-1)
        for m = 2:(M-1)
            u(m,n+1) = u(m,n) - a*lambda*(u(m+1,n) - u(m,n));
        end
        
        if a>0
            u(M,n+1) = u(M - 1,n + 1); % warunek numeryczny
        else
            u(1,n+1) = u(1,n) - a*lambda*(u(1+1,n) - u(1,n)); % ??
        end
    end
elseif strcmp(method, 'forward-time backward-space')
    for n=1:(N-1)
        for m = 2:(M-1)
            u(m,n+1) = u(m,n) - a*lambda*(u(m,n) - u(m-1,n));
        end
        
        if a>0
            u(M,n+1) = u(M,n) - a*lambda*(u(M,n) - u(M-1,n));
        else
            u(1,n+1)=u(2,n+1);
        end
    end
elseif strcmp(method, 'forward-time central-space')
    for n=1:(N-1)
        for m=2:(M-1)
            u(m,n+1) = u(m,n) - a*lambda*(u(m+1,n) - u(m - 1,n))/2;
        end
        
        if a>0
            u(M,n+1) = u(M-1,n+1);
        else
            u(1,n+1)=u(2,n+1);
        end
    end
elseif strcmp(method, 'Lax-Friedrich')
    for n=1:(N-1)
        for m = 2:(M-1)
            u(m,n+1) = (u(m+1,n) - u(m-1,n))/2 - a*lambda*(u(m+1,n)-u(m-1,n))/2;
        end
        
        if a>0
            u(M,n+1)=u(M-1,n);
        else
            u(1,n+1)=u(2,n);
        end
    end
end