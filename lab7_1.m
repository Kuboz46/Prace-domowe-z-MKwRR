%% R�WNANIE TRANSPORTU
% u_t + a u_x = 0,       t0 < t < t1, x0 < x < x1
% u(0,x) = u0(x),
% + warunki brzegowe: je�li a>0 to u(t,x0) = g(t),
% je�li a<0 to u(t,x1) = g(t).

a = [1 -1];
u0 = @(x) (cos(pi*x)).^2.*(abs(x)<=0.5);
g = @(t) zeros(size(t));
x = [-3 3];
t = [0 2.4];

h = [1/10 1/20 1/40];
lambda = [4/5 8/5];
method = 'Lax-Friedrich';

%Kolejne warto�ci parametru a
for i = 1:2
    %Kolejne warto�ci parametru lambda
    for j=1:2
        %Kolejne d�ugo�ci kroku przestrzennego h
        for k = 1:3
[u, X, T, M, N] = NumericSolve(a(i), u0, g, t, x, lambda(j), h(k), method);
%%
for n = 1:N
    plot(X, u(:,n),'g', X, u0(X - a(i)*T(n)), 'r')
    ylim([-0.1 1.6])
    xlabel('x'); ylabel('u(t,x)');
    title(['time: ', num2str(T(n), '%.2f'), ' h= ', num2str(h(k)), ' lambda ', num2str(lambda(j)), ' a = ', num2str(a(i))])
    legend(method,'solution')
    pause(0.2/k)
end

        end
    end
end

%%Wnioski
%Dla a = 1, h ze zbioru {1/10, 1/20, 1/40} i lambda r�wnego 4/5 lub 8/5 rozwi�zanie jest stabilne oraz metoda jest zbie�na. 
%Widzimy, �e je�li k^(-1) h -> 0 (prawostronnie), to Lk,h(u) - L(u) -> 0, gdy k, h ->
%0(prawostronnie) (Lk,h(u(tn, xm)) jest schematem r�nicowym Laxa-Friedrichsa).
%Ale dla a = -1, h ze zbioru {1/10, 1/20, 1/40} i lambda r�wnego 4/5 lub
%8/5 rozwi�zanie nie jest stabilne oraz metoda nie jest zbie�na.
%Zauwa�my, �e je�li k^(-1) h ->0 (prawostronnie), to Lk,h(u) - L(u) nie
%zbiega do zera przy k,h ->0 (prawostronnie).
%Zgadza si� z teori� podan� na wyk�adzie.
