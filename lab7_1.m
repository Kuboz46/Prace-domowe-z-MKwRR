%% RÓWNANIE TRANSPORTU
% u_t + a u_x = 0,       t0 < t < t1, x0 < x < x1
% u(0,x) = u0(x),
% + warunki brzegowe: jeœli a>0 to u(t,x0) = g(t),
% jeœli a<0 to u(t,x1) = g(t).

a = [1 -1];
u0 = @(x) (cos(pi*x)).^2.*(abs(x)<=0.5);
g = @(t) zeros(size(t));
x = [-3 3];
t = [0 2.4];

h = [1/10 1/20 1/40];
lambda = [4/5 8/5];
method = 'Lax-Friedrich';

%Kolejne wartoœci parametru a
for i = 1:2
    %Kolejne wartoœci parametru lambda
    for j=1:2
        %Kolejne d³ugoœci kroku przestrzennego h
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
%Dla a = 1, h ze zbioru {1/10, 1/20, 1/40} i lambda równego 4/5 lub 8/5 rozwi¹zanie jest stabilne oraz metoda jest zbie¿na. 
%Widzimy, ¿e jeœli k^(-1) h -> 0 (prawostronnie), to Lk,h(u) - L(u) -> 0, gdy k, h ->
%0(prawostronnie) (Lk,h(u(tn, xm)) jest schematem ró¿nicowym Laxa-Friedrichsa).
%Ale dla a = -1, h ze zbioru {1/10, 1/20, 1/40} i lambda równego 4/5 lub
%8/5 rozwi¹zanie nie jest stabilne oraz metoda nie jest zbie¿na.
%Zauwa¿my, ¿e jeœli k^(-1) h ->0 (prawostronnie), to Lk,h(u) - L(u) nie
%zbiega do zera przy k,h ->0 (prawostronnie).
%Zgadza siê z teori¹ podan¹ na wyk³adzie.
