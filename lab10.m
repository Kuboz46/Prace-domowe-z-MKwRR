g = 9.81; mu = 0.5;
L = 2; t0 = 10;
v0 = @(x) (x - L).^2 / 4;% 0; % (sin(x* (pi/2)/(L/6))-1)/10.*abs(x<=L/6); % (x - L)/3; %
v1 = @(x)0; % (sin(x*(pi/2)/(L/6)) - 1).*abs(x<=L/6); %
h = 1/200;
lambda = 1/10;

opcja = 5;
zakonczprogram = 4;
while opcja ~= zakonczprogram
    opcja = menu('Podaj wartoœæ oporu powietrza',...
        'nu = 2',...
        'nu = 5',...
        'nu = 15',...
        'koniec');
    if opcja == 1
        nu = 2;
        [X, T, u, M, N] = HangingChain(g, nu, v0, v1, [0,t0], [0,L], lambda, h);
        [uotrzymane, X, T, M, N] = HangingChainBessel(g, nu, v0, v1, t0, L, lambda, h);
        figure(1)
        for n=1:10:N
            plot(uotrzymane(:,n), X,'b-', u(:,n), X,'r-') 
            xlabel('u(x,t)')
            ylabel('x')
            pause(0.01)
        end
    elseif opcja == 2
        %nu = 5
        nu = 5;
        [X, T, u, M, N] = HangingChain(g, nu, v0, v1, [0,t0], [0,L], lambda, h);
        [uotrzymane, X, T, M, N] = HangingChainBessel(g, nu, v0, v1, t0, L, lambda, h);
        figure(1)
        for n=1:10:N
            plot(uotrzymane(:,n), X,'b-', u(:,n), X,'r-') 
            xlabel('u(x,t)')
            ylabel('x')
            pause(0.01)
        end
    elseif opcja == 3
        %nu = 15
        nu = 15;
        [X, T, u, M, N] = HangingChain(g, nu, v0, v1, [0,t0], [0,L], lambda, h);
        [uotrzymane, X, T, M, N] = HangingChainBessel(g, nu, v0, v1, t0, L, lambda, h);
        figure(1)
        for n=1:10:N
            plot(uotrzymane(:,n), X,'b-', u(:,n), X,'r-') 
            xlabel('u(x,t)')
            ylabel('x')
            pause(0.01)
        end
    else
        break;
    end
end

%% Wnioski
%Widzimy, ¿e b³¹d wyraŸnie zale¿y od wspó³czynnika oporu powietrza,