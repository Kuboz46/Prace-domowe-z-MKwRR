L = 10;
x =[0, L];

D = @(x) 2*ones(size(x));
fL = @(t) 0;
fR = @(t) 0;

mu = 1/4; h = 1/20;


%%
% Zbadamy rozwiazania dla wszystkich mozliwych kombinacji typow warunkow brzegowych.

% 1. Z lewej: Dirichlet, prawej: Dirichlet
t = [0, 5];
u0 = @(x) sin(x*pi/L);
[u,X,T,M,N] = DiffusionEquation(D,u0,fL,fR,{'Dirichlet','Dirichlet'},t,x,mu,h);
for n = 1:50:N
    plot( [x(1) - h/2, X, x(2) + h/2], ...
        u(:, n), 'g')
    title('Z lewej: Dirichlet, prawej: Dirichlet');
    xlim([0,L]);
    ylim([0,1]);
    pause(0.1)
end

% 2. Z lewej: Dirichlet, prawej: Neumann
t = [0, 5];
u0 = @(x) sin(x*pi/L - pi/2) + 1;
[u,X,T,M,N] = DiffusionEquation(D,u0,fL,fR,{'Dirichlet','Neumann'},t,x,mu,h);
for n = 1:50:N
    plot( [x(1) - h/2, X, x(2) + h/2], ...
        u(:, n), 'g')
    title('Z lewej: Dirichlet, prawej: Neumann');
    xlim([0,L]);
    ylim([0,1]);
    pause(0.1)
end

% 3. Z lewej: Neumann, prawej: Dirichlet
t = [0, 5];
u0 = @(x) sin(3*pi*x/(4*L))-sin(3*pi/4);
[u,X,T,M,N] = DiffusionEquation(D,u0,fL,fR,{'Neumann','Dirichlet'},t,x,mu,h);
for n = 1:50:N
    plot( [x(1) - h/2, X, x(2) + h/2], ...
        u(:, n), 'g')
    title('Z lewej: Neumann, prawej: Dirichlet');
    xlim([0,L]);
    ylim([0,1]);
    pause(0.1)
end

% 4. Z lewej: Neumann, prawej: Neumann
t = [0, 5];
u0 = @(x) sin(3*pi*x/(4*L));
[u,X,T,M,N] = DiffusionEquation(D,u0,fL,fR,{'Neumann','Neumann'},t,x,mu,h);
for n = 1:50:N
    plot( [x(1) - h/2, X, x(2) + h/2], ...
        u(:, n), 'g')
    title('Z lewej: Neumann, prawej: Neumann');
    xlim([0,L]);
    ylim([0,1]);
    pause(0.1)
end

%% Interpretacja fizyczna warunkow brzegowych
% Rownaniem dyfuzji nazywamy rowniez rownaniem przewodnictwa ciepla.
% 1. War. Dirichleta.
% Warunek Dirichleta na lewym brzegu w punkcie x0 okresla poczatkowa temperature w punkcie x0.
% Warunek Dirichleta na prawym brzegu w punkcie x0 okresla koncowa temperature w punkcie x0.
% 2. War. Neumanna
% Warunek Neumanna na lewym brzegu w punkcie x0 mowi o poczatkowej predkosci
% rozchodzenia sie ciepla w punkcie x0.
% Warunek Neumanna na prawym brzegu w punkcie x0 mowi o koncowej predkosci rozchodzenia
% sie ciepla w punkcie x0.