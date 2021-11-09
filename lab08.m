c = 1;
L = 1.2;
x = [0, L];
t = [0, 200];

v0 = @(x) zeros(size(x));
u0 = @(x) (sin(2*pi*x/L));
fL = @(t) 0;
fR = @(t) 0;


lambda = 0.8;
boundary = 'Dirichlet';
h = 1/50;

k = lambda*h;
X = (x(1)+h/2):h:(x(2)-h/2);
T = t(1):k:t(2);
M = length(X);
N = length(T);
z=zeros(M,N);
u = z;
for l = 2:11
    u(:,:,l) = z;
end
    
for l = 1:11
    k = (1.5 + 0.1*l);
    F=@(t,x) (x>=0.8 && x<=1).*17*(cos((10*x-9)*pi/2))^2*cos(k*c*pi*t/L);
    [u(:,:,l),X,T,M,N]=WaveEquation1D(c,F,u0,v0,fL,fR,boundary,t,x,lambda,h);
end

for n = 1:N
    %k = 1,5
    subplot(6,2,1);
    plot(X,u(:,n,1),'g')
    title(['k=',num2str(1.5 + 0.1*(1-1)), ' t= 0 + ', num2str(n / N)])
    %k = 1,6
    subplot(6,2,2);
    plot(X,u(:,n,2),'g')
    title(['k=',num2str(1.5 + 0.1*(2-1)), ' t= 0 + ', num2str(n / N)])
    %k = 1,7
    subplot(6,2,3);
    plot(X,u(:,n,3),'g')
    title(['k=',num2str(1.5 + 0.1*(3-1)), ' t= 0 + ', num2str(n / N)])
    %k = 1,8
    subplot(6,2,4);
    plot(X,u(:,n,4),'g')
    title(['k=',num2str(1.5 + 0.1*(4-1)), ' t= 0 + ', num2str(n / N)])
    %k = 1,9
    subplot(6,2,5);
    plot(X,u(:,n,5),'g')
    title(['k=',num2str(1.5 + 0.1*(5-1)), ' t= 0 + ', num2str(n / N)])
    %k = 2
    subplot(6,2,6);
    plot(X,u(:,n,6),'g')
    title(['k=',num2str(1.5 + 0.1*(6-1)), ' t= 0 + ', num2str(n / N)])
    %k = 2,1
    subplot(6,2,7);
    plot(X,u(:,n,7),'g')
    title(['k=',num2str(1.5 + 0.1*(7-1)), ' t= 0 + ', num2str(n / N)])
    %k = 2,2
    subplot(6,2,8);
    plot(X,u(:,n,8),'g')
    title(['k=',num2str(1.5 + 0.1*(8-1)), ' t= 0 + ', num2str(n / N)])
    %k = 2,3
    subplot(6,2,9);
    plot(X,u(:,n,9),'g')
    title(['k=',num2str(1.5 + 0.1*(9-1)), ' t= 0 + ', num2str(n / N)])
    %k = 2,4
    subplot(6,2,10);
    plot(X,u(:,n,10),'g')
    title(['k=',num2str(1.5 + 0.1*(10-1)), ' t= 0 + ', num2str(n / N)])
    %k = 2,5
    subplot(6,2,11);
    plot(X,u(:,n,11),'g')
    title(['k=',num2str(1.5 + 0.1*(11-1)), ' t= 0 + ', num2str(n / N)])
    xlim([0,L])
    ylim([-1.1,1.1])
    pause(0.001)
end

% Tu zaobserwujemy zjawisko rezonansu dla k = 2.
for n = 1:N
    plot(X,u(:,n,6),'g')
    title(['k=',num2str(2), ' t= 0 + ', num2str(n / N)])
    xlim([0,L])
    ylim([-1.1,1.1])
    pause(0.001)
end
%% Wnioski
% Biorac przedzial czasowy [0, 200] widzimy z wykresow rozwiazan rownania niejednorodnego dla k = 1.5, 1.6, 
% ..., 2.5, ze okresowosc rozwiazan dla k = 1.5, 1.6, ..., 2.5 "burzy sie" w krotkich odstepach czasu.