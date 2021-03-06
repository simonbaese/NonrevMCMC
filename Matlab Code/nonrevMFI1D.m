function nonrevMFI1D(theta,plotRun,calcIndicators,n)
%% Non-reversible Mean Field Ising model

% -------------------------------------------------------------------------
%% Initialise variables

%Check number of input arguments
if nargin < 1   
    theta = 0.01;
    plotRun = 0;
    calcIndicators = 0;
    n = 1000;
end

if nargin < 2
    plotRun = 0;
    calcIndicators = 0;
    n = 1000;
end

if nargin < 3
    calcIndicators = 0;
    n = 1000;
end

if nargin < 4
    n = 1000;
end

%Number of Samples
%n = 5000;

%Samples (#,value,direction)
samples = zeros(n,1,2);

%Number of Spins (even number)
N = 80;

%Temperatur and Inverse Temperatur
T = 1;
beta = 1/T;

%Exchange Energy
J = 1;

%Theta
%theta = 0.01;


%Magnetization
%M = 0;
%M = randsrc*2*(randi(N + 1) - 1); 
M = randsrc*2*(randi(N/4 + 1) - 1);
dir = randsrc;
samples(1,:,1) = M;
samples(1,:,2) = dir;

%Accepted samples in step A
accepted = 0;

%Prepare true distribution
x = -N:2:N;
yy = zeros(1,N+1);

for i = 1:N+1
    m = x(i);
    factor = exp(gammaln(N+1) - (gammaln(-m/2 + N/2 + 1) + gammaln(m/2 + N/2 + 1)));
    yy(i) = factor*exp(beta*J*m^2/(2*N));
end
z = sum(yy);
yy = yy/z;

%Prepare checkerboard
if M == 0
   A = ones(1,N);
   B = randperm(numel(A));
   A(B(1:N/2)) = -1; 
else
   m = N/2 - sign(M)*M/2;
   A = sign(M)*ones(1,N);
   B = randperm(numel(A));
   A(B(1:m)) = -sign(M);
end

%Prepare figure
close all;
if plotRun == 0 || plotRun == 1
figure;
set(gcf, 'Position', get(0,'Screensize'));
end

% -------------------------------------------------------------------------
%% Calculate Samples and Draw Checkerboard

for i = 2:n  
    
    %Step A ------------------------------------------
    %Read last step
    M = samples(i-1,:,1);
    dir = samples(i-1,:,2);
    proposal = M + 2*dir;

    %Acceptance Propability
    factor = (N - dir*M)/(N + dir*M + 2);
    acceptance = factor*exp(2*beta*J*(dir*M + 1)/N);
    
    %Accept Proposal
    if rand < acceptance
       samples(i,:,1) =  proposal;
       samples(i,:,2) = -dir;
       accepted = accepted + 1;
       acc = 1;
       
    %Deny Proposal
    else
       samples(i,:,:) = samples(i-1,:,:);
       acc = 0;
    end %if
    %end step A 

    %Step B ------------------------------------------
    if rand < 1 - theta
        
        %Keep direction
        samples(i,:,2) = -samples(i,:,2);
        
    end
    %end step B 
    
    %Update checkerboard -----------------------------
    if acc
      
    A = A(1,:); k = 1; u = 1; look = randsrc;
    
    %Read move
    move = (proposal - M)/2;
    
    %Draw move
    while k <= 3
        
        %Propose location
        rr = randi(N);
        
        %Correct proposed location
        while A(rr) == move && u < 3 
            if rr == N && look == 1
               look = -1;
               u = u + 1;
            elseif rr == 1 && look == -1
               look = 1;
               u = u + 1;
            else
               rr = rr + look;
            end
        end %correct
        
        %Make touchy move
        if rr == N || rr == 1
            A(rr) = move;
            break;
        elseif A(rr+1) == move || A(rr-1) == move
            A(rr) = move;
            break;
        elseif k == 3
            A(rr) = move;
            break;
        end
        
        k = k + 1;
        
    end %while
    
    A = repmat(A,6,1);
    
    end %if
    %end update checkerboard
    
    %Plot --------------------------------------------
    if plotRun == 1
    
    %Clear current figure
    clf;
    
    %Plot samples
    subplot(2,1,1);
    title('Samples');
    grid off;
    xlim([-N N]);
    xlabel('Magnetization');
    ylabel('Probability');
    hold('on');   
    y = samples(1:i,:,1);
    a = histc(y,x);
    zz = sum(a);
    bar(x,a/zz,'b','EdgeColor',[0 0 0.6])
    plot(x,yy,'g','LineWidth',3);
    legend('Samples','True Distribution');
    
    %Plot checkerboard
    subplot(2,1,2);
    title('Checkerboard');
    set(gca, 'YTickLabel',[])
    hold('on');
    colormap(gray); 
    imagesc(A);
    axis image;
    
    %Pause
    tic; while toc < 0.0000001; end
    drawnow;
 
    end %if
    %end plot
    
end %for
%end sampling

% -------------------------------------------------------------------------
%% Output after Calculation

%Read samples for plots
y = samples(:,:,1);

%Plot sample histogram and movement of samples
if plotRun == 0
    
    %Plot samples
    subplot(2,1,1);
    title('Samples');
    grid off;
    xlim([-N N]);
    xlabel('Magnetization');
    ylabel('Probability');
    hold('on');   
    a = histc(y,x);
    zz = sum(a);
    bar(x,a/zz,'b','EdgeColor',[0 0 0.6]);
    plot(x,yy,'g','LineWidth',3);
    legend('Samples','True Distribution');
    
    %Movement of magnetization
    subplot(2,1,2);
    title('Movement');
    ylabel('Magnetization');
    xlabel('Samples');
    set(gca, 'XTickLabel',[]);
    hold('on');
    y = y(1:end);
    plot(1:length(y),y,'r','LineWidth',1);
    
end %plot

%Display indicators
if calcIndicators
    
   samples = y;
   L = 100;
   m = mean(samples);
   v = cov(samples);
   autocorrs = zeros(L,1);
   
   %Autocorrelation
   for i = 1:L
       autosum = 0;
       for j = 1:n-i
           autosum = autosum + (samples(j,:)-m)/(2*v)*(samples(j+i,:)-m)';
       end
       autocorrs(i) = autosum/(n-i-1);
    end

   INEFFICIENCY = 1+2*sum(autocorrs)
   ACCEPTANCERATE = accepted/n
    
end %indicators

end %main