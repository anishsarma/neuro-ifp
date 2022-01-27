%% Code to play with IFPs, Instability, and Delays
%% Initialize parameters
n=1; % Number of states (may encounter stability issues for n > 2)
% Initialize state matrices and useful placeholders
dummyI=eye(n);
% Check delays up to this many steps
maxDelay = 20;
aScalarStability=1.1;

% Choose one of these:
problemDef = 'State Fdbk And Act Delay';
% problemDef = 'Full Ctrl And Sense Delay';

%% Sweep through different delays
zeroRows = zeros(n,(maxDelay-2)*n);
if isempty(zeroRows); zeroRows = []; end
allCost = zeros(1,maxDelay);
for delInd = 1:maxDelay
    % Build A-hat
    aStack = [];
    for stackInd = 1:maxDelay
        aStack = blkdiag(aStack,dummyI);
    end
    Ahat = zeros(n*(maxDelay+1));
    Ahat(1:n,1:n) = aScalarStability;
    Ahat(1:(end-n),n+1:end) = aStack;
    % Build B or C, depending on case
    bcStack = [];
    for stackInd = 1:maxDelay
        if stackInd == delInd
            isActive = 1;
        else
            isActive = 0;
        end
        bcStack = blkdiag(bcStack,dummyI*isActive);
    end
    BC = [zeros(n,n*maxDelay);bcStack];
    
    % State Costs / State Disturbance (IFP States are Free/Perfect)
    QW=blkdiag(dummyI,bcStack*0); 
    % Control Costs / Estimator Noise (IFP States are Free/Perfect)
    RN=eps*eye(size(QW)-n);
    % Solve a Riccati equation to get K
    K = dlqr(Ahat,BC,QW,RN);

    switch problemDef
        case 'State Fdbk And Act Delay'
            B = BC;
            C = eye(size(Ahat,1));
            Acl = Ahat-B*K*C;
        case 'Full Ctrl And Sense Delay'
            B= eye(size(Ahat,1));
            Ahat = Ahat';
            C = BC';
            K = K';
            Acl = (Ahat-B*K*C)';
        otherwise
            error
    end

    J = dlyap(Acl,QW);
    
    allCost(delInd) = ones(1,size(Ahat,1))*J*ones(size(Ahat,1),1);
end

%% Make summary figure
figure('Position',[315   279   685   519]);

plot(1:maxDelay, allCost,'ko-','LineWidth',2)
hold on
set(gca,'Ycolor','k')


ylabel('State Cost')
xlabel('Delay Steps')
title(problemDef)

set(gca,'FontSize',16)
axis square

