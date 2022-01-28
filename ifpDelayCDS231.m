%% Code to play with IFPs, Instability, and Delays
%% Initialize parameters
n=1; % Number of states (start with scalar)
% Initialize state matrices and useful placeholders
dummyI=eye(n);
% Check delays up to this many steps
maxDelay = 20;
aScalarStability=1.1; % Max eigenvalue of the discrete-time ring
delaySweep = 0:maxDelay;


% Initialize figure 
hold off; % to allow comparison of different n or a, start off 

% Choose one of these:
problemDef = 'State Fdbk And Act Delay';
% problemDef = 'Full Ctrl And Sense Delay';


%% Sweep through different delays
zeroRows = zeros(n,(maxDelay-2)*n);
if isempty(zeroRows); zeroRows = []; end
allCost = zeros(1,length(delaySweep));
for delInd = 1:length(delaySweep)
    delStep = delaySweep(delInd);
    % Build A-hat
    aStack = [];
    for stackInd = 1:maxDelay
        if delStep == 0
        aStack = blkdiag(aStack,dummyI*0);
        else
        aStack = blkdiag(aStack,dummyI);
        end
    end
    ringTmp = dummyI(:,[2:n 1]);
    ringStructure=(dummyI+ringTmp+ringTmp');
    numNeighbors=sum(ringStructure(:,1));
    ringScaled = ringStructure*aScalarStability/numNeighbors;
    
    Ahat = zeros(n*(maxDelay+1));
    Ahat(1:n,1:n) = ringScaled;
    Ahat(1:(end-n),n+1:end) = aStack;
    % Build B or C, depending on case
    bcStack = [];
    for stackInd = 1:maxDelay
        if stackInd == delStep
            isActive = 1;
        else
            isActive = 0;
        end
        bcStack = blkdiag(bcStack,dummyI*isActive);
    end
    if delStep == 0
        BC = eye(size(Ahat,1));
    else
        BC = [zeros(n,n*maxDelay);bcStack];

    end
    % State Costs / State Disturbance (IFP States are Free/Perfect)
    QW=blkdiag(dummyI,bcStack*0); 
    % Control Costs / Estimator Noise (IFP States are Free/Perfect)
    RN=eps*eye(size(BC,2));
    % Solve a Riccati equation to get K
    K = dlqr(Ahat,BC,QW,RN);

    
    % We now have all the matrices we need.
    % If you wanted to "knock down" IFPs by some scaling factor, 
    % you could do that here.
    
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
    
    % The ones of the impulse vector define the states for which
    % disturbances are allowed to enter the system. If only the first n
    % states are allowed, then the plant can be disturbed while the
    % controller is delayed but internally noiseless.
    impvec = zeros(1,size(Ahat,1));
    impvec(1:n) = 1;
    allCost(delInd) = impvec*J*impvec';
end

%% Make summary figure

% figure; % Might uncomment this depending on how you do your tests

plot(0:(maxDelay), allCost,'ko-','LineWidth',2)
hold on
set(gca,'Ycolor','k')


ylabel('State Cost')
xlabel('Delay Steps')
title(problemDef)

set(gca,'FontSize',16)

