function [xhat] = BiGAMP( opt, PrioriIn, optIn )

%% Setup why define these setting

nitMin  = opt.nitMin;           % minimum number of iterations
step    = opt.step;             % step size
stepMin = opt.stepMin;          % minimum step size
stepMax = opt.stepMax;          % maximum step size
stepFilter = opt.stepFilter;    % step filter setting, <1 for no effect
adaptStep = opt.adaptStep;      % adaptive step size
stepIncr = opt.stepIncr;        % step inc on succesful step
stepDecr = opt.stepDecr;        % step dec on failed step
stepWindow = opt.stepWindow;    % step size check window size
tol = opt.tol;                  % Convergence tolerance
stepTol = opt.stepTol;          % minimum allowed step size
maxBadSteps = opt.maxBadSteps;  % maximum number of allowed bad steps
maxStepDecr = opt.maxStepDecr;  % amount to decrease maxStep after failures
zvarToPvarMax = opt.zvarToPvarMax;  % maximum zvar/pvar ratio



%Get problem dimensions
M = optIn.M;
U = 1 ;
L = optIn.L;
nit = optIn.nit;
% test
Y = PrioriIn.Y;
Sam = PrioriIn.Sam;
QAM = PrioriIn.QAM;
nuz = PrioriIn.ZVar;
nuw = PrioriIn.noiseVar ;
%Assign Avar and xvar mins
xvarMin = opt.xvarMin;
Avarmin = opt.AvarMin;
state = [];

%% Initialization

B = randi([0,1], U * sqrt(QAM), L);
B(:,1) = 1;
xhat = Constell_Modulate(B, QAM);
xvar = ones(U, L);
Ahat = ones(M, U);
Avar = nuz * ones(M, U);

% xhat = zeros(U, L);
% xvar = ones(U, L);
% Ahat = nuz * randn(M, U);
% Avar = nuz * ones(M, U);

xhatBar = 0;
AhatBar = 0;
shat = 0;
svar = 0;
pvarOpt = 0;
zvarOpt = 0 ;

pvarMin = opt.pvarMin;

%Placeholder initializations

valOpt = [];
val = zeros(nit.outer,1);
zhatOpt = 0;


%% Iterations

%Start timing first iteration


%Control variable to end the iterations
stop = false;
it = 0;
failCount = 0;

step1 = 1;
valIn = -inf;    
% Main iteration loop
while ~stop
    
    % Iteration count
    it = it + 1;
    
    % Check for final iteration
    if it >= nit.outer
        stop = true;
    end
    
    Ahat2 = abs(Ahat).^2;
    xhat2 = abs(xhat).^2;
    zvar = Avar*xhat2 + Ahat2*xvar;
    pvar = zvar + Avar*xvar;
       
    %Include pvar step   ?
    pvar = step1*pvar + (1-step1)*pvarOpt;
    zvar = step1*zvar + (1-step1)*zvarOpt;
    
    %Update zhat
    zhat = Ahat * xhat; 
    % Continued output step
    phat = zhat - shat.* zvar;
    
    % Compute log likelihood at the output and add it the total negative   
    % K-L distance at the input.
    wvar_pos = max(1e-20, nuw);
    predErr = - ((abs(Y - phat).^2  ) + pvar./wvar_pos + log(wvar_pos));
    
    valOut = sum(sum( predErr ));   
    val(it) = valOut + valIn;
    
    % Determine if candidate passed
    if ~isempty(valOpt)
        
        %Check against worst value in last stepWindow good steps
        stopInd = length(valOpt);
        startInd = max(1,stopInd - stepWindow);
        
        %Check the step
        pass = (val(it) > min(valOpt(startInd:stopInd))) ||...
            ~adaptStep || (step <= stepMin);
        
    else
        pass = true;
    end

    % If pass, set the optimal values and compute a new target shat and
    % snew.
    if (pass)
        
        %Slightly inrease step size after pass if using adaptive steps
         step = stepIncr*step;

        
        % Set new optimal values
        shatOpt = shat;
        svarOpt = svar;
        xhatBarOpt = xhatBar;
        xhatOpt = xhat;
        AhatBarOpt = AhatBar;
        AhatOpt = Ahat;
        pvarOpt = pvar;
        zvarOpt = zvar;
        
        %Bound pvar
        pvar = max(pvar, pvarMin);
        
        %We keep a record of only the succesful step valOpt values
        valOpt = [valOpt val(it)]; %#ok<AGROW>
      

        % Compute posterior mean and variance
        %wvar = obj.wvar;
        gain = pvar./(pvar + nuw);
        zhat0 = gain.*(Y - phat) + phat;
        zvar0 = nuw.*gain;
 
        %Compute 1/pvar
        pvarInv =  1./pvar;
        
        shatNew = pvarInv.*(zhat0-phat);
        svarNew = pvarInv.*(1-min(zvar0./pvar,zvarToPvarMax));
        
        %Enforce step size bounds
        step = min([max([step stepMin]) stepMax]);
        
    else        
        %Check on failure count
        failCount = failCount + 1;
        if failCount > maxBadSteps
            failCount = 0;
            stepMax = max(stepMin,maxStepDecr*stepMax);
        end
        % Decrease step size
        step = max(stepMin, stepDecr*step);
        
        %Check for minimum step size
        if step < stepTol
            stop = true;
        end
    end
         
    % Check for convergence if step was succesful
    if pass
        if any(isnan(zhat(:))) || any(isinf(zhat(:)))
            stop = true;
        else
            %testVal = norm(zhat(:) - zhatOpt(:)) / norm(zhat(:));
            testVal = norm(zhat(:) - zhatOpt(:)) / norm(zhat(:));
            if (it > 1) && ...
                    (testVal < tol)
                stop = true;
            end
        end

    end
   
    % Create new candidate shat
    if it > 1 
        step1 = step;
        if stepFilter >= 1
            step1 = step1*it/(it+stepFilter);
        end
    end
    shat = (1-step1)*shatOpt + step1*shatNew;
    svar = (1-step1)*svarOpt + step1*svarNew;
    xhatBar = (1-step1)*xhatBarOpt + step1*xhatOpt;
    AhatBar = (1-step1)*AhatBarOpt + step1*AhatOpt;
    
    %Compute rvar and correct for infinite variance
    rvar = 1./((abs(AhatBar).^2)'*svar);  
    rvar(rvar > opt.varThresh) = opt.varThresh;

    %Update rhat
    rGain = (1 - (rvar.*(Avar'*svar)));
    rGain = min(1,max(0,rGain));
    rhat = xhatBar.*rGain + rvar.*(AhatBar'*shat);   
    rvar = max(rvar, xvarMin);
    
    % Input linear step for A
    qvar = 1./(svar*(abs(xhatBar).^2)');
    qvar(qvar > opt.varThresh) = opt.varThresh;
    
    %Update qhat
    qGain = (1 - (qvar.*(svar*xvar')));
    qGain = min(1,max(0,qGain));
    qhat = AhatBar.*qGain + qvar.*(shat*xhatBar');
    qvar = max(qvar,Avarmin);

    % Input nonlinear step
     [xhat,xvar,valInX] = BiGAMP_Xestim(rhat, rvar, Sam);
     [Ahat,Avar,valInA] = BiGAMP_Aestim(qhat, qvar, nuz);

    
    valIn = sum( valInX(:) ) + sum ( valInA(:) );

    
    %Don't stop before minimum iteration count
    if it < nitMin
        stop = false;
    end
    
end


%% Save the options object that was used

%Estimates of the two matrix factors
xhat = xhatOpt;
% estFin.xvar = xvarOpt;
% estFin.Ahat = AhatOpt;
% estFin.Avar = AvarOpt;



