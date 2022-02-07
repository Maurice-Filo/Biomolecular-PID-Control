function [T, X] = SSA(Propensity, Stoichiometry, NetworkParameters, X0, tf, N_MaxRxnEvents)
N_Species = length(X0);
T = zeros(1, N_MaxRxnEvents);
X = [X0, zeros(N_Species, N_MaxRxnEvents)];
RxnCounter = 1;
while T(RxnCounter) < tf
    % Calculate propensity at the current (x(t), t)
    a = Propensity(X(:,RxnCounter), NetworkParameters); 
    % Generate a random time at which the next reaction will occur
    a0 = sum(a);
    r1 = rand(1);
    tau = -log(r1) / a0;  
    % Randomly identify which reaction will occur
    r2 = rand(1);
    j = find((cumsum(a) >= r2*a0), 1);
    % Carry out the reaction and update time
    RxnCounter = RxnCounter + 1;
    T(RxnCounter) = T(RxnCounter - 1) + tau;
    X(:, RxnCounter) = X(:, RxnCounter - 1) + Stoichiometry(:,j);
    if RxnCounter >= N_MaxRxnEvents
        break
    end        
end
T(RxnCounter:end) = [];
X(:, RxnCounter:end) = [];
end

