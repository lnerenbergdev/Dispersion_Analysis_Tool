function [X] = Trajectory(X0,V0,A,t)
    syms X
    X = X0 + V0*t + (1/2)*A*t^2;
end

