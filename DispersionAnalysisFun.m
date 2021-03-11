function [bayTraj, lunaTrajs, out3, out4, lunaImpactTs, lunaImpactPos] = DispersionAnalysisFun(g_,h_,satM_,bayM_, phi_, theta_, K_,DX_,k_,dx_,nCols_,nRows_,r_,RPS_)
    G = [0,0,-g_];
    h = h_;
    
    satM = satM_;
    bayM = bayM_;
    
    phi = phi_*pi/180;
    theta = theta_;
    
    K = K_;
    DX = DX_;
    
    k = k_;
    dx = dx_;
    
    nCols = nCols_;
    nRows = nRows_;
    
    r = r_;
    RPS = RPS_;
    
    liv = (RPS)*2*pi*r;     % Lunasat Linear velocity magnitude
    
    % Calculations
    
    rowM = nCols*satM;        % mass of a row of lunasats. 
    satMTotal = nRows*rowM;     % mass of all lunasats
    M = bayM + satMTotal;       % total mass of lunasat bay (with all lunasats)
    N = nCols * nRows;        % total number of lunasats
    springPE = (1/2)*K*DX^2;            
    biv = sqrt(springPE * 2 / M);   % bay initial velocity calculation
    
    rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
    roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
    rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
    
    % Assume launch position on z axis
    Xb0 = [0,0,h];
    
    %Calculate Bay Initial Velocity Vector
    Vb0 = ([0,0,1]*biv)*roty(phi); % Bay initial velocity vector
    
    ejectRelPos = [0, 1, 0;  1,  0, 0;
                0, -1, 0; -1, 0, 0]*r;
        
    ejectRelVel = [-1, 0, 0; 0, 1, 0;
                1, 0, 0; 0, -1, 0]*liv;
            
        
    syms Xb(t) Ab(t) eject1(t) eject2(t) eject3(t) eject4(t) Theta(t) XL1(et,t) XL2(et,t) XL3(et,t) XL4(et,t)
    
    Theta(t) = RPS * t * 2*pi;
   
    Ab(t) = G;
    
    Xb(t) = Xb0 + Vb0 * t + (1/2)*Ab(t)*t^2;
    
    bayTraj = Xb(t);
    
    % Get Lunasat Ejection position
    
    ejectPos1(t) = (ejectRelPos(1,:)*rotz(Theta(t)))*roty(phi);
    ejectPos2(t) = (ejectRelPos(2,:)*rotz(Theta(t)))*roty(phi);
    ejectPos3(t) = (ejectRelPos(3,:)*rotz(Theta(t)))*roty(phi);
    ejectPos4(t) = (ejectRelPos(4,:)*rotz(Theta(t)))*roty(phi);
    
    ejectVel1(t) = (ejectRelVel(1,:)*rotz(Theta(t)))*roty(phi);
    ejectVel2(t) = (ejectRelVel(2,:)*rotz(Theta(t)))*roty(phi);
    ejectVel3(t) = (ejectRelVel(3,:)*rotz(Theta(t)))*roty(phi);
    ejectVel4(t) = (ejectRelVel(4,:)*rotz(Theta(t)))*roty(phi);

    XL1(et,t) = (ejectPos1(et) + Xb(et)) + (ejectVel1(et) + Vb0) * t + (1/2)*G*t^2;
    XL2(et,t) = (ejectPos2(et) + Xb(et)) + (ejectVel2(et) + Vb0) * t + (1/2)*G*t^2;
    XL3(et,t) = (ejectPos3(et) + Xb(et)) + (ejectVel3(et) + Vb0) * t + (1/2)*G*t^2;
    XL4(et,t) = (ejectPos4(et) + Xb(et)) + (ejectVel4(et) + Vb0) * t + (1/2)*G*t^2;
    
    disp([XL1(0,0),XL2(0,0),XL3(0,0),XL4(0,0)])
    

    syms lunaTraj(et,t) lunaImpactPos;

    XL = [XL1(et,t); XL2(et,t); XL3(et,t); XL4(et,t)];
    
    lunaTraj = [XL1(et,t); XL2(et,t); XL3(et,t); XL4(et,t)];
    lunaTrajs = [XL1(0,t); XL2(0,t); XL3(0,t); XL4(0,t)];
    
    
    lunaImpactTs = zeros(1,4);
        
  
    disp("Test")
 
   %lunaImpactPos = [XL1(0,0);XL2(0,0);XL3(0,0);XL4(0,0)];
   lunaImpactPos = [];
   
   syms traj(et,t) trajT(t) 
   disp('ttt')
   n=1;
   lunaImpactPos = zeros(500,3);
   
   for i = 1:4
       traj(et,t) = lunaTraj(i,:);
       for j = 1:125
           trajT(t) = traj((j-1)*0.003,t);

           trajTZ = trajT(t);

           eqn = trajTZ(3) == 0;
           impact = double(solve(eqn,t));
           impactT = impact(2);
           impactPOS = double(trajT(impactT));

           lunaImpactPos(n,:) = impactPOS(:);
           %disp(impactT)
           disp(n);
           n = n + 1;
       end
   end
       
    % Output ejection port locations at t=0 to check configuration
    out4 = [double(ejectPos1(0));double(ejectPos2(0));double(ejectPos3(0));double(ejectPos4(0))];
     
    % Solve for and extract bay impact time
    bayImpact = bayTraj(3)==0;
    bayImpactT = solve(bayImpact,t,'Real',true);
    bayImpactT = double(bayImpactT);
    bayImpactT = bayImpactT(2);
    
    out3 = double(Xb(bayImpactT));
 
end

