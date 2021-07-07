function [] = dispersionAnalysisFunction_v2(g, h_i, M_bay, I_bay, M_sat, dt, theta, phi, k_bay, dx_bay, k_ej, dx_ej, theta_ej, nCols, nRows, r_in, r_out, rps)

%% variable definitions
%{
g       acceleration due to gravity
h_i     initial height of bottom of lunasat bay
M_bay   mass of luna sat bay
I_bay   moment of inertia of lunasat bay
M_sat   mass of individual luna sat
dt      time step
theta   polar angle of launcher
phi     azimuthal angle of launcher
k_bay   spring coefficient of launcher spring
dx_bay  spring displacement of launcher spring
k_ej    spring coefficient of ejector spring
dx_ej   spring displacement of ejector spring
nCols   number of columns of luna sats
nRows   number of rows of lunasats
r_in    radius of lunasat location from center of bay
r_out   radius of lunasat ejection location from center of bay
rps     revolutions per second of launcher
%}

%% initial calculations
g       = [0, 0, -g];   %turn g to vector
w       = rps*2*pi;     %bay radial velocity
M_total = M_bay + nCols*nRows*M_sat;
PE      = .5*k_bay*dx_bay^2;    %spring PE
V_mag_i = sqrt(PE*2/M_total);   %bay initial velocity
N_sats  = nCols*nRows;

%% launch setup
if nCols == 2
    ejectionSites = [0 1 0; 0 -1 0]*r_out;
    ejectionNormals = [1 0 0; -1 0 0];
elseif nCols == 4
    ejectionSites = [0 1 0; 0 -1 0; 1 0 0; -1 0 0]*r_out;
    ejectionNormals = [1 0 0; -1 0 0; 0 -1 0; 0 1 0];
end

%rotate ejection sites
for i = 1:nCols
    temp = rotation(ejectionSites(i,:)',[0,1,0]',phi);
    ejectionSites(i,:) = rotation(temp,[0,0,1]',theta)';
    temp = rotation(ejectionNormals(i,:)',[0,1,0]',phi);
    ejectionNormals(i,:) = rotation(temp,[0,0,1]',theta)';
end


lunasatTrajectoryStart  = zeros(nCols*nRows, 3);
lunasatInitialVels      = zeros(nCols*nRows, 3);
lunasatLandingSites     = zeros(nCols*nRows, 3);
    
X0  = [0, 0, h_i];  %initial positioning of lunasat launcher
V0  = V_mag_i*[sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)]; %orientation of bay ejection
[ejTimes, ejVels]   = ejectionTimes(k_ej, dx_ej, theta_ej, .0017, M_bay, M_sat, nRows, nCols);
ejVels              = ejVels;

% plot to ensure everything is ejecting in the right direction
% plot3([0,V0(1)/norm(V0)],[0,V0(2)/norm(V0)],[0,V0(3)/norm(V0)]); hold on
% for i=1:nCols
%     plot3([0,ejectionSites(i,1)],[0,ejectionSites(i,2)],[0,ejectionSites(i,3)],'k','linewidth',2)
%     plot3([ejectionSites(i,1),ejectionSites(i,1)+ejectionNormals(i,1)],...
%         [ejectionSites(i,2),ejectionSites(i,2)+ejectionNormals(i,2)],...
%         [ejectionSites(i,3),ejectionSites(i,3)+ejectionNormals(i,3)],'r','linewidth',2)
% end

X   = X0; %position vector
V   = V0; %velocity vector
W   = [w];

i = 1;
n = 0;
k = 0;
bayAngle = 0;
while X(i,3)>0
    X = [X; X(i,:)+V(i,:)*dt];
    V = [V; V(i,:)+g*dt];
    bayAngle = bayAngle+w*i*dt;
    %account for impulse
    if dt*i > ejTimes(n+1) && N_sats > 0
        for j = 1:nCols
            lunasatTrajectoryStart(nCols*k+j,:) = X(i,:)+rotation(ejectionSites(j,:)',V0'/norm(V0),bayAngle)';
            temp = rotation(ejectionNormals(j,:)',V0'/norm(V0),bayAngle)';
            lunasatInitialVels(nCols*k+j,:) = temp/norm(temp)*w*r_out+V(i,:);
        end
        M_new = M_total-nCols*M_sat;
        w       = (I_bay*w+N_sats*M_sat*r_in^2*w-4*M_sat*r_out^2*w)/(I_bay+(N_sats-nCols)*M_sat*r_in^2);
        N_sats  = N_sats-nCols;
        M_total = M_new;
        V(i+1,:) = V(i+1,:) + calcImpulse(M_total, nCols*M_sat, V0, ejVels(n+1), theta_ej);
        
        k = k+1;
        %ejectionSites = 
        if N_sats>0
            n = n+1;
        elseif N_sats==0
            lastEjection = i;
        end
    end
    W(i) = w;
    i = i+1;
end

for i=1:nCols*nRows
    lunaX = ballisticTrajectory(lunasatTrajectoryStart(i,:),lunasatInitialVels(i,:),g,dt,0);
    lunasatLandingSites(i,:) = lunaX(end,:);
end
close all;
figure
scatter(lunasatLandingSites(:,1),lunasatLandingSites(:,2))

figure
ballisticX = ballisticTrajectory(X0, V0, g, dt, 0);
plot(ballisticX(:,1),ballisticX(:,3),'b','linewidth',2); hold on
plot(X(:,1),X(:,3),'k--','linewidth',2)
plot(X(lastEjection,1),X(lastEjection,3),'ro','linewidth',2)
title('LunaSat Bay Trajectories with LunaSat Mass of '+string(M_sat*1000)+'g')
legend('Ballistic Trajectory','Trajectory Accounting for Impulse','Location of last LunaSat Ejection','location','south')


end

%% functions
function velChange = calcImpulse(M_bay,M_sat, V0, ejVel, theta_ej)
    velChange = M_sat/M_bay*ejVel*sind(theta_ej)*V0/norm(V0);
end

function X = ballisticTrajectory(X0, V0, g, dt, makePlot)
    X = X0;
    V = V0;
    i = 1;
    while X(i,3)>0
        X = [X; X(i,:)+V(i,:)*dt];
        V = [V; V(i,:)+g*dt];
        i = i+1;
    end
    if makePlot==1
        plot(X(:,1),X(:,3));
    end
end

function newPos = rotation(pos, axis, theta)
    newPos = (cos(theta)*eye(3) + (1-cos(theta))*(axis*axis')+sin(theta)*makeSkew(axis))*pos;
end

function M = makeSkew(u)
    M = [0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
end
