clear all
% Parameters
IMAX = 301;
C = zeros(IMAX, 2);
X = zeros(IMAX, 1);

% Read input data
fid = fopen('datain.dat', 'r');
C0    = fscanf(fid, '%f', 1);
A0    = fscanf(fid, '%f', 1);
NX    = fscanf(fid, '%d', 1);
KAPPA = fscanf(fid, '%f', 1);
DX    = fscanf(fid, '%f', 1);
TTIME = fscanf(fid, '%f', 1);
fclose(fid);

% Grid generation
for I = 1:NX+1
    X(I) = (I-1) * DX;
end

% Initial condition
for I = 1:NX+1
    C(I,1) = A0;
    C(I,2) = A0;
end
C(1,2) = A0 + C0;

% Time-stepping setup
TIME = 0;
DELT = 0.5 * DX^2 / KAPPA;

% Time loop
while TIME <= TTIME
    TIME = TIME + DELT;
    C = yokaihou(IMAX, C, NX, KAPPA, DELT, DX, A0, C0);
end

% Output
output(C, IMAX, NX);
figure;
plot(linspace(0,NX, IMAX),C(:,2))
grid on
grid minor
ylabel('$C [g/m^3]$','interpreter','latex')
xlabel('$X [m]$','interpreter','latex')
title(['Concentration of contaminant at t=',num2str(TTIME),' seconds.'])

function C_ = yokaihou(IMAX, C, NX, KAPPA, DT, DX, A0, C0)

    % Compute new time level 
    for I = 2:NX
        C(I,2) = C(I,1) + (KAPPA * DT / DX^2) * (C(I+1,1) - 2*C(I,1) + C(I-1,1));
    end

    C(NX+1,2) = C(NX,2);

    % Update the current time level
    for I = 1:NX+1
        C(I,1) = C(I,2);
    end
    C_=C;
end

function output(C, IMAX, NX)
    fileID = fopen('dataout.txt', 'w');
    for I = 1:NX+1
        fprintf(fileID, '%4d %7.2f\n', I, C(I,2));
    end
    fclose(fileID);
end
