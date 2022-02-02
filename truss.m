%clear; clc;
clear; clc;
% input data.
NC = input('Please, Enter the Coordinates of Each Node as a Matrix:\n');
NR = input('Please, Enter the restraints of the nodes as a matrix:\n');
EC = input('Please, Enter the End Connections of the Elements as a Matrix:\n');
p = input('Please, Enter the Loads on Each Node as a Matrix in order:\n');
E = input('Please, Enter " E " as a Matrix:\n');
A = input('Please, Enter " A " as a Matrix:\n');
fprintf('\n\n\n');
% ------------------Specify the degree of freedom D.O.F--------------------
NN = size(NC, 1); % Number of Nodes.
M_R = zeros(NN, 2);
c = 0;
for i=1: NN
   for j=1: 2
      if NR(i,j) == 0
          M_R(i,j) = 0;
      elseif NR(i,j) == 1
          c = c+1;
          M_R(i,j) = c;
      end
   end
end
DOF = max(max(M_R));
% -----------calculate the (L, l, m) Matricies for each element------------
NE = size(EC, 1); % Number of Elements.
L = zeros(NE, 1); % Length of each Element.
l = zeros(NE, 1); % l of each element l=((X2-X1)/L).
m = zeros(NE, 1); % m of each element m=((Y2-Y1)/L).
for i= 1: NE
    for e= 1: NN
        if EC(i,1) == e
            X1 = NC(e, 1);
            Y1 = NC(e, 2);
        elseif EC(i, 2) == e
            X2 = NC(e, 1);
            Y2 = NC(e, 2);
        end
    end
    L(i,1) = sqrt((X2-X1)^2+(Y2-Y1)^2);
    l(i,1) = (X2-X1)/L(i,1);
    m(i,1) = (Y2-Y1)/L(i,1);
end
% ------------------------Calculate the KG Matrix--------------------------
KE = zeros(4, 4); % K of element.
KG = zeros(NE*4); % K matrix of all Elements.
C = 1;
for e = 1: NE
    KE(1, 1) = (E(e,1)*A(e,1)/L(e,1))*(l(e,1)^2);
    KE(1, 2) = (E(e,1)*A(e,1)/L(e,1))*(l(e,1)*m(e,1));
    KE(1, 3) = (E(e,1)*A(e,1)/L(e,1))*(-(l(e,1)^2));
    KE(1, 4) = (E(e,1)*A(e,1)/L(e,1))*(-(l(e,1)*m(e,1)));
    KE(2, 1) = KE(1, 2);
    KE(2, 2) = (E(e,1)*A(e,1)/L(e,1))*(m(e,1)^2);
    KE(2, 3) = KE(1, 4);
    KE(2, 4) = (E(e,1)*A(e,1)/L(e,1))*(-(m(e,1)^2));
    KE(3, 1) = KE(1, 3);
    KE(3, 2) = KE(2, 3);
    KE(3, 3) = KE(1, 1);
    KE(3, 4) = KE(1, 2);
    KE(4, 1) = KE(1, 4);
    KE(4, 2) = KE(2, 4);
    KE(4, 3) = KE(3, 4);
    KE(4, 4) = KE(2, 2);
    for i = 0:3
        for j = 0:3
            KG(C+i, C+j) = KE(i+1, j+1);
            
        end
    end
    C = C+4;
end
% -----------------Create the A Matrix and Calculate the KM----------------
AM = zeros(NE*4, NN*2); % A Matrix.
for i= 1: NE
    for e= 1: NN
        if EC(i,1) == e
            AM(i*4-3, e*2-1) = 1;
            AM(i*4-2, e*2) = 1;
        elseif EC(i, 2) == e
            AM(i*4-1, e*2-1) = 1;
            AM(i*4, e*2) = 1;
        end
    end
end
KM = AM' * KG * AM;
% ---------Create the (P , K) Matricies and Calculate the U Matrix---------
P = zeros(DOF, 1); % Matrix of forces without reaction places.
K = zeros(DOF, DOF); % K Matrix without reaction places.
counter = 0;
for i = 1:NN
    for j = 1:2
        if M_R(i, j) == 0
           counter = counter +1;
        end
        if M_R(i ,j) ~= 0
            P(M_R(i,j), 1) = p(i, j);
            for dc = 1: DOF
                inc = 0;
                plus = 0;
                for x = 1: NN
                    for y = 1: 2
                        inc = inc + 1;
                        if M_R(x, y) == dc
                            plus = inc;
                        end
                    end
                end
                if j == 1
                    K(M_R(i,j), dc) = KM(i*j+(i-1), plus);
                elseif j == 2
                    K(M_R(i,j), dc) = KM(i*j, plus);
                end
            end
        end
    end
end
U = K\P;
% ---------------------Dicplacement in Local Axes--------------------------
KL = zeros(NE*2, 4); % KL = [l m 0 0; 0 0 l m] for all members at once.
UL = zeros(4, NE); % UL = [up vp; uq vq] for all members at once.
UP = zeros(NE, 1); % UP = up telda for all members at once.
UQ = zeros(NE, 1); % UQ = uq telda for all members at once.
% Create the local stiffness matrix of all elements in one matrix.
for i = 1: NE
    KL(i*2-1, 1) = l(i, 1);
    KL(i*2-1, 2) = m(i, 1);
    KL(i*2, 3) = l(i, 1);
    KL(i*2, 4) = m(i, 1);
end
% Create the local displacement matrix of all elements in one matrix.
M_VC = zeros(NE, 4); % VC Table. 
for i=1: 1: NE
    for j=1: 1: 2
        for e=1: 1: NN
            if EC(i,j) == e
                if j == 1
                    M_VC(i, 1) = M_R(e, 1);
                    M_VC(i, 2) = M_R(e, 2);
                else
                    M_VC(i, 3) = M_R(e, 1);
                    M_VC(i, 4) = M_R(e, 2);
                end
            end
        end
    end
end
MVC = M_VC';
for i = 1: NN
   for j = 1: NE
       if MVC(i, j) ~= 0
          UL(i, j) = U(MVC(i, j), 1);
       end
   end
end
% Value of multiplying the (UL, KL) matricies.
EL = zeros(NE*2, 1); % Element displacements {up telda, uq telda}.
for i = 1: NE
    for j = 1: 2
        if j == 1
            EL(i*2-1, 1) = KL(i*2-1, :) * UL(:, i);
        elseif j == 2
            EL(i*j, 1) = KL(i*2, :) * UL(:, i);
        end
    end
end
% Exporting the Forces in Each Member.
EF = zeros(NE, 1); % Forces in each member.
for i = 1: NE
   for j = 1: 2
      if j == 1
         up = EL(i*2-1, 1);
      elseif j == 2
          uq = EL(i*2, 1);
      end
   end
   EF(i, 1) = ((uq - up) / L(i, 1))*(E(i, 1)*A(i, 1));
end
% Printing the values and types of forces in each member.
F = zeros(NE, 1);
X = 1;
for i = 1: NE
    fprintf('---------------------------------------------------------\n');
    if EF(i, 1) > 1
        fprintf('Member %i = %10.2f  Tension Member\n', X, EF(i,1));
    elseif EF(i, 1) == 0
        fprintf('Member %i = %10.2f  ZERO Member\n', X, EF(i,1));
    elseif EF(i, 1) < 0
        fprintf('Member %i = %10.2f  Compression Member\n', X, EF(i,1));
    end
    X = X + 1;
end
fprintf('---------------------------------------------------------\n');
% Student: SAHER ASHRAF RAMADAN MOHAREB.
% THANKS.
M_R;
DOF;
L;
l;
m;
KE;
KG;
AM;
KM;
P;
K;
U;
KL;
UL;
M_VC;
MVC;
EL;
EF;
