%================================== FAILURE ANALYSIS OF LAMINATE ===================================
%===================================================================================================
% NAME: SURAJ KIRAN SHAH                  ROLL NO.: 204103334                       M.Tech
%===================================================================================================

clc
clear all

%% Main Loop for Degradation Method (1:Partial, 2:Complete) ****************************************
Plot = figure(1);
for Degradation_Method = 1:2

%% User Input **************************************************************************************
% Reading from Data Sheet
Data = readmatrix('Data_Input_Sheet.xlsx');
Data(:,1:2) = [];
Engg_Constant = Data(:,[1:4]);
Strength_Parameters = Data(:,[5:9]);
Alpha = Data(:,[10:11]);
Beta = Data(:,[12:13]);
Angle = Data(:, 14);
Thickness = Data(:, 15);
n = size(Data, 1);    % Number of plies

% Loading & Hygrothermal Input ---------------------------------------------------------------------
N = [100000; 0; 0; 0; 0; 0];  % in N/m
delta_T = 50;                 % Change in temperature
delta_C = 0;                  % Moisture absorption in Kg/Kg

% Finding Zk matrix --------------------------------------------------------------------------------
h = sum(Thickness);  % Thickness of laminate
Zk = zeros(n+1,1);
Zk(1,1) = -h/2;
temp = -h/2;
for i = 1:n
    temp = temp + Thickness(i);
    Zk(i+1,1) = temp;
end

% Generating Output Text File ----------------------------------------------------------------------
if Degradation_Method == 1
    file = fopen('Partial_Degradation.txt','w');
else
    file = fopen('Complete_Degradation.txt','w');
end

%% Printing basic information in output file *******************************************************
fprintf(file,'**************** FAILURE ANALYSIS OF LAMINATE *******************\n\n');
fprintf(file,'Failure Theory: Maximum Stress Theory\n');

if Degradation_Method == 1
    fprintf(file,'Degradation Method: Partial Degradation\n\n');
else
    fprintf(file,'Degradation Method: Complete Degradation\n\n');
end

fprintf(file,'*****************************************************************');

%% Main Program ************************************************************************************
zero = zeros(3);
for i = 1:n
    SR_Old_mat(:,:,i) = zero;
    SR_Old_abs_mat(:,:,i) = zero;
end
M = 1; % for Complete Degradation to skip Q and Qbar calculation
k = 1; % for while loop
knee = 1; % for printing in output file
PF_Load = 0; % for first iteration setting it to zero
No_of_ply = n;

%%
while k <= No_of_ply
    
% Finding Q of each ply ----------------------------------------------------------------------------
if M == 1
for i = 1:n
    nu_21 = Engg_Constant(i,3)*Engg_Constant(i,2)/Engg_Constant(i,1);
    Q11 = Engg_Constant(i,1)/(1 - Engg_Constant(i,3)*nu_21);
    Q22 = Engg_Constant(i,2)/(1 - Engg_Constant(i,3)*nu_21);
    Q12 = Engg_Constant(i,3)*Engg_Constant(i,2)/(1 - Engg_Constant(i,3)*nu_21);
    Q66 = Engg_Constant(i,4);
    
    Q = [Q11, Q12, 0;
         Q12, Q22, 0;
          0,   0, Q66];

    Q_mat(:,:,i) = Q;
end

% Finding Qbar of each ply -------------------------------------------------------------------------
for i = 1:n
    x = Angle(i);
    R = [ 1 0 0;
          0 1 0;
          0 0 2];
  
    T = [   cosd(x)^2,       sind(x)^2,      2*sind(x)*cosd(x);
            sind(x)^2,       cosd(x)^2,     -2*sind(x)*cosd(x);
        -sind(x)*cosd(x), sind(x)*cosd(x),  cosd(x)^2-sind(x)^2 ];

    Qbar = inv(T)*Q_mat(:,:,i)*R*T*inv(R);
    Q_bar_mat(:,:,i) = Qbar;
end
end

%% Finding ABBD Matrix *****************************************************************************
% Finding Matrix A
A = zeros(3);
for i = 1:n
    temp = Q_bar_mat(:,:,i)*(Zk(i+1) - Zk(i));
    A = A + temp;
end

% Finding Matrix B ---------------------------------------------------------------------------------
B = zeros(3);
for i = 1:n
    temp = 0.5*Q_bar_mat(:,:,i)*(Zk(i+1)^2 - Zk(i)^2);
    B = B + temp;
end

% Finding Matrix D ---------------------------------------------------------------------------------
D = zeros(3);
for i = 1:n
    temp = (1/3)*Q_bar_mat(:,:,i)*(Zk(i+1)^3 - Zk(i)^3);
    D = D + temp;
end

% Clubbing A, B & D matrix -------------------------------------------------------------------------
ABBD = zeros(6);

ABBD(1:3,1:3) = ABBD(1:3,1:3) + A;
ABBD(1:3,4:6) = ABBD(1:3,4:6) + B;
ABBD(4:6,1:3) = ABBD(4:6,1:3) + B;
ABBD(4:6,4:6) = ABBD(4:6,4:6) + D;

%% Finding Mid Surface In-plane Strain & Curvature *************************************************
Mid_Surface = inv(ABBD)*N;
Mid_Surface_Strain = Mid_Surface(1:3, :);
Mid_Surface_Curvature = Mid_Surface(4:6, :);

%% Finding Strain and Stress in each ply ***********************************************************
% Finding Global Axes(xy) Strain in each ply
for i = 1:n
    for j = 1:3
        if j==1
            strain_xy_top = Mid_Surface_Strain + Zk(i).*Mid_Surface_Curvature;
        
        elseif j==2
            strain_xy_middle = Mid_Surface_Strain + 0.5*(Zk(i)+Zk(i+1)).*Mid_Surface_Curvature;
        
        else
            strain_xy_bottom = Mid_Surface_Strain + Zk(i+1).*Mid_Surface_Curvature;
        end
    end
    strain_xy = [strain_xy_top, strain_xy_middle, strain_xy_bottom];
    strain_xy_mat(:,:,i) = strain_xy;
end
        
% Finding Material Axis(12) Strain in each ply -----------------------------------------------------
for i = 1:n
    x = Angle(i);
    R = [ 1 0 0;
          0 1 0;
          0 0 2];
    T = [   cosd(x)^2,       sind(x)^2,      2*sind(x)*cosd(x);
            sind(x)^2,       cosd(x)^2,     -2*sind(x)*cosd(x);
        -sind(x)*cosd(x), sind(x)*cosd(x),  cosd(x)^2-sind(x)^2 ];

    for j = 1:3
        if j==1
            strain_12_top = R*T*inv(R)*strain_xy_mat(:,1,i);
            
        elseif j==2
            strain_12_middle = R*T*inv(R)*strain_xy_mat(:,2,i);
            
        else
            strain_12_bottom = R*T*inv(R)*strain_xy_mat(:,3,i);
        end
    end
    
    strain_12 = [strain_12_top, strain_12_middle, strain_12_bottom];
    strain_12_mat(:,:,i) = strain_12;
end

% Finding Matrial Axes(12) Stress in eacch ply -----------------------------------------------------
for i = 1:n
    for j = 1:3
        if j==1
            stress_12_top = Q_mat(:,:,i)*strain_12_mat(:,1,i);
            
        elseif j==2
            stress_12_middle = Q_mat(:,:,i)*strain_12_mat(:,2,i);
            
        else
            stress_12_bottom = Q_mat(:,:,i)*strain_12_mat(:,3,i);
        end
    end
    
    stress_12 = [stress_12_top, stress_12_middle, stress_12_bottom];
    stress_12_mat(:,:,i) = stress_12;  
end

%% Hygrothermal Stress in Laminate *****************************************************************
% Finding Alpha of each ply in Global Axes(xy)
Alpha = [Alpha(:,:), zeros(n,1)];
for i = 1:n
    x = Angle(i);
    R = [ 1 0 0;
          0 1 0;
          0 0 2];
    T = [   cosd(x)^2,       sind(x)^2,      2*sind(x)*cosd(x); 
            sind(x)^2,       cosd(x)^2,     -2*sind(x)*cosd(x);
        -sind(x)*cosd(x), sind(x)*cosd(x),  cosd(x)^2-sind(x)^2 ];

    alpha_xy = R*inv(T)*inv(R)*(Alpha(i,:))';
     
    alpha_xy_mat(:,:,i) = alpha_xy;
end
Alpha(:,3) = [];

% Finding Equivalent Thermal Load ------------------------------------------------------------------
N_T = zeros(6,1);
for i = 1:n
    thermal_load = delta_T*Q_bar_mat(:,:,i)*alpha_xy_mat(:,:,i)*(Zk(i+1) - Zk(i));
    thermal_moment = 0.5*delta_T*Q_bar_mat(:,:,i)*alpha_xy_mat(:,:,i)*(Zk(i+1)^2 - Zk(i)^2);
    N_T(1:3,:) = N_T(1:3,:) + thermal_load;
    N_T(4:6,:) = N_T(4:6,:) + thermal_moment;
end

% Finding Mid Surface Strain due to Thermal Load ---------------------------------------------------
Mid_Surface_Thermal = inv(ABBD)*N_T;
Mid_Surface_Thermal_Strain = Mid_Surface_Thermal(1:3, :);
Mid_Surface_Thermal_Curvature = Mid_Surface_Thermal(4:6, :);

% Finding Global Axes(xy) Strain in each ply due to Thermal Load -----------------------------------
for i = 1:n
    strain_xy_thermal = Mid_Surface_Thermal_Strain + 0.5*(Zk(i)+Zk(i+1)).*Mid_Surface_Thermal_Curvature;
    strain_xy_thermal_mat(:,:,i) = strain_xy_thermal;
end

% Finding Free Thermal Strain in each ply ----------------------------------------------------------
for i = 1:n
    strain_xy_thermal_free = delta_T*alpha_xy_mat(:,:,i);
    strain_xy_thermal_free_mat(:,:,i) = strain_xy_thermal_free;
end

% Finding Residual Thermal Strain in each ply ------------------------------------------------------
for i = 1:n
    residual_thermal_strain = strain_xy_thermal_mat(:,:,i) - strain_xy_thermal_free_mat(:,:,i);
    residual_thermal_strain_mat(:,:,i) = residual_thermal_strain;
end

% Finding Residual Thermal Stress in each ply in Global Axes(xy) -----------------------------------
for i = 1:n
    residual_thermal_stress_xy = Q_bar_mat(:,:,i)*residual_thermal_strain_mat(:,:,i);
    residual_thermal_stress_xy_mat(:,:,i) = residual_thermal_stress_xy;
end

% Finding Residual Thermal Stress in each ply in Material Axes(12) ---------------------------------
for i = 1:n
    x = Angle(i);

    T = [   cosd(x)^2,       sind(x)^2,      2*sind(x)*cosd(x);
            sind(x)^2,       cosd(x)^2,     -2*sind(x)*cosd(x);
        -sind(x)*cosd(x), sind(x)*cosd(x),  cosd(x)^2-sind(x)^2 ];

    residual_thermal_stress_12 = T*residual_thermal_stress_xy_mat(:,:,i);   
    residual_thermal_stress_12_mat(:,:,i) = residual_thermal_stress_12;
end

%% Hygroscopic Stress ******************************************************************************
% Finding Beta of each ply in Global Axes(xy)
Beta = [Beta(:,:), zeros(n,1)];
for i = 1:n
    x = Angle(i);
    R = [ 1 0 0;
          0 1 0;
          0 0 2];
    T = [   cosd(x)^2,       sind(x)^2,      2*sind(x)*cosd(x);
            sind(x)^2,       cosd(x)^2,     -2*sind(x)*cosd(x);
        -sind(x)*cosd(x), sind(x)*cosd(x),  cosd(x)^2-sind(x)^2 ];

    beta_xy = R*inv(T)*inv(R)*(Beta(i,:))';    
    beta_xy_mat(:,:,i) = beta_xy;
end
Beta(:,3) = [];

% Finding Equivalent Hygroscopic Load --------------------------------------------------------------
N_H = zeros(6,1);
for i = 1:n
    hygroscopic_load = delta_C*Q_bar_mat(:,:,i)*beta_xy_mat(:,:,i)*(Zk(i+1) - Zk(i));
    hygroscopic_moment = 0.5*delta_C*Q_bar_mat(:,:,i)*beta_xy_mat(:,:,i)*(Zk(i+1)^2 - Zk(i)^2);
    N_H(1:3,:) = N_H(1:3,:) + hygroscopic_load;
    N_H(4:6,:) = N_H(4:6,:) + hygroscopic_moment;
end

% Finding Mid Surface Strain due to Hygroscopic Load -----------------------------------------------
Mid_Surface_Hygroscopic = inv(ABBD)*N_H;
Mid_Surface_Hygroscopic_Strain = Mid_Surface_Hygroscopic(1:3, :);
Mid_Surface_Hygroscopic_Curvature = Mid_Surface_Hygroscopic(4:6, :);

% Finding Global Axes(xy) Strain in each ply due to Hygroscopic Load -------------------------------
for i = 1:n
    strain_xy_hygroscopic = Mid_Surface_Hygroscopic_Strain + 0.5*(Zk(i)+Zk(i+1)).*Mid_Surface_Hygroscopic_Curvature;
    strain_xy_hygroscopic_mat(:,:,i) = strain_xy_hygroscopic;
end

% Finding Free Hygroscopic Strain in each ply ------------------------------------------------------
for i = 1:n
    strain_xy_hygroscopic_free = delta_C*beta_xy_mat(:,:,i);
    strain_xy_hygroscopic_free_mat(:,:,i) = strain_xy_hygroscopic_free;
end

% Finding Residual Hygroscopic Strain in each ply --------------------------------------------------
for i = 1:n
    residual_hygroscopic_strain = strain_xy_hygroscopic_mat(:,:,i) - strain_xy_hygroscopic_free_mat(:,:,i);  
    residual_hygroscopic_strain_mat(:,:,i) = residual_hygroscopic_strain;
end

% Finding Residual Hygroscopic Stress in each ply in Global Axes(xy) -------------------------------
for i = 1:n
    residual_hygroscopic_stress_xy = Q_bar_mat(:,:,i)*residual_hygroscopic_strain_mat(:,:,i);  
    residual_hygroscopic_stress_xy_mat(:,:,i) = residual_hygroscopic_stress_xy;
end

% Finding Residual Hygroscopic Stress in each ply in Material Axes(12) -----------------------------
for i = 1:n
    x = Angle(i);

    T = [   cosd(x)^2,       sind(x)^2,      2*sind(x)*cosd(x);
            sind(x)^2,       cosd(x)^2,     -2*sind(x)*cosd(x);
        -sind(x)*cosd(x), sind(x)*cosd(x),  cosd(x)^2-sind(x)^2 ];

    residual_hygroscopic_stress_12 = T*residual_hygroscopic_stress_xy_mat(:,:,i);    
    residual_hygroscopic_stress_12_mat(:,:,i) = residual_hygroscopic_stress_12;
end

% Finding Residual Hygrothermal Stress in each ply in Material Axes(12) ----------------------------
for i = 1:n
    residual_stress_12 = residual_thermal_stress_12_mat(:,:,i) + residual_hygroscopic_stress_12_mat(:,:,i);
    residual_stress_12_mat(:,:,i) = residual_stress_12;
end

%% Finding Ply Failure Load considering Hygrothermal-Mechanical Loading ****************************
% Finding New Strength Parameters with Residual Stress
New_Strength_Parameters = zeros(8,5);
for i = 1:n
    for j = 1:5
        if j==1 || j==2
        New_Strength_Parameters(i,j) = Strength_Parameters(i,j) - residual_stress_12_mat(1,1,i);
        elseif j==3 || j==4
        New_Strength_Parameters(i,j) = Strength_Parameters(i,j) - residual_stress_12_mat(2,1,i);
        else
        New_Strength_Parameters(i,j) = Strength_Parameters(i,j) - residual_stress_12_mat(3,1,i);
        end
    end
end

% Finding Strength Ratio of each ply ---------------------------------------------------------------
% Strength Ratio in order Longitudnal, Transverse & Shear
for i = 1:n
    for j = 1:3
        if stress_12_mat(1,j,i)>0
            SRL = stress_12_mat(1,j,i)/New_Strength_Parameters(i,1);
        else
            SRL = stress_12_mat(1,j,i)/New_Strength_Parameters(i,2);
        end
        
        if stress_12_mat(2,j,i)>0
            SRT = stress_12_mat(2,j,i)/New_Strength_Parameters(i,3);
        else
            SRT = stress_12_mat(2,j,i)/New_Strength_Parameters(i,4);
        end
        
        SRS = abs(stress_12_mat(3,j,i)/New_Strength_Parameters(i,5));
        
        if j==1
            SR_top = [SRL; SRT; SRS];
            SR_top_abs = [abs(SRL); abs(SRT); SRS];
            
        elseif j==2
            SR_middle = [SRL; SRT; SRS];
            SR_middle_abs = [abs(SRL); abs(SRT); SRS];
            
        else
            SR_bottom = [SRL; SRT; SRS];
            SR_bottom_abs = [abs(SRL); abs(SRT); SRS];
        end
    end
    
    SR = [SR_top, SR_middle, SR_bottom];
    SR_abs = [SR_top_abs, SR_middle_abs, SR_bottom_abs];
    SR_mat(:,:,i) = SR;
    SR_abs_mat(:,:,i) = SR_abs;
end

% Finding New Strength Ratio by considering Total Stress -------------------------------------------
for i = 1:n
    SR_New_mat(:,:,i) = zero;
    SR_New_abs_mat(:,:,i) = zero;
end
 
SR_New_mat = SR_mat + SR_Old_mat;
SR_New_abs_mat = SR_abs_mat + SR_Old_abs_mat;

% Finding Ply Failure Load -------------------------------------------------------------------------
PF_Load = (N(1,1) + PF_Load)/max(SR_New_abs_mat(:));  % Previous PF_Load has to be added as we added previous SR value
Stress(knee) = PF_Load/h;

%% Finding Actual Stress and Strain values when Ply Failure Load is applied ************************
% Defining Ply Failure Load
N_PF = [PF_Load; 0; 0; 0; 0; 0];

% Finding Mid Surface In-plane Strain & Curvature --------------------------------------------------
Mid_Surface = inv(ABBD)*N_PF;
Mid_Surface_Strain = Mid_Surface(1:3, :);
Mid_Surface_Curvature = Mid_Surface(4:6, :);
Strain(knee) = Mid_Surface_Strain(1,1);

% Finding Global Axes(xy) Strain in each ply -------------------------------------------------------
for i = 1:n
    for j = 1:3
        if j==1
            strain_xy_top = Mid_Surface_Strain + Zk(i).*Mid_Surface_Curvature;
        
        elseif j==2
            strain_xy_middle = Mid_Surface_Strain + 0.5*(Zk(i)+Zk(i+1)).*Mid_Surface_Curvature;
        
        else
            strain_xy_bottom = Mid_Surface_Strain + Zk(i+1).*Mid_Surface_Curvature;
        end
    end
    strain_xy = [strain_xy_top, strain_xy_middle, strain_xy_bottom];
    strain_xy_mat(:,:,i) = strain_xy;
end
        
% Finding Material Axis(12) Strain in each ply -----------------------------------------------------
for i = 1:n
    x = Angle(i);
    R = [ 1 0 0;
          0 1 0;
          0 0 2];
    T = [   cosd(x)^2,       sind(x)^2,      2*sind(x)*cosd(x);
            sind(x)^2,       cosd(x)^2,     -2*sind(x)*cosd(x);
        -sind(x)*cosd(x), sind(x)*cosd(x),  cosd(x)^2-sind(x)^2 ];

    for j = 1:3
        if j==1
            strain_12_top = R*T*inv(R)*strain_xy_mat(:,1,i);
            
        elseif j==2
            strain_12_middle = R*T*inv(R)*strain_xy_mat(:,2,i);
            
        else
            strain_12_bottom = R*T*inv(R)*strain_xy_mat(:,3,i);
        end
    end
    
    strain_12_real = [strain_12_top, strain_12_middle, strain_12_bottom];
    strain_12_real_mat(:,:,i) = strain_12_real;
end

% Finding Matrial Axes(12) Stress in eacch ply -----------------------------------------------------
for i = 1:n
    for j = 1:3
        if j==1
            stress_12_top = Q_mat(:,:,i)*strain_12_real_mat(:,1,i);
            
        elseif j==2
            stress_12_middle = Q_mat(:,:,i)*strain_12_real_mat(:,2,i);
            
        else
            stress_12_bottom = Q_mat(:,:,i)*strain_12_real_mat(:,3,i);
        end
    end
    
    stress_12_real = [stress_12_top, stress_12_middle, stress_12_bottom];
    stress_12_real_mat(:,:,i) = stress_12_real;  
end

%% Printing result in output file ******************************************************************
fprintf(file,'\n\n===================== Ply Failure Load: %d =======================\n\n',knee);
Position = [" Top  " "Middle" "Bottom"];

fprintf(file,'1) Material Axes Strain in each ply:\n\n');
fprintf(file,'Ply No.    Position      epsilon_1     epsilon_2       gamma_12\n');
for l = 1:n
    fprintf(file,' %d (%d)\n',l, Angle(l));
    for m = 1:3
        fprintf(file,'            %s       %.3d     %.3d     %.3d\n', Position(m), strain_12_real_mat(1,m,l), strain_12_real_mat(2,m,l), strain_12_real_mat(3,m,l));
    end
end
fprintf(file,'\n-----------------------------------------------------------------\n\n');

fprintf(file,'2) Material Axes Stress(Pa) in each ply:\n\n');
fprintf(file,'Ply No.    Position       sigma_1       sigma_2       tau_12\n');
for l = 1:n
    fprintf(file,' %d (%d)\n',l, Angle(l));
    for m = 1:3
        fprintf(file,'            %s       %.3d     %.3d     %.3d\n', Position(m), stress_12_real_mat(1,m,l), stress_12_real_mat(2,m,l), stress_12_real_mat(3,m,l));
    end
end
fprintf(file,'\n-----------------------------------------------------------------\n\n');

fprintf(file,'3) Strength Ratio of each ply:\n');
fprintf(file,'(+ve SR = Tension)   (-ve SR = Compression)\n\n');
fprintf(file,'Ply No.    Position     Longitudnal     Transverse      Shear\n');
for l = 1:n
    fprintf(file,' %d (%d)\n',l, Angle(l));
    for m = 1:3
        fprintf(file,'            %s        %.4f          %.4f        %.4f\n', Position(m), SR_New_mat(1,m,l), SR_New_mat(2,m,l), SR_New_mat(3,m,l));
    end
end
fprintf(file,'\n-----------------------------------------------------------------\n\n');

fprintf(file,'4) Failed ply with mode of failure:\n');

%% Failed Ply **************************************************************************************
% Finding failed ply
index = find(SR_New_abs_mat(:) == max(SR_New_abs_mat(:)));
for i = 0:n
    Cum_index(i+1) = i*9;
end

for i = 1:size(index, 1)
    for j = 1:n
        if index(i)>Cum_index(j) && index(i)<=Cum_index(j+1)
            temp = j;
        end
    end
    failed_ply(i)=temp;
end

failed_ply_no = unique(failed_ply);
degradation = zeros(size(failed_ply_no));


% Finding Failure Mode of failed ply ---------------------------------------------------------------
temp = 0;
j = 1;
for i = 1:size(index, 1)
    if failed_ply(i)==temp
        continue
    end
        
    if mod(index(i), 3)==1
        if SR_New_mat(index(i)) > 0
            fprintf(file,'Ply %d failed in Longitudnal Tensile mode of failure\n', failed_ply(i));
        else
            fprintf(file,'Ply %d failed in Longitudnal Compressive mode of failure\n', failed_ply(i));
        end
        degradation(j) = 1;
    end
    
    if mod(index(i), 3)==2
        if SR_New_mat(index(i)) > 0
            fprintf(file,'Ply %d failed in Transverse Tensile mode of failure\n', failed_ply(i));
        else
            fprintf(file,'Ply %d failed in Transverse Compressive mode of failure\n', failed_ply(i));
        end
        degradation(j) = 2;
    end
    
    if mod(index(i), 3)==0
        fprintf(file,'Ply %d failed in Shear mode of failure\n', failed_ply(i));
        degradation(j) = 3;
    end
    
    temp = failed_ply(i);
    j = j+1; 
end
failed_ply = [];

%% Stiffness Degradation of Failed Ply *************************************************************
if Degradation_Method == 1
    for i = 1:size(failed_ply_no, 2)
        if degradation(i)==1
            Engg_Constant(failed_ply_no(i),1) = 1e-20;
        else
            Engg_Constant(failed_ply_no(i),2) = 1e-20;
            Engg_Constant(failed_ply_no(i),4) = 1e-20;
        end
    end
    
else
    for i = 1:size(failed_ply_no, 2)
        Q_mat(:,:,failed_ply_no(i)) = 1e-20;
        Q_bar_mat(:,:,failed_ply_no(i)) = 1e-20;
    end
    M = M+1;
end

%% Finding Actual Strength Ratio of each ply for next iteration ************************************
for i = 1:n
    for j = 1:3
        if stress_12_real_mat(1,j,i)>0
            SRL = stress_12_real_mat(1,j,i)/New_Strength_Parameters(i,1);
        else
            SRL = stress_12_real_mat(1,j,i)/New_Strength_Parameters(i,2);
        end
        
        if stress_12_real_mat(2,j,i)>0
            SRT = stress_12_real_mat(2,j,i)/New_Strength_Parameters(i,3);
        else
            SRT = stress_12_real_mat(2,j,i)/New_Strength_Parameters(i,4);
        end
        
        SRS = abs(stress_12_real_mat(3,j,i)/New_Strength_Parameters(i,5));
        
        if j==1
            SR_top = [SRL; SRT; SRS];
            SR_top_abs = [abs(SRL); abs(SRT); SRS];
            
        elseif j==2
            SR_middle = [SRL; SRT; SRS];
            SR_middle_abs = [abs(SRL); abs(SRT); SRS];
            
        else
            SR_bottom = [SRL; SRT; SRS];
            SR_bottom_abs = [abs(SRL); abs(SRT); SRS];
        end
    end
    
    SR = [SR_top, SR_middle, SR_bottom];
    SR_abs = [SR_top_abs, SR_middle_abs, SR_bottom_abs];
    SR_Real_mat(:,:,i) = SR;
    SR_Real_abs_mat(:,:,i) = SR_abs;
end

% Making SR value of failed ply equal to zero ------------------------------------------------------
for i = 1:size(failed_ply_no, 2)
    if degradation(i) == 1
        SR_Real_mat(1,:,failed_ply_no(i)) = 0;
        SR_Real_abs_mat(1,:,failed_ply_no(i)) = 0;
    else
        SR_Real_mat(2,:,failed_ply_no(i)) = 0;
        SR_Real_abs_mat(2,:,failed_ply_no(i)) = 0;
        SR_Real_mat(3,:,failed_ply_no(i)) = 0;
        SR_Real_abs_mat(3,:,failed_ply_no(i)) = 0;
    end
end

%Storing Strength Ratio for next iteration ---------------------------------------------------------
SR_Old_mat = SR_Old_mat + SR_Real_mat;
SR_Old_abs_mat = SR_Old_abs_mat + SR_Real_abs_mat;

%% Printing Ply Failure Load in output file --------------------------------------------------------
fprintf(file,'\n-----------------------------------------------------------------\n');
fprintf(file,'\n5) Ply Failure Load while considering Hygrothermal-Mechanical \nLoading is: %d (N/m)\n', PF_Load);

%% Increment
knee = knee + 1;   % increment in knee
k = k + size(failed_ply_no, 2); % increment in main while loop variable

end
%************************************** End of while loop ******************************************
fclose(file);

%% Plotting of Stress Vs Strain Curve **************************************************************
Stress = [0 Stress];
Strain = [0 Strain];
Cum_Stress = cumsum(Stress)./10e6;  % in MPa
Cum_Strain = cumsum(Strain);

set(gcf, 'Position', get(0,'Screensize'));
if Degradation_Method == 1
   P1 = plot(Cum_Strain, Cum_Stress,'Color','#0072BD','LineWidth', 2);
else
   P2 = plot(Cum_Strain, Cum_Stress,'Color','#D95319','LineWidth', 2);
end
    
title('\bfStress Vs Strain Curve (Partial & Complete Degradation)','FontSize', 18)
xlabel('\bfNormal Strain \epsilon_{x}','FontSize', 15)
ylabel('\bfNormal Stress N_{x}/h (MPa)','FontSize', 15)
hold on

% Marking Ply Failure
for i = 1:(size(Cum_Stress, 2) - 1)
    plot(Cum_Strain(i+1),Cum_Stress(i+1),'o','MarkerEdgeColor','red','MarkerSize',10)
end

for i = 1:(size(Cum_Stress, 2) - 1)
    if i ~= size(Cum_Stress, 2) - 1
        text(Cum_Strain(i+1),Cum_Stress(i+1),'\leftarrow Ply Failure','FontSize', 11);
    else
        text(Cum_Strain(i+1),Cum_Stress(i+1),'\leftarrow Ultimate Laminate Failure','FontSize', 11);
    end
end

for i = 1:(size(Cum_Stress, 2)-1)
    yline(Cum_Stress(i+1),'-.b');
end

Stress = [];   % Setting Stress and Strain to zero for next degradation method
Strain = [];

end
%************************************* END OF MAIN FOR LOOP ****************************************

fprintf("\nBoth Partial Degradation and Complete Degradation solution are store in separate file with respective name");
legend([P1 P2],{'Partial Degradation','Complete Degradation'},'Location','northwest');
saveas(Plot,'Stress Vs Strain (Partial & Complete)','png');

% ======================================= END OF PROGRAM ===========================================