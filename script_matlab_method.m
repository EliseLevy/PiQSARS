%% Calculation of maxima, inflexion points and slopes at the infection points on fluorescence ratio curves
%{
Written by Davy Martin and Elise Lévy, Molecular Virology and Immunology Unit (VIM-UR892), INRAE, Université Paris-Saclay, 78352 Jouy-en-Josas, France 
Use and distribution of this sript is free for academic purposes only, not for commercial use.
For more information see publication: 
Contact address: elise.levy@inrae.fr, davy.martin@inrae.fr, laurence.vernis@cnrs.fr                                               
%}

clear all;
close all;

% --- Variables that you might need to modify

filename = 'donnees_assemblages.xlsx';    % Modify if your own file name is different
AdditionTime=11;                        % It corresponds to the number of datapoints before H2O2 addition
SmoothingParameters=[15,3];             % Order and framelength of the Stavitsky-Golay filtering. 
                                        % The higher these parameters are, the less smooth the curves are

% --- Data importing

[M,Id] = xlsread(filename);
sizeM=size(M);
ncel=sizeM(2)-1;
ntime=sizeM(1);
%Id = cell identifier (1st line of the excel sheet)
% M = spread data exported by R

% ---------------------------------------------------------------------------
% Savitsky-Golay smoothing while separating before and after H2O2 addition --> Ms matrix
% ---------------------------------------------------------------------------
Ms=zeros(sizeM);
for i = 1:ncel
    Ms(:,1) = M(:,1);
    Ms(1:AdditionTime,i+1) = smoothdata(M(1:AdditionTime,i+1),'sgolay',SmoothingParameters(1),'Degree',SmoothingParameters(2));
    Ms(AdditionTime+1:ntime,i+1) = smoothdata(M(AdditionTime+1:ntime,i+1),'sgolay',SmoothingParameters(1),'Degree',SmoothingParameters(2));
end

%{
figure;
for i = 2:ncel
    hold on;
    plot(M(:,1),M(:,i))
end
figure;
for i = 2:ncel
    hold on;
    plot(Ms(:,1),Ms(:,i))
end
%} 

% Median value before addition --> vector M0

M0=zeros(1,ncel);
for i = 1:ncel
    M0(1,i) = median (Ms(1:AdditionTime,i+1));
end

% ------------------------------------
% Maxima calculation --> tmax matrix
% ------------------------------------

tmax=zeros(ncel,2);
for i = 1:ncel
    tmax(i,1) = Ms(1,1); 
    tmax(i,2) = Ms(1,i+1); 
    for k = 2:ntime 
        if Ms(k,i+1) > tmax(i,2) 
            tmax(i,2) = Ms(k,i+1);
            tmax(i,1) = Ms(k,1);
        end
    end
end
%for i = 1:ncel
%    if tmax(i,2)/M0(1,i)<1.5||tmax(i,1)>30
%        tmax(i,1)=missing;
%    end
%end

% --------------------------------
% Smoothed data 1st derivative and smoothed 1st derivative --> N and Ns matrix
% --------------------------------

D1 = diff(Ms,1);
for i = 1:ncel
    D1(:,i+1) = D1(:,i+1) ./ D1(:,1); 
end
N = D1;

for k = 1:ntime-1
    N(k,1) = (Ms(k,1) + Ms(k+1,1)) / 2; 
end

% Stavitzky-Golay smoothing of the 1st derivative

Ns=zeros(size(N));
Ns(:,1) = N(:,1);
for i = 1:ncel
    Ns(1:AdditionTime,i+1) = smoothdata(N(1:AdditionTime,i+1),'sgolay',SmoothingParameters(1),'Degree',SmoothingParameters(2));
    Ns(AdditionTime+1:ntime-1,i+1) = smoothdata(N(AdditionTime+1:ntime-1,i+1),'sgolay',SmoothingParameters(1),'Degree',SmoothingParameters(2));
end

% --------------------------------------------
% Max. inflexion point identification --> tinf
% --------------------------------------------

% tinf = ncel*3 matrix containing in the ith line 3 values for the ith cell
%   - 1st column : time where Ns maximum is reached
%   - 2nd column : derivative value at this inflextion point
%   - 3rd column : to be filled later with the ratio value at the inflexion point
% tinf is filled only for curves of which rmax >= 2*rinit & t(ratio max)<30 min

tinf=zeros(ncel,3);
for i = 1:ncel
%    if tmax(i,2)/M0(1,i)<1.5
 %       tinf(i,1)=missing;
  %      tinf(i,2)=missing;
   %     tinf(i,3)=missing;%}
    %else
        tinf(i,1) = Ns(1,1);
        tinf(i,2) = Ns(1,i+1);
        for k = 2:ntime-1
            if Ns(k,i+1) > tinf(i,2)
                tinf(i,2) = Ns(k,i+1);
                tinf(i,1) = Ns(k,1);
     %       end
        end
    end
end


% -------------------------------------------------------------------------
% lag time calculation --> tlag vector
% -------------------------------------------------------------------------

% Msint = ratio values between the observed times (ie 0.25:64.75:0.5 min)

Msint=zeros(ntime-1,ncel+1);
for i = 1:ncel
    for k = 1:ntime-1
        Msint (k,1) = (Ms(k,1) + Ms(k+1,1)) / 2;
        Msint (k,i+1) = (Ms(k,i+1) + Ms(k+1,i+1)) / 2;
    end
end

% tinf 3rd column filling

for i = 1:ncel
    for k=1:ntime-1
        if Msint(k,1) == tinf(i,1)
            tinf(i,3) = Msint(k,i+1);
        end
    end
end

% tlag = l*ncel matrix

tlag=zeros(1,ncel);
for i = 1:ncel
    tlag(i) = (M0(i) - (tinf(i,3) - (tinf(i,2) .* tinf(i,1)))) ./ tinf(i,2);
end

% -------------------------------------
% Concatenation of the data to be exported with cells identifiers
% -------------------------------------

% Cells identifiers conversion from char to string

for i = 1:ncel
    Idmat(i,1)=string(Id{i+1});
end

% Concatenation with the calculated parameters after conversion into string

R=[Idmat,num2str(tmax(:,1)),num2str(tmax(:,2)),num2str(tinf(:,1)),num2str(tinf(:,2)),num2str(tinf(:,3)),num2str(transpose(tlag))];

titles=["id_cell","tmax","rmax","tinfl","infl","rinfl","tlag"];
Rt=array2table(R);
xlswrite('data_assemblages.xlsx',titles);
writetable(Rt,'data_assemblages.xlsx','WriteVariableNames',false,'Range','A2');