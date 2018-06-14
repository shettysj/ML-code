%------------------ CODE WRITTEN BY AKSHAY LAKKARAJU ---------------------% 
%---------------------- PERCEPTRON ALGORITHM -----------------------------% 
clear; clc; 
%--------------- CREATING TEST & TRAINING SETS----------------------------%
data = xlsread('C:\Users\Akshay\Documents\MATLAB\Home Work 2\Training Set.xlsx'); 
% Data normalization % lmax = max(data(:,1)); lmin = min(data(:,1)); pmax = max(data(:,2)); pmin = min(data(:,2)); lNorm = (data(:,1)-lmin)/(lmax-lmin); pNorm = (data(:,1)-pmin)/(pmax-pmin); 
%data(:,1) = lNorm; data(:,2) = pNorm; 
% Initializing test set %
testingRows = zeros(9,60);
testing = zeros(60,3,9); 
% Initializing training set % 
variable = 1:1:300; 
trainingRows = zeros(9,240); 
training = zeros(240,3,9);
for i=1:9 
    %-------------- TEST SET ------------------% 
    testingRows(i,:) = randperm(300,60);
    testingRows(i,:) = sort(testingRows(i,:));
    for j=1:60 
        testing(j,:,i) = data(testingRows(i,j),:);
    end
    %------------- TRAINING SET ---------------% 
    common = intersect(variable',testingRows(i,:)); 
    trainingRows(i,:) = setxor(variable',common);
    for k=1:240 training(k,:,i) = data(trainingRows(i,k),:); 
    end
end
%-------------------------------------------------------------------------% 
%-------------------------- MAIN ALGORITHM -------------------------------% 
%-------------------------------------------------------------------------% 
% Error Curve % 
errorPercentageTraining = zeros(10,1,9); 
sensitivityTr = zeros(9,1);
sensitivityTe = zeros(9,1); 
specificityTr = zeros(9,1); 
specificityTe = zeros(9,1); 
PPVTr = zeros(9,1); 
PPVTe = zeros(9,1); 
NPVTr = zeros(9,1); 
NPVTe = zeros(9,1); 
fig = 1;
for set=1:9
    %-------------------------- TRAINING -------------------------------------%
    % Initializing weights %
    w1 = 0; w2 = 0; b = 0;
    % Initializing learning rate % 
    a = 0.025; y = 0; 
    % Initializing outputs % 
    outputTraining = zeros(240,1); 
    errorTraining = 0; 
    for iterations=1:10 
        x1 = training(:,1,set); 
        x2 = training(:,2,set); 
        for i=1:240 y = x1(i)*w1 + x2(i)*w2 + b; 
            if(y>=1) outputTraining(i,1) = 1; 
            else outputTraining(i,1) = 0;
            end
            % Adjusting the weights % 
            errorTraining = training(i,3,set) - outputTraining(i,1); 
            w1 = w1 + errorTraining*a*x1(i);
            w2 = w2 + errorTraining*a*x2(i);
            b = b + a*errorTraining; 
        end
        errorPercentageTraining(iterations,1,set) = (nnz(outputTraining-training(:,3,set)))/240; 
    end
    figure(fig)
    plot(1:1:10,errorPercentageTraining(:,1,set)) 
    name = sprintf('Error Curve for Perceptron For Set: %d',set);
    title(name) 
    xlabel('Number of Iterations') 
    ylabel('Error Percentage')
    % True Positives & False Negatives % 
    onesOutput = find(outputTraining); 
    onesActual = find(training(:,3,set)); 
    commonOnes = size(intersect(onesOutput,onesActual)); 
    truePositivesTr = commonOnes(1);
    falseNegativesTr = nnz(outputTraining)-truePositivesTr; 
    % True Negatives & False Positives %
    zerosOutput = find(ones(240,1)-outputTraining);
    zerosActual = find(ones(240,1)-training(:,3,set)); 
    commonZeros = size(intersect(zerosOutput,zerosActual)); 
    trueNegativesTr = commonZeros(1); 
    falsePositivesTr = nnz(outputTraining)-trueNegativesTr; 
    % Testing Data % 
    sensitivityTr(set,1) = truePositivesTr/(truePositivesTr + falseNegativesTr);
    specificityTr(set,1) = trueNegativesTr/(falsePositivesTr + trueNegativesTr);
    PPVTr(set,1) = truePositivesTr/(truePositivesTr + falsePositivesTr);
    NPVTr(set,1) = trueNegativesTr/(trueNegativesTr + falseNegativesTr);
%-------------------------- TESTING --------------------------------------% 
% Initializing outputs % 
outputTesting = zeros(60,1); 
y = 0; x1 = testing(:,1,set);
x2 = testing(:,2,set);
for i=1:60 y = x1(i)*w1 + x2(i)*w2 + b; 
    if(y>=1) outputTesting(i,1) = 1;
    else outputTesting(i,1) = 0; 
    end
end % True Positives & False Negatives % 
onesOutput = find(outputTesting); 
onesActual = find(testing(:,3,set));
commonOnes = size(intersect(onesOutput,onesActual)); 
truePositivesTe = commonOnes(1); 
falseNegativesTe = nnz(outputTesting)-truePositivesTe; 
% True Negatives & False Positives %
zerosOutput = find(ones(60,1)-outputTesting); 
zerosActual = find(ones(60,1)-testing(:,3,set)); 
commonZeros = size(intersect(zerosOutput,zerosActual)); 
trueNegativesTe = commonZeros(1); 
falsePositivesTe = nnz(outputTesting)-trueNegativesTe; 
% Testing Data %
sensitivityTe(set,1) = truePositivesTe/(truePositivesTe + falseNegativesTe); 
specificityTe(set,1) = trueNegativesTe/(falsePositivesTe + trueNegativesTe); 
PPVTe(set,1) = truePositivesTe/(truePositivesTe + falsePositivesTe); 
NPVTe(set,1) = trueNegativesTe/(trueNegativesTe + falseNegativesTe); 
fig = fig+1;
end 
% %---------------------- PLOTTING GRAPHS ----------------------------------% % Trialwise Errors % errorPercentagePlotting = zeros(10,9); for set=1:9 errorPercentagePlotting(:,set) = errorPercentageTraining(:,1,set); end figure(fig) plot(1:1:10,errorPercentagePlotting) name = sprintf('Error Curve for Perceptron'); title(name) xlabel('Number of Iterations') ylabel('Error Percentage') legend('Set 1', 'Set 2', 'Set 3', 'Set 4', 'Set 5', 'Set 6', 'Set 7', 'Set 8', 'Set 9') fig = fig+1; % Mean Training Error % figure(fig) meanErrorPercentagePlotting = mean(errorPercentagePlotting,2); standardDevPercentagePlotting = std(errorPercentagePlotting,0,2); errorbar(1:1:10,meanErrorPercentagePlotting,standardDevPercentagePlotting) name = sprintf('Mean Training Error for Perceptron');
% title(name) xlabel('Number of Iterations') ylabel('Mean Training Error') fig = fig+1; % Sensitivity % figure(fig) bar([sensitivityTr sensitivityTe]) name = sprintf('Sensitivity Training Vs Testing'); title(name) xlabel('Set Number') ylabel('Sensitivity') legend('Training', 'Testing') fig = fig+1; % Specificity % figure(fig) bar([specificityTr specificityTe]) name = sprintf('Specificity Training Vs Testing'); title(name) xlabel('Set Number') ylabel('Specificity') legend('Training', 'Testing') fig = fig+1; % PPV % figure(fig) bar([PPVTr PPVTe]) name = sprintf('PPV Training Vs Testing'); title(name) xlabel('Set Number') ylabel('PPV') legend('Training', 'Testing') fig = fig+1; % NPV % figure(fig) bar([NPVTr NPVTe]) name = sprintf('NPV Training Vs Testing'); title(name) xlabel('Set Number') ylabel('NPV') legend('Training', 'Testing') %--------------- AVERAGE PERFORMANCE METRICS -----------------------------% % Training Mean & Standard Deviation % meanSensitivityTr = mean(sensitivityTr) standardDevSensitivityTr = std(sensitivityTr) meanSpecificityTr = mean(specificityTr) standardDevSpecificityTr = std(specificityTr) meanPPVTr = mean(PPVTr) standardDevPPVTr = std(PPVTr) meanNPVTr = mean(NPVTr) standardDevNPVTr = std(NPVTr) % Testing Mean & Standard Deviation % meanSensitivityTe = mean(sensitivityTe) standardDevSensitivityTe = std(sensitivityTe) meanSpecificityTe = mean(specificityTe) standardDevSpecificityTe = std(specificityTe) meanPPVTe = mean(PPVTe) standardDevPPVTe = std(PPVTe) meanNPVTe = mean(NPVTe) standardDevNPVTe = std(NPVTe)