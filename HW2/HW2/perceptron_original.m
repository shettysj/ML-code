clear all
data = xlsread('C:\Users\sanjana\Documents\intelligent systems\HW2\HW2\data.xlsx');
data1(:,1)=data(:,1);
data1(:,2)=data(:,2);
data2(:,1)=data(:,3);

TrainingSet=zeros(240,2);
TestingSet=zeros(60,2);


normL=[];
normP=[];

 Q11=0;
 Q10=0;
 Q01=0;
 Q00=0;
 count=1;

%normalize data%
minL= min(data1(:,1));
maxL= max(data1(:,1));
minP= min(data1(:,2));
maxP= max(data1(:,2));

rangeL= maxL-minL;
rangeP= maxP-minP;

for i=1:300
    normL(i,1)= (data1(i,1)-minL)/rangeL;
    normP(i,1)= (data1(i,2)-minL)/rangeP;
end



%%%%%%%%%% assign values to training set and test set%%%%%%%

TrainingEpochError=[];
NoEpochChange=[];
%dividing data into training and testing set%
TrainingSet(:,1)=normL(1:240,1); %normL
TrainingSet(:,2)=normP(1:240,1); %normP

TestingSet(:,1)=normL(241:300,1);
TestingSet(:,2)=normP(241:300,1);


TestingError_1=[];
%TrainingOutput=zeros(1:240,1);
nep=5:5:100;
%initializing weights and bias for training%
w1=0;%rand(1,1);
w2=0;%rand(1,1);
b=0;%rand(1,1);

%initializing learning rate%
n=0.01;

%initialising inputs and outputs%

%training the perceptron for 10:100 epochs%
for limit=nep
%     TrainingEpochError=zeros();
y=zeros();
TrainingError=0;
ExpectedOutput=data2(1:240,:);

for epoch=1:limit
    
    x1= TrainingSet(:,1);
    x2= TrainingSet(:,2);
    TrainingOutput=zeros();
    
for i=1:240%length(TrainingSet)
    
    
     y(i)= x1(i)*w1 + x2(i)*w2 + b;
    
    if(y(i)>=1)
        TrainingOutput(i,1)= 1;
    else
        TrainingOutput(i,1)= 0;
    end;
   
  TrainingError=  ExpectedOutput(i,1)-TrainingOutput(i,1) ;
  
  w1= w1 + TrainingError*n*x1(i);
  w2= w2 + TrainingError*n*x2(i);
  b= b + n*(TrainingError);
 end

TrainingEpochError(epoch,trial) = (nnz(ExpectedOutput - TrainingOutput))/240;
end
for r=1:240
    if TrainingOutput(r)==1
        if  TrainingOutput(r)== ExpectedOutput(r)
            Q11=Q11+1;
        else
            Q01=Q01+1;
        end
    else
       if  TrainingOutput(r)== ExpectedOutput(r)
           Q00=Q00+1;
        else
            Q10=Q10+1;
       end
    end
end

   hit_rate_train(limit)=(Q11+Q00)/(Q11+Q10+Q01+Q00); 
% figure();
% scatter(1:1:10,TrainingEpochError);
% xlabel('error');
% ylabel('number');
% 
% figure();
% scatter(1:240,ExpectedOutput(1:240));
% hold on
% scatter(1:240,TrainingOutput(1:240),'g');
% hold off
% xlabel('test case');
% ylabel('output');

% NoEpochChange(end+1) = TrainingEpochError(end);
% end

% figure();
% scatter(10:10:100,NoEpochChange);
% xlabel('error');
% ylabel('number of epochs');

%%%%%%%%%testing the perceptron%%%%%%%%%%%%

    x1= TestingSet(:,1);
    x2= TestingSet(:,2);
    TestingOutput=zeros(1:60,1);
    y=zeros();
    ExpectedOutput=data2(241:300,:);

for i=1:60%length(TrainingSet)
    
    
    y(i)= x1(i)*w1 + x2(i)*w2 + b;
    
    if(y(i)>=1)
        TestingOutput(i,1)= 1;
    else
        TestingOutput(i,1)= 0;
    end;
   
  TestingError_1(end+1)= ExpectedOutput(i)- TestingOutput(i,1);
  
  %check=nnz(TestingError_1);
end
%%TestingErrorPerc = (nnz(ExpectedOutput - TestingOutput))/60;

Q11=0;
 Q10=0;
 Q01=0;
 Q00=0;
 
for r=1:60
    if TestingOutput(r)==1
        if  TestingOutput(r)== ExpectedOutput(r)
            Q11=Q11+1;
        else
            Q01=Q01+1;
        end
    else
       if  TestingOutput(r)== ExpectedOutput(r)
           Q00=Q00+1;
        else
            Q10=Q10+1;
       end
    end
end

   hit_rate_test(limit)=(Q11+Q00)/(Q11+Q10+Q01+Q00);
   
end
%end
hit_rate_test(hit_rate_test==0)=nan;
hit_rate_train(hit_rate_train==0)=nan;

figure();
scatter(1:100,1-hit_rate_test);
hold on
scatter(1:100,1-hit_rate_train,'r');
hold off
ylabel('Error rate');
xlabel('number of epochs');
title('test and training set');


legend('test','training');

