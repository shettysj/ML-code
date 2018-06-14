clear all
data = xlsread('C:\Users\sanjana\Documents\intelligent systems\HW2\HW2\data.xlsx');
data1(:,1)=data(:,1);
data1(:,2)=data(:,2);
data2(:,1)=data(:,3);

TrainingSet=zeros(240,2);
TestingSet=zeros(60,2);
TotalSet=zeros(300,3);
for trial=1:9
for i=1:300
test_1(i)=i;
end
for i=1:300
test_2=test_1(randperm(length(test_1)));
end

for i=1:300
    TotalSet(i,1:2)= data1(test_2(i),:);
    TotalSet(i,3) = data2(test_2(i),:);
end
% for i=1:240
%     TrainingSet(i,1)= data1(test_2(i),:);
%     TrainingSet(i,3)= data2(test_2(i),:);
% end
% for i=241:300
%     TestingSet(count,1)= data1(test_2(i),:);
%     
%     TestSet2(count,1)= data2(test_2(i),:);
%     count=count+1;
% 
% end

%for train=1:9
normL=[];
normP=[];

 Q11=0;
 Q10=0;
 Q01=0;
 Q00=0;
 count=1;

%normalize data%
minL= min(TotalSet(:,1));
maxL= max(TotalSet(:,1));
minP= min(TotalSet(:,2));
maxP= max(TotalSet(:,2));

rangeL= maxL-minL;
rangeP= maxP-minP;

for i=1:300
    normL(i,1)= (TotalSet(i,1)-minL)/rangeL;
    normP(i,1)= (TotalSet(i,2)-minL)/rangeP;
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
nep=10:10:100;
%initializing weights and bias for training%
w1=0;%rand(1,1);
w2=0;%rand(1,1);
b=0;%rand(1,1);

%initializing learning rate%
n=0.01;

%initialising inputs and outputs%
y=zeros();
TrainingError=0;
ExpectedOutput=TotalSet(1:240,3);

%training the perceptron for 10:100 epochs%
% for limit=nep
%     TrainingEpochError=zeros();
for epoch=1:10 %limit
    
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

TrainingEpochError(epoch) = (nnz(ExpectedOutput - TrainingOutput))/240;
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

   hit_rate_train=(Q11+Q00)/(Q11+Q10+Q01+Q00); 
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
    ExpectedOutput=TotalSet(241:300,3);

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
TestingErrorPerc = (nnz(ExpectedOutput - TestingOutput))/60;

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

   hit_rate_test(trial)=(Q11+Q00)/(Q11+Q10+Q01+Q00);
   
end
%end
% figure();
% plot(1:60,ExpectedOutput);
% hold on
% plot(1:60,TestingOutput,'g');
% hold off
% xlabel('test case');
% ylabel('output');