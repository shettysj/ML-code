clear all;
data = xlsread('C:\Users\sanjana\Documents\intelligent systems\HW2\HW2\data.xlsx');
data1(:,1)=data(:,1);
data1(:,2)=data(:,2);
data2(:,1)=data(:,3);


%%%%%%%%%%%%%%%%%%%%%%% arrays needed to store hit rate and error%%%%%%%%%
hit_rate_knn=[];
hit_rate_rad=[];
hit_rate_percep=[];
error_rate_knn=[];
erro_rate_rad=[];
erro_rate_percep=[];

%%%%%%%%%%%%%%%% sensi_tivity, NPV, PPV, specificity%%%%%%%%%%%%%%%
sensi_tivity=[];
specificity=[];
NPV=[];
PPV=[];



for trial=1:9
    
    TotalSet=zeros(300,3);
    TrainingSet=zeros(240,2);
    TestSet=zeros(60,2);
    count=1;
 %%%%%%%%%%% generate array of random numbers%%%%%%% 
 for i=1:300
test_1(i)=i;
end
for i=1:300
test_2=test_1(randperm(length(test_1)));
end

%%%%%%%%%% assign values to training set and test set%%%%%%%

for i=1:300
    TotalSet_percep(i,1:2)= data1(test_2(i),:);
    TotalSet_percep(i,3) = data2(test_2(i),:);
end

for i=1:240
    TrainingSet(i,:)= data1(test_2(i),:);
    TrainingSet2(i,:)= data2(test_2(i),:);
end
for i=241:300
    TestSet(count,:)= data1(test_2(i),:);
    TestSet2(count,:)= data2(test_2(i),:);
    count=count+1;

end
% 
 TotalSet(1:240,1:2)= TrainingSet(1:240,:);
  TotalSet(241:300,1:2)=  TestSet(1:60,:);
  TotalSet(1:240,3)= TrainingSet2(1:240,1);
  TotalSet(241:300,3)=  TestSet2(1:60,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KNN%%%%%%%%%%%%%%%%%%%%%%%%%%%%


j=1; i=1;
eucdist_sq=[];
temp=zeros(2,300);
data3=[];
k=9;                         %%%%%%best value of k
negative=zeros(300,6);
positive=zeros(300,6);

negative_rad=zeros(300);
positive_rad=zeros(300);

final_result=[];

final_result_rad=[];


for j=1:60
for i=1:240
    if j~=i
        eucdist_sq(i,j)=((TrainingSet(i,1)-TestSet(j,1))^2)+((TrainingSet(i,2)-TestSet(j,2))^2);
     
    else
    
        eucdist_sq(i,j)=0; 
    end;
end;
end;




   
    
    for c=1:60          %positive and negatives for each value; c corresponds to the value
         tes=k;
         j=1;
        TrainingSet2(:,2)=eucdist_sq(:,c);
        data3=sortrows(TrainingSet2,2);
    while(tes~=0)
    if(data3(j,2)~=0)
       if(data3(j,1)==1)
        positive(c)=positive(c)+1;
       else
           negative(c)=negative(c)+1;
       end
       tes=tes-1;
    end
    
    j=j+1;
    
     
    
    end
    
    
    end
   % hit_rate(i)= (sum(positive) +sum(negative))/



%%%%%%%%%%%%%%%%%% hit rate%%%%%%%%%%%%%%%%%%%%%%
    Q10=0;
    Q01=0;
    Q11=0;
    Q00=0;
    for j=1:60
        if positive(j)>=negative(j)
            final_result(j)=1;
             if TestSet2(j)==final_result(j)
                   Q11=Q11+1;
             else
                 Q01=Q01+1;
             end
        else
            final_result(j)=0;
            if TestSet2(j)==final_result(j)
                   Q00=Q00+1;
             else
                 Q10=Q10+1;
            end
        end
        
    end
    
    hit_rate_knn(trial)=(Q11+Q00)/(Q10+Q01+Q11+Q00);
    sensi_tivity(trial,1)= (Q11)/(Q11+ Q10);
    specificity(trial,1)= (Q00)/(Q01+ Q00);
    PPV(trial,1)= (Q11)/(Q11+ Q01);
    NPV(trial,1)= (Q00)/(Q00+Q10);
    error_rate_knn(trial)=1-hit_rate_knn(trial);
    





           
Y_mean=mean(mean(eucdist_sq));
Y_min=min(min(eucdist_sq(eucdist_sq>0)));
Y_max=max(max(eucdist_sq));
Y_std= std(std(eucdist_sq));

data2_neighbor=zeros();
data3=zeros();
i=1;

radius=0: 20 : Y_std;

i=1;
neighbor=zeros(300);
inc=100;                                %%%%%%%%%%%%%%%%optimum radius%%%%%%%%%%%%
  
  i=i+1;
%end
   
    
   for c=1:60          %positive and negatives for each value; c corresponds to the value
         tes=k;
         j=1;
         neighbor=0;
        TrainingSet2(:,2)=eucdist_sq(:,c);
        data3=sortrows(TrainingSet2,2);
     while(data3(j,2)<=inc)
          if(data3(j,2)~=0)
              neighbor=neighbor+1;
              if(neighbor>300)
                  break;
              end
             if(data3(j,1)==1)
                  positive_rad(c)=positive_rad(c)+1;
             else
                  negative_rad(c)=negative_rad(c)+1;
             end
          end
    
    j=j+1;
    
     
    
    end
   if neighbor==0
       neighbor=1;
        positive_rad(c)=1;
      
        negative_rad(c)=0;
   end
        
    
end
  


  
    Q10=0;
    Q01=0;
    Q11=0;
    Q00=0;
    for j=1:60
        if positive_rad(j)>=negative_rad(j)
            final_result_rad(j)=1;
             if TestSet2(j)==final_result_rad(j)
                   Q11=Q11+1;
             else
                 Q01=Q01+1;
             end
        else
            final_result_rad(j)=0;
            if TestSet2(j)==final_result_rad(j)
                   Q00=Q00+1;
             else
                 Q10=Q10+1;
            end
        end
    end
           hit_rate_rad(trial)=(Q11+Q00)/(Q10+Q01+Q11+Q00);
           sensi_tivity(trial,2)= (Q11)/(Q11+ Q10);
            specificity(trial,2)= (Q00)/(Q01+ Q00);
            PPV(trial,2)= (Q11)/(Q11+ Q01);
            NPV(trial,2)= (Q00)/(Q00+Q10);
            error_rate_rad(trial)=1-hit_rate_rad(trial);


disp(hit_rate_knn);
disp(hit_rate_rad);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%perceptron%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% data = xlsread('C:\Users\sanjana\Documents\intelligent systems\HW2\HW2\data.xlsx');
% data1(:,1)=data(:,1);
% data1(:,2)=data(:,2);
% data2(:,1)=data(:,3);

TrainingSet=zeros(240,2);
TestingSet=zeros(60,2);
TotalSet=zeros(300,3);
%for trial=1:9
% for i=1:300
% test_1(i)=i;
% end
% for i=1:300
% test_2=test_1(randperm(length(test_1)));
% end


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
minL= min(TotalSet_percep(:,1));
maxL= max(TotalSet_percep(:,1));
minP= min(TotalSet_percep(:,2));
maxP= max(TotalSet_percep(:,2));

rangeL= maxL-minL;
rangeP= maxP-minP;

for i=1:300
    normL(i,1)= (TotalSet_percep(i,1)-minL)/rangeL;
    normP(i,1)= (TotalSet_percep(i,2)-minL)/rangeP;
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
ExpectedOutput=TotalSet_percep(1:240,3);

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

   hit_rate_train_percep(trial)=(Q11+Q00)/(Q11+Q10+Q01+Q00); 
   sensi_tivity(trial,3)= (Q11)/(Q11+ Q10);
    specificity(trial,3)= (Q00)/(Q01+ Q00);
    PPV(trial,3)= (Q11)/(Q11+ Q01);
    NPV(trial,3)= (Q00)/(Q00+Q10);
    error_rate_train_percep(trial)=1-hit_rate_train_percep(trial);
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
    ExpectedOutput=TotalSet_percep(241:300,3);

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

   hit_rate_percep(trial)=(Q11+Q00)/(Q11+Q10+Q01+Q00);
   sensi_tivity(trial,4)= (Q11)/(Q11+ Q10);
    specificity(trial,4)= (Q00)/(Q01+ Q00);
    PPV(trial,4)= (Q11)/(Q11+ Q01);
    NPV(trial,4)= (Q00)/(Q00+Q10);
    error_rate_percep(trial)=1-hit_rate_knn(trial);
   
end

disp(hit_rate_percep);


%%%%%%%%%%%%%%%%%%PLOTS%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%perceptron%%%%%%%%%%%%%%%%%%%%%%%
figure('name','Perceptron Sensitivity after training and after testing');
title('Bar Graph Sensitivity vs Trial Number Perceptron');
 bar([sensi_tivity(1:9,3),sensi_tivity(1:9,4)]);
 xlabel('Trial Number');
ylabel('Sensitivity');
legend('training', 'testing');

figure('name','Perceptron specificity after training and after testing');
title('Bar Graph Specificity vs Trial Number Perceptron');
 bar([specificity(1:9,3),specificity(1:9,4)]);
 xlabel('Trial Number');
ylabel('Specificity');
legend('training', 'testing');

figure('name','Perceptron NPV after training and after testing');
title('Bar Graph NPV vs Trial Number Perceptron');
 bar([NPV(1:9,3),NPV(1:9,4)]);
 xlabel('Trial Number');
ylabel('NPV');
legend('training', 'testing');

figure('name','Perceptron PPV after training and after testing');
title('Bar Graph PPV vs Trial Number Perceptron');
 bar([PPV(1:9,3),PPV(1:9,4)]);
 xlabel('Trial Number');
ylabel('PPV');
legend('training', 'testing');

%%%%%%%%%%%%%%%%%%Nearest neighbor%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('name','Nearest neighbor Sensitivity ');
title('NEC vs Trial Number ');
 bar([sensi_tivity(1:9,2),specificity(1:9,1),NPV(1:9,1),PPV(1:9,1)]);
 xlabel('Trial Number');
legend('sensitivity','specificity','NPV','PPV');


%%%%%%%%%%%%%%%%KNN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
figure('name','KNN Sensitivity ');
title('KNN vs Trial Number ');
 bar([sensi_tivity(1:9,1),specificity(1:9,1),NPV(1:9,1),PPV(1:9,1)]);
 xlabel('Trial Number');
ylabel('Sensitivity');
legend('sensitivity','specificity','NPV','PPV');



%%%%%%%%%%%%%%%mean 
%Calculating the K-Boundary
% minL=min(data1(:,1));
% maxL=max(data1(:,1));
% featureSpaceL=[minL:0.1:maxL];
% 
% minP=min(data1(:,2));
% maxP=max(data1(:,2));
% rangeP=minP-maxP;
% featureSpaceP=[minP:0.1:maxP];
% 
% positive1=0;
% negative1=0;
% boundary=[];
% euc_dist=[];
% for sweepL=featureSpaceL
%     for sweepP=featureSpaceP
%         for l=1:length(TrainingSet)
%             euc_dist(end+1)= ((TrainingSet(i,1)- sweepL)^2)+((TrainingSet(i,2)- sweepP)^2);
%         
%             end;
%             tes=k;
%          j=1;
%          
%         TrainingSet2(:,2)=euc_dist(:,c);
%         data3=sortrows(TrainingSet2,2);
%    
%         while(tes~=0)
%     if(data3(j,2)~=0)
%        if(data3(j,1)==1)
%         positive1=positive1+1;
%        else
%            negative1=negative1+1;
%        end
%        tes=tes-1;
%     end
%     
%     j=j+1;
%         end
%         
%         if(positive==floor(k/2))
%             boundary=(boundary(sweepL,sweepP));
%         end
% end
% end
%  figure('Name','Figure 1: KNN Boundary','NumberTitle','off');
% scatter(boundary(1,:),boundary(2,:));
% title('Pnormalized v/s Lnormalized');
% xlabel('L');
% ylabel('P');
% hold   

    %%%%%%%%%%KNN%%%%%%%%%%%
mean_arr(1,1)=mean(NPV(:,1));
mean_arr(1,2)=mean(sensi_tivity(:,1));
mean_arr(1,3)=mean(specificity(:,1));
mean_arr(1,4)=mean(PPV(:,1));
mean_arr(1,5)=mean(hit_rate_knn);

sd_arr(1,1)=std(NPV(:,1));
sd_arr(1,2)=std(sensi_tivity(:,1));
sd_arr(1,3)=std(specificity(:,1));
sd_arr(1,4)=std(PPV(:,1));
sd_arr(1,5)=std(hit_rate_knn);

%%%%%%%%%%NEC%%%%%%%%%%%%%
mean_arr(2,1)=mean(NPV(:,2));
mean_arr(2,2)=mean(sensi_tivity(:,2));
mean_arr(2,3)=mean(specificity(:,2));
mean_arr(2,4)=mean(PPV(:,2));
mean_arr(2,5)=mean(hit_rate_rad);

sd_arr(2,1)=std(NPV(:,2));
sd_arr(2,2)=std(sensi_tivity(:,2));
sd_arr(2,3)=std(specificity(:,2));
sd_arr(2,4)=std(PPV(:,2));
sd_arr(2,5)=std(hit_rate_rad);

%%%%%%%%%%%%percep%%%%%%%%%%%
mean_arr(3,1)=mean(NPV(:,4));
mean_arr(3,2)=mean(sensi_tivity(:,4));
mean_arr(3,3)=mean(specificity(:,4));
mean_arr(3,4)=mean(PPV(:,4));
mean_arr(3,5)=mean(hit_rate_percep);

sd_arr(3,1)=std(NPV(:,4));
sd_arr(3,2)=std(sensi_tivity(:,4));
sd_arr(3,3)=std(specificity(:,4));
sd_arr(3,4)=std(PPV(:,4));
sd_arr(3,5)=std(hit_rate_percep);



disp(mean_arr);
disp(sd_arr);

xlswrite('tst.xls',sd_arr);