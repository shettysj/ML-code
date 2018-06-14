clear all;
 
count=1;
input_num=784;
data_points=5000;
num_neurons_hidden=100;
num_neurons_output=10;   %%%% same as number of images %%%%
trainSubsetSize=100;     %%%%% number of data points to be used in each epoch%%%
 
testSetSize = 1000;
trainSet_num = data_points-testSetSize;
 
thresh_high=1;%0.75;
thresh_low=0;%0.25;
 
label_d2=input_num+1;
 
weight_min=-1;
weight_max=1;
weight_range=weight_max-weight_min;
 
learningRate_ouput=0.05;
learningRate_hidden=0.05;
 
alpha=0.6;  %%%%%%%%%% momentum %%%%%%%%
delta_weight_ij=zeros();
delta_weight_ij_t_1=zeros();
 
sum_output_hidden=0;
sum_output=0;
 
epoch_num=1;
epoch_store=1;
 
conf_matrix_train=zeros(10,10);
conf_matrix_test=zeros(10,10); 
 
%MNISTnumImages5000
%MNISTnumLabels5000
%%%%%%%%%%% read data from text file %%%%%%%%%%%%
%image intensity%%%%%%%%%%%%%%%%5
fid_image = fopen('MNISTnumImages5000.txt');
formatSpecI = '%f %f';
sizeI = [input_num data_points];
I= fscanf(fid_image,formatSpecI,sizeI);
Data_orig_image=I';
fclose(fid_image);
 
%image label%%%%%%%%%%%%%%%%%%%%
fid_label = fopen('MNISTnumLabels5000.txt');
formatSpecL = '%f %f';
sizeL = [1 data_points];
L= fscanf(fid_label,formatSpecL,sizeL);
Data_orig_label=L';
fclose(fid_label);
 
% % %%%%%%%%%% normalize the inputs %%%%%%%%%%%%%%
% 
% Norm_min= min(Data_orig_image(:,:));
% Norm_max= max(Data_orig_image(:,:));
% 
% Norm_range= Norm_max-Norm_min;
% 
% for i=1:data_points
%     for j=1:input_num
% Norm_data(i,j)= (Data_orig_image(i,j)-Norm_min)/Norm_range;%(Data_orig_image(i,j)*2)-1;
%     end
% end
 
% %%%%%%%%%%%%%%%%%% select random samples for testing and training%%%%%%%%%
% 
% % %%%%%%%%%%% generate array of random numbers%%%%%%%
for i=1:data_points
test_1(i)=i;
end
for i=1:data_points
test_2=test_1(randperm(length(test_1)));
end
%%%%%%%%%% assign values to training set and test set%%%%%%%
for i=1:data_points
TotalSet_image(i,1:input_num)= Data_orig_image(test_2(i),:);
 
TotalSet_image(i,(label_d2))= Data_orig_label(test_2(i),:);    %%%% storing the corresponding labels in the totalset %%%%%
 
end
 
for i=1: trainSet_num
TrainingSet(i,:)= TotalSet_image(i,:);
 
end
for i=trainSet_num+1:data_points
TestSet(count,:)= TotalSet_image(i,:);
count=count+1;
end
 
%%%%%%%%%%% random weight generation %%%%%%%%%%%%%%%%%%
 
weight_input_hidden_jk= ((weight_range)*(rand(num_neurons_hidden,input_num)))+weight_min;
weight_hidden_output_ij= ((weight_range)*(rand(num_neurons_output,num_neurons_hidden)))+weight_min;
 
%%%%%%%%%%%%%%%% setting label vector %%%%%%%%%%
%%%%%%% it sets a 1 to the point which corresponds to the label %%%%%%%%
%%%%%%% example if label is 7, 1 will be set to the 7th element of the
%%%%%%% array
%%%%%%% this is then treated as expected output from the 10 neurons and
%%%%%%% used for error correction
for c=1: data_points
    temp_label_store= TotalSet_image(c,label_d2);
    for i=1:num_neurons_output
        if i == (temp_label_store+1)
            expected_label_array(c,i) = 1;
        else
            expected_label_array(c,i) = 0;
        end
    end
end
 
expected_label_train = expected_label_array(1:trainSet_num,:);
expected_label_test = expected_label_array((trainSet_num+1):data_points,:);
 
% %%%%%%% choose a random start and end value for subsets on which to train
% %%%%%%% each epoch  %%%%%%%%%%%%%%
% start_trainSubset=randi((trainSet_num-trainSubsetSize));
% end_trainSubset= start_trainSubset+trainSubsetSize;
 
%%%%%%%%%%%%%% compute output %%%%%%%%%%%%%%%%%%%%%%
for e=1:epoch_num
    
%%%%%%% choose a random start and end value for subsets on which to train
%%%%%%% each epoch  %%%%%%%%%%%%%%
start_trainSubset=randi((trainSet_num-trainSubsetSize));
end_trainSubset= start_trainSubset+trainSubsetSize;
 
    
 d2=1;  %%%%%%%%% used as d2 for output and hidden neuron%%% 
  Q11=0;
%%%%%%%% loop for data_points %%%%%%%%%
for d= start_trainSubset:end_trainSubset
    
    
     
    %%%%%%%%%% Forward pass %%%%%%%%%%%%                
   %%%% output layer neuron i %%%% 
  
    for i=1: num_neurons_output
        
   %%%% hidden layer neuron j %%%%
       sum_output=0;
       for j=1:num_neurons_hidden
           sum_output_hidden=0;
          for k= 1:input_num 
             sum_output_hidden= sum_output_hidden+(weight_input_hidden_jk(j,k) * TrainingSet(d,k));
          end
         hidden_output_fj(d2,j)=1/(1+exp(-sum_output_hidden));
       
         sum_output= sum_output + (weight_hidden_output_ij(i,j) * hidden_output_fj(d2,j));
       end
       actual_output_fi(d2,i)=1/(1+exp(-sum_output));
       
         %%%%%%%%%%%%%%%%%%% creating output labelling %%%%%%%%%%%%%%%%
     if actual_output_fi(d2,i) >=thresh_high        %%%%% set threshold; anything greater than 0.75 is considered as high                     
         output_label(d2,i)=1;
         
     elseif actual_output_fi(d2,i) <=thresh_low   %%%%%%% anything lesser -0.75 is considered as low (0)
         output_label(d2,i)=0;
         
     else
         output_label(d2,i)=actual_output_fi(d2,i);
     end
    end
    
    %%%%%%%%%%%%%% start loops for backpropogation %%%%%%%%%%%
    
     %%%%% below part for momentum %%%%%%%
    delta_weight_ij_t_1= delta_weight_ij;
    %%%%%% above part for momentum %%%%%
    
       for j=1:num_neurons_hidden  %%%%% loop for change in weights %%%%%%%%%%
           
           sum_delta_i=0;
           
           for i2 = 1:num_neurons_output
     %%%%%%% find error by comparing expected output with acquired output %%
             error_output(d2,i2) = expected_label_train(d,i2)-output_label(d2,i2);
             
             
             
    
    %%%%%%%%%% find delta i %%%%%%%%%%%%%%%%%
             delta_output_i (d2,i2) = (1-actual_output_fi(d2,i2))*actual_output_fi(d2,i2)*(error_output(d2,i2));
    
   
    %%%%%%%%%%%%%%% find weight change required for weights between output
    %%%%%%%%%%%%%%% layer and hidden layer %%%%%%%%%%%
    
   
    %%% without momentum: delta_weight_ij(i2,j) = learningRate_ouput * delta_output_i (d2,i2) * hidden_output_fj(d2,j);
    %%%%% with momentum: %%%%%%%%%%%%%%%%%%
    if d2==1      
        delta_weight_ij(i2,j) = learningRate_ouput * delta_output_i (d2,i2) * hidden_output_fj(d2,j);
    else
        delta_weight_ij(i2,j) = (learningRate_ouput * delta_output_i (d2,i2) * hidden_output_fj(d2,j)) +  (alpha * delta_weight_ij_t_1(i2,j));
    end
            
    %%%%%%%% summation old weights * delta i %%%%%%%%%%%%%%%%%%
            
    sum_delta_i = sum_delta_i + (weight_hidden_output_ij(i2,j) * delta_output_i (d2,i2));
             
    %%%%%%%%%%%% make weight changes to output to hidden layer weights %%%%%%%%%%%
             
              
    weight_hidden_output_ij(i2,j)= weight_hidden_output_ij(i2,j)+ delta_weight_ij(i2,j);
             
                 
 
            
           end
       
            
            
       delta_hidden_j (d2,j) = ((1- hidden_output_fj(d2,j)) * hidden_output_fj(d2,j)) * sum_delta_i;
       %%%%%%%%%%%%% loop inputs for weight change corresponding to hidden
       %%%%%%%%%%%%% neuron- input (jk)
       for k = 1:input_num
       delta_weight_jk (j,k) = learningRate_hidden * delta_hidden_j (d2,j) * TrainingSet(d,k);
      
       
        %%%%%%%%%%%% make weight changes to hidden to input layer weights %%%%%%%%%%%
        weight_input_hidden_jk(j,k)= weight_input_hidden_jk(j,k)+ delta_weight_jk (j,k);
         end
       end
      
       %%%%%%%%%%%%% classification %%%%%%%%%%%%%
       %%%%%%% check if highest value at the ouput label corresponds to the
       %%%%%%% position of 1 at the expected label
       [m1,idx_output]=max(output_label(d2,:));
       [m2,idx_expected]=max(expected_label_train(d,:));
       if idx_output == idx_expected
          Q11=Q11+1;      
       end
       
       %%%%%%%%%%%%%%% confusionn matrix for training %%%%%%%%%%%%%%%
       if e == epoch_num
     conf_matrix_train(idx_output,idx_expected)= conf_matrix_train(idx_output,idx_expected)+1;
       end
       
  d2=d2+1;
 
end
  %%%%%%%%%%% hit rate %%%%%%%%%%%%%%
  hit_rate_temp(e)=Q11/trainSubsetSize;
  %if (epoch_store*11)<=e
  if e==1 || mod(e,10)==1
       hit_rate(epoch_store)=Q11/trainSubsetSize;
       epoch_store = epoch_store+1;
 % end
  end
end
 






%%%%%%% testing %%%%%%%%%%%%%%%
   Q11_test=0;
%%%%%%%% loop for data_points %%%%%%%%%
for d= 1:testSetSize
    
   %%%%%%%%% Forward pass %%%%%%%%%%%%                
   %%%% output layer neuron i %%%% 
  
    for i=1: num_neurons_output
        
   %%%% hidden layer neuron j %%%%
       sum_output=0;
       for j=1:num_neurons_hidden
           sum_output_hidden=0;
          for k= 1:input_num 
             sum_output_hidden= sum_output_hidden+(weight_input_hidden_jk(j,k) * TestSet(d,k));
          end
         hidden_output_fj_test(d,j)=1/(1+exp(-sum_output_hidden));
       
         sum_output= sum_output + (weight_hidden_output_ij(i,j) * hidden_output_fj_test(d,j));
       end
       actual_output_fi_test(d,i)=1/(1+exp(-sum_output));
       
         
    %%%%%%%%%%%% labelling %%%%%%%%%%%%%%%%%%
     if actual_output_fi_test(d,i) >=thresh_high        %%%%% set threshold; anything greater than 0.75 is considered as high                     
         output_label_test(d,i)=1;
         
     elseif actual_output_fi_test(d,i) <=thresh_low    %%%%%%% anything lesser -0.75 is considered as low (0)
         output_label_test(d,i)=0;
         
     else
         output_label_test(d,i)=actual_output_fi_test(d,i);
     end
    end
    
     %%%%%%%%%%%%% classification %%%%%%%%%%%%%
       %%%%%%% check if highest value at the ouput label corresponds to the
       %%%%%%% position of 1 at the expected label
       [m1,idx_output_test]=max(output_label_test(d,:));
       [m2,idx_expected_test]=max(expected_label_test(d,:));
       
       if idx_output_test == idx_expected_test
          Q11_test = Q11_test + 1;
       end
       
       %%%%%%%%%%%%% confusion matrix %%%%%%%%%%%%%%%
        conf_matrix_test(idx_output_test,idx_expected_test)= conf_matrix_test(idx_output_test,idx_expected_test)+1;
end
hit_rate_test = Q11_test/testSetSize;

figure();
plot(1:10:epoch_num,hit_rate);
xlabel('epoch number');
ylabel('hit rate');
title('training plot');
legend('hit rate');

hold on
plot(1:10:epoch_num,1-hit_rate,'r');

legend('1-hit rate');



