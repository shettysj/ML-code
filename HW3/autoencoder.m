clear all;
 
count=1;
input_num=784;
data_points=5000;
num_neurons_hidden=100;
num_neurons_output=784;   %%%% same as number of images %%%%
trainSubsetSize=100;     %%%%% number of data points to be used in each epoch%%%
 
testSetSize = 1000;
trainSet_num = data_points-testSetSize;
 
thresh_high=1;%0.75;
thresh_low=0;%0.25;
 
label_d2=input_num+1;
 
weight_min=-1;
weight_max=1;
weight_range=weight_max-weight_min;
 
learningRate_ouput=0.1;%06;
learningRate_hidden=0.1;%06;
 
alpha=0.6;  %%%%%%%%%% momentum %%%%%%%%
delta_weight_ij=zeros();
delta_weight_ij_t_1=zeros();
 
sum_output_hidden=0;
sum_output=0;
 
epoch_num=400;
epoch_store=1;
 
sum_error_square=0;
sum_error_square_test=0;
sum_error_square_train=0;
loss_func_e=zeros();
loss_func_q=zeros();
counter=1; 
sum_loss_func_train=0;
loss_func_label_train= zeros(1,10);
loss_func_label_test = zeros(1,10);
 
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
sum_loss_func_e=0;
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
       
         
    end
    %%%%%%%%%%%%%% start loops for backpropogation %%%%%%%%%%%
    
     %%%%% below part for momentum %%%%%%%
    delta_weight_ij_t_1= delta_weight_ij;
    %%%%%% above part for momentum %%%%%
    
       for j=1:num_neurons_hidden  %%%%% loop for change in weights %%%%%%%%%%
           
           
           sum_error_square=0;
           sum_delta_i=0;
        
           
           for i2 = 1:num_neurons_output
               
                %%%%%%% find error by comparing expected output which is the input with acquired output %%
             error_output(d2,i2) = TrainingSet(d,i2)-actual_output_fi(d2,i2);
             
             
           
             
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
       %%%%%%% check if the input va;ues correspond to the ouptut
        for i2 = 1:num_neurons_output
            %%%%%%%%%%% calculate sum of square of errors
             sum_error_square = (sum_error_square + ((error_output(d2,i2))^2));
       if actual_output_fi(d2,i2)== TrainingSet(d,i2)
          Q11=Q11+1;      
       end
        end
        
        loss_func_q(e,d2)= 0.5 * sum_error_square;    %%%%%% loss function of each data point %%%%%%
         %%%%%% loss function of each epoch %%%%%%
  sum_loss_func_e = sum_loss_func_e + loss_func_q(e,d2);
  
  d2=d2+1;
 
 
end
loss_func_e(e)=sum_loss_func_e;
  %%%%%%%%%%% hit rate %%%%%%%%%%%%%%
  hit_rate_temp(e)=Q11/trainSubsetSize;
  %if (epoch_store*11)<=e
  if e==1 || mod(e,11)==0
       loss_func_tenth(epoch_store)=loss_func_e(e);
       epoch_store = epoch_store+1;
 % end
  end
end

for d= 1:1000
    
    
     
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
         hidden_output_fj(d,j)=1/(1+exp(-sum_output_hidden));
       
         sum_output= sum_output + (weight_hidden_output_ij(i,j) * hidden_output_fj(d,j));
       end
       actual_output_fi(d,i)=1/(1+exp(-sum_output));
       
         
    end
    %%%%%%%%%%%%%% start loops for backpropogation %%%%%%%%%%%
    
     %%%%% below part for momentum %%%%%%%
    delta_weight_ij_t_1= delta_weight_ij;
    %%%%%% above part for momentum %%%%%
    
       for j=1:num_neurons_hidden  %%%%% loop for change in weights %%%%%%%%%%
           
           
           sum_error_square_train=0;
           sum_delta_i=0;
           
           
        
           
           for i2 = 1:num_neurons_output
               
                %%%%%%% find error by comparing expected output which is the input with acquired output %%
             error_output(d,i2) = TrainingSet(d,i2)-actual_output_fi(d,i2);
             
             
           
             
       %%%%%%%%%% find delta i %%%%%%%%%%%%%%%%%
             delta_output_i (d,i2) = (1-actual_output_fi(d,i2))*actual_output_fi(d,i2)*(error_output(d,i2));
    
   
    %%%%%%%%%%%%%%% find weight change required for weights between output
    %%%%%%%%%%%%%%% layer and hidden layer %%%%%%%%%%%
    
   
    %%% without momentum: delta_weight_ij(i2,j) = learningRate_ouput * delta_output_i (d2,i2) * hidden_output_fj(d2,j);
    %%%%% with momentum: %%%%%%%%%%%%%%%%%%
    if d==1      
        delta_weight_ij(i2,j) = learningRate_ouput * delta_output_i (d,i2) * hidden_output_fj(d,j);
    else
        delta_weight_ij(i2,j) = (learningRate_ouput * delta_output_i (d,i2) * hidden_output_fj(d,j)) +  (alpha * delta_weight_ij_t_1(i2,j));
    end
            
    %%%%%%%% summation old weights * delta i %%%%%%%%%%%%%%%%%%
            
    sum_delta_i = sum_delta_i + (weight_hidden_output_ij(i2,j) * delta_output_i (d,i2));
             
    %%%%%%%%%%%% make weight changes to output to hidden layer weights %%%%%%%%%%%
             
              
    weight_hidden_output_ij(i2,j)= weight_hidden_output_ij(i2,j)+ delta_weight_ij(i2,j);
             
                 
 
            
           end
       
            
            
       delta_hidden_j (d,j) = ((1- hidden_output_fj(d,j)) * hidden_output_fj(d,j)) * sum_delta_i;
       %%%%%%%%%%%%% loop inputs for weight change corresponding to hidden
       %%%%%%%%%%%%% neuron- input (jk)
       for k = 1:input_num
       delta_weight_jk (j,k) = learningRate_hidden * delta_hidden_j (d,j) * TrainingSet(d,k);
      
       
        %%%%%%%%%%%% make weight changes to hidden to input layer weights %%%%%%%%%%%
        weight_input_hidden_jk(j,k)= weight_input_hidden_jk(j,k)+ delta_weight_jk (j,k);
         end
       end
      
        %%%%%%%%%%%%% classification %%%%%%%%%%%%%
       %%%%%%% check if the input va;ues correspond to the ouptut
        for i2 = 1:num_neurons_output
            %%%%%%%%%%% calculate sum of square of errors
             sum_error_square_train = (sum_error_square_train + ((error_output(d,i2))^2));
       if actual_output_fi(d,i2)== TrainingSet(d,i2)
          Q11=Q11+1;      
       end
        end
        loss_func_q_train(d)= 0.5 * sum_error_square_train;    %%%%%% loss function of each data point %%%%%%
     label_loss_train = TrainingSet(d,785) + 1;
     loss_func_label_train(label_loss_train)= loss_func_label_train(label_loss_train) + loss_func_q_train(d);
      
 
  
  d2=d2+1;
 
 
end
loss_func_train=sum(loss_func_q_train);


%%%%%%% testing %%%%%%%%%%%%%%%
   Q11_test=0;
%%%%%%%% loop for data_points %%%%%%%%%
for d= 1:1000 %testSetSize
   
    loss_func_q_test(d)= 0.5 * sum_error_square_test;
    sum_error_square_test=0;
    label_loss_test = TestSet(d,785) + 1;
    loss_func_label_test(label_loss_test)= loss_func_label_test(label_loss_test) + loss_func_q_test(d);
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
       
         %%%%%%% find error by comparing expected output which is the input with acquired output %%
             error_output_test(d,i) = TestSet(d,i)-actual_output_fi_test(d,i);
             
             %%%%%%%%%%% calculate sum of square of errors
             sum_error_square_test = (sum_error_square_test + ((error_output_test(d,i))^2));  
    end
    
     %%%%%%%%%%%%% hit rate %%%%%%%%%%%%%%%%%%
       
      for i2 = 1:num_neurons_output
       if actual_output_fi_test(d,i)== TestSet(d,i)
          Q11_test = Q11_test + 1;
       end
      end
       
    
end
hit_rate_test = Q11_test/testSetSize;
loss_func_e_test = sum(loss_func_q_test);

   
%%%%%%%%%%%%%%%%% gray scale image %%%%%%%%%%%%%%%%%%%%%
figure();
while(counter<=100)
    total_pts=1;
    for g_r=1:28
        for g_c=1:28
            grey_image(g_r,g_c)= weight_hidden_output_ij(total_pts,counter);
            total_pts=total_pts+1;
        end
    end
    ABCD=mat2gray(grey_image);
   subplot(10,10,counter), imshow(ABCD); 
   counter=counter+1;
end 

figure();
c = categorical({'training','testing'});
title('loss function of training set vs test set');
bar(c,[loss_func_train , loss_func_e_test]);

trans_label_train=transpose(loss_func_label_train); 
trans_label_test=transpose(loss_func_label_test);

figure();
bar([trans_label_train(1:10,1), trans_label_test(1:10,1)]);
xlabel('labels')
ylabel('loss function')
title('loss function of training set vs test set for corresponding labels');
legend('training','testing');

figure();
plot(1:10:epoch_num,loss_func_tenth);
xlabel('epoch number');
ylabel('loss_func');
title('training plot');
 
