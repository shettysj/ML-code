clear all
data = xlsread('C:\Users\sanjana\Documents\intelligent systems\HW2\HW2\data.xlsx');
%data=[5,3,1;2,4,0;4,1,1];
data1(:,1)=data(:,1);
data1(:,2)=data(:,2);
data2(:,1)=data(:,3);
j=1; i=1;
eucdist_sq=[];
temp=zeros(2,300);
data3=[];
k=[1,3,5,7,9,11];%,3,5,7,9,11};
negative=zeros(300,6);
positive=zeros(300,6);
final_result=[];
hit_rate=[];
%{
disp (data1(1,1));
disp (data1(1,2));
disp(data2(1));
%}
for j=1:300
for i=1:300
    if j~=i
        eucdist_sq(i,j)=((data1(i,1)-data1(j,1))^2)+((data1(i,2)-data1(j,2))^2);
     
    else
    
        eucdist_sq(i,j)=0; 
    end;
end;
end;



for i=1:length(k)        %corresponds to the k value 
   
    
    for c=1:300           %positive and negatives for each value; c corresponds to the value
         tes=k(i);
         j=1;
        data2(:,2)=eucdist_sq(:,c);
        data3=sortrows(data2,2);
    while(tes~=0)
    if(data3(j,2)~=0)
       if(data3(j,1)==1)
        positive(c,i)=positive(c,i)+1;
       else
           negative(c,i)=negative(c,i)+1;
       end
       tes=tes-1;
    end
    
    j=j+1;
    
     
    
    end
    
    
    end
   % hit_rate(i)= (sum(positive) +sum(negative))/
end
for i=1:length(k)
    false_negative=0;
    false_positive=0;
    true_positive=0;
    true_negative=0;
    for j=1:300
        if positive(j,i)>=negative(j,i)
            final_result(j,i)=1;
             if data(j,3)==final_result(j,i)
                   true_positive=true_positive+1;
             else
                 false_positive=false_positive+1;
             end
        else
            final_result(j,i)=0;
            if data(j,3)==final_result(j,i)
                   true_negative=true_negative+1;
             else
                 false_negative=false_negative+1;
            end
        end
    end
           hit_rate(i)=(true_positive+true_negative)/(false_negative+false_positive+true_positive+true_negative);

end

figure();
scatter(hit_rate,k);
xlabel('hit rate');
ylabel('k');


            



