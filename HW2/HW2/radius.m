%ta=[5,3,1;2,4,0;2,4,1];
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

negative_rad=zeros(300);
positive_rad=zeros(300);

final_result=[];
hit_rate=[];
final_result_rad=[];
hit_rate_rad=[];


%%%%%%%%%%%%%%radius%%%%%%%%%%%%%%%%%%%%
            
Y_mean=mean(mean(eucdist_sq));
Y_min=min(min(eucdist_sq(eucdist_sq>0)));
Y_max=max(max(eucdist_sq));

data2_neighbor=zeros();
data3=zeros();
i=1;

radius=0: 0.05: 1;
% 
% for inc=radius
%     if(eucdist_sq<=inc)
%         neighbour(end+1)=eudist_sq;
%     end
% end
c=1;
for inc=radius
 
 

   
    
    for c=1:300           %positive and negatives for each value; c corresponds to the value
         
      j=1;                %j corresponds to euclidian distance neigbor number
       neighbor=0;
     data2(:,2)=eucdist_sq(:,c);
    
     data3=sortrows(data2,2);
     while(data3(j,2)<=inc)
          if(data3(j,2)~=0)
              neighbor=neighbor+1;
             if(data3(j,1)==1)
                  positive_rad(c,neighbor)=positive_rad(c,neighbor)+1;
             else
                  negative_rad(c,neighbor)=negative_rad(c,neighbor)+1;
             end
          end
    
    j=j+1;
    
     
    
    end
   if neighbor==0
       neighbor=1;
        positive_rad(c,neighbor)=1;
        
        negative_rad(c,neighbor)=0;
   end
        
    
end
   % hit_rate(i)= (sum(positive) +sum(negative))/


end  
  for i=1:21
    Q10=0;
    Q01=0;
    Q11=0;
    Q00=0;
    for j=1:300
        if positive_rad(j,i)>=negative_rad(j,i)
            final_result_rad(j,i)=1;
             if data(j,3)==final_result_rad(j,i)
                   Q11=Q11+1;
             else
                 Q01=Q01+1;
             end
        else
            final_result_rad(j,i)=0;
            if data(j,3)==final_result_rad(j,i)
                   Q00=Q00+1;
             else
                 Q10=Q10+1;
            end
        end
    end
           hit_rate_rad(i)=(Q11+Q00)/(Q10+Q01+Q11+Q00);

end  
% finalHit(end+1)=


 figure();
scatter(hit_rate_rad,neighbor);
xlabel('hit rate');
ylabel('neighbor');