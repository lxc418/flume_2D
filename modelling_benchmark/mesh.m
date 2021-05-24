clear
%% input

x0_1 = 0;
x0_2 = 1.2;
nx_1 = 60;


y(1,1)=0;
y(1,2)=0;

ny(1)=10;

y(2,1)=0.2;
y(2,2)=0.1;

ny(2)=20;

y(3,1)=0.3;
y(3,2)=0.2;

ny(3)=30;

y(4,1)=0.4;
y(4,2)=0.3;

%% process
x_array = x0_1:(x0_2-x0_1)/nx_1:x0_2;
ny_sum=sum(ny);
y_keyline = size (y,1); %control lines in height

for i=1:y_keyline

	y_key_interval=(y(i,2)-y(i,1))/nx_1;

	if y_key_interval==0
		y_array(i,1:nx_1+1) = y(i,1);
	else
		y_array(i,1:nx_1+1) = y(i,1):y_key_interval:y(i,2);
	end

end

for i=1:nx_1+1

        x_nod_mtx ( (i-1)*(ny_sum+1)+1: i*(ny_sum+1) )= x_array(i);

	for j=1:y_keyline-1
        
        y_interval (j,i) = (y_array(j+1,i)-y_array(j,i))./ny(j);      
        y_location (1)   = 0;
        y_location (j+1) = sum(ny(1:j));
                
        y_nod_mtx ( ny_sum*(i-1)+i + y_location(j) :  ny_sum*(i-1)+i + y_location (j+1) )= y_array(j,i): y_interval(j,i) : y_array(j+1,i) ; 

     end
    
	
end