ney     =sum(ney_section);
nx      =nex+1;
ny      =ney+1;

if BLOCK == 1
    extra_nodes = 0;
else
    node_end = length (x_nod_mtx);
    extra_nodes = node_end-(ney+1);
end

x_array = x1:(x2-x1)/nex:x2;
y_keyline = size (y,1); %control lines in height


for i=1:y_keyline

	y_key_interval=(y(i,2)-y(i,1))/nex;

	if y_key_interval==0
		y_array(i,1:nx) = y(i,1);
	else
		y_array(i,1:nx) = y(i,1):y_key_interval:y(i,2);
	end

end

for i=1:nx

        x_nod_mtx ( (i-1)*ny+1 + extra_nodes : i*ny + extra_nodes )= x_array(i);

	for j=1:y_keyline-1
        
        y_interval (j,i) = (y_array(j+1,i)-y_array(j,i))./ney_section(j);      
        y_location (1)   = 0;
        y_location (j+1) = sum(ney_section(1:j));
                
        y_nod_mtx ( ney*(i-1)+i + y_location(j)+ extra_nodes : ney*(i-1)+i + y_location (j+1)+ extra_nodes )=...
            y_array(j,i): y_interval(j,i) : y_array(j+1,i) ; 
    end
    
	
end