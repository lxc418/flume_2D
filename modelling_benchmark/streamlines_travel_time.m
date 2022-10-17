%%% CALCULATE TRAVEL TIME OF THE STREAMLINES, MUST RUN BASED ON THE NOD AND
%%% ELE MTX AND AFTER GENERATE STREAMLINES

%%%%%%% TRAVEL TIME OF STREAMLINES
for l = 1:smline_num
    x_line = smline(l).XData;
    y_line = smline(l).YData;
    x_line = rmmissing(x_line);
    y_line = rmmissing(y_line);
    
    
    for k = 1:length(x_line)
    [xminValue, xclosestIndex] = min(abs(x_line(k) - x_ele_matrix(1,:)));
    [yminValue, yclosestIndex] = min(abs(y_line(k) - y_ele_matrix(:,xclosestIndex)));
    searchU(k) = vx_matrix(yclosestIndex,xclosestIndex);
    searchV(k) = vy_matrix(yclosestIndex,xclosestIndex);
    x_coord(k) = x_line(k);
    y_coord(k) = y_line(k);
    end
%     searchUfileName = sprintf('U_of_line_%d.mat',l);
%     searchVfileName = sprintf('V_of_line_%d.mat',l);
%     save(searchUfileName,'searchU');
%     save(searchVfileName,'searchV');
    
    for p = 1:length(x_line)-1
        x_dis(p) = x_line(p)-x_line(p+1);
        y_dis(p) = y_line(p)-y_line(p+1);
        travel_dis(p) = sqrt(x_dis(p)^2+y_dis(p)^2);
        speed(p) = sqrt(searchU(p)^2+searchV(p)^2);
    end
%     disfileName = sprintf('dis_of_line_%d.mat',l);
%     speedfileName = sprintf('speed_of_line_%d.mat',l);
%     save(disfileName,'travel_dis');
%     save(speedfileName,'speed');
    
    for q = 1:length(speed)
        time(q) = travel_dis(q)./speed(q);
    end
    total_time(l) = sum(time)/3600/24;
    
%     timeArrayfileName = sprintf('timearray_of_line_%d.mat',l);
%     totaltimefileName = sprintf('totaltime_of_line_%d.mat',l);
%     save(timeArrayfileName,'time');
%     save(totaltimefileName,'total_time');
    clear searchU searchV x_coord y_coord x_dis y_dis travel_dis speed time
end 
writematrix(total_time,'Streamline_traveltime_days.txt','Delimiter','tab');
