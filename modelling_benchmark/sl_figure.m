% clear
% load plot.mat
fclose('all');
c=ConstantObj();

time_step = length(bcof);
time_day  = [bcof.tout]/3600/24;%second to day
time_nod_day= arrayfun(@(y) y.tout,nod) * c.dayPsec;

x_matrix = reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);%inp.nn2 is number of nodes in y direction 
y_matrix = reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);
x_ele_matrix = reshape(ele(1).terms{xele_idx},[inp.nn1-1,inp.nn2-1]);
y_ele_matrix = reshape(ele(1).terms{yele_idx},[inp.nn1-1,inp.nn2-1]);

%% evaporation data
evapo_kgs = zeros(inp.nn2,time_step);
for i=1:inp.nn2
    
if i<inp.nn2   
    area1(1:i)    = (x_matrix(1,i+1)-x_matrix(1,i)); %evaporation area 
else
    area1(1:i)    = (x_matrix(1,i)-x_matrix(1,i-1)); %the right end node
end 

    evapo_kgs(i,:)  = -arrayfun(@(y) y.qin(i),bcof);
    evapo_mmday     = evapo_kgs/area1(i)*86400; %evaporation rate of every surface node
    
end

evapo_mmday(1,:)    =   evapo_mmday(1,:)*2;
evapo_mmday(end,:)  =   evapo_mmday(end,:)*2;%avoid boundary effect by double the two nodes
total_evapo_mmday(1:time_step)   =  sum (evapo_mmday(:,1:time_step)); %total evaporation rate

cumulative_evapo_mm =   zeros(1,time_step); %cumulative evaporation amount
for i=2:time_step
    cumulative_evapo_mm(1) =  total_evapo_mmday(1)*inp.scalt*inp.nbcfpr*c.dayPsec;
    cumulative_evapo_mm(i) =  total_evapo_mmday(i)*inp.scalt*inp.nbcfpr*c.dayPsec + cumulative_evapo_mm(i-1);
end


%% plot control
fig_pos.left   = 0.05;
fig_pos.bottom = 0.7;
fig_pos.length = 0.35;
fig_pos.height = 0.26;

%nt=10;
%a.fig=figure;
a.fs = 10;
a.lw = 2; %line width
a.cz = 8; %the size of the marker
fs   = 2; % sampling frequency
% Creating the movie : quality = 100%, no compression is used and the
% timing is adjusted to the sampling frequency of the time interval
qt = 100;
%A number from 0 through 100. Higher quality numbers result in higher video quality and larger file sizes
a.fig = figure;
set (gcf,'Position',[0,0,1600,1200]);
% set(gcf,'Units','normalized', 'OuterPosition',[0 0 1 1]);  % maximize the plotting figure
mov           =  VideoWriter('linux.avi');% avifile('pvc1.avi','quality',qt,'compression','indeo5','fps',fs);
mov.FrameRate = 5;
mov.Quality=qt;
open(mov);

for nt=1:round(time_step/50):time_step
    %% -------------  sub 1 ET over time  --------------
    a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom,...
          fig_pos.length-0.035,fig_pos.height]);
yyaxis left
    a.plot1=plot(time_day(1:nt),total_evapo_mmday(1:nt),...
             'k-','linewidth',a.lw);hold on
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',10);
    ylabel('Evt (mm/day)','FontSize',a.fs);
    axis([-0.1 time_day(end) -0.5 inf])
yyaxis right
    a.plot1=plot(time_day(1:nt),cumulative_evapo_mm(1:nt),...
             'b--','linewidth',a.lw);hold off
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',10);
    txt=sprintf('Result at day %.2f',nod(nt).tout*c.dayPsec);
    title(txt);
    xlabel('Time (day)','FontSize',a.fs);
    ylabel('Cumulative Evt (mm)','FontSize',a.fs);
    axis([-0.1 time_day(end) 0 inf])

    %title('ead profile');

    %% -------------  sub 2 ET for all surface nodes  --------------
    a.sub2=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom,...
          fig_pos.length-0.035,fig_pos.height]);
    a.plot2=plot(x_matrix(1,:), evapo_mmday(:,nt),...
             'k-','linewidth',a.lw);hold off
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',10);
    xlabel('x','FontSize',a.fs);
    ylabel('Evt (mm/day)','FontSize',a.fs);
    axis([0 1.2 0 inf])
    
    %% -------------  sub 3 velocity vector for each node  --------------
    a.sub3=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom-0.64,...
          fig_pos.length-0.035,fig_pos.height]);
    vx_matrix = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1]);
    vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);
    a.plot3=quiver(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix);hold off
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',10);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Velocity')
    axis([0 1.2 0 0.4])
   
%% -------- contour plot on Saturation ---------
    a.sub4=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-0.32,...
          fig_pos.length,fig_pos.height]);		  
    % write pressure and conc in matrix form.
    s_matrix  = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2]);
    a.plot4=contourf(x_matrix,y_matrix,s_matrix,'EdgeColor','none');hold off

%    scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    color = jet;
    color = flipud(color);
    colormap(gca,color);
	caxis([0 1])
    cbsal = colorbar;
    cbsal.Label.String = 'Saturation (-)';
    get(gca,'xtick');
    set(gca,'fontsize',10);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Saturation (-)')


%% -------- contour plot on concentration ---------
    a.sub5=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom-0.32,...
          fig_pos.length,fig_pos.height]);
    
    % write pressure and conc in matrix form.
    c_matrix = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2]);    
    a.plot5=contourf(x_matrix,y_matrix,c_matrix,'EdgeColor','none');hold off
	
%	scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    color = jet;
    colormap(gca,color);
	caxis([0 0.264])
    cbsal = colorbar;
    cbsal.Label.String = 'Concentration (-)';
    get(gca,'xtick');
    set(gca,'fontsize',10);
    title('Concentration (kg/kg)');
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    %axis([10, 40,9,10])
    
    %% -------- contour plot on temperature ---------
    a.sub6=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-0.64,...
          fig_pos.length,fig_pos.height]);
    
    % write pressure and conc in matrix form.
    temp_matrix = reshape(nod(nt).terms{temp_idx},[inp.nn1,inp.nn2]);    
    a.plot6=contourf(x_matrix,y_matrix,temp_matrix-273.15,'EdgeColor','none');hold off
	
%	scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    color = jet;
    colormap(gca,color);
	caxis([20 40])
    cbsal = colorbar;
    cbsal.Label.String = 'Temperature (°C)';
    get(gca,'xtick');
    set(gca,'fontsize',10);
    title('Temperature (°C)');
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    %axis([10, 40,9,10])
    
    %% ------------- concentration profile at given x --------------
    a.sub7=subplot('position'...
         ,[fig_pos.left+0.77,fig_pos.bottom-0.32,...
          0.15,fig_pos.height]);
    c_profile = c_matrix(:,11:20:51);
    a.plot7=plot(c_profile,y_matrix(:,11:20:51),'-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',10);
    set(gca,'yaxislocation','right');
    xlabel('Concentration (-)','FontSize',a.fs);
    legend('x=0.2 m','x=0.6 m','x=1.0 m','Location','southeast')
    axis([0 0.3 0.2 0.4])
    %% ------------- zoom in velocity --------------
    a.sub7=subplot('position'...
         ,[fig_pos.left+0.77,fig_pos.bottom-0.64,...
          0.15,fig_pos.height]);
    vx_matrix = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1]);
    vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);
    a.plot7=quiver(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix);hold off
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',10);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Velocity')
    axis([1. 1.2 0.2 0.35])
    
    
    
	F = getframe(gcf); % save the current figure
    writeVideo(mov,F);% add it as the next frame of the movie

end

close(mov);
close(a.fig);
%figure
% [csal,hcsal] = contourf(x_matrix,y_matrix,p_mtx);
    % hold on
    % set(hcsal,'EdgeColor','none');
    % color = jet;
    % colormap(gca,color);
    % cbsal = colorbar;
    % cbsal.Label.String = 'Pressure (-)';
% scatter (a1(1,(1:f3(1)),1),a1(2,(1:f3(1)),1),2,'filled','w')
% savefig('pressure_contour.fig')						 

figure
[csal,hcsal] = contourf(x_matrix,y_matrix,c_matrix);
    hold on
    set(hcsal,'EdgeColor','none');
    color = jet;
    colormap(gca,color);
    cbsal = colorbar;
	caxis([0 0.246])
    cbsal.Label.String = 'Concentration (-)';
scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
savefig('concentration_contour.fig')						 

figure
[csal,hcsal] = contourf(x_matrix,y_matrix,s_matrix);
    hold on
    set(hcsal,'EdgeColor','none');
    color = jet;
    color = flipud(color);
    colormap(gca,color);
    cbsal = colorbar;
	caxis([0 1])
    cbsal.Label.String = 'Saturation (-)';
scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
savefig('saturation_contour.fig')	
