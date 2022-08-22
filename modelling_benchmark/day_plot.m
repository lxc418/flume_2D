%% this file is used for plot transient results at a given time(day) 
% clear
% load plot.mat
% fclose('all');
c=ConstantObj();
day_output = [1,6,11,30];%set days need to be plotted			 																  					 																  		

time_step = length(et);
time_day  = [bcof.tout]/3600/24;%second to day
time_nod_day = arrayfun(@(y) y.tout,nod) * c.dayPsec;
							
x_matrix = reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);%inp.nn2 is number of nodes in x direction 
y_matrix = reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);
x_ele_matrix = reshape(ele(1).terms{xele_idx},[inp.nn1-1,inp.nn2-1]);
y_ele_matrix = reshape(ele(1).terms{yele_idx},[inp.nn1-1,inp.nn2-1]);
x_area_m2    = reshape(mesh.area_xz,[inp.nn2,inp.nn1]);
x_area_m2    = x_area_m2';
y_area_m2    = reshape(mesh.area_yz,[inp.nn2,inp.nn1]);
y_area_m2	 = y_area_m2';

%% boundary condition 							 
if abs(inp.ipbc(1)-inp.ipbc(2))==1
boundary_condition = 2;%1/bottom boundary   2/right side boundary 
flux_plot_level_y  = y_ele_matrix(2,end);%set the level for streamline plotting (the bottom ele cannot used to plot streamline)
else
boundary_condition = 1;
end		 
%boundary area
top_boundary_x_area_m2   = x_area_m2(1,:);
if boundary_condition == 1
	bottom_boundary_x_area_m2= x_area_m2(1,:);
elseif boundary_condition == 2
	right_boundary_y_area_m2 = y_area_m2(inp.nn1-length(inp.ipbc)+1:end,end);%from top to bottom
	right_boundary_y_plot    = flip(y_matrix(1:length(inp.ipbc),end));%from top to bottom
end

%% spatial_related plot
%salt and solution flux of side (from bcop) 
if boundary_condition == 1
	for i= 1:length(inp.ipbc) 
		salt_kgs(:,i)      = arrayfun(@(y) y.qpu(i),bcop);
		salt_kgsm2(:,i)    = salt_kgs(:,i)./bottom_boundary_x_area_m2(i);% for bottom boundary condition
	end
	for i= 1:length(inp.ipbc)
    solution_kgs(:,i)    = arrayfun(@(y) y.qpl(i),bcop);
    solution_ms(:,i)     = solution_kgs(:,i)/(c.rhow_pure_water+700*0.035)/bottom_boundary_x_area_m2(i);
	end
elseif boundary_condition == 2
	for i= 1:length(inp.ipbc) %from top to bottom
		salt_kgs(:,i)      = arrayfun(@(y) y.qpu(i),bcop);
		salt_kgsm2(:,i)    = salt_kgs(:,i)./right_boundary_y_area_m2(i);
	end

	for i= 1:length(inp.ipbc)
    solution_kgs(:,i)    = arrayfun(@(y) y.qpl(i),bcop);
    solution_ms(:,i)     = solution_kgs(:,i)/(c.rhow_pure_water+700*0.035)/right_boundary_y_area_m2(i);
	end
end	

%find nodes at a given level (e.g. groundwater table) 
if boundary_condition == 1
		x_flux_level = x_ele_matrix(1,:);
		y_flux_level = y_ele_matrix(2,:);
	elseif boundary_condition == 2
		for i = 1:inp.nn2-1
			[yminValue(i), yclosestIndex_plot(i)]=min(abs(flux_plot_level_y - y_ele_matrix(:,i))); %get the node closest to water table by comparing y 
			x_flux_level(i) = x_ele_matrix(yclosestIndex_plot(i),i);
			y_flux_level(i) = y_ele_matrix(yclosestIndex_plot(i),i);
		end
end

%% temporal_related plot
%evaporation rate from QV.DAT with the vapor flow accounted
for i=2:time_step %The first time step starts from the initial condition of not reaching equilibrium, so start from step 2

	evapo_mmday(1,:)  		 =  reshape(et(1).terms{et_idx},[1,inp.nn2])*c.ms2mmday;
	total_evapo_mmday(1)     =  sum (evapo_mmday(1,:))./inp.nn2; %the evp rate from the whole surface
	evapo_mmday(i,:)  		 =  reshape(et(i).terms{et_idx},[1,inp.nn2])*c.ms2mmday;
	total_evapo_mmday(i)     =  sum (evapo_mmday(i,:))./inp.nn2;%cannot use if the mesh along x is not even
	
    cumulative_evapo_mm(1)   =  total_evapo_mmday(1)*inp.scalt*inp.nbcfpr*c.dayPsec;
    cumulative_evapo_mm(i)   =  total_evapo_mmday(i)*inp.scalt*inp.nbcfpr*c.dayPsec + cumulative_evapo_mm(i-1);

end

%% plot control
fig_pos.left   = 0.05;
fig_pos.bottom = 0.7;
fig_pos.length = 0.35;
fig_pos.height = 0.26;

a.fs = 10;
a.lw = 1.75; %line width
a.cz = 8; %the size of the marker
fs   = 2; % sampling frequency
% Creating the movie : quality = 100%, no compression is used and the
% timing is adjusted to the sampling frequency of the time interval
qt = 100;
%A number from 0 through 100. Higher quality numbers result in higher video quality and larger file sizes
a.fig = figure;
set (gcf,'Position',[0,0,1920,1080]); %resolution 1080p
% set(gcf,'Units','normalized', 'OuterPosition',[0 0 1 1]);  % maximize the plotting figure

%calculate time_step for output
timestep_output = round(day_output*c.secPday/inp.nprint/inp.scalt)+1;
% timestep_output = time_step; %just plot the end of simulation

for nt = timestep_output
    s_matrix  = reshape(nod(nt+1).terms{s_idx},[inp.nn1,inp.nn2]);	 %write in matrix form.
    s_surface_matrix = s_matrix(inp.nn1,:);	   
	c_matrix = reshape(nod(nt+1).terms{c_idx},[inp.nn1,inp.nn2])*1000;%change unit to ppt
	c_surface_matrix = c_matrix(inp.nn1,:);
    vx_matrix = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1]);
    vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);
%% -------------  sub 1 sat & evt  ---------------------
    a.sub1=subplot('position'...
         ,[fig_pos.left+0.35,fig_pos.bottom+fig_pos.height/2+0.005,...
          fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left
    a.plot1=plot(x_matrix(1,:), s_surface_matrix,...
             '-','linewidth',a.lw,'color',[0.4660 0.6740 0.1880]);hold off		  			 
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
	axis([0 x_matrix(1,end) 0 1])
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
	set(gca,'YColor',[0.4660 0.6740 0.1880]);
    yticks([0.5,1])
	ylabel('Saturation (-)','FontSize',a.fs);		  
yyaxis right		   
    a.plot1=plot(x_matrix(1,:), evapo_mmday(nt,:),...
             '-','linewidth',a.lw,'color',[0 0.4470 0.7410]);hold off
	grid on
	% grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
    % xlabel('x','FontSize',a.fs);
	set(gca,'YColor',[0 0.4470 0.7410]);
	yticks([5,10])
    ylabel('Evt (mm/day)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 10])
		  
%% -------------  sub 2 Salinity & solid salt  ---------------------
    a.sub2=subplot('position'...
         ,[fig_pos.left+0.35,fig_pos.bottom,...
          fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left	  
    a.plot2=plot(x_matrix(1,:), c_surface_matrix,...
             '-','linewidth',a.lw,'color',[0.4940 0.1840 0.5560]);hold off
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
	axis([0 x_matrix(1,end) 0 300])
    yticks([0,100,200,300])
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
	set(gca,'YColor',[0.4940 0.1840 0.5560]);
	ylabel('Salinity (ppt)','FontSize',a.fs);
yyaxis right
    solidmass_matrix_kg = reshape(nod(nt+1).terms{sm_idx},[inp.nn1,inp.nn2]);
    solidmass_surface_kg(1:inp.nn2) = solidmass_matrix_kg(inp.nn1,:);
    solidmass_thickness_mm(1:inp.nn2) = solidmass_surface_kg./c.density_solid_nacl_kgPm3./top_boundary_x_area_m2*c.m2mm;
    a.plot2 = plot(x_matrix(1,:), solidmass_thickness_mm(1:inp.nn2),...
             'r-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
    ylabel('Solid salt (mm)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 0.12])
    yticks([0,0.04,0.08,0.12])
	set(gca,'YColor','r');
	grid on
	% grid minor
   
%% -------- sub 3 concentration contour and streamline ---------
    a.sub3=subplot('position'...
         ,[fig_pos.left+0.35,fig_pos.bottom-0.205,...
          fig_pos.length-0.015,0.2]);
% plot log contour
   a.plot3=contourf(x_matrix,y_matrix,log10(c_matrix+1-inp.ubc(1)*1000),'LevelStep',0.005,'LineStyle','None');hold on
    color = jet;
    colormap(gca,color);
	caxis([0.05 2.35])
	set(gca,'ColorScale','log')%log scale of colormap
	cbsal = colorbar;
	% cbsal.Location = 'southoutside';
    cbsal.Label.String = 'Salinity (ppt)';
	cbsal.TicksMode    = 'manual';
	cbsal.Ticks        = [0.06,0.3,0.78,1.42,2.35];
	cbsal.TickLabels   = [35,36,40,60,260];		
% plot normal contour
    % a.plot3=contourf(x_matrix,y_matrix,c_matrix,'EdgeColor','none');hold on
    % color = jet;
    % colormap(gca,color);
	% caxis([35 264])
    % cbsal = colorbar;
    % cbsal.Label.String = 'Concentration (-)';
%plot vapor flow vector
   % qvx_mtx = qv(nt).qvx; % vector of vapor is acquired from the concentration difference bewteen two nodes,
    % qvy_mtx = qv(nt).qvy; % which is used to calculate the vector in element center.
    % qvx_plot_mtx = zeros(inp.nn1-1,inp.nn2-1);  
    % qvy_plot_mtx = zeros(inp.nn1-1,inp.nn2-1);
    % for i = 1:inp.nn2-1 %along x
        % for j = 1:inp.nn1-1 %along y
            % qvx_plot_mtx(j,i) = (qvx_mtx(j,i)+qvx_mtx(j+1,i))/2;
            % qvy_plot_mtx(j,i) = (qvy_mtx(j,i)+qvy_mtx(j,i+1))/2;
        % end
    % end
    % a.plot3=quiver(x_ele_matrix,y_ele_matrix,qvx_plot_mtx,qvy_plot_mtx,'k-'); %plot vapor flow vector
	% hold on	
    a.plot3   = quiver(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix);
	hold on

    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    % xlabel('x (m)','FontSize',a.fs);
    ylabel('z (m)','FontSize',a.fs);
	%plot stream line
	startx    = x_flux_level(round(inp.nn2/12):round(inp.nn2/12):inp.nn2-round(inp.nn2/12));
	starty    = y_flux_level(round(inp.nn2/12):round(inp.nn2/12):inp.nn2-round(inp.nn2/12));
	% startx    = [x_ele_matrix(round(inp.nn1/12):round(inp.nn1/12):inp.nn1,160);x_ele_matrix(round(inp.nn1/12):round(inp.nn1/12):inp.nn1,200)];
	% starty    = [y_ele_matrix(round(inp.nn1/12):round(inp.nn1/12):inp.nn1,160);y_ele_matrix(round(inp.nn1/12):round(inp.nn1/12):inp.nn1,200)];
	a.plot3   = scatter(startx,starty,'filled','w');												 
	a.smline   = streamline(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix,startx,starty);
	smline_num = length (a.smline);
	domain_x = [x_matrix(1,1),x_matrix(1,end),x_matrix(end,end),x_matrix(end,1),x_matrix(1,1)];%define the area with line 
	domain_y = [y_matrix(1,1),y_matrix(1,end),y_matrix(end,end),y_matrix(end,1),y_matrix(1,1)];
	for i = 1:smline_num %extract the coordinates of streamline 
		x_smline = a.smline(i).XData;
		y_smline = a.smline(i).YData;
		indomain = inpolygon(x_smline,y_smline,domain_x,domain_y); % to determine if the points are inside of the domain polygon
		x_smline(~indomain) = NaN;
		y_smline(~indomain) = NaN;
		% fname = sprintf('Line_%d_x_y.csv',i); %save the data of streamline
		% writematrix([x_line',y_line'],fname);
		idx_firstNan = find(isnan(x_smline),i); % find the first NaN location, then change all value after that to be NaN
		x_smline(idx_firstNan:end) = NaN;
		y_smline(idx_firstNan:end) = NaN;
		smline(i) = plot(x_smline,y_smline,'linewidth',a.lw,'Color',[1 1 1]);
		hold on;
	end
	delete(a.smline);
	hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
	% xlabel('x (m)','FontSize',a.fs);
    ylabel('z (m)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 y_matrix(end,1)+0.01])	
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
	
%calculate stream line travel time 	
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
	txt_name=sprintf('day_%.2f.txt',nod(nt+1).tout*c.dayPsec);
	writematrix(total_time',txt_name,'Delimiter','tab');

%% -------- sub 4 horizontal flux (from .ele for side boundary and from .bcop for bot. boundary) ---------
if boundary_condition == 1 
    a.sub4=subplot('position'...
         ,[fig_pos.left+0.35,fig_pos.bottom-0.34,...
          fig_pos.length-0.05,0.1]);
yyaxis left										   
    a.plot4=plot(x_matrix(1,:), solution_ms(nt,:),...
             '-','color',[0.9290 0.6940 0.1250]	,'linewidth',a.lw);hold off
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
	set(gca,'YColor',[0.9290 0.6940 0.1250]);
    set(gca,'fontsize',a.fs);
    ylabel(['bottom water flux' sprintf('\n') '(m/s)'],'FontSize',a.fs);
    axis([0 x_matrix(1,end) -5e-5 5e-5])
yyaxis right		   		  
    a.plot4=plot(x_matrix(1,:), salt_kgsm2(nt,:),...
             'k-','linewidth',a.lw);hold off	
    get(gca,'xtick');
	set(gca,'YColor','k');
    set(gca,'fontsize',a.fs);
    xlabel('x (m)','FontSize',a.fs);
    ylabel(['bottom salt flux' sprintf('\n') '(kg/s/m^2)'],'FontSize',a.fs);
    axis([0 x_matrix(1,end) -5e-3 5e-3])
	% legend('water','salt','FontSize',a.fs,'Location','northwest' )	
end

if boundary_condition == 2 
%plot flux at a given horizontal level
    a.sub4=subplot('position'...
         ,[fig_pos.left+0.35,fig_pos.bottom-1.5*fig_pos.height,...
          fig_pos.length-0.05,fig_pos.height/2]);
	x_flux_matrix_ms  = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1])*inp.por(1);
	y_flux_matrix_ms  = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1])*inp.por(1);
	for i=1:inp.nn2-1
	x_flux_ms(i) = x_flux_matrix_ms(yclosestIndex_plot(i),i);
	y_flux_ms(i) = y_flux_matrix_ms(yclosestIndex_plot(i),i);
	end
	flux_level_ms = y_flux_ms./abs(y_flux_ms).*(x_flux_ms.^2+y_flux_ms.^2).^0.5;
	
% yyaxis left										   
    a.plot4=plot(x_ele_matrix(1,:), flux_level_ms,...
             '-','color',[0.9290 0.6940 0.1250]	,'linewidth',a.lw);hold off
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('water flux (g/day)','FontSize',a.fs);
    axis([0 x_matrix(1,end) -0.002 0.002])
% yyaxis right		   		  
    % a.plot4=plot(x_matrix(1,:), salt_kgsm2(nt,:),...
             % 'k-','linewidth',a.lw);hold off	
    % get(gca,'xtick');
	% set(gca,'YColor','k');
    % set(gca,'fontsize',a.fs);
    % xlabel('x (m)','FontSize',a.fs);
    % ylabel('bottom salt flux (g/day)','FontSize',a.fs);
    % axis([0 x_matrix(1,end) -5e-3 5e-3])
	% legend('water','salt','FontSize',a.fs,'Location','northwest' )	
	
%plot side boundary flux
    a.sub4=subplot('position'...
         ,[fig_pos.left+0.32+fig_pos.length,fig_pos.bottom-0.9*fig_pos.height,...
          fig_pos.length-0.25,0.75*fig_pos.height]);
% the flux plot shares the same y axis
	ax1=gca;
    ax1_pos = ax1.Position;
    a.plot8=plot( solution_ms(nt,:),right_boundary_y_plot',...
        '-','color',[0.9290 0.6940 0.1250]	,'linewidth',a.lw);hold off
	grid on
	grid minor
	ax1.GridAlpha = 0.4;
	ax1.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
	set(gca,'XColor',[0.9290 0.6940 0.1250]);
    set(gca,'fontsize',a.fs);
    ylabel(['boundary water flux' sprintf('\n') '(m/s)'],'FontSize',a.fs);
    axis([-1e-2 1e-2 0 y_matrix(end,end)])

	ax2 = axes('Position', ax1_pos, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
	a.plot4=plot( salt_kgsm2(nt,:),right_boundary_y_plot',...
             'k-','linewidth',a.lw);hold off
	ax2.XAxisLocation = 'top';
	ax2.YAxisLocation = 'right';
 	ax2.Color = 'none';		
    get(gca,'xtick');
	set(gca,'XColor','k');
    set(gca,'fontsize',a.fs);
	xlabel('x (m)','FontSize',a.fs);
    ylabel(['boundary salt flux' sprintf('\n') '(kg/s/m^2)'],'FontSize',a.fs);
    axis([-0.5 0.5 0 y_matrix(end,end)])
end																

	figure_name=sprintf('day_%.2f.fig',nod(nt+1).tout*c.dayPsec);
	saveas(a.fig,figure_name)

end


