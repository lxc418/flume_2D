% clear
% load plot.mat
% fclose('all');
c=ConstantObj();

time_step = length(et);
time_day  = [bcof.tout]/3600/24;%second to day
time_nod_day = arrayfun(@(y) y.tout,nod) * c.dayPsec;
water_table  = inp.pbc/9800;

x_matrix = reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);%inp.nn2 is number of nodes in x direction 
y_matrix = reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);
x_ele_matrix = reshape(ele(1).terms{xele_idx},[inp.nn1-1,inp.nn2-1]);
y_ele_matrix = reshape(ele(1).terms{yele_idx},[inp.nn1-1,inp.nn2-1]);
area1_m2    = (x_matrix(1,2)-x_matrix(1,1))*inp.z(1);

%% evaporation data (from bcof without the vapor contribution)
% evapo_kgs = zeros(time_step,inp.nn2);
% for i=1:inp.nn2
    
% if i<inp.nn2   
    % area1_m2(1:i)    = (x_matrix(1,i+1)-x_matrix(1,i))*inp.z(1); %evaporation area 
% else
    % area1_m2(1:i)    = (x_matrix(1,i)-x_matrix(1,i-1))*inp.z(1); %the right end node
% end 

    % evapo_kgs(i,:)  = -arrayfun(@(y) y.qin(i),bcof);
    % evapo_mmday     = evapo_kgs/area1_m2(i)*86400; %evaporation rate of every surface node
    
% end

%% spatial_related plot
%solute and solution inflow rates of bottom (from bcop)
for i= 1:inp.nn2
    solute_kgs(i,:)  = arrayfun(@(y) y.qpu(i),bcop);
    
end
	solute_gday= solute_kgs'.*c.kg2g*c.secPday;

for i= 1:inp.nn2
    solution_kgs(i,:)  = arrayfun(@(y) y.qpl(i),bcop);
    
end
	solution_gday= solution_kgs'.*c.kg2g*c.secPday;

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
day_output = [1,5,10];
timestep_output = round(day_output*c.secPday/inp.nprint/inp.scalt);

for nt = timestep_output
    s_matrix  = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2]);	 %write in matrix form.
    s_surface_matrix = s_matrix(inp.nn1,:);	   
	c_matrix = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2]);
	c_surface_matrix = c_matrix(inp.nn1,:);
    vx_matrix = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1]);
    vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);
%% -------------  sub 1 sat & evt  ---------------------
    a.sub1=subplot('position'...
         ,[fig_pos.left+0.35,fig_pos.bottom+fig_pos.height/2,...
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
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
    xlabel('x','FontSize',a.fs);
	set(gca,'YColor',[0 0.4470 0.7410]);
	yticks([5,10])
    ylabel('Evt (mm/day)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 10])
		  
%% -------------  sub 2 concentration & solid salt  ---------------------
    a.sub2=subplot('position'...
         ,[fig_pos.left+0.35,fig_pos.bottom,...
          fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left	  
    a.plot2=plot(x_matrix(1,:), c_surface_matrix,...
             '-','linewidth',a.lw,'color',[0.4940 0.1840 0.5560]);hold off
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
	axis([0 x_matrix(1,end) 0 0.3])
    yticks([0,0.1,0.2,0.3])
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
	set(gca,'YColor',[0.4940 0.1840 0.5560]);
	ylabel('concentration (-)','FontSize',a.fs);
yyaxis right
    solidmass_matrix_kg = reshape(nod(nt).terms{sm_idx},[inp.nn1,inp.nn2]);
    solidmass_surface_kg(1:inp.nn2) = solidmass_matrix_kg(inp.nn1,:);
    solidmass_thickness_mm(1:inp.nn2) = solidmass_surface_kg./c.density_solid_nacl_kgPm3./area1_m2*c.m2mm;
    a.plot2 = plot(x_matrix(1,:), solidmass_thickness_mm(1:inp.nn2),...
             'r-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs,'XTickLabel',{[]});
    ylabel('Solid salt (mm)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 0.12])
    yticks([0,0.04,0.08,0.12])
	grid on
	grid minor
	ax = gca;
	ax.YAxis(2).Color = 'r';	
	ax.XAxis.Visible = 'off';	
   
%% -------- sub 3 concentration contour and streamline ---------
    a.sub3=subplot('position'...
         ,[fig_pos.left+0.35,fig_pos.bottom-fig_pos.height,...
          fig_pos.length-0.015,fig_pos.height]);
    a.plot3=contourf(x_matrix,y_matrix,c_matrix,'EdgeColor','none','LevelStep',0.02);hold on

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
    % a.plot4=quiver(x_ele_matrix,y_ele_matrix,qvx_plot_mtx,qvy_plot_mtx,'k-'); %plot vapor flow vector
	% hold off
    color = jet;
    colormap(gca,color);
	caxis([0.035 0.264])
    cbsal = colorbar;
    cbsal.Label.String = 'Concentration (-)';
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
	%plot stream line
	startx    = x_ele_matrix(1,round(inp.nn2/6):round(inp.nn2/6):inp.nn2-round(inp.nn2/6));
	starty    = y_ele_matrix(2,round(inp.nn2/6):round(inp.nn2/6):inp.nn2-round(inp.nn2/6));%cannot plot when choose the bottom elements y as starty
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
    xlabel('x (m)','FontSize',a.fs);
    ylabel('z (m)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 y_matrix(end,1)+0.01])	
	

	figure_name=sprintf('day_%.2f.fig',nod(nt).tout*c.dayPsec);
	saveas(a.fig,figure_name)
end


