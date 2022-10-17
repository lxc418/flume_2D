% clear
% load plot.mat
% fclose('all');
c=ConstantObj();

time_step = length(et);
time_day  = [bcof.tout]/3600/24;%second to day
time_nod_day = arrayfun(@(y) y.tout,nod) * c.dayPsec;
water_table  = inp.pbc/(c.rhow_pure_water+700*0.035);							
							
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
flux_plot_level_y  = y_ele_matrix(2,end);%set the level for streamline plotting when side boundary applied (e.g. water table)										 
										 %(the bottom ele cannot used to plot streamline)
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
%boundary salt and solution flux(from bcop) 
if boundary_condition == 1 %bottom boundary
	for i= 1:length(inp.ipbc) 
		salt_kgs(:,i)      = arrayfun(@(y) y.qpu(i),bcop);
		salt_kgsm2(:,i)    = salt_kgs(:,i)./bottom_boundary_x_area_m2(i);% for bottom boundary condition
	end
	for i= 1:length(inp.ipbc)
    solution_kgs(:,i)    = arrayfun(@(y) y.qpl(i),bcop);
    solution_ms(:,i)     = solution_kgs(:,i)/(c.rhow_pure_water+700*0.035)/bottom_boundary_x_area_m2(i);
	end
elseif boundary_condition == 2 %right side boundary
	for i= 1:length(inp.ipbc) %from top to bottom
		salt_kgs(:,i)      = arrayfun(@(y) y.qpu(i),bcop);
		salt_kgsm2(:,i)    = salt_kgs(:,i)./right_boundary_y_area_m2(i);
	end
	for i= 1:length(inp.ipbc)
    solution_kgs(:,i)    = arrayfun(@(y) y.qpl(i),bcop);
    solution_ms(:,i)     = solution_kgs(:,i)/(c.rhow_pure_water+700*0.035)/right_boundary_y_area_m2(i);
	end
end	

%find nodes at a given level (e.g. bottom boundary or groundwater table for streamline plot) 
if boundary_condition == 1
		x_flux_level = x_ele_matrix(1,:);
		y_flux_level = y_ele_matrix(2,:); %bottom boundary as the start of stream lines
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
	evapo_mmday(1,(evapo_mmday(1,:)==0)) = max(evapo_mmday(1,:)); %make the surface below the water table has the same evp as full saturated surface 
															% specifically for case with inundation
	total_evapo_mmday(1)     =  sum (evapo_mmday(1,:))./inp.nn2; %the evp rate from the whole surface
	evapo_mmday(i,:)  		 =  reshape(et(i).terms{et_idx},[1,inp.nn2])*c.ms2mmday;
	evapo_mmday(i,(evapo_mmday(i,:)==0)) = max(evapo_mmday(i,:));
	total_evapo_mmday(i)     =  sum (evapo_mmday(i,:))./inp.nn2;%cannot use if the mesh along x is not even
	
    cumulative_evapo_mm(1)   =  total_evapo_mmday(1)*inp.scalt*inp.nbcfpr*c.dayPsec;
    cumulative_evapo_mm(i)   =  total_evapo_mmday(i)*inp.scalt*inp.nbcfpr*c.dayPsec + cumulative_evapo_mm(i-1);

end
% total solution/salt mass and surface solution/salt mass (from nod)
for i=2:time_step
	total_solutionmass_kgm(i)  = sum(nod(i).terms{m_solution_idx})/inp.z(1);%times 100 is because the z direction length is 0.01m
	total_saltmass_kgm(i)      = sum(nod(i).terms{m_solute_idx})/inp.z(1);

    solutionmass_surface_kgm(i) = sum(nod(i).terms{m_solution_idx}(inp.nn1:inp.nn1:inp.nn))/inp.z(1);
    saltmass_surface_kgm(i)     = sum(nod(i).terms{m_solute_idx}(inp.nn1:inp.nn1:inp.nn))/inp.z(1);
	
end
% cumulative_bottom or side solution/salt mass (from bcop)
if boundary_condition == 1
	for i=2:time_step
		boundary_salt_kgsm2(i)  	   = sum(salt_kgs(i,:))./sum(bottom_boundary_x_area_m2); 
		boundary_solution_ms(i)        = sum(solution_kgs(i,:))./(c.rhow_pure_water+700*0.035)./sum(bottom_boundary_x_area_m2);
		cumu_boundary_salt_kgm2(1)	   = 0;
		cumu_boundary_solution_m(1)    = 0;
		cumu_boundary_salt_kgm2(i)	   = boundary_salt_kgsm2(i)*inp.scalt*inp.nbcfpr + cumu_boundary_salt_kgm2(i-1);
		cumu_boundary_solution_m(i)    = boundary_solution_ms(i)*inp.scalt*inp.nbcfpr + cumu_boundary_solution_m(i-1);
	end
end

if boundary_condition == 2
	for i=2:time_step
		boundary_salt_kgsm2(i)  	   = sum(salt_kgs(i,:))./sum(right_boundary_y_area_m2); 
		boundary_solution_ms(i)        = sum(solution_kgs(i,:))./(c.rhow_pure_water+700*0.035)./sum(right_boundary_y_area_m2);
		cumu_boundary_salt_kgm2(1)	   = 0;
		cumu_boundary_solution_m(1)    = 0;
		cumu_boundary_salt_kgm2(i)	   = boundary_salt_kgsm2(i)*inp.scalt*inp.nbcfpr + cumu_boundary_salt_kgm2(i-1);
		cumu_boundary_solution_m(i)    = boundary_solution_ms(i)*inp.scalt*inp.nbcfpr + cumu_boundary_solution_m(i-1);
	end
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

mov           =  VideoWriter('linux.avi');% avifile('pvc1.avi','quality',qt,'compression','indeo5','fps',fs);
mov.FrameRate = 5;
mov.Quality   = qt;
open(mov);

for nt=2:round(time_step/40):time_step														 
    s_matrix  = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2]);	 %write in matrix form.
    s_surface_matrix = s_matrix(inp.nn1,:);	   
	c_matrix = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2])*1000;%change unit to ppt
	c_surface_matrix = c_matrix(inp.nn1,:);
    vx_matrix = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1]);
    vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);																																	   																																																														 
%% -------------  sub 1_1 ET over time  --------------
    a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom+fig_pos.height/2,...
          fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left
    a.plot1=plot(time_day(2:nt),total_evapo_mmday(2:nt),...
             'k-','linewidth',a.lw);hold off
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('avg Evp (mm/day)','FontSize',a.fs);
    axis([0 time_day(end) -0.1 inf])	
yyaxis right
    a.plot1=plot(time_day(2:nt),cumulative_evapo_mm(2:nt),...
             'b--','linewidth',a.lw);hold off
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
	set(gca,'Xticklabel',[])
    txt=sprintf('Result at day %.2f & avg evp is %.2f mm/day',nod(nt).tout*c.dayPsec, cumulative_evapo_mm(nt)/(nod(nt).tout*c.dayPsec));
    title(txt);
    ylabel('Cumulative Evt (mm)','FontSize',a.fs);
    axis([0 time_day(end) 0 inf])
	ax = gca;
	ax.YAxis(1).Color = 'k';
	ax.YAxis(2).Color = 'b';
%% -------------  sub 1_2 surface water/salt mass over time --------------	
	a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom,...
          fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left		   		  
    a.plot1=plot(time_day(2:nt),solutionmass_surface_kgm(2:nt),...
             'g-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('water(kg/m)','FontSize',a.fs);
    xlim([0 time_day(end)])
yyaxis right
    a.plot1=plot(time_day(2:nt),saltmass_surface_kgm(2:nt),...
             'r-','linewidth',a.lw);hold off
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
	set(gca,'Xticklabel',[])
    ylabel('salt(kg/m)','FontSize',a.fs);
	xlim([0 time_day(end)])
	legend('surface water mass','surface salt mass','FontSize',a.fs,'Location','northwest' )	
	ax.YAxis(1).Color = 'g';
	ax.YAxis(2).Color = 'r';	
%% -------------  sub 1_3 total water/salt mass over time --------------	
	a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-fig_pos.height/2,...
          fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left		   		  
    a.plot1=plot(time_day(2:nt),total_solutionmass_kgm(2:nt),...
             'k-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('water(kg/m)','FontSize',a.fs);
    xlim([0 time_day(end)])
yyaxis right
    a.plot1=plot(time_day(2:nt),total_saltmass_kgm(2:nt),...
             'b-','linewidth',a.lw);hold off
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
	set(gca,'Xticklabel',[])
    ylabel('salt(kg/m)','FontSize',a.fs);
	xlim([0 time_day(end)])
	legend('total water mass','total salt mass','FontSize',a.fs,'Location','northwest' )
	ax.YAxis(1).Color = 'k';
	ax.YAxis(2).Color = 'b';
%% -------------  sub 1_4 total boundary water flux --------------	
	a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-fig_pos.height,...
          fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left		   		  
    a.plot1=plot(time_day(2:nt),boundary_solution_ms(2:nt),...
             'g-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('water flux (m/s)','FontSize',a.fs);
    xlim([0 time_day(end)])
yyaxis right
    a.plot1=plot(time_day(2:nt),cumu_boundary_solution_m(2:nt),...
             'r-','linewidth',a.lw);hold off
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
	set(gca,'Xticklabel',[])
    ylabel('cum. water flux (m)','FontSize',a.fs);
	xlim([0 time_day(end)])
	legend('water flux','cum. water flux','FontSize',a.fs,'Location','northwest' )
	ax.YAxis(1).Color = 'g';
	ax.YAxis(2).Color = 'r';
%% -------------  sub 1_5 total boundary salt flux --------------	
	a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-fig_pos.height-fig_pos.height/2,...
          fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left		   		  
    a.plot1=plot(time_day(2:nt),boundary_salt_kgsm2(2:nt),...
             'k-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('salt flux (kg/s/m^2)','FontSize',a.fs);
    xlim([0 time_day(end)])
yyaxis right
    a.plot1=plot(time_day(2:nt),cumu_boundary_salt_kgm2(2:nt),...
             'b-','linewidth',a.lw);hold off
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
	xlabel('Time (day)','FontSize',a.fs);
    ylabel('cum. salt flux (kg/m^2)','FontSize',a.fs);
	xlim([0 time_day(end)])
	legend('salt flux','cum. salt flux','FontSize',a.fs,'Location','northwest' )
	ax.YAxis(1).Color = 'k';
	ax.YAxis(2).Color = 'b';
	
	
%% -------------  sub 2_1 horizontal flux (from .ele for side boundary and from .bcop for bot. boundary) ---------------------
    a.sub2=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom+fig_pos.height/2,...
          fig_pos.length-0.05,fig_pos.height/2]);
% yyaxis left
% plot flux at a given level for side boundary condition 
	if boundary_condition == 1
yyaxis left										   
    a.plot2=plot(x_matrix(1,:), solution_ms(nt,:),...
             '-','color',[0.4940 0.1840 0.5560],'linewidth',a.lw);hold off		
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
	set(gca,'YColor',[0.9290 0.6940 0.1250]);
    set(gca,'fontsize',a.fs);
    ylabel('bottom water flux (m/s)','FontSize',a.fs);
    axis([0 x_matrix(1,end) -5e-5 5e-5])
yyaxis right		   		  
    a.plot2=plot(x_matrix(1,:), salt_kgsm2(nt,:),...
             'k-','linewidth',a.lw);hold off	
    get(gca,'xtick');
	set(gca,'YColor','k');
    set(gca,'fontsize',a.fs);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('bottom salt flux (kg/s/m^2)','FontSize',a.fs);
    axis([0 x_matrix(1,end) -5e-3 5e-3])
	legend('water','salt','FontSize',a.fs,'Location','northwest' )
%plot flux at a given level	
	elseif boundary_condition == 2
	x_flux_matrix_ms  = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1])*inp.por(1);
	y_flux_matrix_ms  = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1])*inp.por(1);
		for i=1:inp.nn2-1
			x_flux_ms(i) = x_flux_matrix_ms(yclosestIndex_plot(i),i);
			y_flux_ms(i) = y_flux_matrix_ms(yclosestIndex_plot(i),i);
		end
	flux_level_ms = y_flux_ms./abs(y_flux_ms).*(x_flux_ms.^2+y_flux_ms.^2).^0.5;	
    a.plot2=plot(x_ele_matrix(1,:), flux_level_ms,...
             '-','color',[0.4940 0.1840 0.5560],'linewidth',a.lw);hold off
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('flux (m/s)','FontSize',a.fs);
    axis([0 x_matrix(1,end) -0.002 0.002])
end
% yyaxis right		   		  
    % a.plot2=plot(x_matrix(1,:), salt_gday(nt,:),...
             % '-','color',[0.8500 0.3250 0.0980],'linewidth',a.lw);hold off	
    % get(gca,'xtick');
    % set(gca,'fontsize',a.fs);
	% set(gca,'Xticklabel',[])
    % ylabel('bottom salt flux (g/day)','FontSize',a.fs);
    % axis([0 x_matrix(1,end) -5 5])
	% legend('solution','salt','FontSize',a.fs,'Location','northwest' )
	% ax = gca;
%% -------------  sub 2_2 ET for all surface nodes  --------------
    a.sub2=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom,...
          fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left		   
    a.plot2=plot(x_matrix(1,:), evapo_mmday(nt,:),...
             'k-','linewidth',a.lw);hold off
	grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    xlabel('x','FontSize',a.fs);
    ylabel('Evt (mm/day)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 25])
yyaxis right
    solidmass_matrix_kg = reshape(nod(nt).terms{sm_idx},[inp.nn1,inp.nn2]);
    solidmass_surface_kg(1:inp.nn2) = solidmass_matrix_kg(inp.nn1,:);
    solidmass_thickness_mm(1:inp.nn2) = solidmass_surface_kg./c.density_solid_nacl_kgPm3./top_boundary_x_area_m2*c.m2mm;
    a.plot2 = plot(x_matrix(1,:), solidmass_thickness_mm(1:inp.nn2),...
             'r-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('Solid salt (mm)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 0.2])
	ax = gca;
	ax.YAxis(1).Color = 'k';
	ax.YAxis(2).Color = 'r';	
   
%% -------- sub 3 contour plot on Saturation --------- %the location overlaps with sub1-3 to 1-5
    % a.sub3=subplot('position'...
         % ,[fig_pos.left,fig_pos.bottom-0.32,...
          % fig_pos.length,fig_pos.height]);		  
    % % write pressure and conc in matrix form.
% yyaxis left
    vapori_plane = qv(nt).vapori_plane_y;
    % a.plot3=contourf(x_matrix,y_matrix,s_matrix,'EdgeColor','none');hold on
    % a.plot3=plot(x_matrix(1,:), vapori_plane,...
             	% '-','color',[0.72,0.27,1.00],'linewidth',a.lw);hold off
% %    scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    % color = jet;
    % color = flipud(color);
    % colormap(gca,color);
	% caxis([0 1])
    % cbsal = colorbar;
    % cbsal.Label.String = 'Saturation (-)';
    % get(gca,'xtick');
    % set(gca,'fontsize',a.fs);
    % xlabel('x (m)','FontSize',a.fs);
    % ylabel('Elevation (m)','FontSize',a.fs);
    % title('Saturation (-)')
% yyaxis right
    % a.plot3=plot(x_matrix(1,:), s_surface_matrix,...
             % 'w-','linewidth',a.lw);hold off
    % % ylabel('surface saturation (-)','FontSize',a.fs);
    % axis([0 x_matrix(1,end) -0.1 2])
    % yticks([0,0.5,1])

%% -------- sub 4 contour plot on concentration and quiver of vapor flow ---------
    a.sub4=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom-0.32,...
          fig_pos.length,fig_pos.height]);
% yyaxis left
%    a.plot4=contourf(x_matrix,y_matrix,c_matrix,'EdgeColor','none');hold on
%plot log colorbar
	a.plot4=contourf(x_matrix,y_matrix,log10(c_matrix+1-inp.ubc(1)*1000),'LevelStep',0.005,'LineStyle','None');hold on
    color = jet;
    colormap(gca,color);
	caxis([0.04 2.35])
	set(gca,'ColorScale','log')%log scale of colormap
	cbsal = colorbar;
    cbsal.Label.String = 'Salinity (ppt)';
	cbsal.TicksMode    = 'manual';
	cbsal.Ticks        = [0.04,0.3,0.78,1.2,1.82,2.35];
	cbsal.TickLabels   = [35,36,40,50,100,260];	
	
    a.plot4=plot(x_matrix(1,:), vapori_plane,...
             	'-','color',[0.72,0.27,1.00],'linewidth',a.lw);hold on
    qvx_mtx = qv(nt).qvx; % vector of vapor is acquired from the concentration difference bewteen two nodes,
    qvy_mtx = qv(nt).qvy; % which is used to calculate the vector in element center.
    qvx_plot_mtx = zeros(inp.nn1-1,inp.nn2-1);  
    qvy_plot_mtx = zeros(inp.nn1-1,inp.nn2-1);
    for i = 1:inp.nn2-1 %along x
        for j = 1:inp.nn1-1 %along y
            qvx_plot_mtx(j,i) = (qvx_mtx(j,i)+qvx_mtx(j+1,i))/2;
            qvy_plot_mtx(j,i) = (qvy_mtx(j,i)+qvy_mtx(j,i+1))/2;
        end
    end
%    a.plot4=quiver(x_ele_matrix,y_ele_matrix,qvx_plot_mtx,qvy_plot_mtx,'k-'); %plot vapor flow vector
	hold off
%	scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w'); %plot node location
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    title('Concentration (kg/kg)');
    xlabel('x (m)','FontSize',a.fs);
    ylabel('z (m)','FontSize',a.fs);
% yyaxis right
    % a.plot4=plot(x_matrix(1,:), c_surface_matrix,...
             % 'w-','linewidth',a.lw);hold off
    % ylabel('surface concentration (-)','FontSize',a.fs);
	% yticks([0,0.1,0.2,0.3])
    axis([0 x_matrix(1,end) 0 y_matrix(end,1)])
    
%% -------- sub 5 plot on surface sat vs concentration & solid salt vs evp  ---------
    a.sub5=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-0.64,...
           fig_pos.length-0.05,fig_pos.height/2]);
yyaxis left
    a.plot5=plot(x_matrix(1,:), s_surface_matrix,...
             '-','linewidth',a.lw,'color',[0.4660 0.6740 0.1880]);hold on		  
    a.plot5=plot(x_matrix(1,:), c_surface_matrix,...
             '.','linewidth',a.lw,'color',[0.4940 0.1840 0.5560]);hold off
	grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
	axis([0 x_matrix(1,end) 0 1.05])
    set(gca,'fontsize',a.fs);
	set(gca,'YColor',[0.4660 0.6740 0.1880]);
	ylabel('con or sat (-)','FontSize',a.fs);

yyaxis right	
    a.plot6=plot(x_matrix(1,:), evapo_mmday(nt,:),...
             '-','linewidth',a.lw,'color',[0 0.4470 0.7410]);hold off
	axis([0 x_matrix(1,end) 0 25])
    set(gca,'fontsize',a.fs);
	set(gca,'YColor',[0 0.4470 0.7410]);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Evt (mm/day)','FontSize',a.fs);
	legend('saturation','concentration','evaporation','FontSize',a.fs,'Location','northwest' )

% temperature contour plot    
    % write pressure and conc in matrix form.
    % temp_matrix = reshape(nod(nt).terms{temp_idx},[inp.nn1,inp.nn2]);    
    % a.plot6=contourf(x_matrix,y_matrix,temp_matrix-273.15,'EdgeColor','none');hold off
	
% %	scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    % color = jet;
    % colormap(gca,color);
	% caxis([20 40])
    % cbsal = colorbar;
    % cbsal.Label.String = 'Temperature (°C)';
    % get(gca,'xtick');
    % set(gca,'fontsize',a.fs);
    % title('Temperature (°C)');
    % xlabel('x (m)','FontSize',a.fs);
    % ylabel('Elevation (m)','FontSize',a.fs);
	
%% -------------  sub 6 velocity vector for each node  --------------
    a.sub6=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom-0.64,...
          fig_pos.length-0.05,fig_pos.height]);
    a.plot6   = quiver(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix);
	hold on
	startx    = x_flux_level(round(inp.nn2/24):round(inp.nn2/24):inp.nn2-round(inp.nn2/24));
	starty    = y_flux_level(round(inp.nn2/24):round(inp.nn2/24):inp.nn2-round(inp.nn2/24));
	% startx    = [x_ele_matrix(round(inp.nn1/12):round(inp.nn1/12):inp.nn1,160);x_ele_matrix(round(inp.nn1/12):round(inp.nn1/12):inp.nn1,200)];%vertical
	% starty    = [y_ele_matrix(round(inp.nn1/12):round(inp.nn1/12):inp.nn1,160);y_ele_matrix(round(inp.nn1/12):round(inp.nn1/12):inp.nn1,200)];
	a.plot6   = streamline(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix,startx,starty);
	hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Velocity')
    axis([0 x_matrix(1,end) 0 y_matrix(end,1)])	
	
%% ------------- sub 7 total solid mass of salt or zoom in vapor vector --------------
    a.sub7=subplot('position'...
         ,[fig_pos.left+0.75,fig_pos.bottom,...
          0.17,fig_pos.height]);
		  
% 	solidmass_column_g(1:inp.nn2) = sum (solidmass_matrix_kg(:,1:inp.nn2))*1000;
%     solidmass_total_g = zeros(1,time_step);
% 	solidmass_total_g(nt)   = sum(solidmass_column_g);
%     a.plot7 = scatter(time_day(nt),solidmass_total_g(nt),10,'bo');hold on
%     get(gca,'xtick');
%     set(gca,'fontsize',a.fs);
%     set(gca,'yaxislocation','right');
%     xlabel('time (day)','FontSize',a.fs);
% 	ylabel('Solid salt(g)','FontSize',a.fs);
%     axis([-0.1 time_day(end) 0 inf])	

    a.plot7=contourf(x_matrix,y_matrix,c_matrix,'EdgeColor','none');hold on
    a.plot7=plot(x_matrix(1,:), vapori_plane,...
             	'-','color',[0.72,0.27,1.00],'linewidth',a.lw);hold on
    a.plot7=quiver(x_ele_matrix,y_ele_matrix,qvx_plot_mtx,qvy_plot_mtx,...
                'color',[0.00,0.45,0.74]);hold off   
    color = jet;
    colormap(gca,color);
	caxis([0.035 0.264])
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
%    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Vapor flow')
    axis([0.75 0.85 0.31 0.345])    
    
%% ------------- sub 8 side boundary flux or porosity at given x --------------
    a.sub8=subplot('position'...
         ,[fig_pos.left+0.77,fig_pos.bottom-0.32,...
          0.15,fig_pos.height],'Color', 'none');
 %concentration profile at given x			
%     c_profile = c_matrix(:,(inp.nn2-1)*0.25+1:(inp.nn2-1)*0.25:inp.nn2);% integer is needed for colon operator (should be modified based on the x length)
%     a.plot8=plot(c_profile,y_matrix(:,(inp.nn2-1)*0.25+1:(inp.nn2-1)*0.25:inp.nn2),'-','linewidth',a.lw);hold off
%     get(gca,'xtick');
%     set(gca,'fontsize',a.fs);
%     set(gca,'yaxislocation','right');
%     xlabel('Concentration (-)','FontSize',a.fs);
%     legend('x=0.3 m','x=0.6 m','x=0.9 m','x=1.2 m','Location','southeast')
%     axis([0 0.3 0.2 0.4])
 
 %porosity at given x 
    % poro_matrix = reshape(nod(nt).terms{poro_idx},[inp.nn1,inp.nn2]);
    % poro_profile = poro_matrix(:,(inp.nn2-1)*0.25+1:(inp.nn2-1)*0.25:inp.nn2);          
    % a.plot8=plot(poro_profile, y_matrix(:,(inp.nn2-1)*0.25+1:(inp.nn2-1)*0.25:inp.nn2),'o-','linewidth',a.lw);hold off
    % get(gca,'xtick');
    % set(gca,'fontsize',a.fs);
    % set(gca,'yaxislocation','right');
    % xlabel('prorsity (-)','FontSize',a.fs);
    % axis([inp.por(1)-0.1 inp.por(1)+0.02 0.2 0.4])
    % legend('x=0.3 m','x=0.6 m','x=0.9 m','x=1.2 m','Location','southeast')
    % set(gca, 'XDir','reverse')

 %flux at side boundary
 % the flux plot shares the same y axis
 if boundary_condition ==2
	ax1=gca;
    ax1_pos = ax1.Position;
    a.plot8 = plot( solution_ms(nt,:),right_boundary_y_plot',...
        '-','color',[0.9290 0.6940 0.1250]	,'linewidth',a.lw);hold on
	grid on
	grid minor
	ax1.GridAlpha = 0.4;
	ax1.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
	set(gca,'XColor',[0.9290 0.6940 0.1250]);
    set(gca,'fontsize',a.fs);
    xlabel('boundary water flux (m/s)','FontSize',a.fs);
    axis([-10e-3 10e-3 0 y_matrix(end,end)])

	ax2 = axes('Position', ax1_pos, 'XAxisLocation', 'top', 'YAxisLocation', 'right');
	a.plot4=plot( salt_kgsm2(nt,:),right_boundary_y_plot',...
             'k-','linewidth',a.lw);hold off
	ax2.XAxisLocation = 'top';
	ax2.YAxisLocation = 'right';
 	ax2.Color = 'none';		
    get(gca,'xtick');
	set(gca,'XColor','k');
    set(gca,'fontsize',a.fs);
    xlabel('boundary salt flux (kg/s/m^2)','FontSize',a.fs);
    ylabel('z (m)','FontSize',a.fs);
    axis([-0.5 0.5 0 y_matrix(end,end)])
end
%% ------------- sub 9 zoom in velocity or flux profile along middle-x--------------
    a.sub9=subplot('position'...
         ,[fig_pos.left+0.77,fig_pos.bottom-0.64,...
          0.15,fig_pos.height]);
    vx_matrix = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1]);
    vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);
	% a.plot9=contourf(x_matrix,y_matrix,c_matrix,'EdgeColor','none');hold on
    % a.plot9=quiver(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix,20,'y');hold off

    % get(gca,'xtick');
    % set(gca,'fontsize',a.fs);
    % xlabel('x (m)','FontSize',a.fs);
    % ylabel('Elevation (m)','FontSize',a.fs);
    % title('Velocity')
    % axis([0.75 0.85 0.31 0.345])    
    
  	x_flux_matrix_ms  = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1])*inp.por(1);
	y_flux_matrix_ms  = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1])*inp.por(1);
	% locate the middle x
		x_middle_no = round((inp.nn2+1)/2);
		x_flux_middle_ms = x_flux_matrix_ms(:,x_middle_no);
		y_flux_middle_ms = y_flux_matrix_ms(:,x_middle_no);
		
	flux_middle_ms = (x_flux_middle_ms.^2+y_flux_middle_ms.^2).^0.5;	
    a.plot2=plot(flux_middle_ms,y_ele_matrix(:,x_middle_no), ...
             '-','color',[0.4940 0.1840 0.5560],'linewidth',a.lw);hold off
    grid on
	grid minor
	ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;	
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    xlabel('flux (m/s)','FontSize',a.fs);
    % axis([0 x_matrix(1,end) -0.002 0.002])					
    
	F = getframe(gcf); % save the current figure
    writeVideo(mov,F);% add it as the next frame of the movie

end
saveas(a.fig,'a','fig')
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
	caxis([0.035 0.264])
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
