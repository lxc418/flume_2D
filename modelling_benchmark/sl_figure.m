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

%% solute inflow from bottom (from bcof without the vapor contribution)

for i= 1:inp.nn2

    solute_kgs(i,:)  = -arrayfun(@(y) y.qpu(i),bcop);
    
end
solute_gday= solute_kgs'.*c.kg2g*c.secPday;

area1_m2    = (x_matrix(1,2)-x_matrix(1,1))*inp.z(1);
for i=2:time_step

	evapo_mmday(1,:)  		 =  reshape(et(1).terms{et_idx},[1,inp.nn2])*c.ms2mmday;
	total_evapo_mmday(1)     =  sum (evapo_mmday(1,:))./inp.nn2; %the evp rate from the whole surface
	evapo_mmday(i,:)  		 =  reshape(et(i).terms{et_idx},[1,inp.nn2])*c.ms2mmday;
	total_evapo_mmday(i)     =  sum (evapo_mmday(i,:))./inp.nn2;
	
    cumulative_evapo_mm(1)   =  total_evapo_mmday(1)*inp.scalt*inp.nbcfpr*c.dayPsec;
    cumulative_evapo_mm(i)   =  total_evapo_mmday(i)*inp.scalt*inp.nbcfpr*c.dayPsec + cumulative_evapo_mm(i-1);

end


%% plot control
fig_pos.left   = 0.05;
fig_pos.bottom = 0.7;
fig_pos.length = 0.35;
fig_pos.height = 0.26;

%nt=10;
a.fs = 15;
a.lw = 2; %line width
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

for nt=1:round(time_step/40):time_step
%% -------------  sub 1 ET over time  --------------
    a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom,...
          fig_pos.length-0.05,fig_pos.height]);
yyaxis left
    a.plot1=plot(time_day(1:nt),total_evapo_mmday(1:nt),...
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
    axis([-0.1 time_day(end) -0.5 inf])	
yyaxis right
    a.plot1=plot(time_day(1:nt),cumulative_evapo_mm(1:nt),...
             'b--','linewidth',a.lw);hold off
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    txt=sprintf('Result at day %.2f',nod(nt).tout*c.dayPsec);
    title(txt);
    xlabel('Time (day)','FontSize',a.fs);
    ylabel('Cumulative Evt (mm)','FontSize',a.fs);
    axis([-0.1 time_day(end) 0 inf])

	ax = gca;
	ax.YAxis(1).Color = 'k';
	ax.YAxis(2).Color = 'b';
    %title('ead profile');
%% -------------  sub 2_1 bottom boundary solute flux  --------------
    a.sub2=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom+fig_pos.height/2,...
          fig_pos.length-0.05,fig_pos.height/2]);	   
    a.plot2=plot(x_matrix(1,:), solute_gday(nt,:),...
             '-','color',[0.8500 0.3250 0.0980],'linewidth',a.lw);hold off
	grid on
	grid minor
		ax = gca;
	ax.GridAlpha = 0.4;
	ax.MinorGridAlpha = 0.5;
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('bottom solute flux (g/day)','FontSize',a.fs);
    axis([0 x_matrix(1,end) -5 5])

	ax = gca;	
	ax.XAxis.Visible = 'off';
	ax.YAxis.Color = [0.8500 0.3250 0.0980];
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
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    xlabel('x','FontSize',a.fs);
    ylabel('Evt (mm/day)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 25])
yyaxis right
    solidmass_matrix_kg = reshape(nod(nt).terms{sm_idx},[inp.nn1,inp.nn2]);
    solidmass_surface_kg(1:inp.nn2) = solidmass_matrix_kg(inp.nn1,:);
    solidmass_thickness_mm(1:inp.nn2) = solidmass_surface_kg./c.density_solid_nacl_kgPm3./area1_m2*c.m2mm;
    a.plot2   =  scatter(x_matrix(1,:), solidmass_thickness_mm(1:inp.nn2),10,...
             'r*');hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    ylabel('Solid salt (mm)','FontSize',a.fs);
    axis([0 x_matrix(1,end) 0 6])
	ax = gca;
	ax.YAxis(1).Color = 'k';
	ax.YAxis(2).Color = 'r';	
   
%% -------- sub 3 contour plot on Saturation ---------
    a.sub4=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-0.32,...
          fig_pos.length,fig_pos.height]);		  
    % write pressure and conc in matrix form.
    s_matrix  = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2]);
yyaxis left
    vapori_plane = qv(nt).vapori_plane_y;
    a.plot4=contourf(x_matrix,y_matrix,s_matrix,'EdgeColor','none');hold on
    a.plot4=plot(x_matrix(1,:), vapori_plane,...
             	'-','color',[0.72,0.27,1.00],'linewidth',a.lw);hold off
%    scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    color = jet;
    color = flipud(color);
    colormap(gca,color);
	caxis([0 1])
    cbsal = colorbar;
    cbsal.Label.String = 'Saturation (-)';
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Saturation (-)')
yyaxis right
    s_surface_matrix = s_matrix(inp.nn1,:);
    a.plot4=plot(x_matrix(1,:), s_surface_matrix,...
             'w-','linewidth',a.lw);hold off
    % ylabel('surface saturation (-)','FontSize',a.fs);
    axis([0 x_matrix(1,end) -0.1 2])
    yticks([0,0.5,1])

%% -------- sub 4 contour plot on concentration ---------
    a.sub5=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom-0.32,...
          fig_pos.length,fig_pos.height]);
    % write pressure and conc in matrix form.
    c_matrix = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2]);
yyaxis left
    a.plot5=contourf(x_matrix,y_matrix,c_matrix,'EdgeColor','none');hold on
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
    a.plot5=quiver(x_ele_matrix,y_ele_matrix,qvx_plot_mtx,qvy_plot_mtx,'k-');hold off
%	scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    color = jet;
    colormap(gca,color);
	caxis([0.035 0.264])
    cbsal = colorbar;
    cbsal.Label.String = 'Concentration (-)';
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    title('Concentration (kg/kg)');
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    %axis([10, 40,9,10])
yyaxis right
    c_surface_matrix = c_matrix(inp.nn1,:);
    a.plot5=plot(x_matrix(1,:), c_surface_matrix,...
             'w-','linewidth',a.lw);hold off
    % ylabel('surface concentration (-)','FontSize',a.fs);
    axis([0 x_matrix(1,end) -0.05 0.6])
    yticks([0,0.1,0.2,0.3])
    
%% -------- sub 5 plot on surface sat vs concentration & solid salt vs evp  ---------
    a.sub6=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-0.64,...
           fig_pos.length-0.05,fig_pos.height]);
yyaxis left
    a.plot6=plot(x_matrix(1,:), s_surface_matrix,...
             '-','linewidth',a.lw,'color',[0.4660 0.6740 0.1880]);hold on		  
    a.plot6=plot(x_matrix(1,:), c_surface_matrix,...
             '.','linewidth',a.lw,'color',[0.4660 0.6740 0.1880]);hold off
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
    a.sub3=subplot('position'...
         ,[fig_pos.left+0.4,fig_pos.bottom-0.64,...
          fig_pos.length-0.05,fig_pos.height]);
    vx_matrix = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1]);
    vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);
    a.plot3   = quiver(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix);
	hold on
	startx    = 0.105:0.1:1.105;
	starty    = y_ele_matrix(1,11:10:111);
	a.plot3   = streamline(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix,startx,starty);
	hold off
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Velocity')
    axis([0 x_matrix(1,end) 0 0.4])	
	
%% ------------- total solid mass of salt or zoom in vapor vector --------------
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
    axis([0.6 1 0.3 0.38])    
    
%% ------------- concentration profile or porosity at given x --------------
    a.sub8=subplot('position'...
         ,[fig_pos.left+0.77,fig_pos.bottom-0.32,...
          0.15,fig_pos.height],'Color', 'none');
%     c_profile = c_matrix(:,(inp.nn2-1)*0.25+1:(inp.nn2-1)*0.25:inp.nn2);% integer is needed for colon operator (should be modified based on the x length)
%     a.plot8=plot(c_profile,y_matrix(:,(inp.nn2-1)*0.25+1:(inp.nn2-1)*0.25:inp.nn2),'-','linewidth',a.lw);hold off
%     get(gca,'xtick');
%     set(gca,'fontsize',a.fs);
%     set(gca,'yaxislocation','right');
%     xlabel('Concentration (-)','FontSize',a.fs);
%     legend('x=0.3 m','x=0.6 m','x=0.9 m','x=1.2 m','Location','southeast')
%     axis([0 0.3 0.2 0.4])
    
    poro_matrix = reshape(nod(nt).terms{poro_idx},[inp.nn1,inp.nn2]);
    poro_profile = poro_matrix(:,(inp.nn2-1)*0.25+1:(inp.nn2-1)*0.25:inp.nn2);          
    a.plot8=plot(poro_profile, y_matrix(:,(inp.nn2-1)*0.25+1:(inp.nn2-1)*0.25:inp.nn2),'o-','linewidth',a.lw);hold off
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    set(gca,'yaxislocation','right');
    xlabel('prorsity (-)','FontSize',a.fs);
    axis([inp.por(1)-0.1 inp.por(1)+0.02 0.2 0.4])
    legend('x=0.3 m','x=0.6 m','x=0.9 m','x=1.2 m','Location','southeast')
    set(gca, 'XDir','reverse')

%% ------------- zoom in velocity --------------
    a.sub9=subplot('position'...
         ,[fig_pos.left+0.77,fig_pos.bottom-0.64,...
          0.15,fig_pos.height]);
    vx_matrix = reshape(ele(nt).terms{vx_idx},[inp.nn1-1,inp.nn2-1]);
    vy_matrix = reshape(ele(nt).terms{vy_idx},[inp.nn1-1,inp.nn2-1]);
    a.plot9=quiver(x_ele_matrix,y_ele_matrix,vx_matrix,vy_matrix);hold off
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',a.fs);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Velocity')
    axis([x_matrix(1,end)-0.3 x_matrix(1,end) 0.2 0.3])
    
    
    
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
,