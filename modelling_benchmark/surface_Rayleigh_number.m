%calculate and plor Rayleigh number along soil surface
%can be used for evenly mesh Y of surface layer 
c=ConstantObj();

layer_thickness_m = 0.01;%the thickness of the surface layer

density_variation = 700; %for sodium chloride d(tho_water)/d(concentration)
porosity          = inp.por(1);
dynamic_viscosity = 1e-3;
diffusivity       = inp.sigmaw;
permeability      = inp.pmax(1);

day_output = [1,5,10];%set days need to be plotted			 																  					 																  		
timestep_output = round(day_output*c.secPday/inp.nprint/inp.scalt)+1;
timestep_output (timestep_output>time_step)=time_step;

for nt = timestep_output
time_step = nt;
time_day  = [bcof.tout]/3600/24;%second to day
x_matrix = reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);%inp.nn2 is number of nodes in x direction 
y_matrix = reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);
x_ele_matrix = reshape(ele(1).terms{xele_idx},[inp.nn1-1,inp.nn2-1]);
y_ele_matrix = reshape(ele(1).terms{yele_idx},[inp.nn1-1,inp.nn2-1]);

kr_matrix = reshape(ele(time_step).terms{kr_idx},[inp.nn1-1,inp.nn2-1]);
c_matrix  = reshape(nod(time_step+1).terms{c_idx},[inp.nn1,inp.nn2]);
y_layer_bottom = y_matrix(end,:)-layer_thickness_m;
[y_diff,layer_nod_index] = min(abs(y_matrix-y_layer_bottom));%find the location of the layer bottom in matrix

surface_kr= kr_matrix(layer_nod_index(1):end,:);
average_surface_kr = mean(surface_kr);
real_surface_permeability = average_surface_kr*permeability;
concentration_differnece  = c_matrix(end,:)-c_matrix(layer_nod_index(1),:);%calculate concentration difference between the surface and bottom of the layer

Rayleigh_number = density_variation*c.g*real_surface_permeability.*concentration_differnece(1:end-1)*layer_thickness_m./...
    (porosity*dynamic_viscosity*diffusivity);

%plot
a.fs = 10;
a.lw = 1.75; %line width
a.cz = 8; %the size of the marker

figure
tiledlayout(3,1);
ax1 = nexttile;
semilogy (x_ele_matrix(1,:),real_surface_permeability,'k-','linewidth',a.lw);hold off
    grid on
	grid minor
% 	ax = gca;
	ax1.GridAlpha = 0.4;
	ax1.MinorGridAlpha = 0.5;	
    ylabel('permeability (mm^2)','FontSize',a.fs);
    xlim ([0 x_matrix(1,end)])
	yticks([1e-19 1e-16 1e-13 1e-10])
	yticklabels({'1e-19', '1e-16', '1e-13', '1e-10'})
ax2 = nexttile;
plot (x_matrix(1,:),concentration_differnece,'b-','linewidth',a.lw);hold off
    grid on
	grid minor
% 	ax = gca;
	ax2.GridAlpha = 0.4;
	ax2.MinorGridAlpha = 0.5;	
    ylabel('\Delta C (-)','FontSize',a.fs);

ax3 = nexttile;
plot (x_ele_matrix(1,:),Rayleigh_number,'r-','linewidth',a.lw);hold off
    grid on
	grid minor
% 	ax = gca;
	ax3.GridAlpha = 0.4;
	ax3.MinorGridAlpha = 0.5;	
    ylabel('Ra (-)','FontSize',a.fs);
    xlabel('x (m)','FontSize',a.fs);

	figure_name=sprintf('Rayleigh_day_%.2f.fig',nod(nt+1).tout*c.dayPsec);
	saveas(gcf,figure_name)

end


