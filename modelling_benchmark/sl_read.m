clear all
close all
% project name
name='COLUMN';

% read input files
%inp=i0('desaline');
fil  =  readFIL;
inp  = vapinpObj(fil.basename,'block_reading','yes');
nod  = readNOD(fil.basename);
ele  = readELE( fil.basename);
bcop = readBCOP(fil.basename);
bcof = readBCOF(fil.basename);

%index for p c and s
x_idx  = strcmp(nod(1).label,'X');
y_idx  = strcmp(nod(1).label,'Y');
p_idx  = strcmp(nod(1).label,'Pressure');
c_idx  = strcmp(nod(1).label,'Concentration');
s_idx  = strcmp(nod(1).label,'Saturation');
sm_idx = strcmp(nod(1).label,'SOLIDMASS(KG)');
poro_idx  = strcmp(nod(1).label,'POROSITY');
temp_idx  = strcmp(nod(1).label,'TEMP');
vy_idx    = strcmp(ele(1).label,'Y velocity');
vx_idx    = strcmp(ele(1).label,'X velocity');
xele_idx  = strcmp(ele(1).label,'X origin');
yele_idx  = strcmp(ele(1).label,'Y origin');
kr_idx    = strcmp(ele(1).label,'RELK');

%read vapor transport
qv   = readQV(inp,ele,nod);
%readlabdata;
save('plot.mat','-v7.3') 