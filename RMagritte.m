function RMagritte
p=get_parameter('init');
prepare_mainfig(p);
add_log('Rene Magritte. v1.3 by SHJK. (last update FEB 28, 2020)');

function prepare_mainfig(p)
fhd=findobj('Tag','Rene_main');
if isempty(fhd)
    fhd=figure;
    fhd.MenuBar='none';
    fhd.Name='Rene Magritte';
    fhd.Tag='Rene_main';
    fhd.NumberTitle='off';
    fhd.Position=[500 100 610 400];
    fhd.Resize='off';
end

figure(fhd);clf

g_pos1=[20 340 0 0];
hp = uipanel('Title','File path','FontSize',8,...
    'Unit','pixel',...
    'Position',g_pos1+[-5 -5 478 50]);

uicontrol('Style','text',    'Tag','RM_pth',   'String',p.pth,...
    'HorizontalAlignment','Left',...
    'Position',g_pos1+[10 -3 460 30]);

g_pos1=[500 390 80 30];

uicontrol('Style','pushbutton',    'String','Load data',...
    'Position',g_pos1+[0 -40 0 0],'Callback',@load_data);

uicontrol('Style','pushbutton',    'String','Remove junks',...
    'Position',g_pos1+[0 -70 0 0],'Callback',@rm_junk);

uicontrol('Style','pushbutton',    'String','Find peaks',...
    'Position',g_pos1+[0 -100 0 0],'Callback',@find_binding_events);

uicontrol('Style','pushbutton',    'String','Analyze peaks',...
    'Position',g_pos1+[0 -130 0 0],'Callback',@ana_dwells);

uicontrol('Style','checkbox','Tag','RM_flg2go_ana_dwells',...
    'Position',g_pos1+[82 -125 0 -15]);

uicontrol('Style','pushbutton',    'String','Show traces',...
    'Position',g_pos1+[0 -200 0 0],'Callback',@show_idv_mols);
uicontrol('Style','checkbox','Tag','RM_flg2go_show_idv_mols',...
    'Enable','off',...
    'Position',g_pos1+[82 -195 0 -15]);

uicontrol('Style','pushbutton',    'String','Show dwells',...
    'Position',g_pos1+[0 -230 0 0],'Callback',@show_dwells);

uicontrol('Style','checkbox','Tag','RM_flg2go_show_dwells',...
    'Position',g_pos1+[82 -225 0 -15]);

uicontrol('Style','pushbutton',    'String','Heterogeneity',...
    'Position',g_pos1+[0 -260 0 0],'Callback',@show_hetero);

uicontrol('Style','checkbox','Tag','RM_flg2go_show_hetero',...
    'Position',g_pos1+[82 -255 0 -15]);

uicontrol('Style','pushbutton',    'String','Save results',...
    'Position',g_pos1+[0 -300 0 0],'Callback',@save_results);

uicontrol('Style','checkbox','Tag','RM_flg2go_save_results',...
    'Position',g_pos1+[82 -295 0 -15]);

uicontrol('Style','pushbutton',    'String','Load results',...
    'Position',g_pos1+[0 -330 0 0],'Callback',@load_results);

uicontrol('Style','pushbutton',    'String','Save traces',...
    'Position',g_pos1+[0 -360 0 0],'Callback',@save_traces);

g_pos1=[10 312 95 15];
g_pos2=[110 310 60 20];

uicontrol('Style','text',    'String','Time unit (s)','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 0 0 0]);
uicontrol('Style','edit','Tag','RM_time_unit','String',p.time_unit,...
    'Position',g_pos2+[0 0 0 0]);

uicontrol('Style','text',    'String','Leakage','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -25 0 0]);
uicontrol('Style','edit','Tag','RM_leakage','String',p.leakage,...
    'Position',g_pos2+[0 -25 0 0]);

uicontrol('Style','text',    'String','Gamma','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -50 0 0]);
uicontrol('Style','edit','Tag','RM_gamma','String',p.gamma,...
    'Position',g_pos2+[0 -50 0 0]);

uicontrol('Style','text',    'String','BG donor','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -75 0 0]);
uicontrol('Style','edit','Tag','RM_BG_d','String',p.BG_d,...
    'Position',g_pos2+[0 -75 0 0]);

uicontrol('Style','text',    'String','BG acceptor','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -100 0 0]);
uicontrol('Style','edit','Tag','RM_BG_a','String',p.BG_a,...
    'Position',g_pos2+[0 -100 0 0]);

uicontrol('Style','text',    'String','auto BG correction',...
    'Position',g_pos1+[40 -125 0 0]);
uicontrol('Style','checkbox','Tag','RM_auto_BG','Value',p.auto_BG,...
    'Position',g_pos2+[40 -125 0 0]);


g_pos1=[10 152 95 15];
g_pos2=[110 150 60 20];

uicontrol('Style','text',    'String','Int max.','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 0 0 0]);
uicontrol('Style','edit','Tag','RM_y_max','String',p.y_max,...
    'Position',g_pos2+[0 0 0 0]);

uicontrol('Style','text',    'String','Int min.','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -25 0 0]);
uicontrol('Style','edit','Tag','RM_y_min','String',p.y_min,...
    'Position',g_pos2+[0 -25 0 0]);

uicontrol('Style','text',    'String','Time max (s)','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -50 0 0]);
uicontrol('Style','edit','Tag','RM_xmax','String',p.xmax,...
    'Position',g_pos2+[0 -50 0 0]);

uicontrol('Style','text',    'String','update plots','HorizontalAlignment','Right',...
    'Position',g_pos1+[40 -75 0 0]);
uicontrol('Style','checkbox','Tag','RM_update_figures','Value',1,...
    'Position',g_pos2+[40 -75 0 0]);


g_pos1=[170 312 95 15];
g_pos2=[270 310 60 20];

uicontrol('Style','text',    'String','Keyword','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 0 0 0]);
uicontrol('Style','edit','Tag','RM_user_keyword','String',p.user_keyword,...
    'Position',g_pos2+[0 0 0 0]);

uicontrol('Style','text',    'String','# binning','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -25 0 0]);
uicontrol('Style','edit','Tag','RM_N_bin','String',p.N_bin,...
    'Position',g_pos2+[0 -25 0 0]);

uicontrol('Style','text',    'String','Threshold','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -50 0 0]);
uicontrol('Style','edit','Tag','RM_auto_th','String',p.auto_th,...
    'Position',g_pos2+[0 -50 0 0]);

uicontrol('Style','text',    'String','min peak len (s)','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -75 0 0]);
uicontrol('Style','edit','Tag','RM_min_dw_len_sel','String',p.min_dw_len_sel,...
    'Position',g_pos2+[0 -75 0 0]);


uicontrol('Style','text',    'String','min peak int.','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -100 0 0]);
uicontrol('Style','edit','Tag','RM_min_peak_int','String',p.peak_int(1),...
    'Position',g_pos2+[0 -100 0 0]);

uicontrol('Style','text',    'String','max peak int.','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -125 0 0]);
uicontrol('Style','edit','Tag','RM_max_peak_int','String',p.peak_int(2),...
    'Position',g_pos2+[0 -125 0 0]);



g_pos1=[335 312 85 15];
g_pos2=[425 310 60 20];
uicontrol('Style','text',    'String','E tolerance','HorizontalAlignment','Right',...
    'Position',g_pos1+[0 -75 0 0]);
uicontrol('Style','edit','Tag','RM_E_tolerance','String',p.E_tolerance,...
    'Position',g_pos2+[0 -75 0 0]);


g_pos1=[317 237 125 15];
g_pos2=[447 235 40 20];
uicontrol('Style','text',    'String','min # peaks per state','HorizontalAlignment','Right',...
    'Position',g_pos1+[40 -30 0 0]);
uicontrol('Style','edit','Tag','RM_min_N_peak_in_state','String',p.min_N_peak_in_state,...
    'Position',g_pos2+[0 -50 0 0]);

hdl_state_find_algorithm=uibuttongroup('Title','','Visible','on','Tag','RM_state_find_algorithm',...
    'Unit','pixel','Position',[390 268 95 60]);
r1 = uicontrol(hdl_state_find_algorithm,'Style','radiobutton','Tag','RM_peak_by_threshold','String','Threshold',...
    'Unit','pixel','Position',[7 10 70 15]);
r2 = uicontrol(hdl_state_find_algorithm,'Style','radiobutton','Tag','RM_peak_by_k-means','String','k-means','Value',1,...
    'Unit','pixel','Position',[7 32 70 15]);


t=uitable('Tag','RM_states_tbl','Data',{true 'A'  0 0;true 'B'   0 0; true 'C'   0 0},...
    'Position',[265 72 222 100],'ColumnWidth',{30,40,60,60});
t.ColumnName = {'USE','Name','Position','Width'};
t.ColumnEditable = true;

uicontrol('Style','text',    'Tag','RM_log_text',   'String','',...
    'HorizontalAlignment','Left',...
    'Position',[50 30 450 30]);

g_pos1=[50 20 430 10];
hdl_prog=axes('Unit','pixel','Position',g_pos1,'Tag','RM_prog_bar');
axis([0 1 0 1]);axis off;
set_prog(hdl_prog,0)

function p=get_parameter(r_mode,p)
if nargin<2
    % case of 'init' or 'fetch'
end
switch r_mode
    case 'init'
        mfile_path_ful=mfilename('fullpath');
        machine_type=computer;
        if strfind(machine_type,'MAC')
            tokens=strfind(mfile_path_ful,'/');
        else
            tokens=strfind(mfile_path_ful,'\');
        end
        mfile_path=mfile_path_ful(1:tokens(end)-1);
        mfile_name=mfile_path_ful(tokens(end)+1:end);
        if exist(fullfile(mfile_path, 'pgdata',['pgdata_' mfile_name '.mat']),'file')
            load(fullfile(mfile_path, 'pgdata',['pgdata_' mfile_name '.mat']));
        end
        if ~exist('pth','var') || ~isstr(pth)
            if strfind(machine_type,'MAC')
                pth='/Users/sunghyunkim/Documents/MATLAB';
            else
                pth='c:\';
            end
        end
        
        p.pth=pth;
        p.y_max=1000;
        p.y_min=-p.y_max/10;
        p.xmax=300;
        p.user_keyword='blink';
        p.N_bin=0;
        p.auto_th=2000;
        p.update_figures=1;
        p.auto_BG=1;
        p.min_dw_len_sel=0.3;
        p.state_find_algorithm='k-means';
        p.min_N_peak_in_state=5;
%         p.N_FRET_state=4;
        p.time_unit=0.1;
        p.leakage=0;
        p.gamma=1;
        p.BG_d=0;
        p.BG_a=0;
        p.peak_int=[0 Inf];
        p.E_tolerance=0.5;
        
    case 'fetch'
        tmphdl=findobj('Tag','RM_pth');        p.pth=tmphdl.String;
        tmphdl=findobj('Tag','RM_y_max');        p.y_max=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_y_min');        p.y_min=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_xmax');        p.xmax=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_user_keyword');        p.user_keyword=tmphdl.String;
        tmphdl=findobj('Tag','RM_N_bin');        p.N_bin=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_auto_th');        p.auto_th=str2num(tmphdl.String);
        
        tmphdl=findobj('Tag','RM_min_peak_int');        p.peak_int(1)=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_max_peak_int');        p.peak_int(2)=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_E_tolerance');        p.E_tolerance=str2num(tmphdl.String);
        
        
        tmphdl=findobj('Tag','RM_update_figures');        p.update_figures=tmphdl.Value;
        tmphdl=findobj('Tag','RM_auto_BG');        p.auto_BG=tmphdl.Value;
        tmphdl=findobj('Tag','RM_min_dw_len_sel');        p.min_dw_len_sel=str2num(tmphdl.String);
        
        tmphdl=findobj('Tag','RM_state_find_algorithm');        p.state_find_algorithm=tmphdl.SelectedObject.String;
        tmphdl=findobj('Tag','RM_min_N_peak_in_state');        p.min_N_peak_in_state=str2num(tmphdl.String);
%         tmphdl=findobj('Tag','RM_N_FRET_state');        p.N_FRET_state=str2num(tmphdl.String);
        
        tmphdl=findobj('Tag','RM_time_unit');        p.time_unit=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_leakage');        p.leakage=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_gamma');        p.gamma=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_BG_d');        p.BG_d=str2num(tmphdl.String);
        tmphdl=findobj('Tag','RM_BG_a');        p.BG_a=str2num(tmphdl.String);
        p.hdl_prog=findobj('Tag','RM_prog_bar');
        
    case 'set'
        tmphdl=findobj('Tag','RM_pth');        tmphdl.String=p.pth;
        tmphdl=findobj('Tag','RM_y_max');        tmphdl.String=p.y_max;
        tmphdl=findobj('Tag','RM_y_min');        tmphdl.String=p.y_min;
        tmphdl=findobj('Tag','RM_xmax');        tmphdl.String=p.xmax;
        tmphdl=findobj('Tag','RM_user_keyword');       tmphdl.String=p.user_keyword;
        tmphdl=findobj('Tag','RM_N_bin');       tmphdl.String=p.N_bin;
        tmphdl=findobj('Tag','RM_auto_th');       tmphdl.String=p.auto_th;
        
        tmphdl=findobj('Tag','RM_min_peak_int');        tmphdl.String=p.peak_int(1);
        tmphdl=findobj('Tag','RM_max_peak_int');        tmphdl.String=p.peak_int(2);
        tmphdl=findobj('Tag','RM_E_tolerance');        tmphdl.String=p.E_tolerance;
        
        tmphdl=findobj('Tag','RM_update_figures');       tmphdl.Value=p.update_figures;
        tmphdl=findobj('Tag','RM_auto_BG');        tmphdl.Value=p.auto_BG;
        tmphdl=findobj('Tag','RM_min_dw_len_sel');       tmphdl.String=p.min_dw_len_sel;
        
        tmphdl=findobj('Tag','RM_state_find_algorithm');        %tmphdl.SelectedObject.String=p.state_find_algorithm;
        switch   p.state_find_algorithm
            case 'HMM'
                tmp_selected=findobj('Tag','RM_peak_by_HMM');
            case 'Threshold'
                tmp_selected=findobj('Tag','RM_peak_by_threshold');
        end
        
        tmphdl=findobj('Tag','RM_min_N_peak_in_state');        tmphdl.String=p.min_N_peak_in_state;
%         tmphdl=findobj('Tag','RM_N_FRET_state');        tmphdl.String=p.N_FRET_state;
        
        tmphdl=findobj('Tag','RM_time_unit');        tmphdl.String=p.time_unit;
        tmphdl=findobj('Tag','RM_leakage');        tmphdl.String=p.leakage;
        tmphdl=findobj('Tag','RM_gamma');        tmphdl.String=p.gamma;
        tmphdl=findobj('Tag','RM_BG_d');        tmphdl.String=p.BG_d;
        tmphdl=findobj('Tag','RM_BG_a');        tmphdl.String=p.BG_a;
        p.hdl_prog=findobj('Tag','RM_prog_bar');
        drawnow;
end

function load_data(~, ~)
% set_path();
p=get_parameter('fetch');
set_prog(p.hdl_prog,0);
hdl_log_text=findobj('Tag','RM_log_text');

%%
choice = questdlg('Choose data format', ...
    'Data format', ...
    'ASCII','Binary','Binary');

switch choice
    case 'ASCII'
        %% get file infor       
        pth=uigetdir(p.pth);

        %% load data
%         u.filelist=dir(fullfile(pth,['\*' p.user_keyword '*']));
        u.filelist=dir(fullfile(pth,['*' p.user_keyword '*']));
        [u.N_mol, ~]=size(u.filelist);
        
        for mid=1:u.N_mol
            add_log(['loading "' u.filelist(mid).name '"  (' num2str(mid) '/' num2str(u.N_mol) ')'],hdl_log_text);
            %% load data
            rawdata=load([pth '\' u.filelist(mid).name]);
            if ~isempty(rawdata)
                %% data read
                u.time{mid}=rawdata(:,1);
                u.len(mid)=length(u.time{mid});
                u.donor{mid}=rawdata(:,2);
                u.acceptor{mid}=rawdata(:,3);
            end
            if rem(mid,floor(u.N_mol/20))==0
                set_prog(p.hdl_prog,mid/u.N_mol);
            end
        end
        p.time_unit=u.time{1}(2)-u.time{1}(1);
        
    case 'Binary'
        %% get file info
        [fname,pth,~]=uigetfile(fullfile(p.pth, '*.traces*'),'Chooose traces files');
%          [fname,pth,~]=uigetfile([p.pth '\*.traces*'],'Chooose traces files');     
        %% get molecule position info
        tokens=strfind(fname,'.');
        fname_base=fname(1:tokens(end)-1);
        
        if strcmp(fname(tokens(end)+1:end),'traces2')
            binary_format_def='32-32-16';
        else
            binary_format_def='32-16-16';
        end
        if exist(fullfile(pth, [fname_base '.pos']),'file')
            mol_pos=load(fullfile(pth, [fname_base '.pos']));
        else
            mol_pos=[];
        end
        
        %% read data
        switch binary_format_def
            case '32-16-16'   
                disp('trace1 loaded');
                fid=fopen([pth fname],'r');
                len=fread(fid,1,'int32');
                Ntraces=fread(fid,1,'int16');
                u.N_mol=Ntraces/2;
                raw=fread(fid,Ntraces*len,'int16');
                fclose(fid);
            case '32-32-16'   
                disp('trace2 loaded');
                fid=fopen([pth fname],'r');
                len=fread(fid,1,'int32');
                Ntraces=fread(fid,1,'int32');
                u.N_mol=Ntraces/2;
                raw=fread(fid,Ntraces*len,'int16');
                fclose(fid);
        end
        %% convert into donor and acceptor traces
        Data=reshape(raw,Ntraces,len);

%         fr2ana=[1 len-50];
        fr2ana=[1 len];
        frid_selected=fr2ana(1):fr2ana(2); 
        msgbox(['frames' num2str(fr2ana(1)) ' to ' num2str(fr2ana(2)) ' selected']);
        len=fr2ana(2)-fr2ana(1)+1;
        donor_raw=Data(1:2:end,fr2ana(1):fr2ana(2));
        acceptor_raw=Data(2:2:end,fr2ana(1):fr2ana(2));
        
        u.donor=mat2cell(donor_raw',len,ones(u.N_mol,1));
        u.acceptor=mat2cell(acceptor_raw',len,ones(u.N_mol,1));
        u.time(1:u.N_mol)={(1:len)'*p.time_unit};
        u.len(1:u.N_mol)=len;
        
        %% make dummy file list name
        for mid=1:u.N_mol
            u.filelist(mid).name=[fname_base ' tr' num2str(mid)];
        end
        
        
end

if u.N_mol>0
    u.time_org=u.time;
    u.donor_org=u.donor;
    u.acceptor_org=u.acceptor;
    u.N_mol_org=u.N_mol;
end

fhd=findobj('Tag','Rene_main');
fhd.UserData.u=u;
fhd.UserData.dw={};

%% save file path for next use
mfile_path_ful=mfilename('fullpath');
% system_type=computer;
if contains(computer,'MAC')
    tokens=strfind(mfile_path_ful,'/');
else
    tokens=strfind(mfile_path_ful,'\');
end
mfile_path=mfile_path_ful(1:tokens(end)-1);
mfile_name=mfile_path_ful(tokens(end)+1:end);
if ~exist([mfile_path '\pgdata'],'dir'),    mkdir([mfile_path '\pgdata']);  end
save([mfile_path '\pgdata\pgdata_' mfile_name '.mat'],'pth');

%% update gui data
p.xmax=max(u.len)*(u.time{1}(2)-u.time{1}(1));

p.pth=pth;
p=get_parameter('set',p);

add_log([num2str(u.N_mol_org) ' traces loaded']);
set_prog(p.hdl_prog,1);

%% plot averaged traces
plot_avr_trace(u,p)
rm_junk(1,1);

function rm_junk(~,~)
p=get_parameter('fetch');
add_log('getting info for junks...');
set_prog(p.hdl_prog,0);

hdl_Rene_main=findobj('Tag','Rene_main');
u=hdl_Rene_main.UserData.u;

fr_max=min(u.len);

for mid=1:u.N_mol
    tot_i(mid,:)=u.donor_org{mid}(1:fr_max)+u.acceptor_org{mid}(1:fr_max);
end

% tot_i=( cell2mat(u.donor_org)+cell2mat(u.acceptor_org) )';
allan_corr_ana(tot_i);

function allan_corr_ana(totali)
p=get_parameter('fetch');
add_log('Calculating Allan deviation...');
set_prog(p.hdl_prog,0.1);

%% allan deviation
allan_tau=1:2:20;
N_allan_tau=length(allan_tau);

xbin=0:10:500;
hist_allani_sel=[];
for i=1:N_allan_tau
    allan_i_sel(:,i)=sqrt(.5*mean((totali(:,1:end-allan_tau(i))-totali(:,allan_tau(i)+1:end)).^2,2));
    [hist_allani_sel(:,i),hist_allani_sel_bin]=hist(allan_i_sel(:,i),xbin);
end

%% correlation series
add_log('Calculating time correlation...');
set_prog(p.hdl_prog,0.3);

corr_tau=[2 4 8 16 32 64 128 256 512 1024];
N_corr_tau=length(corr_tau);
hist_corr_i_sel_bin=-0.5:0.02:2;
g_tau_sel=[];
corr_i_sel=[];
hist_corr_i_sel=[];

if 1 % by no binning
    tmp=mean((totali(:,1:end-1).*totali(:,2:end)),2);
    normfactor=tmp./(mean(totali(:,1:end-1),2).*mean(totali(:,2:end),2))-1;
    for i=1:N_corr_tau
        tmp=mean((totali(:,1:end-corr_tau(i)).*totali(:,corr_tau(i)+1:end)),2);
        corr_i_sel(:,i)=tmp./(mean(totali(:,1:end-corr_tau(i)),2).*mean(totali(:,corr_tau(i)+1:end),2))-1;

        g_tau_sel(i)=mean(corr_i_sel(:,i));
        corr_i_sel(:,i)=corr_i_sel(:,i)./normfactor;

        [tmp1,~]=hist(corr_i_sel(:,i),hist_corr_i_sel_bin);
        hist_corr_i_sel(:,i)=tmp1;
        set_prog(p.hdl_prog,0.3+i/N_corr_tau*.6);
    end
else % by coarsening (NOTE: it turns out that coarsening does not improve the result.
    for i=1:N_corr_tau
        tmp_tr=[];
        for mid=1:size(totali,1)
            tmp_tr(mid,:)=binning(totali(mid,:),corr_tau(i),'coarsening');
        end
        
        tmp=mean(tmp_tr.^2,2);
        normfactor=tmp./(mean(tmp_tr,2).^2)-1;
        
        tmp=mean((tmp_tr(:,1:end-1).*tmp_tr(:,2:end)),2);
        corr_i_sel(:,i)=tmp./(mean(tmp_tr(:,1:end-1),2).*mean(tmp_tr(:,2:end),2))-1;

        g_tau_sel(i)=mean(corr_i_sel(:,i));
        corr_i_sel(:,i)=corr_i_sel(:,i)./normfactor;

        [tmp1,~]=hist(corr_i_sel(:,i),hist_corr_i_sel_bin);
        hist_corr_i_sel(:,i)=tmp1;
        set_prog(p.hdl_prog,0.3+i/N_corr_tau*.6);
    end
end

y_lim_set=[0 max(hist_corr_i_sel(:))];
g_tau_sel=g_tau_sel./g_tau_sel(1);

% plot the series
fhd_corr_allan=findobj('Tag','figure_corr_series');
if ~isempty(fhd_corr_allan)
    figure(fhd_corr_allan);clf;
else
    fhd_corr_allan=figure('Tag','figure_corr_series','Name','Corr. series','NumberTitle','off');
end


for i=1:N_allan_tau
    subplot(N_corr_tau,2,i*2);
    bar(hist_allani_sel_bin,hist_allani_sel(:,i),1,'FaceColor',[.7 0 0],'LineStyle','none');
    set(gca,'XLim',[xbin(1) xbin(end)],'XTickLabel',[],'YTickLabel',[]);    grid on;
    if i==1, title('Allan'); end
    text(xbin(end)-(xbin(end)-xbin(1))*0.25,max(hist_allani_sel(:,i))*0.5,['tau' num2str(allan_tau(i))]);
end
xlabel('Allan of total intensity');

for i=1:N_corr_tau
    subplot(N_allan_tau,2,i*2-1);
    
    bar(hist_corr_i_sel_bin,hist_corr_i_sel(:,i),1,'FaceColor',[1 0 0],'LineStyle','none');hold on;
    %     set(gca,'XLim',[0 3],'YLim',y_lim_set);
    set(gca,'XLim',[-0.25 1.5],'XTickLabel',[],'YTickLabel',[]);grid on;
    if i==1, title('t-corr.'); end
    %     if i~=N_corr_tau,        set(gca,'XTick',[]);    end
    text(1.2,max(hist_corr_i_sel(:,i))/2,['tau' num2str(corr_tau(i))]);
end
xlabel('Corr of total intensity');

%% show full time correlation function
if 0
    figure(221);clf
    subplot(1,2,1)
    semilogx(corr_tau,g_tau_sel,'r.-','MarkerSize',12);
    ylim([0 max(g_tau_sel*1.1)]);
    ylabel('Time correlation');
    xlabel('Tau (frame)');
    
    subplot(1,2,2)
    semilogx(corr_tau,corr_i_sel','r-');
    ylim([-1 100]);
    ylabel('Time correlation');
    xlabel('Tau (frame)');
end

%% show allan vs. T-corr
% plot the series
fhd_corr_allan=findobj('Tag','figure_corr_allan');
if ~isempty(fhd_corr_allan)
    figure(fhd_corr_allan);clf;
else
    fhd_corr_allan=figure('Tag','figure_corr_allan','Name','T-corr. vs. Allan','NumberTitle','off');
end
fhd_corr_allan.Units='Normalized';

for ci=1:N_corr_tau,  corr_tau_name{ci}=num2str(corr_tau(ci)); end
for ci=1:N_corr_tau,  allan_tau_name{ci}=num2str(allan_tau(ci)); end

ax_corr=axes('Units','Normalized','Position', [0.1 0.6 0.35 0.33]);
ax_allan=axes('Units','Normalized','Position', [0.6 0.6 0.35 0.33]);
ax_cvsa=axes('Units','Normalized','Position', [0.1 0.1 0.8 0.35]);

uicontrol('Style', 'text','String', 't-corr tau',...
    'Units','Normalized','Position', [0.2 0.95 0.1 0.04]);
uicontrol('Style', 'text','String', 'allan tau',...
    'Units','Normalized','Position', [0.7 0.95 0.1 0.04]);

uicontrol('Style', 'popup','String', corr_tau_name,'Tag','corr_tau_sel','Value',2,...
    'Units','Normalized','Position', [0.3 0.95 0.1 0.05],...
    'Callback', {@show_corr_vs_allan,fhd_corr_allan,ax_corr,ax_allan,ax_cvsa,corr_i_sel,allan_i_sel,hist_corr_i_sel,hist_allani_sel,hist_corr_i_sel_bin,hist_allani_sel_bin,corr_tau,allan_tau});
uicontrol('Style', 'popup','String', allan_tau_name,'Tag','allan_tau_sel',...
    'Units','Normalized','Position', [0.8 0.95 0.1 0.05],...
    'Callback', {@show_corr_vs_allan,fhd_corr_allan,ax_corr,ax_allan,ax_cvsa,corr_i_sel,allan_i_sel,hist_corr_i_sel,hist_allani_sel,hist_corr_i_sel_bin,hist_allani_sel_bin,corr_tau,allan_tau});
uicontrol('Style','pushbutton','String','Limits',...
    'Units','Normalized','Position',[0.905 0.3 0.09 0.08],...
    'Callback',{@show_corr_vs_allan,fhd_corr_allan,ax_corr,ax_allan,ax_cvsa,corr_i_sel,allan_i_sel,hist_corr_i_sel,hist_allani_sel,hist_corr_i_sel_bin,hist_allani_sel_bin,corr_tau,allan_tau,1})

show_corr_vs_allan('','',fhd_corr_allan,ax_corr,ax_allan,ax_cvsa,corr_i_sel,allan_i_sel,hist_corr_i_sel,hist_allani_sel,hist_corr_i_sel_bin,hist_allani_sel_bin,corr_tau,allan_tau);

add_log('job done!')
set_prog(p.hdl_prog,1);

function show_corr_vs_allan(~,~,~,ax_corr,ax_allan,ax_cvsa,corr_i_sel,allan_i_sel,hist_corr_i_sel,hist_allani_sel,hist_corr_i_sel_bin,hist_allani_sel_bin,corr_tau,allan_tau,flg)
if nargin<15
    flg=0;
end
%% fetch the tau selection info
tmp_hdl=findobj('Tag','corr_tau_sel');
selected_tau_corr_id=tmp_hdl.Value;
rep_corr_sel=corr_i_sel(:,selected_tau_corr_id);

tmp_hdl=findobj('Tag','allan_tau_sel');
selected_tau_allan_id=tmp_hdl.Value;
rep_allan_sel=allan_i_sel(:,selected_tau_allan_id);

%% plot
% axis(ax_cvsa);
bar(ax_corr,hist_corr_i_sel_bin,hist_corr_i_sel(:,selected_tau_corr_id),1,'r');
ax_corr.XLim=[-0.2 1.5];
ax_corr.XLabel.String=['G (tau=' num2str(corr_tau(selected_tau_corr_id)) ')'];
ax_corr.XGrid='on';
ax_corr.YGrid='on';

bar(ax_allan,hist_allani_sel_bin,hist_allani_sel(:,selected_tau_allan_id),1,'k');
ax_allan.XLim=[0 5000];
ax_allan.XLabel.String='Allan deviation';
ax_allan.XGrid='on';
ax_allan.YGrid='on';

%% choose rep_corr_sel
tmp_hdl=findobj('Tag','corr_tau_sel');
selected_tau_corr_id=tmp_hdl.Value;
rep_corr_sel=corr_i_sel(:,selected_tau_corr_id);

tmp_hdl=findobj('Tag','allan_tau_sel');
selected_tau_allan_id=tmp_hdl.Value;
rep_allan_sel=allan_i_sel(:,selected_tau_allan_id);

%%
if ~flg
    corr_range=[-Inf Inf];
    allan_range=[-Inf Inf];
    
    vid1=rep_corr_sel > corr_range(1) & rep_corr_sel < corr_range(2);
    vid2=rep_allan_sel > allan_range(1) & rep_allan_sel < allan_range(2);
    vid=vid1&vid2;
    
    corr_cert= rep_corr_sel(vid);
    corr_fall=rep_corr_sel(~vid);
    allan_cert= rep_allan_sel(vid);
    allan_fall=rep_allan_sel(~vid);
    
    axis(ax_cvsa);
    plot(ax_cvsa,corr_cert,allan_cert,'r.','MarkerSize',12);hold on;
    plot(ax_cvsa,corr_fall,allan_fall,'.','Color',[.6 .6 .6],'MarkerSize',12);hold off;
    
    xlabel(['Correlation (tau=' num2str(corr_tau(selected_tau_corr_id)) ')']);
    ylabel(['Allan (tau=' num2str(allan_tau(selected_tau_allan_id)) ')']);
    xlim([-.2 1.2]);
    grid on;
    text(ax_cvsa.XLim(2)*0.7,ax_cvsa.YLim(2)*0.9,['t-corr (' num2str(corr_range(1),'%.2f') ', ' num2str(corr_range(2),'%.2f') ')']);
    text(ax_cvsa.XLim(2)*0.7,ax_cvsa.YLim(2)*0.8,['Allan (' num2str(allan_range(1),'%.0f') ', ' num2str(allan_range(2),'%.0f') ')']);
    
else
    corr_range=[-Inf Inf];
    allan_range=[-Inf Inf];
    while 1
        vid1=rep_corr_sel > corr_range(1) & rep_corr_sel < corr_range(2);
        vid2=rep_allan_sel > allan_range(1) & rep_allan_sel < allan_range(2);
        vid=vid1&vid2;
        
        corr_cert= rep_corr_sel(vid);
        corr_fall=rep_corr_sel(~vid);
        allan_cert= rep_allan_sel(vid);
        allan_fall=rep_allan_sel(~vid);
        
        axis(ax_cvsa);cla;
        plot(ax_cvsa,corr_cert,allan_cert,'r.','MarkerSize',12);hold on;
        plot(ax_cvsa,corr_fall,allan_fall,'.','Color',[.6 .6 .6],'MarkerSize',12);hold off;
        xlim([-.2 1.2]);
        xlabel(['Correlation (tau=' num2str(corr_tau(selected_tau_corr_id)) ')']);
        ylabel(['Allan (tau=' num2str(allan_tau(selected_tau_allan_id)) ')']);
        grid on;
        text(ax_cvsa.XLim(2)-(ax_cvsa.XLim(2)-ax_cvsa.XLim(1))*0.2,ax_cvsa.YLim(2)-(ax_cvsa.YLim(2)-ax_cvsa.YLim(1))*0.1,...
            ['t-corr (' num2str(corr_range(1),'%.2f') ', ' num2str(corr_range(2),'%.2f') ')']);
        text(ax_cvsa.XLim(2)-(ax_cvsa.XLim(2)-ax_cvsa.XLim(1))*0.2,ax_cvsa.YLim(2)-(ax_cvsa.YLim(2)-ax_cvsa.YLim(1))*0.2,...
            ['Allan (' num2str(allan_range(1),'%.0f') ', ' num2str(allan_range(2),'%.0f') ')']);
        text(ax_cvsa.XLim(2)-(ax_cvsa.XLim(2)-ax_cvsa.XLim(1))*0.2,ax_cvsa.YLim(2)-(ax_cvsa.YLim(2)-ax_cvsa.YLim(1))*0.3,...
            [num2str(sum(vid)) ' molecules selected']);
        
        [mx,my,mbutton]=ginput(1);
        if mx > ax_cvsa.XLim(2)
            mx=Inf;
        elseif mx < ax_cvsa.XLim(1)
            mx=-Inf;
        end
        if my > ax_cvsa.YLim(2)
            my=Inf;
        elseif my < ax_cvsa.YLim(1)
            my=-Inf;
        end
        switch mbutton
            case 1
                allan_range(1)=my;
                corr_range(1)=mx;
            case 2
                break;
            case 3
                allan_range(2)=my;
                corr_range(2)=mx;
        end
    end
    
    %% find junks
    hdl_Rene_main=findobj('Tag','Rene_main');
    u=hdl_Rene_main.UserData.u;
    
    u.allan_range=[allan_range(1) allan_range(2)];
    u.allan_tau=allan_tau(selected_tau_allan_id);
    
    u.t_corr_range=[corr_range(1) corr_range(2)];
    u.t_corr_tau=corr_tau(selected_tau_corr_id);
    
    for mid=1:u.N_mol_org
        % reject by Allan
        totali=u.donor_org{mid}+u.acceptor_org{mid};
        allan_dev(mid)=sqrt(.5*mean((totali(1:end-u.allan_tau)-totali(u.allan_tau+1:end)).^2));
        
        % reject by T-corr
        tmp=mean((totali(1:end-u.t_corr_tau).*totali(u.t_corr_tau+1:end)));
        t_corr_amp(mid)=tmp./(mean(totali(1:end-u.t_corr_tau)).*mean(totali(u.t_corr_tau+1:end)))-1;
        
        tmp=mean((totali(1:end-1).*totali(2:end)));
        normfactor=tmp./(mean(totali(1:end-1)).*mean(totali(2:end)))-1;
        
        t_corr(mid)=t_corr_amp(mid)./normfactor;
    end
    valid_id(1,:) = allan_dev > u.allan_range(1) & allan_dev < u.allan_range(2);
    valid_id(2,:) = t_corr > u.t_corr_range(1) & t_corr < u.t_corr_range(2);
    u.isjunk= ~(valid_id(1,:) & valid_id(2,:));
    
    tmp=1:u.N_mol_org;
    u.mol_id_of_org=tmp(~u.isjunk);
    u.donor=u.donor_org(~u.isjunk);
    u.acceptor=u.acceptor_org(~u.isjunk);
    u.time=u.time_org(~u.isjunk);
    u.N_mol=sum(~u.isjunk);
    
    %% Update parameters
    hdl_Rene_main.UserData.u=u;
    
end

function tr_bin=binning(tr_org,N_bin,binning_mode)
binlen=floor(length(tr_org)/N_bin);
tr_bin=zeros(binlen,1);

switch binning_mode
    case 'coarsening'
        for m=1:binlen
            for mm=0:N_bin-1
                tr_bin(m)=tr_bin(m)+tr_org(N_bin*m-mm);
            end
            tr_bin(m)=tr_bin(m)/N_bin;
        end
    case 'median'
        for m=1:binlen
            tr_bin(m)=median(tr_org(N_bin*(m-1)+1:N_bin*m));
        end
end

function find_binding_events(hObject,eventdata)
p=get_parameter('fetch');

add_log('working on... ');
set_prog(p.hdl_prog,0);

tic
%% fetch the data
main_fhd=findobj('Tag','Rene_main');
if ~isfield(main_fhd.UserData,'u')
    load_data(hObject, eventdata)
end
u=main_fhd.UserData.u;
u.N_mol=u.N_mol_org;

%% select non-junk molecules
add_log('sort junks... ');
set_prog(p.hdl_prog,0.05);

if ~isfield(u,'isjunk')
    u.isjunk=zeros(1,u.N_mol);
end
u.donor=u.donor_org(~u.isjunk);
u.acceptor=u.acceptor_org(~u.isjunk);
u.time=u.time_org(~u.isjunk);
u.N_mol=sum(~u.isjunk);

%% make time vector
u.t_tot_obs=0;
for mid=1:u.N_mol
    timeunit=u.time{mid}(2)-u.time{mid}(1);
    u.t_tot_obs=u.t_tot_obs+u.len(mid)*timeunit;
end

%% binning
add_log('binning... ');
set_prog(p.hdl_prog,0.1);

if p.N_bin ~= 0
    for mid=1:u.N_mol
        set_prog(p.hdl_prog,0.1+mid/u.N_mol*.9);
        u.binD{mid}=binning(u.donor{mid},p.N_bin,'median');
        u.binA{mid}=binning(u.acceptor{mid},p.N_bin,'median');
        timeunit=u.time{mid}(2)-u.time{mid}(1);
        u.bintime{mid}=(1:floor(u.len(mid)/p.N_bin))*timeunit*p.N_bin;
    end
end

%% find state
add_log('finding binding events... ');
set_prog(p.hdl_prog,0.2);
u=find_state(u,p);

%% remove molecules having too small or many peaks 
% This option is reserved for future use
rm_by_loudeness=0;
if rm_by_loudeness
    min_bind_event=10;
    x_sigma=3;
    
    for mid=1:u.N_mol
        N_bind(mid)=sum(diff(u.stateI{mid})==1);
        N_idss(mid)=sum(diff(u.stateI{mid})==1);
    end
    
    tmp_gfit=fitgmdist(N_bind',1);
    vid=N_bind > min_bind_event & N_bind < tmp_gfit.mu+sqrt(tmp_gfit.Sigma)*x_sigma;
    
    u.N_mol=sum(vid);
    u.time=u.time(vid);
    u.donor=u.donor(vid);
    u.acceptor=u.acceptor(vid);
    u.stateI=u.stateI(vid);
    
    if p.N_bin ~= 0
        u.binD=u.binD(vid);
        u.binA=u.binA(vid);
        u.bintime=u.bintime(vid);
    end
end

%% BG_correction
if p.auto_BG
    add_log('Background correction ... ');
    set_prog(p.hdl_prog,0.9);
    for mid=1:u.N_mol    
        %% for the binned data
        if p.N_bin ~= 0
            BG_d(mid)=median(u.binD{mid}(~u.stateI{mid}));
            BG_a(mid)=median(u.binA{mid}(~u.stateI{mid}));
            u.binD{mid}=u.binD{mid}-BG_d(mid);
            u.binA{mid}=u.binA{mid}-BG_a(mid);
        else
            BG_d(mid)=median(u.donor{mid}(~u.stateI{mid}));
            BG_a(mid)=median(u.acceptor{mid}(~u.stateI{mid}));

        end
        u.donor{mid}=u.donor{mid}-BG_d(mid);
        u.acceptor{mid}=u.acceptor{mid}-BG_a(mid);
    end
    u.BG_d_ind=BG_d;
    u.BG_a_ind=BG_a;
else
    for mid=1:u.N_mol
        u.donor{mid}=u.donor{mid}-p.BG_d;
        u.acceptor{mid}=u.acceptor{mid}-p.BG_a;
        %% for the binned data
        if p.N_bin ~= 0
            u.binD{mid}=u.binD{mid}-p.BG_d;
            u.binA{mid}=u.binA{mid}-p.BG_a;
        end
    end
end

%% Leakage & gamma correction
for mid=1:u.N_mol
    u.acceptor{mid}=u.acceptor{mid}*p.gamma;

    u.donor{mid}=u.donor{mid}+u.donor{mid}*p.leakage;
    u.acceptor{mid}=u.acceptor{mid} - u.donor{mid}*p.leakage;   
    
    %% for the binned data
    if p.N_bin ~= 0
        u.binA{mid}=u.binA{mid}*p.gamma;
        
        u.binD{mid}=u.binD{mid} + u.binD{mid} *p.leakage;
        u.binA{mid}=u.binA{mid} - u.binD{mid} *p.leakage;
    end
end

%% calc FRET
add_log('Calculating FRET efficiency... ');
set_prog(p.hdl_prog,0.95);
for mid=1:u.N_mol
    %% calc FRET
    u.totali{mid}=u.donor{mid}+u.acceptor{mid};
    u.FRETe{mid}=u.acceptor{mid}./u.totali{mid};
    
    %% binning
    if p.N_bin ~= 0
        u.binI{mid}=u.binD{mid}+u.binA{mid};
        u.binE{mid}=u.binA{mid}./u.binI{mid};
    end
end


%% split peaks by FRET change
if 1,     u=split_peaks_by_FRET(p,u);  end


%% update gui data
main_fhd.UserData.u=u;

add_log([num2str(u.N_mol) ' traces processed in ' num2str(toc,'%.2f') 's']);
set_prog(p.hdl_prog,1);

%% auto run next steps
tmp_hdl=findobj('Tag','RM_flg2go_ana_dwells');
if tmp_hdl.Value
    ana_dwells(1,1);
end

function u=split_peaks_by_FRET(p,u)
%%
% p.E_tolerance=0.1;
for mid=1:u.N_mol
    
    %%
    if isfield(u,'binD')
%         tmp_d=u.binD{mid};
%         tmp_a=u.binA{mid};
%         tmp_i=u.binI{mid};
        tmp_e=u.binE{mid};
    else
%         tmp_d=u.donor{mid};
%         tmp_a=u.acceptor{mid};
%         tmp_i=u.totali{mid};
        tmp_e=u.FRETe{mid};
    end
    tmp_e_sel=tmp_e.*u.stateI{mid};
    trans_pos=abs(diff(tmp_e_sel))>p.E_tolerance;
    
    new_stateI{mid}=u.stateI{mid};
    new_stateI{mid}(trans_pos)=false;
    
end
u.stateI=new_stateI;

function u=find_state(u,p)
for mid=1:u.N_mol
    if p.N_bin == 0
        temp_tr{mid}=u.donor{mid}+u.acceptor{mid};
    else
        temp_tr{mid}=u.binD{mid}+u.binA{mid};
    end
end
switch p.state_find_algorithm
    case 'Threshold'
        %% give the threshold
        for mid=1:u.N_mol,        u.stateI{mid}=temp_tr{mid}>p.auto_th;    end        
    case 'k-means'
%         t1=tic;
%         telap=0;
        for mid=1:u.N_mol
%             t_est=telap/(mid-1)*(u.N_mol-mid);
%             add_log(['k-means on int. trace of mol ' num2str(mid) '/' num2str(u.N_mol) ' (time elapsed: '  num2str(floor(telap/60)) 'm ' num2str(rem(telap,60),'%.0f') 's, remaining: ' num2str(floor(t_est/60)) 'm ' num2str(rem(t_est,60),'%.0f') 's)']);
            set_prog(p.hdl_prog,.1+( mid/u.N_mol*.8));
            
            u.stateI{mid}=find_kmeans_state(temp_tr{mid},p);
%             telap=toc(t1);
        end
end

function stateI=find_kmeans_state(seq,p)
%%
N_state=2;
[k_mean_state,k_mean_rep,sumd,D] = kmeans(seq,N_state);

[~,min_state_id]=min(k_mean_rep);
stateI=k_mean_state~=min_state_id;

%% draw k-mean result
hdl_update_figures=findobj('Tag','RM_update_figures');
if hdl_update_figures.Value==1
    for sti=1:N_state
        vid=k_mean_state==sti;
        k_mean_tr(vid)=k_mean_rep(sti);
        k_mean_tr_ind{sti}=seq(vid);
    end

    fhd_kmean_tmp=findobj('Tag','RM_HMM');
    if isempty(fhd_kmean_tmp)
        fhd_kmean_tmp=figure('Tag','RM_HMM','Numbertitle','off');
    end
    set(groot,'CurrentFigure',fhd_kmean_tmp);clf;
    fhd_kmean_tmp.Name='k-means';
    
    subplot(3,2,1:2);box on;hold on;

    color_scheme=[.7 .7 .7; 1 .6 0];
    xmin=min(seq);
	xmax=max(seq);
    xbin=xmin:(xmax-xmin)/50:xmax;
    
    for sti=1:N_state
        [tmp_hist,xbin]=hist(k_mean_tr_ind{sti},xbin);
        bar(xbin,tmp_hist,1,'FaceColor',color_scheme(sti,:));
        ymax=max(tmp_hist)*1.2;
        plot(k_mean_rep(sti),[ymax ymax],'r*','MarkerSize',12,'LineWidth',2);
    end

    
    subplot(3,2,3:6);hold on;
    plot(seq,'Color',[.7 .7 .7]);
    plot(k_mean_tr,'Color',[0 0 0]);
    ylim([p.y_min p.y_max]);
    drawnow
end

function ana_dwells(~,~)
p=get_parameter('fetch');

add_log('getting dwell info... ');
set_prog(p.hdl_prog,0);

%% fetch the FRET data
main_fhd=findobj('Tag','Rene_main');
u=main_fhd.UserData.u;

%% collect dwell times
[u,dw]=collect_dwell(u,p);
dw=select_longer_dw(dw,p);

%% make dwell and FRET hist
dw=make_hist(dw,p,u);

%% update gui data
main_fhd.UserData.u=u;
main_fhd.UserData.dw=dw;

add_log([num2str(sum(dw.N_up_survived_ind)) ' peaks selected from ' num2str(sum(dw.N_up_ind)) ' peaks']);
set_prog(p.hdl_prog,1);

tmp_hdl=findobj('Tag','RM_flg2go_show_dwells');
if tmp_hdl.Value
    show_dwells(1,1);
end

function [u,dw]=collect_dwell(u,p)
%%
hdl_log_text=findobj('Tag','RM_log_text');
for mid=1:u.N_mol
    add_log('Collecting dwell time info... ',hdl_log_text);
    if rem(mid,10)==0
        set_prog(p.hdl_prog,mid/u.N_mol*.7);
    end
    if p.N_bin==0
        E_inuse=u.FRETe{mid};
        I_inuse=u.totali{mid};
        timeunit_inuse=u.time{mid}(2)-u.time{mid}(1);
    else
        E_inuse=u.binE{mid};
        I_inuse=u.binI{mid};
        timeunit_inuse=u.bintime{mid}(2)-u.bintime{mid}(1);
    end
    stateI=u.stateI{mid};
    
    %% determine dwell time
    tmp_dwell=1;
    tmpE=E_inuse(1);
    tmpI=I_inuse(1);
    tmp_dw=[];
    tmp_up=[];
    N_peak_per_mol=0;
    Peak_E_raw={};
    Peak_E_avr=[];
    Peak_I_raw={};
    Peak_I_avr=[];
    Peak_t=[];
    
    t_stamp=1;
    for tti=2:length(stateI)
        if stateI(tti)==stateI(tti-1)
            tmp_dwell=tmp_dwell+1;
            tmpE=[tmpE E_inuse(tti)];
            tmpI=[tmpI I_inuse(tti)];
        else
            tmp_dwell=tmp_dwell*timeunit_inuse;
            %                     if N_bin ~= 0,    tmp_dwell=tmp_dwell*N_bin;    end
            if stateI(tti-1)==0
                tmp_dw=[tmp_dw tmp_dwell];
                t_stamp=tti;
            else
                tmp_up=[tmp_up tmp_dwell];
                N_peak_per_mol=N_peak_per_mol+1;
                Peak_E_raw{N_peak_per_mol}=tmpE;
                Peak_E_avr(N_peak_per_mol)=mean(tmpE);
                Peak_I_raw{N_peak_per_mol}=tmpI;
                Peak_I_avr(N_peak_per_mol)=mean(tmpI);
                Peak_t(N_peak_per_mol,1:2)=[t_stamp tti-1]*timeunit_inuse;
                t_stamp=tti;
            end
            tmpI=I_inuse(tti);
            tmpE=E_inuse(tti);
            tmp_dwell=1;
        end
    end
    
    dw.dw_ind{mid}=tmp_dw;
    dw.up_ind{mid}=tmp_up;
    dw.N_dw_ind(mid)=length(tmp_dw);
    dw.N_up_ind(mid)=length(tmp_up);

    dw.up_peak_t{mid}=Peak_t;
    
    dw.up_peak_E_raw{mid}=Peak_E_raw;
    dw.up_peak_E_avr{mid}=Peak_E_avr;
    
    dw.up_peak_I_raw{mid}=Peak_I_raw;
    dw.up_peak_I_avr{mid}=Peak_I_avr;
end

function dw=select_longer_dw(dw,p)
%% collect dwells longer than N data points
N_mol=length(dw.up_ind);

dw.up_survived_ind={};
dw.N_up_survived_ind=[];
dw.up_peak_E_raw_survived_ind={};
dw.up_peak_I_raw_survived_ind={};

for mid=1:N_mol
    add_log('Removing short dwells');
    if rem(mid,25)==0
        set_prog(p.hdl_prog,mid/N_mol*.28+.7);
    end
    
    % remove peaks with too short
    vid_up1=dw.up_ind{mid}>=p.min_dw_len_sel;
    
	% remove peaks by intensity    
    vid_up2=[];
    for pki=1:length(dw.up_peak_I_raw{mid})
        c_int=mean(dw.up_peak_I_raw{mid}{pki});
        if c_int > p.peak_int(1) && c_int < p.peak_int(2)
            vid_up2(pki) = 1 ;
        else
            vid_up2(pki) = 0;
        end
    end
    
    % remove peaks with too much std
    vid_up3=[];
    for pki=1:length(dw.up_peak_E_raw{mid})
        if length(dw.up_peak_E_raw{mid}{pki})>1
            vid_up3(pki) = std(dw.up_peak_E_raw{mid}{pki}) < 0.1;
        else
            vid_up3(pki) = 1;
        end
    end

    % remove peaks with strange intensity
    vid_up4=[];
    for pki=1:length(dw.up_peak_I_raw{mid})
        if length(dw.up_peak_I_raw{mid}{pki})>1
            tmp_I_allan=sqrt(.5*median((dw.up_peak_I_raw{mid}{pki}(:,1:end-1)-dw.up_peak_I_raw{mid}{pki}(:,2:end)).^2,2));
            max_diff=max(dw.up_peak_I_raw{mid}{pki})-min(dw.up_peak_I_raw{mid}{pki});
            vid_up4(pki) = max_diff < tmp_I_allan*12;
        else
            vid_up4(pki) = 1;
        end
    end

    vid_up=vid_up1 & vid_up2;
    
    dw.up_survived_ind{mid}=dw.up_ind{mid}(vid_up);
    dw.N_up_survived_ind(mid)=sum(vid_up);
    dw.up_peak_E_raw_survived_ind{mid}=dw.up_peak_E_raw{mid}(vid_up);
    dw.up_peak_I_raw_survived_ind{mid}=dw.up_peak_I_raw{mid}(vid_up);
    dw.up_peak_t_survived{mid}=dw.up_peak_t{mid}(vid_up,:);
end

function dw=make_hist(dw,p,u)
%% build FRET and I hist
E_raw_selected_all=[];
E_avr_selected_all=[];
I_raw_selected_all=[];
I_avr_selected_all=[];
D_selected_all=[];

for mid=1:u.N_mol
    tmp_E_raw=[];
    tmp_E_avr=[];
    tmp_I_raw=[];
    tmp_I_avr=[];
    tmp_D=[];

    for pki=1:length(dw.up_peak_E_raw_survived_ind{mid})
        tmp_E_raw=[tmp_E_raw dw.up_peak_E_raw_survived_ind{mid}{pki}];
        tmp_I_raw=[tmp_I_raw dw.up_peak_I_raw_survived_ind{mid}{pki}];
        
        tmp_E_avr=[tmp_E_avr mean(dw.up_peak_E_raw_survived_ind{mid}{pki})];
        tmp_I_avr=[tmp_I_avr mean(dw.up_peak_I_raw_survived_ind{mid}{pki})];
        
        tmp_D=[tmp_D length(dw.up_peak_I_raw_survived_ind{mid}{pki})];
    end
    
    E_raw_selected_all=[E_raw_selected_all tmp_E_raw];
    I_raw_selected_all=[I_raw_selected_all tmp_I_raw];
    E_avr_selected_all=[E_avr_selected_all tmp_E_avr];
    I_avr_selected_all=[I_avr_selected_all tmp_I_avr];
    
    D_selected_all=[D_selected_all tmp_D];
    
end

% build E hist
E_bin_fine=-0.1:0.01:1.1;
[Ehist1,~]=hist(E_raw_selected_all,E_bin_fine);
[Ehist_peak_avr,~]=hist(E_avr_selected_all,E_bin_fine);

dw.FRET_hist(:,1)=E_bin_fine;
dw.FRET_hist(:,2)=Ehist1;
dw.FRET_hist(:,3)=Ehist_peak_avr;

%% build I-E hist
dw.edges{1}=floor(min(I_avr_selected_all)*0.9):20:floor(max(I_avr_selected_all)*1.1);
dw.edges{2}=-0.5:0.02:1.5;
dw.IE_2Dhist=hist3([I_avr_selected_all; E_avr_selected_all]','Edges',dw.edges);
dw.IE_2Dhist_Ionly=hist(I_avr_selected_all,dw.edges{1});
dw.IE_2Dhist_Eonly=hist(E_avr_selected_all,dw.edges{2});

%% build D-E hist
log_D=log(D_selected_all);

dw.DE_edges{1}=0:log(2):floor(max(log_D)*1.1);
dw.DE_edges{2}=-0.5:0.02:1.5;
dw.DE_2Dhist=hist3([log_D; E_avr_selected_all]','Edges',dw.DE_edges);
dw.DE_2Dhist_Donly=hist(log_D,dw.DE_edges{1});
dw.DE_2Dhist_Eonly=hist(E_avr_selected_all,dw.DE_edges{2});

%% build the dwell time histogram
dw_all=[];
up_all=[];
for mid=1:u.N_mol
    dw_all=[dw_all dw.dw_ind{mid}];
    up_all=[up_all dw.up_survived_ind{mid}];
end
time_unit=u.time{1}(2)-u.time{1}(1);

% calc hist
bin_step=floor(mean(dw_all)/2/time_unit)*time_unit;
xbin=min(dw_all):bin_step:mean(dw_all)*10;
[down_hist, downX]=hist(dw_all,xbin);
dw.down_hist=down_hist(2:end-1);
dw.downX=downX(2:end-1);

bin_step=floor(mean(up_all)/4/time_unit)*time_unit;
xbin=min(up_all):bin_step:mean(up_all)*6;
[up_hist, upX]=hist(up_all,20);
dw.up_hist=up_hist(2:end-1);
dw.upX=upX(2:end-1);

% single exp fit
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0],...
    'Upper',[Inf Inf],...
    'StartPoint',[dw.down_hist(1),mean(dw_all)]);
ft = fittype('a*exp(-x/tau)','options',fo);
[dw.down_fit,~] = fit(dw.downX',dw.down_hist',ft);

fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0],...
    'Upper',[Inf Inf],...
    'StartPoint',[dw.up_hist(1),mean(up_all)]);
[dw.up_fit,~] = fit(dw.upX',dw.up_hist',ft);

function show_idv_mols(~,~)
%% fetch the FRET data
main_fhd=findobj('Tag','Rene_main');
u=main_fhd.UserData.u;
p=get_parameter('fetch');

%% check if peaks are detected
if isfield(main_fhd.UserData,'dw')
    dw=main_fhd.UserData.dw;
else
    dw={};
end

%% show averaged properties
plot_avr_trace(u,p);

%% run over individual traces
mid=1;
while 1
    %% check and update parameters
    p=get_parameter('fetch');
    
    %% check if mid is  valid
    if mid<1
        mid=1;
    elseif mid>u.N_mol
        choice = questdlg('There is no more molecule.', ...
            'Choose',  'Finish','Go to the first','Go to the last', 'Finish');
        switch choice
            case 'Finish'
                break;
            case 'Go to the first'
                mid=1;
            case 'Go to the last'
                mid=u.N_mol;
        end
    end
    
    %% plot traces
    plot_trace(u,p,mid,dw);
    u_ans=input('[p: previous, g: go to, e: end] ','s');
    switch u_ans
        case 'p'
            mid=mid-2;
        case 'e'
            break;
        case 'g'
            mid2go=input('type molecule id... ');
            mid=mid2go-1;
        case 'r'
            u.isjunk(u.mol_id_of_org(mid))=1;

            vid=true(1,u.N_mol);
            vid(mid)=false;
            
            u.donor=u.donor(vid);
            u.acceptor=u.acceptor(vid);
            u.time=u.time(vid);
            u.mol_id_of_org=u.mol_id_of_org(vid);
            if isfield(u,'stateI')
                u.stateI=u.stateI(vid);
            end
            if isfield(u,'FRETe')
                u.totali=u.totali(vid);
                u.FRETe=u.FRETe(vid);
            end
            if isfield(u,'BG_d_ind')
                u.BG_d_ind=u.BG_d_ind(vid);
                u.BG_a_ind=u.BG_a_ind(vid);
            end            
            
            if isfield(dw,'N_up_survived_ind')
                dw.N_up_survived_ind=dw.N_up_survived_ind(vid);
                dw.up_peak_E_raw_survived_ind=dw.up_peak_E_raw_survived_ind(vid);
                dw.up_peak_I_raw_survived_ind=dw.up_peak_I_raw_survived_ind(vid);
                dw.up_peak_t_survived=dw.up_peak_t_survived(vid);
            end
            
            u.N_mol=u.N_mol-1;
            main_fhd.UserData.u=u;
            
            mid=mid-1;
        case ''
            % do nothing for next mol (DO NOT REMOVE THIS CASE!!!)
        otherwise
            mid=mid-1;
    end
    mid=mid+1;
end

function plot_avr_trace(u,p)
%% make averaged traces
fr_max=min(u.len);
D_avr=zeros(fr_max,1);
A_avr=zeros(fr_max,1);
for mid=1:u.N_mol
    D_avr=D_avr+u.donor{mid}(1:fr_max);
    A_avr=A_avr+u.acceptor{mid}(1:fr_max);
end
%%
D_avr=D_avr./u.N_mol;
A_avr=A_avr./u.N_mol;
T_avr=D_avr+A_avr;
E_avr=A_avr./T_avr;

%% plot
fhd_trace=findobj('Tag','figure_trace_avr');
if ~isempty(fhd_trace)
    figure(fhd_trace);clf;
else
    fhd_trace=figure('Tag','figure_trace_avr','Name','Avr. trace','NumberTitle','off');
end

I_max=max(T_avr);
I_min=min([min(T_avr) min(D_avr) min(A_avr)]);
I_range=[I_min-(I_max-I_min)/10 I_max*1.1];

tmp_tvct=(1:fr_max)*p.time_unit;

% total int trace
subplot('Position',[.15 .2 .73 .70]);
plot(tmp_tvct,T_avr,'Color',[.5 .5 .5]);hold on;
plot(tmp_tvct,D_avr,'Color',[.5 1 .5]);
plot(tmp_tvct,A_avr,'Color',[1 .5 .5]);hold off;
set(gca,'FontSize',12,'XLim',[0 p.xmax],'YLim',I_range,'Layer','top');
ylabel('Avr. intensity (a.u.)');
xlabel('Time (s)');
title(p.pth);
grid on;

% total int hist
subplot('Position',[.89 .2 .09 .78]);
[tmp_hist,tmp_hist_bin]=hist(T_avr,50);
bar(tmp_hist_bin,tmp_hist,1,'FaceColor',[.4 .4 .4],'LineStyle','none');hold on;
[tmp_hist,tmp_hist_bin]=hist(D_avr,50);
bar(tmp_hist_bin,tmp_hist,1,'FaceColor',[.4 1 .4],'LineStyle','none');
[tmp_hist,tmp_hist_bin]=hist(A_avr,50);
bar(tmp_hist_bin,tmp_hist,1,'FaceColor',[1 .4 .4],'LineStyle','none');

set(gca,'FontSize',12,'YTickLabel',[],'XTickLabel',[],'XLim',I_range,'Layer','top');
view([90 270]);
grid on;

function plot_trace(u,p,mid,dw)

if p.N_bin==0
    del_t=p.time_unit;
else
	del_t=p.time_unit*p.N_bin;
end

% check if peaks are found
if isfield(u,'totali')
    peaks_found=1;
else
    peaks_found=0;
end

% check if dwell time analyzed
if isempty(dw)
    peaks_selected=0;
else
    peaks_selected=1;   
end

%% make selected dwell traces
if peaks_selected
    selected_dw_tr_I=dw.up_peak_I_raw_survived_ind{mid};
    selected_dw_tr_E=dw.up_peak_E_raw_survived_ind{mid};
    selected_dw_tr_t=dw.up_peak_t_survived{mid};
end


%% get focus on the trace figure
fhd_trace=findobj('Tag','figure_trace');
if ~isempty(fhd_trace)
    figure(fhd_trace);clf;
else
    fhd_trace=figure('Tag','figure_trace','Name','Idv. trace','NumberTitle','off');
end

%% draw total intensity
subplot('Position',[0.11 0.68 0.75 0.27]);
if ~peaks_found
    plot(u.time{mid},u.donor{mid}+u.acceptor{mid},'Color',[.5 .5 .5]);hold on;
else
    plot(u.time{mid},u.totali{mid},'Color',[.9 .9 .9]);hold on;
    if p.N_bin ==0
        tmp_out=u.totali{mid};
        tmp_out(~u.stateI{mid})=nan;
        plot(u.time{mid},tmp_out,'.-','Color',[0 0 0]);
    else
        plot(u.bintime{mid},u.binI{mid}, 'Color',[.7 .7 .7] ,'LineWidth',1);
        tmp_out=u.binI{mid};
        tmp_out(~u.stateI{mid})=nan;
        plot(u.bintime{mid},tmp_out,'.-', 'Color',[.2 .2 .2] ,'LineWidth',1);
    end

    if strcmp(p.state_find_algorithm,'Threshold')
        if isfield(u,'BG_d_ind')
            tmp_BG_d=u.BG_d_ind(mid);
            tmp_BG_a=u.BG_a_ind(mid);
        else
            tmp_BG_d=0;
            tmp_BG_a=0;
        end
        plot(u.time{mid},u.time{mid}*0+p.auto_th-tmp_BG_d-tmp_BG_a,'r');
    end
end
if peaks_selected
    for pki=1:length(selected_dw_tr_I)
        plot((selected_dw_tr_t(pki,1):del_t:selected_dw_tr_t(pki,2))-p.time_unit,selected_dw_tr_I{pki},'r');
    end
end
hold off;
ylabel('Intensity (a.u.)');
set(gca,'FontSize',12,'XTickLabel',[],'XLim',[0 p.xmax],'YLim',[p.y_min p.y_max],'Layer','top');
title([u.filelist(mid).name '/' num2str(u.N_mol)]);
grid on;

%% draw total int hist
subplot('Position',[0.88 0.68 0.1 0.27]);
if peaks_found
    [tmp_hist,tmp_hist_bin]=hist(u.totali{mid},50);
else
    [tmp_hist,tmp_hist_bin]=hist(u.donor{mid}+u.acceptor{mid},50);
end
bar(tmp_hist_bin,tmp_hist,1,'FaceColor',[.4 .4 .4],'LineStyle','none');
view([90 270]);
set(gca,'FontSize',12,'XTickLabel',[],'YTickLabel',[],'XLim',[p.y_min p.y_max],'Layer','top');
grid on;

%% donor and acceptor intensity trace
subplot('Position',[0.11 0.38 0.75 0.27]);
if ~peaks_found
    plot(u.time{mid},u.donor{mid},'Color',[.3 1 .3]);hold on;
    plot(u.time{mid},u.acceptor{mid},'Color',[1 .3 0.3]);
else
    plot(u.time{mid},u.donor{mid},'Color',[.8 1 .8]);hold on;
    plot(u.time{mid},u.acceptor{mid},'Color',[1 .8 0.8]);
    if p.N_bin ==0
        tmp_out=u.donor{mid};
        tmp_out(~u.stateI{mid})=nan;
        plot(u.time{mid},tmp_out,'Color',[0 .7 0]);
        tmp_out=u.acceptor{mid};
        tmp_out(~u.stateI{mid})=nan;
        plot(u.time{mid},tmp_out,'Color',[.8 0 0]);
    else
        plot(u.bintime{mid},u.binD{mid},'Color',[0.5 1 0.5]);
        tmp_out=u.binD{mid};
        tmp_out(~u.stateI{mid})=nan;
        plot(u.bintime{mid},tmp_out,'.-','Color',[0 .5 0],'LineWidth',1);

        plot(u.bintime{mid},u.binA{mid},'Color',[1 .5 .5]);
        tmp_out=u.binA{mid};
        tmp_out(~u.stateI{mid})=nan;
        plot(u.bintime{mid},tmp_out,'.-','Color',[.7 0 0],'LineWidth',1);

    end
end

grid on;
ylabel('Intensity (a.u.)');
set(gca,'FontSize',12,'XTickLabel',[],'XLim',[0 p.xmax],'YLim',[p.y_min p.y_max],'Layer','top');

%% draw Donor and Acceptor hist
subplot('Position',[0.88 0.38 0.1 0.27]);
[tmp_hist1,tmp_hist_bin1]=hist(u.donor{mid},50);
[tmp_hist2,tmp_hist_bin2]=hist(u.acceptor{mid},50);
bar(tmp_hist_bin1,tmp_hist1,1,'FaceColor',[0 1 0],'LineStyle','none');hold on;
bar(tmp_hist_bin2,tmp_hist2,1,'FaceColor',[1 0 0],'LineStyle','none');hold on;

view([90 270]);
set(gca,'FontSize',12,'XTickLabel',[],'YTickLabel',[],'XLim',[p.y_min p.y_max],'Layer','top');
grid on;

%% FRET trace
subplot('Position',[0.11 0.08 0.75 0.27]);
if ~peaks_found
    tmp_out=u.acceptor{mid}./(u.donor{mid}+u.acceptor{mid});
    plot(u.time{mid},tmp_out,'Color',[0.3 0.3 .7]);hold on;
else
	plot(u.time{mid},u.FRETe{mid},'Color',[.95 .95 1]);hold on;
    if p.N_bin == 0
        tmp_out=u.FRETe{mid};
        tmp_out(~u.stateI{mid})=nan;
        plot(u.time{mid},tmp_out,'.-','Color',[0.7 0.7 .7]);
    else
        plot(u.bintime{mid},u.binE{mid},'Color',[0.85 0.85 1],'LineWidth',1);
        tmp_out=u.binE{mid};
        tmp_out(~u.stateI{mid})=nan;
        plot(u.bintime{mid},tmp_out,'.-','Color',[0 0 .5]);
        %                 plot(bintime,stateI,'y.-');
    end
end
if peaks_selected
    tmp_out=[];
    tmp_out2=[];
    for pki=1:length(selected_dw_tr_E)
        plot(selected_dw_tr_t(pki,1):del_t:selected_dw_tr_t(pki,2),selected_dw_tr_E{pki},'Color',[0 0 .5]);
        tmp_out=[tmp_out selected_dw_tr_E{pki}];
        tmp_out2=[tmp_out2 mean(selected_dw_tr_E{pki})];
    end
end
grid on;
ylabel('FRET');
xlabel('Time (s)');
set(gca,'XLim',[0 p.xmax],'YLim',[-0.1 1.1],'FontSize',12,'Layer','top');

%% draw FRET hist
xbin=-0.15:0.01:1.15;
subplot('Position',[0.88 0.08 0.1 0.27]);
[tmp_hist,tmp_hist_bin]=hist(tmp_out,xbin);
bar(tmp_hist_bin,tmp_hist,1,'FaceColor',[0.1 0.1 .5],'LineStyle','none');hold on;
if peaks_selected
    [tmp_hist2,tmp_hist_bin]=hist(tmp_out2,xbin);
    bar(tmp_hist_bin,tmp_hist2,1,'FaceColor',[1 0 0],'LineStyle','none');
end
view([90 270]);
set(gca,'FontSize',12,'XTickLabel',[],'YTickLabel',[],'XLim',[-0.1 1.1],'Layer','top');
grid on;

function show_dwells(~,~)
%% fetch the dwell data
main_fhd=findobj('Tag','Rene_main');
dw=main_fhd.UserData.dw;
u=main_fhd.UserData.u;

p=get_parameter('fetch');

%% plot
plot_FRET_hist(dw,p);
plot_dwelltime_hist(dw,u);
plot_peak_int_dist(dw,p);

function plot_peak_int_dist(dw,p)
%% FRET vs int.
fhd_FRET_hist=findobj('Tag','figure_peak_int_dist');
if ~isempty(fhd_FRET_hist)
    figure(fhd_FRET_hist);clf;
else
    fhd_FRET_hist=figure('Tag','figure_peak_int_dist','Name','int vs. FRET','NumberTitle','off');
end

% plot FRET vs int hist
X=meshgrid(dw.edges{2},dw.edges{1});
Y=meshgrid(dw.edges{1},dw.edges{2});

ax1=axes('Position',[.15 .15 .6 .6]);
pcolor(ax1,X,Y',dw.IE_2Dhist);axis tight;
colormap('jet');
shading flat;
grid on;
box on;
ylim([dw.edges{1}(1) dw.edges{1}(end)])
xlim([dw.edges{2}(1) dw.edges{2}(end)])
xlabel('FRET');
ylabel('Peak intensity (a.u.)');
% FRET hist
ax2=axes('Position',[.15 .77 .6 .2]);
bar(dw.edges{2},dw.IE_2Dhist_Eonly,1,'FaceColor',[.4 .4 .7],'LineStyle','none')
xlim([dw.edges{2}(1) dw.edges{2}(end)])
ylim([0 max(dw.IE_2Dhist_Eonly)*1.1]);
ax2.XTickLabel=[];
grid on;
% Int hist
ax3=axes('Position',[.77 .15 .2 .6]);
bar(dw.edges{1},dw.IE_2Dhist_Ionly,1,'FaceColor',[.7 .7 .4],'LineStyle','none')
view(90,270);
xlim([dw.edges{1}(1) dw.edges{1}(end)])
ylim([0 max(dw.IE_2Dhist_Ionly)*1.1]);
ax3.XTickLabel=[];
grid on;
% set(ax1,'ButtonDownFcn',@set_int_FRET_limit);

%% FRET vs dwell time.
fhd_FRET_hist=findobj('Tag','figure_FRET_vs_dwell');
if ~isempty(fhd_FRET_hist)
    figure(fhd_FRET_hist);clf;
else
    fhd_FRET_hist=figure('Tag','figure_FRET_vs_dwell','Name','Dwell vs. FRET','NumberTitle','off');
end

% plot FRET vs int hist
X=meshgrid(dw.DE_edges{2},dw.DE_edges{1});
Y=meshgrid(dw.DE_edges{1},dw.DE_edges{2});

ax1=axes('Position',[.15 .15 .6 .6]);
pcolor(ax1,X,Y',dw.DE_2Dhist);axis tight;
colormap('jet');
shading flat;
grid on;
box on;
ylim([dw.DE_edges{1}(1) dw.DE_edges{1}(end)])
xlim([dw.DE_edges{2}(1) dw.DE_edges{2}(end)])
xlabel('FRET');
ylabel('Dwell time (s)');
% FRET hist
ax2=axes('Position',[.15 .77 .6 .2]);
bar(dw.DE_edges{2},dw.DE_2Dhist_Eonly,1,'FaceColor',[.4 .4 .7],'LineStyle','none')
xlim([dw.DE_edges{2}(1) dw.DE_edges{2}(end)])
ylim([0 max(dw.DE_2Dhist_Eonly)*1.1]);
ax2.XTickLabel=[];
grid on;
% Int hist
ax3=axes('Position',[.77 .15 .2 .6]);
bar(dw.DE_edges{1},dw.DE_2Dhist_Donly,1,'FaceColor',[.7 .7 .4],'LineStyle','none')
view(90,270);
xlim([dw.DE_edges{1}(1) dw.DE_edges{1}(end)])
ylim([0 max(dw.DE_2Dhist_Donly)*1.1]);
ax3.XTickLabel=[];
grid on;

function plot_FRET_hist(dw,p)
%% get focus on the FRET hist figure
fhd_FRET_hist=findobj('Tag','figure_FRET_hist');
if ~isempty(fhd_FRET_hist)
    figure(fhd_FRET_hist);clf;
else
    fhd_FRET_hist=figure('Tag','figure_FRET_hist','Name','FRET hist','NumberTitle','off');
    fhd_FRET_hist.Position(2)=fhd_FRET_hist.Position(2)-fhd_FRET_hist.Position(4)/2;
    fhd_FRET_hist.Position(4)=fhd_FRET_hist.Position(4)*1.5;
end

subplot(3,1,1);
bar(dw.FRET_hist(:,1),dw.FRET_hist(:,2),1,'FaceColor',[.8 .8 1],'EdgeColor',[0 0 .5]);
xlabel('FRET');
xlim([-0.05 1.1]);
ylabel('# data pts');
set(gca,'FontSize',15);
grid on


subplot(3,1,2);
bar(dw.FRET_hist(:,1),dw.FRET_hist(:,3),1,'FaceColor',[.8 .8 1],'EdgeColor',[0 0 .5]);
xlabel('FRET');
xlim([-0.05 1.1]);
ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
set(gca,'FontSize',15);
grid on

uicontrol('Style','pushbutton',    'String','Fit',...
    'Unit','normalized',...
    'Position',[.92 .35 .06 .08],'Callback',{@fit_fret_hist,dw,fhd_FRET_hist,p});

uicontrol('Style','pushbutton',    'String','Copy',...
    'Unit','normalized',...
    'Position',[.92 .25 .06 .08],'Callback',{@copy_fret_hist_fit,fhd_FRET_hist});

uicontrol('Style','pushbutton',    'String','Apply',...
    'Unit','normalized',...
    'Position',[.92 .15 .06 .08],'Callback',{@Apply_fret_hist_fit,fhd_FRET_hist});

[~,N_hist]=size(dw.FRET_hist);
if N_hist==4
    subplot(3,1,3);
    bar(dw.FRET_hist(:,1),dw.FRET_hist(:,4),1,'FaceColor',[.8 .8 1],'EdgeColor',[0 0 .5]);
    xlabel('FRET');
    xlim([-0.05 1.1]);
    % ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
    set(gca,'FontSize',15);
    grid on
end

function Apply_fret_hist_fit(~,~,~)

main_fhd=findobj('Tag','Rene_main');
dw=main_fhd.UserData.dw;
% u=main_fhd.UserData.u;
% p=get_parameter('fetch');
fit_para=main_fhd.UserData.dw.FRET_hist_fit_res;

N_state=length(fit_para.a);
tmp_Names=cell(N_state,1);
for sti=1:N_state
    out_put(sti,1:4)={true, char(64+sti), fit_para.xc(sti), fit_para.w(sti)};
end
hdl_states_tbl=findobj('Tag','RM_states_tbl');
hdl_states_tbl.Data=out_put;

function fit_fret_hist(~,~,dw,fhd_FRET_hist,p)
x=dw.FRET_hist(:,1)';
[~,flg]=size(dw.FRET_hist);
if flg==3   % no HMM case
    y=dw.FRET_hist(:,3)';
else
    y=dw.FRET_hist(:,4)';
end

%% clean up the histogram
figure(fhd_FRET_hist);
if flg==3
    subplot(3,1,2);
else
    subplot(3,1,3);
end
% bar(dw.FRET_hist(:,1),dw.FRET_hist(:,4),1,'FaceColor',[.8 .8 1]);hold on;
bar(x,y,1,'FaceColor',[.8 .8 1]);hold on;
xlabel('FRET');
xlim([-0.05 1.1]);
ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
set(gca,'FontSize',15);
grid on

%% fit the histograms
[bx,~,~]=ginput();
N_peak=length(bx);

fn_user='';
init_a=[];
init_xc=[];
init_w=[];

for pki=1:N_peak
    appx=num2str(pki);
    fn_user=[fn_user 'a' appx '*exp( -((x-xc' appx ')/w' appx ')^2/2 ) +'];
    
    vid= x > bx(pki)- 0.01 & x < bx(pki)+0.01;
    init_a=[init_a max(y(vid))];
    init_xc=[init_xc bx(pki)];
    init_w=[init_w 0.01];
end

fn_user=[fn_user 'y0'];
init_Lower=[zeros(1,N_peak*3) 0];
init_Upper=[(1:N_peak)*0+Inf (1:N_peak)*0+0.1 (1:N_peak)*0+1.2 Inf];
init_val=[init_a init_w init_xc 0];

fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',init_Lower,'Upper',init_Upper,...
    'StartPoint',init_val);
ft = fittype(fn_user,'options',fo);

[g_fit,~] = fit(x',y',ft);



%% calc idv peaks
x_fc=-0.1:0.0001:1.1;
fit_curv=zeros(length(x_fc),N_peak+1);
fit_summary='';
bin_size=x_fc(2)-x_fc(1);
for pki=1:N_peak
    fit_para.a(pki)=getfield(g_fit, ['a' num2str(pki)]);
    fit_para.xc(pki)=getfield(g_fit, ['xc' num2str(pki)]);
    fit_para.w(pki)= getfield(g_fit, ['w' num2str(pki)]);
    fit_para.popul(pki)=fit_para.a(pki)/( 1/( fit_para.w(pki) *sqrt(2*pi)   )  )/bin_size;
	fit_para.se(pki)=fit_para.w(pki)/sqrt(fit_para.popul(pki));

    fit_curv(:,pki)=fit_para.a(pki)*...
        exp( -(      (x_fc-fit_para.xc(pki)) /fit_para.w(pki)).^2/2);
    
	fit_curv_se(:,pki)=fit_para.a(pki)*...
        exp( -(      (x_fc-fit_para.xc(pki)) /fit_para.se(pki)).^2/2);
    fit_summary=[fit_summary 'a' num2str(pki) ': ' num2str(fit_para.a(pki)) newline];
    fit_summary=[fit_summary 'xc' num2str(pki) ': ' num2str(fit_para.xc(pki)) newline];
    fit_summary=[fit_summary 'w' num2str(pki) ': ' num2str(fit_para.w(pki)) newline];
end
fit_summary=[fit_summary 'y0: ' num2str(g_fit.y0)];
fit_curv(:,pki+1)=sum(fit_curv,2)+g_fit.y0;
clipboard('copy',fit_summary);
fhd_FRET_hist.UserData.fit_summary=fit_summary;

% update guidata

main_fhd=findobj('Tag','Rene_main');
main_fhd.UserData.dw.FRET_hist_fit_res=fit_para;

%% plot
figure(fhd_FRET_hist);

if flg==3
    subplot(3,1,2);
else
    subplot(3,1,3);
end
plot(x_fc,fit_curv(:,1:N_peak),'LineWidth',.5)
plot(x_fc,fit_curv(:,N_peak+1),'LineWidth',3);hold off;

xlabel('FRET');
xlim([-0.05 1.1]);
ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
set(gca,'FontSize',15);
grid on

for pki=1:N_peak
    text(fit_para.xc(pki)-0.1,fit_para.a(pki)*1.3,...
        [num2str(fit_para.xc(pki),'%.2f') '\pm' num2str(fit_para.w(pki),'%.2f' )],'FontSize',9);
end
xlabel('FRET');
ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
legend off;


%% plot extra
%get focus on the FRET hist figure
fhd_FRET_hist=findobj('Tag','figure_FRET_hist_disp');
if ~isempty(fhd_FRET_hist)
    figure(fhd_FRET_hist);clf;
else
    fhd_FRET_hist=figure('Tag','figure_FRET_hist_disp','Name','FRET hist_disp','NumberTitle','off');
	fhd_FRET_hist.Position(3)=400;
    fhd_FRET_hist.Position(4)=300;
    fhd_FRET_hist.Position(2)=fhd_FRET_hist.Position(2)-fhd_FRET_hist.Position(4)/2;
end

axes(fhd_FRET_hist,'Position',[0.25 0.2 0.72 0.75]);
bar(x,y,1,'FaceColor',[0.2 0.3 .4],'LineStyle','none');hold on;
xlabel('FRET');
xlim([-0.05 1.1]);
ylim([0 max(y(:))*1.25]);
ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
set(gca,'FontSize',15);
% grid on

plot(x_fc,fit_curv(:,1:N_peak),'LineWidth',.5)
plot(x_fc,fit_curv(:,N_peak+1),'LineWidth',2,'Color',[1 0 0]);hold off;

xlabel('FRET');
xlim([-0.05 1.05]);
ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
set(gca,'FontSize',15);
set(gca,'Xtick',[0 .2 .4 .6 .8 1]);
% grid on

for pki=1:N_peak
    text(fit_para.xc(pki)-0.1,fit_para.a(pki)*1.2,...
        [num2str(fit_para.xc(pki),'%.3f') '\pm' num2str(fit_para.w(pki),'%.3f' )],'FontSize',9);
end
xlabel('FRET');
ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
legend off;
grid on;


%% plot se
%get focus on the FRET hist figure
fhd_FRET_hist=findobj('Tag','figure_FRET_hist_disp_se');
if ~isempty(fhd_FRET_hist)
    figure(fhd_FRET_hist);clf;
else
    fhd_FRET_hist=figure('Tag','figure_FRET_hist_disp_se','Name','FRET hist_disp_se','NumberTitle','off');
	fhd_FRET_hist.Position(3)=400;
    fhd_FRET_hist.Position(4)=300;
    fhd_FRET_hist.Position(2)=fhd_FRET_hist.Position(2)-fhd_FRET_hist.Position(4)/2;
end

plot(x_fc,fit_curv_se(:,1:N_peak),'LineWidth',2)

xlabel('FRET');
xlim([-0.05 1.05]);
ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
set(gca,'FontSize',15);
set(gca,'Xtick',[0 .2 .4 .6 .8 1]);
% grid on

for pki=1:N_peak
    text(fit_para.xc(pki)-0.1,fit_para.a(pki)*1.2,...
        [num2str(fit_para.xc(pki),'%.3f') '\pm' num2str(fit_para.se(pki),'%.5f' )],'FontSize',9);
end
xlabel('FRET');
ylabel(['# peaks (>' num2str(p.min_dw_len_sel) 's)']);
legend off;
grid on;


%% save histogram
output=[x;y]';
save([p.pth 'E_hist.dat'], 'output','-ascii');
output=[x_fc', fit_curv];
save([p.pth 'E_hist_fit_std.dat'], 'output','-ascii');
output=[x_fc', fit_curv_se];
save([p.pth 'E_hist_fit_sem.dat'], 'output','-ascii');


function copy_fret_hist_fit(~,~,fhd_FRET_hist)
fit_summary=fhd_FRET_hist.UserData.fit_summary;
clipboard('copy',fit_summary);

function plot_dwelltime_hist(dw,u)
%% get focus on the FRET hist figure
fhd_dwelltime_hist=findobj('Tag','figure_dwelltime_hist');
if ~isempty(fhd_dwelltime_hist)
    figure(fhd_dwelltime_hist);clf;
else
    fhd_dwelltime_hist=figure('Tag','figure_dwelltime_hist','Name','Dwell time','NumberTitle','off');
end

% draw
figure(fhd_dwelltime_hist);clf
fhd12.Position(3)=500;
fhd12.Position(4)=200;

subplot('Position',[.15 .3 .34 .65]);
bar(dw.downX,dw.down_hist,1,'FaceColor',[.7 .7 .7]);hold on;
ph1=plot(dw.down_fit,'k');hold off;
set(ph1,'LineWidth',3);
xlabel('Down dwell time (s)');
ylabel('Population');
legend off;
N_dw=sum(dw.N_dw_ind);
text(dw.down_fit.tau*2,dw.down_hist(1)*3/4,...
    [num2str(N_dw) ' events (' num2str(N_dw/u.t_tot_obs,'%.2f') 'Hz)' ],'FontSize',12);
text(dw.down_fit.tau*2,dw.down_hist(1)/2,...
    ['tau = ' num2str(dw.down_fit.tau,'%.2f') 's' ],'FontSize',12);
set(gca,'Fontsize',12);
xlim([0 dw.downX(end)]);
ylim([0 max(dw.down_hist*1.1)]);

subplot('Position',[.65 .3 .34 .65]);
bar(dw.upX,dw.up_hist,1,'FaceColor',[1 .7 0]);hold on;
ph2=plot(dw.up_fit);hold off;
set(ph2,'LineWidth',3);
xlabel('Up dwell time (s)');
ylabel('');
legend off;
N_up=sum(dw.N_up_survived_ind);
text(dw.up_fit.tau*2,dw.up_hist(1)*3/4,...
    [num2str(N_up) ' events (' num2str(N_up/u.t_tot_obs,'%.2f') 'Hz)' ],'FontSize',12);
text(dw.up_fit.tau*2,dw.up_hist(1)/2,...
    ['tau = ' num2str(dw.up_fit.tau,'%.2f') 's' ],'FontSize',12);
set(gca,'Fontsize',12);
xlim([0 dw.upX(end)]);
ylim([0 max(dw.up_hist*1.1)]);
drawnow;

function load_results(~, ~)
add_log('Loading results...');

p=get_parameter('fetch');

%% fetch the data
main_fhd=findobj('Tag','Rene_main');

%% get file info
set_prog(p.hdl_prog,0.05);
[fname,fpath]=uigetfile([p.pth '*.mat']);
clear('p');
%% load results
load([fpath fname]);

%% update data
main_fhd.UserData.u=u;
main_fhd.UserData.p=p;
main_fhd.UserData.dw=dw;
p=get_parameter('set',p);

set_prog(p.hdl_prog,1);
add_log('Results loaded!');

function save_results(~, ~)
%% fetch the data
main_fhd=findobj('Tag','Rene_main');
dw=main_fhd.UserData.dw;
u=main_fhd.UserData.u;
p=get_parameter('fetch');

%% make a folder
add_log('Saving results...');
set_prog(p.hdl_prog,0.05);

timetag=clock;
timetag=[num2str(timetag(1)) num2str(timetag(2),'%02d') num2str(timetag(3),'%02d') ...
    num2str(timetag(4),'%02d') num2str(timetag(5),'%02d') num2str(floor(timetag(6)),'%02d') ];
subfolder_name=['RM_Result.' timetag ];
mkdir(p.pth,subfolder_name);

set_prog(p.hdl_prog,.1);

%% save
% dwell info
save(fullfile(p.pth, subfolder_name, 'All_dwells.mat'),'dw','u','p');
set_prog(p.hdl_prog,.8);

output=dw.FRET_hist;
save(fullfile(p.pth, subfolder_name, 'FRET_hist.dat'),'output','-ascii');
set_prog(p.hdl_prog,.9);

fid=fopen(fullfile(p.pth, subfolder_name, 'rates.txt'),'wt');

N_dw_all=0;
N_up_all=0;
for mid=1:u.N_mol
    N_dw_all=N_dw_all+dw.N_dw_ind(mid);
    N_up_all=N_up_all+dw.N_up_ind(mid);
end
N_down2up=N_dw_all;
N_up2down=N_up_all;
k_down2up=N_dw_all/u.t_tot_obs;
k_up2down=N_up_all/u.t_tot_obs;
fprintf(fid,'N molecules: %g\n',u.N_mol);
fprintf(fid,'down to up events: %g\n',N_down2up);
fprintf(fid,'up to down events: %g\n',N_up2down);
fprintf(fid,'down to up rate: %g Hz\n',k_down2up);
fprintf(fid,'up to down rate: %g Hz\n',k_up2down);
fprintf(fid,'total obs. time: %g s\n',u.t_tot_obs);
fclose(fid);
set_prog(p.hdl_prog,1);
add_log('Saved!');

function show_hetero(~, ~)
add_log('working on heteroneniety...');

%% fetch the FRET data
main_fhd=findobj('Tag','Rene_main');
u=main_fhd.UserData.u;
dw=main_fhd.UserData.dw;
p=get_parameter('fetch');
set_prog(p.hdl_prog,0.05);

%% Classify the molecules
S=classify_idv_mol(dw,p);

%% update data
main_fhd.UserData.S=S;
set_prog(p.hdl_prog,.3);

%% plot heterogeniety
plot_hetero_mol(S,dw);
add_log('Task finished.');
set_prog(p.hdl_prog,1);

function S=classify_idv_mol(dw,p)
% fetch the state info
hdl_state_tbl=findobj('Tag','RM_states_tbl');
S.state_inuse=cell2mat(hdl_state_tbl.Data(:,1));
S.state_name=hdl_state_tbl.Data(:,2);
S.state_name=S.state_name(S.state_inuse);
S.center=cell2mat(hdl_state_tbl.Data(:,3));
S.center=S.center(S.state_inuse);
S.halfwidth=cell2mat(hdl_state_tbl.Data(:,4));
S.halfwidth=S.halfwidth(S.state_inuse)/2;

%% Classify the peaks
% S.halfwidth=0.1;
S.N_state=sum(S.state_inuse);
boundary=[S.center'-S.halfwidth'; S.center'+S.halfwidth']';
N_mol=length(dw.up_ind);

for sid=1:S.N_state
    for mid=1:N_mol
        E_per_peak_ind=[];
        for pid=1:length(dw.up_peak_E_raw_survived_ind{mid})
            E_per_peak_ind(pid)=mean(dw.up_peak_E_raw_survived_ind{mid}{pid});
        end
        tmpvid = boundary(sid,1) < E_per_peak_ind &  E_per_peak_ind < boundary(sid,2);
        S.peak{mid}{sid}=E_per_peak_ind(tmpvid);
        S.N_peak(mid,sid)=length(S.peak{mid}{sid});
        S.existance_bool(mid,sid)=S.N_peak(mid,sid) > p.min_N_peak_in_state;
    end
end
for mid=1:N_mol
    % assining text for the final results
    S.existance_text{mid}=['' S.state_name{S.existance_bool(mid,:)}];
    if isempty(S.existance_text{mid})
        S.existance_text{mid}='--';
    end
end

%% calc N molecules per state
S.N_comp=2^S.N_state;

for cmi=1:S.N_comp
    tmp=dec2bin(cmi-1,S.N_state);
    for tti=1:S.N_state
        S.comp_id(cmi,tti)=str2double(tmp(tti));
    end
    S.comp_name{cmi}=['' S.state_name{boolean(S.comp_id(cmi,:))}];
    if isempty(S.comp_name{cmi})
        S.comp_name{cmi}='--';
    end
    tmpvid=S.existance_bool==S.comp_id(cmi,:);
    vid(:,cmi)=sum(~tmpvid,2)==0;
end

S.N_mol_per_comp=sum(vid);

function plot_hetero_mol(S,dw)
N_mol=length(dw.up_ind);

%%
fhdl=figure(81);
fhdl.Position(2)=0;
fhdl.Position(3)=300;
fhdl.Position(4)=800;
bar(S.N_mol_per_comp,'FaceColor',[.7 0 0]);
set(gca,'XTick',1:S.N_comp,'XTickLabel',S.comp_name,'FontSize',9);
xlim([0.25 S.N_comp+0.75]);
view([90 90])

%% show FRET hist and Int hist of individual molecule
E_bin=-0.05:0.02:1;
I_bin=0:100:8000;
for tmid=1:N_mol
    %% collect all FRET data points for each molecule
%     FRET_raw{tmid}=[];
%     I_raw{tmid}=[];
%     for pki=1:dw.N_up_survived_ind(tmid)
%         FRET_raw{tmid}=[FRET_raw{tmid} dw.up_peak_E_raw_survived_ind{tmid}{pki}];
%         I_raw{tmid}=[I_raw{tmid} dw.up_peak_I_raw_survived_ind{tmid}{pki}];
%     end
%     hist_FRET_raw{tmid}=hist(FRET_raw{tmid},E_bin);
%     hist_I_raw{tmid}=hist(I_raw{tmid},I_bin);
    %% coleect avr FRET for each mol
    E_per_peak_ind=[];
    for pid=1:length(dw.up_peak_E_raw_survived_ind{tmid})
        E_per_peak_ind(pid)=mean(dw.up_peak_E_raw_survived_ind{tmid}{pid});
    end
    hist_FRET_avr{tmid}=hist(E_per_peak_ind,E_bin);
%     hist_I_avr{tmid}=hist(E_per_peak_ind,I_bin);
end

%%
N_mol_per_page=40;
N_page=ceil(N_mol/N_mol_per_page);
N_col_plot=8;
N_raw_plot=ceil(N_mol_per_page/N_col_plot);
% N_raw_plot=10;
x_span=1/N_col_plot;
y_span=1/N_raw_plot;
%%
for pgi=1:N_page
% for pgi=1:1
    %% FRET hist
    figure(170+pgi);clf
    if pgi~=N_page
        mid2disp=(pgi-1)*N_mol_per_page+1:pgi*N_mol_per_page;
    else
        mid2disp=(pgi-1)*N_mol_per_page+1:(pgi-1)*N_mol_per_page + N_mol-(N_page-1)*N_mol_per_page;
    end

    for tmid=1:length(mid2disp)
        c_col_id=floor((tmid-1)/N_raw_plot)+1;
        c_row_id=rem(tmid-1,N_raw_plot)+1;
        subplot('Position',[(x_span*(c_col_id-1)+x_span*0.05)*0.95+0.05,(y_span*(c_row_id-1)+y_span*.1)*0.95+0.05,...
            x_span*0.9*0.95,y_span*0.85*0.95]);
        bar(E_bin,hist_FRET_avr{mid2disp(tmid)},1,'LineStyle','none','FaceColor',[0.7 0 .2]);hold off;
        [X_tick_pos,X_tick_order]=sort(S.center);
        set(gca,'XTick',X_tick_pos,'YTick',[0 10 20],...
            'XLim',[0 1],'YLim',[0 20],...
            'Box','on','FontSize',12);
        if c_col_id~=1
            set(gca,'YTickLabel',[]);
        end
        if c_row_id~=1
            set(gca,'XTickLabel',[]);
        else
            set(gca,'XTickLabel',S.state_name(X_tick_order));
        end
        grid on;

        ahdl=gca;
        text(ahdl.XLim(2)*0.1,ahdl.YLim(2)*0.9,S.existance_text{mid2disp(tmid)},'FontSize',12);
    end

end


if 0
    
    %% Int hist
    figure(72);clf
    for mid=1:N_mol_per_page
        c_col_id=floor((mid-1)/N_raw_plot)+1;
        c_row_id=rem(mid-1,N_raw_plot)+1;
        subplot('Position',[(x_span*(c_col_id-1)+x_span*0.05)*0.95+0.05,(y_span*(c_row_id-1)+y_span*.1)*0.95+0.05,...
            x_span*0.9*0.95,y_span*0.85*0.95]);
        %     bar(E_bin,hist_FRET_raw{mid},1,'FaceColor',[.7 .7 .8],'LineStyle','none'); hold on;
        bar(I_bin,hist_I_avr{mid},1,'FaceColor',[0.1 0.2 0.4]);hold off;
        set(gca,'XTick',0:2000:8000,'XTickLabel',[],'XLim',[I_bin(1) I_bin(end)],'YLim',[0 20]);
        grid on;
    end
    %% FRET vs Int
    figure(73);clf
    for mid=1:N_mol_per_page
        c_col_id=floor((mid-1)/N_raw_plot)+1;
        c_row_id=rem(mid-1,N_raw_plot)+1;
        subplot('Position',[(x_span*(c_col_id-1)+x_span*0.05)*0.95+0.05,(y_span*(c_row_id-1)+y_span*.1)*0.95+0.05,...
            x_span*0.9*0.95,y_span*0.85*0.95]);
        hold on;
        for tpi=1:length(dw.up_peak_E_raw_survived_ind{mid})
            plot(dw.up_peak_E_raw_survived_ind{mid}{tpi},dw.up_peak_I_raw_survived_ind{mid}{tpi},'.','Color',[.7 .7 .7]);
        end
        plot(dw.up_peak_E_avr_survived_ind{mid},dw.up_peak_I_avr_survived_ind{mid},...
            '.','MarkerSize',9,'Color',[.8 0 0]);hold off;
        set(gca,'XTick',S.center,'YTick',[],...
            'Box','on',...
            'XLim',[0.2 E_bin(end)],'YLim',[1000 I_bin(end)]);
        if c_col_id~=1
            set(gca,'YTickLabel',[]);
        end
        if c_row_id~=1
            set(gca,'XTickLabel',[]);
        else
            set(gca,'XTickLabel',S.state_name);
        end
        grid on;
        ahdl=gca;
        text(ahdl.XLim(2)*0.7,ahdl.YLim(2)*0.9,S.existance_text{mid});
        
        %     for sdi=1:S.N_state
        % %         plot(boundary(sdi,1), boundary(sdi,2),);
        %     end
    end
    
end

function add_log(str,hdl_log_text)
if nargin<2
    hdl_log_text=findobj('Tag','RM_log_text');
end

hdl_log_text.String=str;
drawnow

function set_prog(hdl_prog,p_rate)
% p_rate=0.1;
cla(hdl_prog);
rectangle(hdl_prog,'Position',[0 0 1 1],'EdgeColor',[.7 .7 .7],'FaceColor',[0.9 0.9 .9],'Curvature',[.03 1]);
rectangle(hdl_prog,'Position',[0 0 p_rate 1],'EdgeColor',[.7 .7 .7],'FaceColor',[0.6 0.6 .9],'Curvature',[.03 1]);
if p_rate~=0
    a=text(hdl_prog,p_rate/2,0.5,[num2str(p_rate*100,'%.1f') '%']);
    a.FontSize=7;
end
drawnow

function save_traces(~,~)

%% fetch the data
main_fhd=findobj('Tag','Rene_main');
% dw=main_fhd.UserData.dw;
u=main_fhd.UserData.u;
p=get_parameter('fetch');

%% make a folder
add_log('Saving traces...');
set_prog(p.hdl_prog,0.05);

timetag=clock;
timetag=[num2str(timetag(1)) num2str(timetag(2),'%02d') num2str(timetag(3),'%02d') ...
    num2str(timetag(4),'%02d') num2str(timetag(5),'%02d') num2str(floor(timetag(6)),'%02d') ];
subfolder_name=['RM_selected_tr.' timetag];
mkdir(p.pth,subfolder_name);



%% save traces
for mid=1:u.N_mol
    output=[u.time{mid} u.donor{mid} u.acceptor{mid}];
    tmp_fname=fullfile(p.pth, subfolder_name, [p.user_keyword '_tr' num2str(u.mol_id_of_org(mid)) '.dat']);
    save(tmp_fname,'output','-ascii');
%     save([p.pth '\' subfolder_name '\' p.user_keyword '_tr' num2str(u.mol_id_of_org(mid)) '.dat'],'output','-ascii');
end

add_log('Traces saved!');
set_prog(p.hdl_prog,1);


% tmp_hdl=findobj('Tag','RM_flg2go_show_hetero');
% if tmp_hdl.Value
%     show_hetero(1,1);
% end
% tmp_hdl=findobj('Tag','RM_flg2go_save_results');
% if tmp_hdl.Value
%     save_results(1,1);
% end


