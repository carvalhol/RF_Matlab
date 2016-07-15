%BY LUCIO DE ABREU CORREA (MSSMat - CentraleSupelec)

clear all
close all
clc

PLOT_Temp = 0;
    
for pp =1
    
    if pp==1
        dire = './MOMENT_TEST_CTE';
        info = read_SEM_files(dire);
        f=info.freq;
        vs = info.Vs;
        vp = info.Vp;
    end
    
    cd(dire)
    fid = fopen(info.capteurs_file, 'r' );
    src_pos = fscanf(fid,'%g',[3 Inf]);
    line_sz = (size(src_pos,2)-1)/4;
    fclose(fid);
                
    font_size = 22;
    root = pwd;
    cd('traces')
    
    TXT = strcmp(info.traces_format, 'text');
    N_COLUMNS = sum(info.outVar_Int) + int32(1);
    DISP_START = sum(info.outVar_Int(1:4)) + int32(2);
    
    if TXT
        nome = dir('CAP*.txt');
        for i=1:size(nome)
            %fid = fopen(['CAP_0',num2str(i-1,'%3.3i'),'.txt']);
            fid = fopen(nome(i).name);
            Out_aux = fscanf(fid,'%g',[N_COLUMNS Inf])';
            fclose(fid);
            u(i,:,:) = Out_aux(1:1:end,DISP_START:DISP_START+2)'; %Why hte :2: ?
            t =  Out_aux(1:1:end,1);
            clear Out_aux
        end
    else    
            nome = dir('capteurs*');
            count = 1;
            for i = 1 : size(nome,1)
                info = hdf5info(nome(i).name);
                for j = 1 : size(info.GroupHierarchy.Datasets,2)
                    aux = h5read(nome(i).name, info.GroupHierarchy.Datasets(j).Name);
                    D(count,:,:)=aux(DISP_START:DISP_START+2,:);
                    id(count) = str2double(info.GroupHierarchy.Datasets(j).Name(6:end))+1;
                    count = count+1;
                    t = aux(1,:);
                    clear aux
                end
            end
            [id, SortIndex] = sort(id);
            u = D(SortIndex,:,:);
            mean(diff(id));
            
            clear D
            save([dire,num2str(f),'bruto.mat'],'u','id','src_pos','t','-v7.3')
            % %     return
        
        
    end
    
    cd(root)
    %%
    
    if PLOT_Temp
        plot(t,squeeze(u(1,:,:)))
    end
    
    %%
    %Good resolution
    %pad_x = 512*4;
    %pad_t = 2048*4;
    
    %Normal resolution
    pad_x = 512;
    pad_t = 2048;
    
    % Sources
    sources = 1;
    src = sum(squeeze(u(1:sources,:,:)),1);
    for i = 1 : sources
        ene(i) = rms(src(i,:));
    end
    energy = sum(ene);
    
    %%    
    subs_t = 1;
    subs_x = 1;
    
    n_lines = (size(src_pos,2)-1)/line_sz;
    
    for i = 1 : n_lines
        %d{i} = flipdim(squeeze((u((i-1)*line_sz+sources+1:subs_x:(i)*line_sz+sources,3,1:subs_t:end)))/energy,1);
        d{i} = flip(squeeze(sum(u((i-1)*line_sz+sources+1:subs_x:(i)*line_sz+sources,:,1:subs_t:end),2))/energy,1);
        %d{i} = d{i}.*(ones(50,1)*linspace(0,50^2,4343))  ;
    end
    
    t = t(1:subs_t:end);
    
    t_fa = 1/mean(diff(t));
    x_fa = 1/mean(diff(src_pos(3,end-line_sz+1:subs_x:end)));
    
    w = 0:t_fa/pad_t:t_fa-t_fa/pad_t;
    k = 0:x_fa/pad_x:x_fa-x_fa/pad_x;



    %% plot temporal antes
%     for i =1:size(d,2)
%         a=figure;
%         set(gcf,'renderer','zbuffer');
%         surf(t,1:line_sz,d{i}(:,:))
%         shading flat
%         view(2)
%         axis square
%         colorbar
%         %saveas(a,['/mssmat2/home/abreuc/CILAMCE/figu/',num2str(i),'.jpg'])
%     end
   %return
    %% Filter
    
     band = [2/7 2]*f; %For rickers source
%     
%     for i = 1:size(d,2)    
%         d{i} = OG_filter(d{i}',t_fa,band   ,2,2,'bandpass',0)';
%         d{i} = OG_filter(d{i} ,x_fa,band/vs,2,1,'bandpass',0);
%     end
%     
     %% janelamento
%     
%     for i = 1:size(d,2)
%         d{i} = OG_Janelamento(d{i}(:,:),1,'tukey',0,.2);
%         d{i} = OG_Janelamento(d{i}(:,:),2,'tukey',0,.2);
%     end
    
    %% plot tempora l
    
    %for i =1:size(d,2)
    %    figure;
    %    set(gcf,'renderer','zbuffer');
    %    surf(d{i})
    %    shading flat
    %    view(2)
    %    axis square
    %end
    
    % return
    %% fft2d
    cd(root)
    cd ..
    tau=1/f; %TODO Read tau no input.spec
    for i = 1:size(d,2)
        corte = 1000;
        RIC = fft(src,pad_t)/numel(src);
        %RIC = ricker(t,f,tau);
        %RIC = fft(RIC,pad_t)/numel(RIC);
        Ri = RIC'*ones(1,pad_x);
        au{i} = fft(d{i}(:,:),pad_x,1);
        au{i} = fft(au{i}(:,1:end),pad_t,2)./Ri';
        %save([dire,'fft2d.mat'],'d','w','k','f','t','posicao','Ri','-v7.3')
    end
    
    %% Dispersion
    
    for i = 1:size(d,2)

        lim = logical((w<=band(2)).*(w>=band(1)));
        
        a=figure(i);hold on
        set(gca,'FontSize',font_size)
        set(gcf,'renderer','zbuffer');

        %surf(k,w,abs( au{i}(:,:)'/max(max(au{i}(:,:))) ))
                
        %surf(k(:),w(lim),abs( au{i}(:,lim)'/max(max(au{i}(:,lim))) ))
        
        emdb = 20*log(abs( au{i}(1:end/2,lim)'/max(max(au{i}(1:end/2,lim)))));
        ind = emdb<=-60;
        emdb(ind) = NaN;
        
        surf(k(1:end/2),w(lim), emdb )
                
        shading flat
        xlabel(['Real Wave-number [1/m]'])
        ylabel(['Frequency [Hz]'])
        %ylim([0 400])
        %xlim([0 3])
        view(2)
        axis square
        colorbar

        x = k;
        y = vp*x;
        z = 1*ones(numel(x));
        plot3(x,y,z,'k--','LineWidth',1)
        
        x= k;
        y = vs*x;
        z = 1*ones(numel(x));
        plot3(x,y,z,'k--','LineWidth',1)
        
       
%         saveas(a,[dire,num2str(i),'teste.png']);
%         saveas(a,[dire,num2str(i),'teste.jpg']);

    end
    
end