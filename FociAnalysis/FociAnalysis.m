function FociAnalysis(dirname,paramFit)



%param_fit=50; %% Minimo de puntos para considerar el fit exponencial valido. 
%list_1 = dir('*xy*');% todas las xy carpetas

clearvars 
contents = dir([dirname,'xy*']); %List all xy folders
num_dir_tmp = numel(contents);
nxy = [];
num_xy = 0;

%Creates a struct for all fluor1 directories 

dirnamelist=cell(1,num_dir_tmp);
for i = 1:num_dir_tmp
    if (contents(i).isdir) && (numel(contents(i).name) > 2)
    num_xy = num_xy+1;
    nxy = [nxy, str2double(contents(i).name(3:end))];
    dirnamelist{i} = [dirname,contents(i).name,filesep,'fluor1'];%Set fluor1 channel as the green channel
    end
end

 

for p= 1:num_xy
    
    name_root=list_1(p).name; 
    name_folder=strcat(name_root,'/morphometrics');
    cd(name_folder)
    list_2 = dir('*pill_MESH.mat'); %pill mesh file
    name_morph=list_2(1).name;
    load(name_morph);

    %%%%% Encuentra el frame inicial y final por celula con la
    %%%%% mayor longitud in-interrumpida de frames. Genera la matrix
    %%%%% limites. Cada fila es una celula. Esto podria ser una
    %%%%% funcion?

    %%%%% Inicia aca

    ind=zeros(Nframe,Ncell);
    cont=zeros(Ncell, Nframe);
    limites=zeros(Ncell,2);

    for N=1:1:Ncell
        for fr=1:1:Nframe
            allCN = [frame(fr).object.cellID]; % comma separated list expansion 
            ind = find(allCN == N);
            aux = isempty(ind);
            if aux==1
                cont(N,fr)=-10;    
            end
        end
    end  



    for j=1:1:Ncell
        
         maximo=Nframe;
         min=1;
         temp=find(cont(j,:)==0);
         if isempty(temp)==1
             continue
         else
         loc=find(diff(temp)~=1);
         
         if isempty(loc)==1
            limites(j,1)=temp(1);
            limites(j,2)=temp(end);
            continue
         else
            delta=[temp(1),temp(loc),temp(loc(end)+1),temp(end)];
            dif_d=diff(delta);
            pos=find(dif_d==max(dif_d));
            limites(j,1)=delta(pos(1));
            limites(j,2)=delta(pos(1)+1);
          end
         end
    end

            %%%% Termina aca. 
            
            
    salida=fun_anal4N(Ncell,frame,name_morph,limites,paramFit);

            %fun_ph_growth(Ncell,name,size(tsStack,2),limites,param_fit);
            
    cd ../..
            
            %cd('cell_cycle') 
    name=name_morph(1:end-4); 
    arch_nombre=strcat(name,'_OUTPUT_FOCI.mat');
    arch_nombre=strcat('cell_cycle/',arch_nombre);
    save(arch_nombre,'salida');

    %cd ..
%            
    clear frame
     
 end       
               

end

