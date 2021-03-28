%% INDIVIDUAL SUBJECTS 
% by Michelle Tran

%% set up directory
myDir = uigetdir;
color_folder = dir(myDir);
color_folder([1,2],:) = [];
size(color_folder);

%% define variables
exp = ["Covert 1","Covert 2","Covert 3", "Overt 1","Overt 2","Overt 3"];
fileNames = ["covert1","covert2","covert3", "overt1","overt2","overt3"];
plotLines = ["-", "-", "-", "-.", "-.", "-."];
Color = { (1/255)*[204 121 167], (1/255)*[0 158 115], (1/255)*[86 180 233],(1/255)*[213 94 0], (1/255)*[230 159 0] ,(1/255)*[0 114 178]};
Markers = ['.', '.', '.', '*', '*', '*'];
Subjects = ['A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'K'];
unique_eccen= [-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40];

%% -----------sorting data-------------------------------------------------------------
 % define whether blue or red experiment
if strcmp(color_folder(1).name,'Covert blue') || strcmp(color_folder(1).name,'Overt blue') 
    COLOR='Blue';
    colNum={'[0, 0, 255]'};
    folderName='blue_subjects';  
end

if strcmp(color_folder(1).name,'Covert red') || strcmp(color_folder(1).name,'Overt red') 
    COLOR='Red';
    colNum={'[255, 0, 0]'};
    folderName='red_subjects';  
end

% make new folder to store results/ figures
new_folder=folderName;
mkdir(new_folder)
listing = dir(new_folder);
fpath = listing(1).folder;

% store file directory for covert and overt protocols
covert_folder=dir(fullfile(color_folder(1).folder,color_folder(1).name));
for i = 3: size(covert_folder)
   covert(i-2).name= fullfile(covert_folder(i).folder, covert_folder(i).name);
end

overt_folder=dir(fullfile(color_folder(2).folder,color_folder(2).name));
for i = 3: size(overt_folder)
   overt(i-2).name= fullfile(overt_folder(i).folder, overt_folder(i).name);
end

% sort out different conditions (1,2,3) for overt and covert
for c=1:length(covert)
    filename= char(covert(c).name);
    raw_data = readtable(filename);
    
    trial = raw_data(:,1);
    eccen = raw_data(:,2);%
    RT = raw_data(:,3); %
    dir= raw_data(:,4); %left or right
    color= raw_data(:,5);
    correctness= raw_data(:,6);
    color123= raw_data(:,7);
    name=raw_data(:,8);
    attention=raw_data(:,9);
    
    
    for i=1:size(trial) %number of corrected trials
        if color123{i,1}== 1 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==0
            c1(c).eccentricities{i,1}= (eccen{i,1})*-1;
            c1(c).RT{i,1}= (RT{i,1})*1000; 
        elseif color123{i,1}== 1 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==2
            c1(c).eccentricities{i,1}= eccen{i,1};
            c1(c).RT{i,1}= (RT{i,1})*1000;  
            
        elseif color123{i,1}== 2 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==0
            c2(c).eccentricities{i,1}= (eccen{i,1})*-1;
            c2(c).RT{i,1}= (RT{i,1})*1000; 
        elseif color123{i,1}== 2 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==2
            c2(c).eccentricities{i,1}= eccen{i,1};
            c2(c).RT{i,1}= (RT{i,1})*1000;  
            
        elseif color123{i,1}== 3 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==0
            c3(c).eccentricities{i,1}= (eccen{i,1})*-1;
            c3(c).RT{i,1}= (RT{i,1})*1000; 
        elseif color123{i,1}== 3 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==2
            c3(c).eccentricities{i,1}= eccen{i,1};
            c3(c).RT{i,1}= (RT{i,1})*1000; 
 
        end
    end   
end

for o=1:length(overt)
    filename= char(overt(o).name);
    raw_data = readtable(filename);
    
    trial = raw_data(:,1);
    eccen = raw_data(:,2);%
    RT = raw_data(:,3); %
    dir= raw_data(:,4); %left or right
    color= raw_data(:,5);
    correctness= raw_data(:,6);
    color123= raw_data(:,7);
    name=raw_data(:,8);
    attention=raw_data(:,9);
    
    for i=1:size(trial) %number of corrected trials
        if color123{i,1}== 1 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==0
            o1(o).eccentricities{i,1}= (eccen{i,1})*-1;
            o1(o).RT{i,1}= (RT{i,1})*1000; 
        elseif color123{i,1}== 1 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==2
            o1(o).eccentricities{i,1}= eccen{i,1};
            o1(o).RT{i,1}= (RT{i,1})*1000;  
            
        elseif color123{i,1}== 2 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==0
            o2(o).eccentricities{i,1}= (eccen{i,1})*-1;
            o2(o).RT{i,1}= (RT{i,1})*1000; 
        elseif color123{i,1}== 2 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==2
            o2(o).eccentricities{i,1}= eccen{i,1};
            o2(o).RT{i,1}= (RT{i,1})*1000;  
            
        elseif color123{i,1}== 3 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==0
            o3(o).eccentricities{i,1}= (eccen{i,1})*-1;
            o3(o).RT{i,1}= (RT{i,1})*1000; 
        elseif color123{i,1}== 3 && correctness{i,1}==1 && strcmp(color{i,1},colNum)==1 && dir{i,1}==2
            o3(o).eccentricities{i,1}= eccen{i,1};
            o3(o).RT{i,1}= (RT{i,1})*1000; 
        end
    end  
end

% convert table to matrix
for i=1:length(Subjects)
    cov1(i).eccentricities=cell2mat(c1(i).eccentricities);
    cov2(i).eccentricities=cell2mat(c2(i).eccentricities);
    cov3(i).eccentricities=cell2mat(c3(i).eccentricities);
    ov1(i).eccentricities=cell2mat(o1(i).eccentricities);
    ov2(i).eccentricities=cell2mat(o2(i).eccentricities);
    ov3(i).eccentricities=cell2mat(o3(i).eccentricities);
    cov1(i).RT=cell2mat(c1(i).RT);
    cov2(i).RT=cell2mat(c2(i).RT);
    cov3(i).RT=cell2mat(c3(i).RT);
    ov1(i).RT=cell2mat(o1(i).RT);
    ov2(i).RT=cell2mat(o2(i).RT);
    ov3(i).RT=cell2mat(o3(i).RT);
end
% add all experiments to list  
folders={cov1,cov2,cov3,ov1,ov2,ov3};
%----------------------------- end of sorting data----------------------------------------------------

%% statistical analysis
% loop through each protocol
 for j = 1:length(folders) % cov1,cov2,cov3,ov1,ov2, or ov3
    folder = folders{j}; % defien single protocol
    % Make empty arrays
    Right_slope = [];
    Right_slope_sd = [];
    Right_intercept = [];
    Right_intercept_sd = [];
    Chi2_right = [];
    Red_Chi2_right = [];
    R_Right = [];
    Left_slope = [];
    Left_slope_sd = [];
    Left_intercept = [];
    Left_intercept_sd = [];
    Chi2_left = [];
    Red_Chi2_left = [];
    R_Left = [];
    subject=[];
    protocol=[];
    
    % define aesthetics 
    C = Color{j};
    LS= plotLines(j);
    MS= Markers(j);
    Label = exp(j);

    % loop through each subject 
    for k = 1: length(folder) 
        figure(k)
        % define subject data
        RT=folder(k).RT;
        eccen=folder(k).eccentricities;
        
        % apply outlier function
        r = find((prctile(RT, 25) - 1.5*iqr(RT)) < RT & RT < (prctile(RT, 75) + 1.5*iqr(RT)));
        RT_correct= RT(r, :);
        eccen_correct= eccen(r, :);
        
        % make empty arrays 
        eccen_RT = []; 
        meanRT = [];
        SD= [];
        SE = [];
        
        % group all the same eccentrcities and find parameters
        for g = 1:length(unique_eccen) % 17 degrees
            for v = 1:length(eccen_correct) % total eccentricities for subject
                if eccen_correct(v) == unique_eccen(g) 
                    eccen_RT = [eccen_RT, RT_correct(v)];
                end
            end
            
            % calculate mean, SD, SE
            meanRT = [meanRT, mean(eccen_RT)]; 
            SD = [SD, std(eccen_RT)]; 
            SE = [SE, std(eccen_RT)/sqrt(length(eccen_RT))];
            eccen_RT = []; 
        end 
        % separate negative an positive ecentrcities 
           %negative ecentricities 
        info(1).SE = SE(unique_eccen<=0);
        info(1).meanRT = meanRT(unique_eccen<=0);
        info(1).eccen = unique_eccen(unique_eccen<=0);
        
           %positive ecentiricities
        info(2).SE = SE(unique_eccen>=0);
        info(2).meanRT = meanRT(unique_eccen>=0);
        info(2).eccen = unique_eccen(unique_eccen>=0);
       
 %% apply chi square function and Plot
        grid on
        chiSquareFit = ChiSquarePlot(info, LS, MS, C);
        
        h = plot(0, 3000, 'Color', C, 'MarkerFaceColor', C, 'MarkerEdgeColor', C, ...
            'LineStyle', LS, 'Marker', MS, 'DisplayName', Label);
        
        % legend/ tile
       
        pbaspect([2 3 1])
        xlabel('Retinal Eccentricity (°)');
        ylabel('Reaction Time (ms)');
        title(['Subject ', Subjects(k), ': Reaction Time vs. Eccentricity (',COLOR,')']);
        xlim([min(unique_eccen)-5 max(unique_eccen)+5]);
        xticks(min(unique_eccen):5:max(unique_eccen));
        xtickangle(0);
        xline(0, 'Color', '#C9C9C9', 'HandleVisibility', 'off');
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',8)
        ylim([50 600]);
        legend('show', 'Location', 'northeastoutside');
        
        % save figures 
        saveas(gcf, fullfile(fpath, [Subjects(k), '.png']))
        hold on
       
        % make parameter table and save
        params = struct2table(chiSquareFit);
        param_info = table2array(params);
        
        Right_slope = [Right_slope; param_info(2,1)];
        Right_slope_sd = [Right_slope_sd; param_info(2,4)];
        Right_intercept = [Right_intercept; param_info(2,2)];
        Right_intercept_sd = [Right_intercept_sd; param_info(2,5)];
        Chi2_right = [Chi2_right; param_info(2,3)];
        Red_Chi2_right = [Red_Chi2_right; param_info(2,6)];
        R_Right = [R_Right; param_info(2,7)];
        Left_slope = [Left_slope; param_info(1,1)];
        Left_slope_sd = [Left_slope_sd; param_info(1,4)];
        Left_intercept = [Left_intercept; param_info(1,2)];
        Left_intercept_sd = [Left_intercept_sd; param_info(1,5)];
        Chi2_left = [Chi2_left; param_info(1,3)];
        Red_Chi2_left = [Red_Chi2_left; param_info(1,6)];
        R_Left = [R_Left; param_info(1,7)];
        subject=[subject; Subjects(k)];
        protocol=[protocol; exp(j)];
        
    end
    
    param_table = table(subject, protocol, Right_slope, Right_slope_sd, Left_slope, Left_slope_sd, ...
    Right_intercept, Right_intercept_sd, Left_intercept, Left_intercept_sd,...
    Chi2_right, Red_Chi2_right, Chi2_left, Red_Chi2_left, R_Right, R_Left);
    
    writetable(param_table, fullfile(fpath, fileNames(j) + '_parameters.xlsx'));
        
 end
 
 % define chi square function
 function [chiSquareFit] = ChiSquarePlot(info, LS, MS, Color)
        for vi = 1:2 
            weights = (1./info(vi).SE).^2;
            f = @(x, xPoints, yPoints, w)sum(w.*((yPoints- ((xPoints.*x(1))+x(2))).^2));
            optFun = @(x)f(x, info(vi).eccen, info(vi).meanRT, weights);
            ms = MultiStart;
            OLSFit = polyfit(info(vi).eccen, info(vi).meanRT, 1);
            guessParams = [OLSFit(1), OLSFit(2)];
            problem = createOptimProblem('fmincon', 'x0', guessParams, ...
            'objective', optFun, 'lb', [OLSFit(1)-50, OLSFit(2)-500], 'ub', [OLSFit(1)+50, OLSFit(2)+500]);
            params = run(ms, problem, 25);
            slope = params(1);
            intercept = params(2);
            chiSquareFit(vi).slope = slope;
            chiSquareFit(vi).intercept = intercept;
            chi2Val = optFun(params);
            chiSquareFit(vi).chi2Val = chi2Val;
            syms sErr;
            slopeErr = solve(f([sErr, intercept], info(vi).eccen, info(vi).meanRT, weights)==...
                chi2Val + 1, sErr);
            chiSquareFit(vi).slopeErr = double(slopeErr(2) - slope);
            syms iErr;
            intErr = solve(f([slope, iErr], info(vi).eccen, info(vi).meanRT, weights) == chi2Val+1, iErr);
            chiSquareFit(vi).interceptErr = double(intErr(2)-intercept);
            chiSquareFit(vi).redChiSquare = chi2Val/ (length(info(vi).eccen) - 2);
            R = corrcoef(info(vi).eccen, info(vi).meanRT);
            chiSquareFit(vi).R = R(1,2);

            scatter(info(vi).eccen, info(vi).meanRT, 10, MS, 'MarkerEdgeColor', Color, 'MarkerFaceColor', Color, 'HandleVisibility', 'off');
            hold on;
            plot(info(vi).eccen, polyval(params,info(vi).eccen), 'color', Color, 'linestyle', LS, 'HandleVisibility', 'off');
            hold on;
            errorbar(info(vi).eccen, info(vi).meanRT, info(vi).SE, '.' , 'color' , Color, 'CapSize', 0, 'HandleVisibility',...
                'off');
            hold on;
        end

end
