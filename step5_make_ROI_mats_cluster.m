%%%% make brain RDMs
%% this script extracts values from a beta and creates an ROI.mat
%% take the info from that trial for each voxel in an ROI, re-structure it, and export

%%%% HEY! Check the for loop and make sure you're doing all 246 brainnetome ROIs

clear all
addpath('/mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/spm12')
addpath /mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/
% betas needs to be a 4x120 cell. 4 days and 120 individual trials/objects per day
% these betas (and stimIDs) need to be in stimID order.

% so, I need to load all the renamed betas and sort by day and then by
% stimID. After that I have to put them into a cell.
subjects = {'5007'};
%subjects = {'5002','5017','5019','5020'};%,'5021','5022'};
% next do 1:246 for 5021 and 5022
% 5007, 5011, 5012, 5014 already run and are on desktop
% '5001', '5002'

for subj = 1:length(subjects)

    tic
    disp('subject loop is running')
    subject = subjects{subj};

    %all_betas = dir(sprintf('/Users/matthewslayton/Desktop/ST_localOnly/renamed_betas/%s/*.nii',subject));
    all_betas = dir(sprintf('/mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/renamedBetas/%s/*.nii',subject));
    
    beta_names = cell(length(all_betas),2); %full name, stimID
    for row = 1:length(all_betas)
    
        currBeta = all_betas(row).name;
        currID = str2double(cell2mat(extractBetween(all_betas(row).name,'stimID','.')));
        beta_names{row,1} = currBeta;
        beta_names{row,2} = currID;
    end
    
    % there are 480 per subject, 120 per day (remember 20 trials per 3 runs,
    % two objects per trial is 40 per run for 3 runs = 120)

    % something for if we're missing betas
%     if length(beta_names) ~= 480
%         disp('There are not 480 betas like there should be')
%         continue
%     end
    
    % make four cells, one for each day
    beta_names_day1 = beta_names(1:120,:);
    beta_names_day2 = beta_names(121:240,:);
    beta_names_day3 = beta_names(241:360,:);
    beta_names_day4 = beta_names(361:480,:);
    
    % sort by col 2 which has the stimIDs
    beta_names_day1_sorted = sortrows(beta_names_day1,2);
    beta_names_day2_sorted = sortrows(beta_names_day2,2);
    beta_names_day3_sorted = sortrows(beta_names_day3,2);
    beta_names_day4_sorted = sortrows(beta_names_day4,2);
    
    % the final form needs to be a 4x120 cell of the names only
    betas = vertcat(beta_names_day1_sorted(:,1)',beta_names_day2_sorted(:,1)',...
        beta_names_day3_sorted(:,1)',beta_names_day4_sorted(:,1)');
    
    % final form of stimIDs is a 4x120 double of the stimIDs only
    stimIDs_cell = vertcat(beta_names_day1_sorted(:,2)',beta_names_day2_sorted(:,2)',...
        beta_names_day3_sorted(:,2)',beta_names_day4_sorted(:,2)');
    
    stimIDs = cell2mat(stimIDs_cell);
    
    for X = 1:246 %1:246 %246 ROIs from brainnetome
        %tic %i feel like this is in the wrong spot. who cares if it takes 30ish seconds per ROI? 
        %% i want to know how many subjects have run
        %ROI_file = sprintf('/Users/simonwdavis/Library/CloudStorage/Dropbox/ppp/ROIs/BN_ROI_%03d.nii', X);
        %ROI_file = sprintf('/Users/matthewslayton/Documents/Duke/Simon_Lab/RSA_practice/ModelX_ROIs/ROIs/BN_ROI_%03d.nii.gz', X);
        ROI_file = sprintf('/mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/ROIs/BN_ROI_%03d.nii.gz', X);
        currFileUnzip = gunzip(ROI_file); %loads as a cell
        
        %disp(sprintf('Working on ROI %d'),X) %would this work?
    
        ROI = spm_read_vols(spm_vol(currFileUnzip{:})); %used to load ROI_file directly, but it can't read compressed
        ROI(ROI==0) = NaN;
        
        for j = 1:4  % day loop
            try
                parfor i = 1:size(betas, 2)  % stim loop
                    %BETA = spm_read_vols(spm_vol(sprintf('/Users/simonwdavis/Library/CloudStorage/Dropbox/ppp/betas/%s', betas{j,i}) ));
                    %BETA = spm_read_vols(spm_vol(sprintf('/Users/matthewslayton/Documents/Duke/Simon_Lab/RSA_practice/ModelX_ROIs/betas/%s', betas{j,i}) ));
                    %BETA = spm_read_vols(spm_vol(sprintf('/Users/matthewslayton/Desktop/ST_localOnly/renamed_betas/%s/%s',subject,betas{j,i})));
                    BETA = spm_read_vols(spm_vol(sprintf('/mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/renamedBetas/%s/%s',subject,betas{j,i})));
                    Xvals = ROI .* BETA;
                    vals=(Xvals(~isnan(Xvals)));
                    ROIvals(:,i) = reshape(vals, size(vals,1), 1);
                end
                
                R = corrcoef(ROIvals, 'rows', 'pairwise');
                stimID = stimIDs(j,:);

                % make output folder if it's not already there
                %outputFolder = sprintf('/Users/matthewslayton/Desktop/ST_localOnly/RSA_mats/%s/',subject);
                outputFolder = '/mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/RSA_mats/';
                if ~exist(outputFolder,'dir'); mkdir(outputFolder); end
                subjOutputFolder = sprintf('/mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/RSA_mats/%s/',subject);
                if ~exist(subjOutputFolder,'dir'); mkdir(subjOutputFolder); end
                %filename = sprintf('/Users/matthewslayton/Desktop/ST_localOnly/RSA_mats/%s/%s_ROI%03d_Day%d.mat', subject, subject, X, j);
                filename = sprintf('/mnt/munin2/Simon/NetTMS.01/Analysis/SingleTrialModels/June_2023_LSS/RSA_mats/%s/%s_ROI%03d_Day%d.mat',subject,subject,X,j);
                save(filename, "ROIvals", "R", "stimID"); % save 1 file for each ROI / day
                clear ROIvals R stimID;
            catch
                fprintf('error on day %d, ROI %d\n', j, X);
            end
        end
        %toc
        fprintf('\n Finished ROI %d \n', X);
    end
    toc
end %subj loop
%%% stimIDs = [13	21	29	33	36	42	46	61	63	69	93	107	109	116	118	129	135	140	141	145	165	167	169	170	171	189	193	196	206	207	215	222	230	234	246	250	253	256	262	283	291	303	305	323	366	375	384	399	406	408	418	428	432	433	442	459	485	490	496	497	511	518	530	540	542	560	571	573	586	606	608	609	618	629	634	649	659	668	677	680	697	719	740	756	780	784	785	794	810	816	817	827	829	832	834	851	853	856	870	879	886	888	891	900	933	934	936	939	945	946	958	961	968	969	978	979	982	992	995	1015; 2	3	15	22	32	34	48	77	82	85	86	87	90	92	95	98	119	124	159	164	173	179	181	191	212	218	235	239	245	259	298	307	312	327	328	330	355	378	379	388	391	409	421	423	427	429	434	435	448	450	454	456	458	465	473	476	488	508	509	514	525	534	559	574	582	589	592	611	614	621	623	626	636	637	646	656	662	683	687	693	699	715	727	731	732	735	752	754	762	776	782	795	801	808	831	833	838	842	845	854	871	880	899	902	907	912	914	927	931	940	943	947	957	960	975	983	986	989	991	1017; 7	11	26	40	47	62	64	72	75	80	81	96	100	121	125	131	144	146	149	162	175	176	197	199	202	204	211	220	223	225	238	241	244	248	257	282	294	297	301	302	322	329	358	393	396	397	405	441	443	447	452	466	469	479	486	500	515	517	520	521	523	531	536	537	547	548	550	553	567	576	579	580	593	595	596	600	605	610	632	641	666	676	701	713	717	718	765	766	772	774	793	819	836	857	858	859	862	863	885	893	897	905	906	916	917	919	921	924	935	944	953	954	956	963	971	972	985	987	1009	1016; 5	31	38	50	51	56	58	70	71	83	99	101	111	114	122	134	151	168	180	190	194	203	219	227	242	261	265	272	276	278	288	289	292	299	314	319	335	344	351	357	360	368	377	383	385	387	390	400	401	407	411	412	436	438	445	463	481	483	489	493	498	502	519	526	529	539	543	549	552	557	566	577	578	581	587	594	620	622	630	644	645	652	660	667	672	675	681	692	696	708	734	738	745	760	787	788	798	800	802	830	843	848	849	850	878	884	889	901	913	928	950	955	967	984	993	994	1007	1008	1013	1019];


