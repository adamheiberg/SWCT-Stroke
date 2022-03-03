%% Export HRF and analyze in MATLAB

%% Export

% Data-format: dcAvg( tid , hæmoglobiner , kanaler , stimuli )

% Kanalconfiguration:   (letter=source, number = receiver)
% 1: A1     7: C5
% 2: A2     8: C6
% 3: A3     9: C7
% 4: B2    10: D6
% 5: B3    11: D7
% 6: B4    12: D8

% Stimuli:
% 1: Congruent       2: Incongruent     3: Neutral

load('groupResults.mat')

%% Establish parameters

% Patient-ID
[file,path] = uigetfile('*.nirs');
subj = strsplit(file,'_');
subj=subj{1};

subj_row = strsplit(file,'.');
subj_row = ['Subj_' subj_row{1}] ;
subj_row = find(strcmp({group.subjs.name}, subj_row)==1) ;

% Establish timeline from HRF
t = group.subjs(subj_row).procStream.output.dcAvg.GetTime();
t0=find(t==0);

% Loading test performance
load('Behavior_Control.mat')
tslut_NEU = Behavior_Control.RT_NEUTRAL(subj_row)*5 + 15 ;
[~,~,idx]=unique(abs(t - tslut_NEU)); 
tslut_NEU=find(t==t(idx==1)); 

tslut_CON = Behavior_Control.RT_CONGRUENT(subj_row)*5 + 15 ;
[~,~,idx]=unique(abs(t - tslut_CON)); 
tslut_CON=find(t==t(idx==1)); 

tslut_INC = Behavior_Control.RT_INCONGRUENT(subj_row)*5 + 15 ;
[~,~,idx]=unique(abs(t - tslut_INC)); 
tslut_INC=find(t==t(idx==1)); 

% Import HRF and rename dcAvg from HRF
group.subjs(subj_row).Load();
dcAvg = group.subjs(subj_row).procStream.output.dcAvg.GetDataTimeSeries('reshape');

dcSD = abs(group.subjs(subj_row).procStream.output.dcAvgStd.GetDataTimeSeries('reshape'));

%% NEUTRAL OXY
% Finding peak, time-to-peak (TTP), return to baseline (low) and
% time-to-low (TTL)

% Max. oxy during NEUTRAL
max_oxy_NEU_R_latUp = max(dcAvg(t0:tslut_NEU,1,2,3))-dcAvg(t0,1,2,3); 
max_oxy_NEU_R_latDown = max(dcAvg(t0:tslut_NEU,1,3,3))-dcAvg(t0,1,3,3);
max_oxy_NEU_R_medUp = max(dcAvg(t0:tslut_NEU,1,4,3))-dcAvg(t0,1,4,3);
max_oxy_NEU_R_medDown = max(dcAvg(t0:tslut_NEU,1,5,3))-dcAvg(t0,1,5,3); 
max_oxy_NEU_L_medUp = max(dcAvg(t0:tslut_NEU,1,8,3))-dcAvg(t0,1,8,3);
max_oxy_NEU_L_medDown = max(dcAvg(t0:tslut_NEU,1,9,3))-dcAvg(t0,1,9,3);
max_oxy_NEU_L_latUp = max(dcAvg(t0:tslut_NEU,1,10,3))-dcAvg(t0,1,10,3);
max_oxy_NEU_L_latDown = max(dcAvg(t0:tslut_NEU,1,11,3))-dcAvg(t0,1,11,3);

% Min. oxy during NEUTRAL
min_oxy_NEU_R_latUp = min(dcAvg(t0:tslut_NEU,1,2,3))-dcAvg(t0,1,2,3);
min_oxy_NEU_R_latDown = min(dcAvg(t0:tslut_NEU,1,3,3))-dcAvg(t0,1,3,3);
min_oxy_NEU_R_medUp = min(dcAvg(t0:tslut_NEU,1,4,3))-dcAvg(t0,1,4,3);
min_oxy_NEU_R_medDown = min(dcAvg(t0:tslut_NEU,1,5,3))-dcAvg(t0,1,5,3);
min_oxy_NEU_L_medUp = min(dcAvg(t0:tslut_NEU,1,8,3))-dcAvg(t0,1,8,3);
min_oxy_NEU_L_medDown = min(dcAvg(t0:tslut_NEU,1,9,3))-dcAvg(t0,1,9,3);
min_oxy_NEU_L_latUp = min(dcAvg(t0:tslut_NEU,1,10,3))-dcAvg(t0,1,10,3);
min_oxy_NEU_L_latDown = min(dcAvg(t0:tslut_NEU,1,11,3))-dcAvg(t0,1,11,3);

% Finding peak (min. or max.) and time-to-peak during NEUTRAL 
% If identical, then use the first one
if abs(max_oxy_NEU_R_latUp) > abs(min_oxy_NEU_R_latUp)
   peak_oxy_NEU_R_latUp = max_oxy_NEU_R_latUp ;
   ttp_oxy_NEU_R_latUp = t(dcAvg (:,1,2,3) == max(dcAvg(t0:tslut_NEU,1,2,3))) ;
elseif abs(max_oxy_NEU_R_latUp) < abs(min_oxy_NEU_R_latUp)
   peak_oxy_NEU_R_latUp = min_oxy_NEU_R_latUp ;
   ttp_oxy_NEU_R_latUp = t( dcAvg (:,1,2,3) == min(dcAvg(t0:tslut_NEU,1,2,3))) ;
elseif isnan(max_oxy_NEU_R_latUp)
   ttp_oxy_NEU_R_latUp = NaN;
   peak_oxy_NEU_R_latUp = NaN ;
elseif abs(max_oxy_NEU_R_latUp) == abs(min_oxy_NEU_R_latUp)
   ttp_oxy_NEU_R_latUp_max = t( dcAvg (:,1,2,3) == max(dcAvg(t0:tslut_NEU,1,2,3))) ;
   ttp_oxy_NEU_R_latUp_min = t( dcAvg (:,1,2,3) == min(dcAvg(t0:tslut_NEU,1,2,3))) ;
   if ttp_oxy_NEU_R_latUp_max < ttp_oxy_NEU_R_latUp_min
       ttp_oxy_NEU_R_latUp = ttp_oxy_NEU_R_latUp_max;
       peak_oxy_NEU_R_latUp = max_oxy_NEU_R_latUp ; 
   elseif ttp_oxy_NEU_R_latUp_max > ttp_oxy_NEU_R_latUp_min
       ttp_oxy_NEU_R_latUp = ttp_oxy_NEU_R_latUp_min;
       peak_oxy_NEU_R_latUp = min_oxy_NEU_R_latUp ;
   end
end

ttp_row = find(t==ttp_oxy_NEU_R_latUp);
if isnan(peak_oxy_NEU_R_latUp)
    ttl_oxy_NEU_R_latUp = NaN ;
    low_oxy_NEU_R_latUp = NaN ;
elseif peak_oxy_NEU_R_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,2,3) < dcAvg(t0,1,2,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,2,3),[ttp_oxy_NEU_R_latUp 25],[dcAvg(t0,1,2,3) dcAvg(t0,1,2,3)]) ;
        ttl_oxy_NEU_R_latUp = min(xlin);
        low_oxy_NEU_R_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,2,3) < dcAvg(t0,1,2,3)))
        ttl_oxy_NEU_R_latUp = t( find(dcAvg(ttp_row:tslut_NEU,1,2,3)==min(dcAvg(ttp_row:tslut_NEU,1,2,3)))+ttp_row-1);
        low_oxy_NEU_R_latUp = min(dcAvg(ttp_row:tslut_NEU,1,2,3))-dcAvg(t0,1,2,3);
    end
elseif peak_oxy_NEU_R_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,2,3) > dcAvg(t0,1,2,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,2,3),[ttp_oxy_NEU_R_latUp 25],[dcAvg(t0,1,2,3) dcAvg(t0,1,2,3)]) ;
        ttl_oxy_NEU_R_latUp = min(xlin);
        low_oxy_NEU_R_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,2,3) > dcAvg(t0,1,2,3)))
        ttl_oxy_NEU_R_latUp = t( find(dcAvg(ttp_row:tslut_NEU,1,2,3)==max(dcAvg(ttp_row:tslut_NEU,1,2,3)))+ttp_row-1);
        low_oxy_NEU_R_latUp = max(dcAvg(ttp_row:tslut_NEU,1,2,3))-dcAvg(t0,1,2,3);
    end
end

if abs(max_oxy_NEU_R_latDown) > abs(min_oxy_NEU_R_latDown)
   peak_oxy_NEU_R_latDown = max_oxy_NEU_R_latDown ;
   ttp_oxy_NEU_R_latDown = t( dcAvg (:,1,3,3) == max(dcAvg(t0:tslut_NEU,1,3,3))) ;
elseif abs(max_oxy_NEU_R_latDown) < abs(min_oxy_NEU_R_latDown)
   peak_oxy_NEU_R_latDown = min_oxy_NEU_R_latDown ;
   ttp_oxy_NEU_R_latDown = t( dcAvg (:,1,3,3) == min(dcAvg(t0:tslut_NEU,1,3,3))) ;
elseif isnan(max_oxy_NEU_R_latDown)
   peak_oxy_NEU_R_latDown = NaN;
   ttp_oxy_NEU_R_latDown = NaN;
elseif abs(max_oxy_NEU_R_latDown) == abs(min_oxy_NEU_R_latDown)
   ttp_oxy_NEU_R_latDown_max = t( dcAvg (:,1,3,3) == max(dcAvg(t0:tslut_NEU,1,3,3))) ;
   ttp_oxy_NEU_R_latDown_min = t( dcAvg (:,1,3,3) == min(dcAvg(t0:tslut_NEU,1,3,3))) ;
    if ttp_oxy_NEU_R_latDown_max < ttp_oxy_NEU_R_latDown_min
       ttp_oxy_NEU_R_latDown = ttp_oxy_NEU_R_latDown_max;
       peak_oxy_NEU_R_latDown = max_oxy_NEU_R_latDown ; 
    elseif ttp_oxy_NEU_R_latDown_max > ttp_oxy_NEU_R_latDown_min
       ttp_oxy_NEU_R_latDown = ttp_oxy_NEU_R_latDown_min;
       peak_oxy_NEU_R_latDown = min_oxy_NEU_R_latDown ;
    end
end

ttp_row = find(t==ttp_oxy_NEU_R_latDown);
if isnan(peak_oxy_NEU_R_latDown)
    ttl_oxy_NEU_R_latDown = NaN ;
    low_oxy_NEU_R_latDown = NaN ;
elseif peak_oxy_NEU_R_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,3,3) < dcAvg(t0,1,3,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,3,3),[ttp_oxy_NEU_R_latDown 25],[dcAvg(t0,1,3,3) dcAvg(t0,1,3,3)]) ;
        ttl_oxy_NEU_R_latDown = min(xlin);
        low_oxy_NEU_R_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,3,3) < dcAvg(t0,1,3,3)))
        ttl_oxy_NEU_R_latDown = t( find(dcAvg(ttp_row:tslut_NEU,1,3,3)==min(dcAvg(ttp_row:tslut_NEU,1,3,3)))+ttp_row-1);
        low_oxy_NEU_R_latDown = min(dcAvg(ttp_row:tslut_NEU,1,3,3))-dcAvg(t0,1,3,3);
    end
elseif peak_oxy_NEU_R_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,3,3) > dcAvg(t0,1,3,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,3,3),[ttp_oxy_NEU_R_latDown 25],[dcAvg(t0,1,3,3) dcAvg(t0,1,3,3)]) ;
        ttl_oxy_NEU_R_latDown = min(xlin);
        low_oxy_NEU_R_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,3,3) > dcAvg(t0,1,3,3)))
        ttl_oxy_NEU_R_latDown = t( find(dcAvg(ttp_row:tslut_NEU,1,3,3)==max(dcAvg(ttp_row:tslut_NEU,1,3,3)))+ttp_row-1);
        low_oxy_NEU_R_latDown = max(dcAvg(ttp_row:tslut_NEU,1,3,3))-dcAvg(t0,1,3,3);
    end
end

if abs(max_oxy_NEU_R_medUp) > abs(min_oxy_NEU_R_medUp)
   peak_oxy_NEU_R_medUp = max_oxy_NEU_R_medUp ;
   ttp_oxy_NEU_R_medUp = t( dcAvg (:,1,4,3) == max(dcAvg(t0:tslut_NEU,1,4,3))) ;
elseif abs(max_oxy_NEU_R_medUp) < abs(min_oxy_NEU_R_medUp)
   peak_oxy_NEU_R_medUp = min_oxy_NEU_R_medUp ;
   ttp_oxy_NEU_R_medUp = t( dcAvg (:,1,4,3) == min(dcAvg(t0:tslut_NEU,1,4,3))) ;
elseif isnan(max_oxy_NEU_R_medUp)
   ttp_oxy_NEU_R_medUp = NaN;
   peak_oxy_NEU_R_medUp = NaN ;
elseif abs(max_oxy_NEU_R_medUp) == abs(min_oxy_NEU_R_medUp)
   ttp_oxy_NEU_R_medUp_max = t( dcAvg (:,1,4,3) == max(dcAvg(t0:tslut_NEU,1,4,3))) ;
   ttp_oxy_NEU_R_medUp_min = t( dcAvg (:,1,4,3) == min(dcAvg(t0:tslut_NEU,1,4,3))) ;
   if ttp_oxy_NEU_R_medUp_max < ttp_oxy_NEU_R_medUp_min
       ttp_oxy_NEU_R_medUp = ttp_oxy_NEU_R_medUp_max;
       peak_oxy_NEU_R_medUp = max_oxy_NEU_R_medUp ; 
   elseif ttp_oxy_NEU_R_medUp_max > ttp_oxy_NEU_R_medUp_min
       ttp_oxy_NEU_R_medUp = ttp_oxy_NEU_R_medUp_min;
       peak_oxy_NEU_R_medUp = min_oxy_NEU_R_medUp ;
   end
end

ttp_row = find(t==ttp_oxy_NEU_R_medUp);
if isnan(peak_oxy_NEU_R_medUp)
    ttl_oxy_NEU_R_medUp = NaN ;
    low_oxy_NEU_R_medUp = NaN ;
elseif peak_oxy_NEU_R_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,4,3) < dcAvg(t0,1,4,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,4,3),[ttp_oxy_NEU_R_medUp 25],[dcAvg(t0,1,4,3) dcAvg(t0,1,4,3)]) ;
        ttl_oxy_NEU_R_medUp = min(xlin);
        low_oxy_NEU_R_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,4,3) < dcAvg(t0,1,4,3)))
        ttl_oxy_NEU_R_medUp = t( find(dcAvg(ttp_row:tslut_NEU,1,4,3)==min(dcAvg(ttp_row:tslut_NEU,1,4,3)))+ttp_row-1);
        low_oxy_NEU_R_medUp = min(dcAvg(ttp_row:tslut_NEU,1,4,3))-dcAvg(t0,1,4,3);
    end
elseif peak_oxy_NEU_R_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,4,3) > dcAvg(t0,1,4,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,4,3),[ttp_oxy_NEU_R_medUp 25],[dcAvg(t0,1,4,3) dcAvg(t0,1,4,3)]) ;
        ttl_oxy_NEU_R_medUp = min(xlin);
        low_oxy_NEU_R_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,4,3) > dcAvg(t0,1,4,3)))
        ttl_oxy_NEU_R_medUp = t( find(dcAvg(ttp_row:tslut_NEU,1,4,3)==max(dcAvg(ttp_row:tslut_NEU,1,4,3)))+ttp_row-1);
        low_oxy_NEU_R_medUp = max(dcAvg(ttp_row:tslut_NEU,1,4,3))-dcAvg(t0,1,4,3);
    end
end

if abs(max_oxy_NEU_R_medDown) > abs(min_oxy_NEU_R_medDown)
   peak_oxy_NEU_R_medDown = max_oxy_NEU_R_medDown ;
   ttp_oxy_NEU_R_medDown = t( dcAvg (:,1,5,3) == max(dcAvg(t0:tslut_NEU,1,5,3))) ;
elseif abs(max_oxy_NEU_R_medDown) < abs(min_oxy_NEU_R_medDown)
   peak_oxy_NEU_R_medDown = min_oxy_NEU_R_medDown ;
   ttp_oxy_NEU_R_medDown = t( dcAvg (:,1,5,3) == min(dcAvg(t0:tslut_NEU,1,5,3))) ;
elseif isnan(max_oxy_NEU_R_medDown)
   ttp_oxy_NEU_R_medDown = NaN;
   peak_oxy_NEU_R_medDown = NaN ;
elseif abs(max_oxy_NEU_R_medDown) == abs(min_oxy_NEU_R_medDown)
   ttp_oxy_NEU_R_medDown_max = t( dcAvg (:,1,5,3) == max(dcAvg(t0:tslut_NEU,1,5,3))) ;
   ttp_oxy_NEU_R_medDown_min = t( dcAvg (:,1,5,3) == min(dcAvg(t0:tslut_NEU,1,5,3))) ;
   if ttp_oxy_NEU_R_medDown_max < ttp_oxy_NEU_R_medDown_min
       ttp_oxy_NEU_R_medDown = ttp_oxy_NEU_R_medDown_max;
       peak_oxy_NEU_R_medDown = max_oxy_NEU_R_medDown ; 
   elseif ttp_oxy_NEU_R_medDown_max > ttp_oxy_NEU_R_medDown_min
       ttp_oxy_NEU_R_medDown = ttp_oxy_NEU_R_medDown_min;
       peak_oxy_NEU_R_medDown = min_oxy_NEU_R_medDown ;
   end
end

ttp_row = find(t==ttp_oxy_NEU_R_medDown);
if isnan(peak_oxy_NEU_R_medDown)
    ttl_oxy_NEU_R_medDown = NaN ;
    low_oxy_NEU_R_medDown = NaN ;
elseif peak_oxy_NEU_R_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,5,3) < dcAvg(t0,1,5,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,5,3),[ttp_oxy_NEU_R_medDown 25],[dcAvg(t0,1,5,3) dcAvg(t0,1,5,3)]) ;
        ttl_oxy_NEU_R_medDown = min(xlin);
        low_oxy_NEU_R_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,5,3) < dcAvg(t0,1,5,3)))
        ttl_oxy_NEU_R_medDown = t( find(dcAvg(ttp_row:tslut_NEU,1,5,3)==min(dcAvg(ttp_row:tslut_NEU,1,5,3)))+ttp_row-1);
        low_oxy_NEU_R_medDown = min(dcAvg(ttp_row:tslut_NEU,1,5,3))-dcAvg(t0,1,5,3);
    end
elseif peak_oxy_NEU_R_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,5,3) > dcAvg(t0,1,5,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,5,3),[ttp_oxy_NEU_R_medDown 25],[dcAvg(t0,1,5,3) dcAvg(t0,1,5,3)]) ;
        ttl_oxy_NEU_R_medDown = min(xlin);
        low_oxy_NEU_R_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,5,3) > dcAvg(t0,1,5,3)))
        ttl_oxy_NEU_R_medDown = t( find(dcAvg(ttp_row:tslut_NEU,1,5,3)==max(dcAvg(ttp_row:tslut_NEU,1,5,3)))+ttp_row-1);
        low_oxy_NEU_R_medDown = max(dcAvg(ttp_row:tslut_NEU,1,5,3))-dcAvg(t0,1,5,3);
    end
end

if abs(max_oxy_NEU_L_medUp) > abs(min_oxy_NEU_L_medUp)
   peak_oxy_NEU_L_medUp = max_oxy_NEU_L_medUp ;
   ttp_oxy_NEU_L_medUp = t( dcAvg (:,1,8,3) == max(dcAvg(t0:tslut_NEU,1,8,3))) ;
elseif abs(max_oxy_NEU_L_medUp) < abs(min_oxy_NEU_L_medUp)
   peak_oxy_NEU_L_medUp = min_oxy_NEU_L_medUp ;
   ttp_oxy_NEU_L_medUp = t( dcAvg (:,1,8,3) == min(dcAvg(t0:tslut_NEU,1,8,3))) ;
elseif isnan(max_oxy_NEU_L_medUp)
   ttp_oxy_NEU_L_medUp = NaN;
   peak_oxy_NEU_L_medUp = NaN ;
elseif abs(max_oxy_NEU_L_medUp) == abs(min_oxy_NEU_L_medUp)
   ttp_oxy_NEU_L_medUp_max = t( dcAvg (:,1,8,3) == max(dcAvg(t0:tslut_NEU,1,8,3))) ;
   ttp_oxy_NEU_L_medUp_min = t( dcAvg (:,1,8,3) == min(dcAvg(t0:tslut_NEU,1,8,3))) ;
   if ttp_oxy_NEU_L_medUp_max < ttp_oxy_NEU_L_medUp_min
       ttp_oxy_NEU_L_medUp = ttp_oxy_NEU_L_medUp_max;
       peak_oxy_NEU_L_medUp = max_oxy_NEU_L_medUp ; 
   elseif ttp_oxy_NEU_L_medUp_max > ttp_oxy_NEU_L_medUp_min
       ttp_oxy_NEU_L_medUp = ttp_oxy_NEU_L_medUp_min;
       peak_oxy_NEU_L_medUp = min_oxy_NEU_L_medUp ;
   end
end

ttp_row = find(t==ttp_oxy_NEU_L_medUp);
if isnan(peak_oxy_NEU_L_medUp)
    ttl_oxy_NEU_L_medUp = NaN ;
    low_oxy_NEU_L_medUp = NaN ;
elseif peak_oxy_NEU_L_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,8,3) < dcAvg(t0,1,8,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,8,3),[ttp_oxy_NEU_L_medUp 25],[dcAvg(t0,1,8,3) dcAvg(t0,1,8,3)]) ;
        ttl_oxy_NEU_L_medUp = min(xlin);
        low_oxy_NEU_L_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,8,3) < dcAvg(t0,1,8,3)))
        ttl_oxy_NEU_L_medUp = t( find(dcAvg(ttp_row:tslut_NEU,1,8,3)==min(dcAvg(ttp_row:tslut_NEU,1,8,3)))+ttp_row-1);
        low_oxy_NEU_L_medUp = min(dcAvg(ttp_row:tslut_NEU,1,8,3))-dcAvg(t0,1,8,3);
    end
elseif peak_oxy_NEU_L_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,8,3) > dcAvg(t0,1,8,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,8,3),[ttp_oxy_NEU_L_medUp 25],[dcAvg(t0,1,8,3) dcAvg(t0,1,8,3)]) ;
        ttl_oxy_NEU_L_medUp = min(xlin);
        low_oxy_NEU_L_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,8,3) > dcAvg(t0,1,8,3)))
        ttl_oxy_NEU_L_medUp = t( find(dcAvg(ttp_row:tslut_NEU,1,8,3)==max(dcAvg(ttp_row:tslut_NEU,1,8,3)))+ttp_row-1);
        low_oxy_NEU_L_medUp = max(dcAvg(ttp_row:tslut_NEU,1,8,3))-dcAvg(t0,1,8,3);
    end
end

if abs(max_oxy_NEU_L_medDown) > abs(min_oxy_NEU_L_medDown)
   peak_oxy_NEU_L_medDown = max_oxy_NEU_L_medDown ;
   ttp_oxy_NEU_L_medDown = t( dcAvg (:,1,9,3) == max(dcAvg(t0:tslut_NEU,1,9,3))) ;
elseif abs(max_oxy_NEU_L_medDown) < abs(min_oxy_NEU_L_medDown)
   peak_oxy_NEU_L_medDown = min_oxy_NEU_L_medDown ;
   ttp_oxy_NEU_L_medDown = t( dcAvg (:,1,9,3) == min(dcAvg(t0:tslut_NEU,1,9,3))) ;
elseif isnan(max_oxy_NEU_L_medDown)
   ttp_oxy_NEU_L_medDown = NaN;
   peak_oxy_NEU_L_medDown = NaN ;
elseif abs(max_oxy_NEU_L_medDown) == abs(min_oxy_NEU_L_medDown)
   ttp_oxy_NEU_L_medDown_max = t( dcAvg (:,1,9,3) == max(dcAvg(t0:tslut_NEU,1,9,3))) ;
   ttp_oxy_NEU_L_medDown_min = t( dcAvg (:,1,9,3) == min(dcAvg(t0:tslut_NEU,1,9,3))) ;
   if ttp_oxy_NEU_L_medDown_max < ttp_oxy_NEU_L_medDown_min
       ttp_oxy_NEU_L_medDown = ttp_oxy_NEU_L_medDown_max;
       peak_oxy_NEU_L_medDown = max_oxy_NEU_L_medDown ; 
   elseif ttp_oxy_NEU_L_medDown_max > ttp_oxy_NEU_L_medDown_min
       ttp_oxy_NEU_L_medDown = ttp_oxy_NEU_L_medDown_min;
       peak_oxy_NEU_L_medDown = min_oxy_NEU_L_medDown ;
   end
end

ttp_row = find(t==ttp_oxy_NEU_L_medDown);
if isnan(peak_oxy_NEU_L_medDown)
    ttl_oxy_NEU_L_medDown = NaN ;
    low_oxy_NEU_L_medDown = NaN ;
elseif peak_oxy_NEU_L_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,9,3) < dcAvg(t0,1,9,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,9,3),[ttp_oxy_NEU_L_medDown 25],[dcAvg(t0,1,9,3) dcAvg(t0,1,9,3)]) ;
        ttl_oxy_NEU_L_medDown = min(xlin);
        low_oxy_NEU_L_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,9,3) < dcAvg(t0,1,9,3)))
        ttl_oxy_NEU_L_medDown = t( find(dcAvg(ttp_row:tslut_NEU,1,9,3)==min(dcAvg(ttp_row:tslut_NEU,1,9,3)))+ttp_row-1);
        low_oxy_NEU_L_medDown = min(dcAvg(ttp_row:tslut_NEU,1,9,3))-dcAvg(t0,1,9,3);
    end
elseif peak_oxy_NEU_L_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,9,3) > dcAvg(t0,1,9,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,9,3),[ttp_oxy_NEU_L_medDown 25],[dcAvg(t0,1,9,3) dcAvg(t0,1,9,3)]) ;
        ttl_oxy_NEU_L_medDown = min(xlin);
        low_oxy_NEU_L_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,9,3) > dcAvg(t0,1,9,3)))
        ttl_oxy_NEU_L_medDown = t( find(dcAvg(ttp_row:tslut_NEU,1,9,3)==max(dcAvg(ttp_row:tslut_NEU,1,9,3)))+ttp_row-1);
        low_oxy_NEU_L_medDown = max(dcAvg(ttp_row:tslut_NEU,1,9,3))-dcAvg(t0,1,9,3);
    end
end

if abs(max_oxy_NEU_L_latUp) > abs(min_oxy_NEU_L_latUp)
   peak_oxy_NEU_L_latUp = max_oxy_NEU_L_latUp ;
   ttp_oxy_NEU_L_latUp = t( dcAvg (:,1,10,3) == max(dcAvg(t0:tslut_NEU,1,10,3))) ;
elseif abs(max_oxy_NEU_L_latUp) < abs(min_oxy_NEU_L_latUp)
   peak_oxy_NEU_L_latUp = min_oxy_NEU_L_latUp ;
   ttp_oxy_NEU_L_latUp = t( dcAvg (:,1,10,3) == min(dcAvg(t0:tslut_NEU,1,10,3))) ;
elseif isnan(max_oxy_NEU_L_latUp)
   ttp_oxy_NEU_L_latUp = NaN;
   peak_oxy_NEU_L_latUp = NaN ;
elseif abs(max_oxy_NEU_L_latUp) == abs(min_oxy_NEU_L_latUp)
   ttp_oxy_NEU_L_latUp_max = t( dcAvg (:,1,10,3) == max(dcAvg(t0:tslut_NEU,1,10,3))) ;
   ttp_oxy_NEU_L_latUp_min = t( dcAvg (:,1,10,3) == min(dcAvg(t0:tslut_NEU,1,10,3))) ;
   if ttp_oxy_NEU_L_latUp_max < ttp_oxy_NEU_L_latUp_min
       ttp_oxy_NEU_L_latUp = ttp_oxy_NEU_L_latUp_max;
       peak_oxy_NEU_L_latUp = max_oxy_NEU_L_latUp ; 
   elseif ttp_oxy_NEU_L_latUp_max > ttp_oxy_NEU_L_latUp_min
       ttp_oxy_NEU_L_latUp = ttp_oxy_NEU_L_latUp_min;
       peak_oxy_NEU_L_latUp = min_oxy_NEU_L_latUp ;
   end
end
 
ttp_row = find(t==ttp_oxy_NEU_L_latUp);
if isnan(peak_oxy_NEU_L_latUp)
    ttl_oxy_NEU_L_latUp = NaN ;
    low_oxy_NEU_L_latUp = NaN ;
elseif peak_oxy_NEU_L_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,10,3) < dcAvg(t0,1,10,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,10,3),[ttp_oxy_NEU_L_latUp 25],[dcAvg(t0,1,10,3) dcAvg(t0,1,10,3)]) ;
        ttl_oxy_NEU_L_latUp = min(xlin);
        low_oxy_NEU_L_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,10,3) < dcAvg(t0,1,10,3)))
        ttl_oxy_NEU_L_latUp = t( find(dcAvg(ttp_row:tslut_NEU,1,10,3)==min(dcAvg(ttp_row:tslut_NEU,1,10,3)))+ttp_row-1);
        low_oxy_NEU_L_latUp = min(dcAvg(ttp_row:tslut_NEU,1,10,3))-dcAvg(t0,1,10,3);
    end
elseif peak_oxy_NEU_L_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,10,3) > dcAvg(t0,1,10,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,10,3),[ttp_oxy_NEU_L_latUp 25],[dcAvg(t0,1,10,3) dcAvg(t0,1,10,3)]) ;
        ttl_oxy_NEU_L_latUp = min(xlin);
        low_oxy_NEU_L_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,10,3) > dcAvg(t0,1,10,3)))
        ttl_oxy_NEU_L_latUp = t( find(dcAvg(ttp_row:tslut_NEU,1,10,3)==max(dcAvg(ttp_row:tslut_NEU,1,10,3)))+ttp_row-1);
        low_oxy_NEU_L_latUp = max(dcAvg(ttp_row:tslut_NEU,1,10,3))-dcAvg(t0,1,10,3);
    end
end

if abs(max_oxy_NEU_L_latDown) > abs(min_oxy_NEU_L_latDown)
   peak_oxy_NEU_L_latDown = max_oxy_NEU_L_latDown ;
   ttp_oxy_NEU_L_latDown = t( dcAvg (:,1,11,3) == max(dcAvg(t0:tslut_NEU,1,11,3))) ;
elseif abs(max_oxy_NEU_L_latDown) < abs(min_oxy_NEU_L_latDown)
   peak_oxy_NEU_L_latDown = min_oxy_NEU_L_latDown ;
   ttp_oxy_NEU_L_latDown = t( dcAvg (:,1,11,3) == min(dcAvg(t0:tslut_NEU,1,11,3))) ;
elseif isnan(max_oxy_NEU_L_latDown)
   ttp_oxy_NEU_L_latDown = NaN;
   peak_oxy_NEU_L_latDown = NaN ;
elseif abs(max_oxy_NEU_L_latDown) == abs(min_oxy_NEU_L_latDown)
   ttp_oxy_NEU_L_latDown_max = t( dcAvg (:,1,11,3) == max(dcAvg(t0:tslut_NEU,1,11,3))) ;
   ttp_oxy_NEU_L_latDown_min = t( dcAvg (:,1,11,3) == min(dcAvg(t0:tslut_NEU,1,11,3))) ;
   if ttp_oxy_NEU_L_latDown_max < ttp_oxy_NEU_L_latDown_min
       ttp_oxy_NEU_L_latDown = ttp_oxy_NEU_L_latDown_max;
       peak_oxy_NEU_L_latDown = max_oxy_NEU_L_latDown ; 
   elseif ttp_oxy_NEU_L_latDown_max > ttp_oxy_NEU_L_latDown_min
       ttp_oxy_NEU_L_latDown = ttp_oxy_NEU_L_latDown_min;
       peak_oxy_NEU_L_latDown = min_oxy_NEU_L_latDown ;
   end
end

ttp_row = find(t==ttp_oxy_NEU_L_latDown);
if isnan(peak_oxy_NEU_L_latDown)
    ttl_oxy_NEU_L_latDown = NaN ;
    low_oxy_NEU_L_latDown = NaN ;
elseif peak_oxy_NEU_L_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,11,3) < dcAvg(t0,1,11,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,11,3),[ttp_oxy_NEU_L_latDown 25],[dcAvg(t0,1,11,3) dcAvg(t0,1,11,3)]) ;
        ttl_oxy_NEU_L_latDown = min(xlin);
        low_oxy_NEU_L_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,11,3) < dcAvg(t0,1,11,3)))
        ttl_oxy_NEU_L_latDown = t( find(dcAvg(ttp_row:tslut_NEU,1,11,3)==min(dcAvg(ttp_row:tslut_NEU,1,11,3)))+ttp_row-1);
        low_oxy_NEU_L_latDown = min(dcAvg(ttp_row:tslut_NEU,1,11,3))-dcAvg(t0,1,11,3);
    end
elseif peak_oxy_NEU_L_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,1,11,3) > dcAvg(t0,1,11,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,1,11,3),[ttp_oxy_NEU_L_latDown 25],[dcAvg(t0,1,11,3) dcAvg(t0,1,11,3)]) ;
        ttl_oxy_NEU_L_latDown = min(xlin);
        low_oxy_NEU_L_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,1,11,3) > dcAvg(t0,1,11,3)))
        ttl_oxy_NEU_L_latDown = t( find(dcAvg(ttp_row:tslut_NEU,1,11,3)==max(dcAvg(ttp_row:tslut_NEU,1,11,3)))+ttp_row-1);
        low_oxy_NEU_L_latDown = max(dcAvg(ttp_row:tslut_NEU,1,11,3))-dcAvg(t0,1,11,3);
    end
end

%% NEUTRAL DEOXY
% Finding peak, time-to-peak (TTP), return to baseline (low) and
% time-to-low (TTL)

% Max. deoxy during NEUTRAL
max_deoxy_NEU_R_latUp = max(dcAvg(t0:tslut_NEU,2,2,3))-dcAvg(t0,2,2,3); 
max_deoxy_NEU_R_latDown = max(dcAvg(t0:tslut_NEU,2,3,3))-dcAvg(t0,2,3,3);
max_deoxy_NEU_R_medUp = max(dcAvg(t0:tslut_NEU,2,4,3))-dcAvg(t0,2,4,3);
max_deoxy_NEU_R_medDown = max(dcAvg(t0:tslut_NEU,2,5,3))-dcAvg(t0,2,5,3); 
max_deoxy_NEU_L_medUp = max(dcAvg(t0:tslut_NEU,2,8,3))-dcAvg(t0,2,8,3);
max_deoxy_NEU_L_medDown = max(dcAvg(t0:tslut_NEU,2,9,3))-dcAvg(t0,2,9,3);
max_deoxy_NEU_L_latUp = max(dcAvg(t0:tslut_NEU,2,10,3))-dcAvg(t0,2,10,3);
max_deoxy_NEU_L_latDown = max(dcAvg(t0:tslut_NEU,2,11,3))-dcAvg(t0,2,11,3);

% Min. deoxy during NEUTRAL
min_deoxy_NEU_R_latUp = min(dcAvg(t0:tslut_NEU,2,2,3))-dcAvg(t0,2,2,3);
min_deoxy_NEU_R_latDown = min(dcAvg(t0:tslut_NEU,2,3,3))-dcAvg(t0,2,3,3);
min_deoxy_NEU_R_medUp = min(dcAvg(t0:tslut_NEU,2,4,3))-dcAvg(t0,2,4,3);
min_deoxy_NEU_R_medDown = min(dcAvg(t0:tslut_NEU,2,5,3))-dcAvg(t0,2,5,3);
min_deoxy_NEU_L_medUp = min(dcAvg(t0:tslut_NEU,2,8,3))-dcAvg(t0,2,8,3);
min_deoxy_NEU_L_medDown = min(dcAvg(t0:tslut_NEU,2,9,3))-dcAvg(t0,2,9,3);
min_deoxy_NEU_L_latUp = min(dcAvg(t0:tslut_NEU,2,10,3))-dcAvg(t0,2,10,3);
min_deoxy_NEU_L_latDown = min(dcAvg(t0:tslut_NEU,2,11,3))-dcAvg(t0,2,11,3);

% Finding peak (min. or max.) and time-to-peak during NEUTRAL 
% If identical, then use the first one
if abs(max_deoxy_NEU_R_latUp) > abs(min_deoxy_NEU_R_latUp)
   peak_deoxy_NEU_R_latUp = max_deoxy_NEU_R_latUp ;
   ttp_deoxy_NEU_R_latUp = t( dcAvg (:,2,2,3) == max(dcAvg(t0:tslut_NEU,2,2,3))) ;
elseif abs(max_deoxy_NEU_R_latUp) < abs(min_deoxy_NEU_R_latUp)
   peak_deoxy_NEU_R_latUp = min_deoxy_NEU_R_latUp ;
   ttp_deoxy_NEU_R_latUp = t( dcAvg (:,2,2,3) == min(dcAvg(t0:tslut_NEU,2,2,3))) ;
elseif isnan(max_deoxy_NEU_R_latUp)
   ttp_deoxy_NEU_R_latUp = NaN;
   peak_deoxy_NEU_R_latUp = NaN ;
elseif abs(max_deoxy_NEU_R_latUp) == abs(min_deoxy_NEU_R_latUp)
   ttp_deoxy_NEU_R_latUp_max = t( dcAvg (:,2,2,3) == max(dcAvg(t0:tslut_NEU,2,2,3))) ;
   ttp_deoxy_NEU_R_latUp_min = t( dcAvg (:,2,2,3) == min(dcAvg(t0:tslut_NEU,2,2,3))) ;
   if ttp_deoxy_NEU_R_latUp_max < ttp_deoxy_NEU_R_latUp_min
       ttp_deoxy_NEU_R_latUp = ttp_deoxy_NEU_R_latUp_max;
       peak_deoxy_NEU_R_latUp = max_deoxy_NEU_R_latUp ; 
   elseif ttp_deoxy_NEU_R_latUp_max > ttp_deoxy_NEU_R_latUp_min
       ttp_deoxy_NEU_R_latUp = ttp_deoxy_NEU_R_latUp_min;
       peak_deoxy_NEU_R_latUp = min_deoxy_NEU_R_latUp ;
   end
end
  
ttp_row = find(t==ttp_deoxy_NEU_R_latUp);
if isnan(peak_deoxy_NEU_R_latUp)
    ttl_deoxy_NEU_R_latUp = NaN ;
    low_deoxy_NEU_R_latUp = NaN ;
elseif peak_deoxy_NEU_R_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,2,3) < dcAvg(t0,2,2,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,2,3),[ttp_deoxy_NEU_R_latUp 25],[dcAvg(t0,2,2,3) dcAvg(t0,2,2,3)]) ;
        ttl_deoxy_NEU_R_latUp = min(xlin);
        low_deoxy_NEU_R_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,2,3) < dcAvg(t0,2,2,3)))
        ttl_deoxy_NEU_R_latUp = t( find(dcAvg(ttp_row:tslut_NEU,2,2,3)==min(dcAvg(ttp_row:tslut_NEU,2,2,3)))+ttp_row-1);
        low_deoxy_NEU_R_latUp = min(dcAvg(ttp_row:tslut_NEU,2,2,3))-dcAvg(t0,2,2,3);
    end
elseif peak_deoxy_NEU_R_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,2,3) > dcAvg(t0,2,2,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,2,3),[ttp_deoxy_NEU_R_latUp 25],[dcAvg(t0,2,2,3) dcAvg(t0,2,2,3)]) ;
        ttl_deoxy_NEU_R_latUp = min(xlin);
        low_deoxy_NEU_R_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,2,3) > dcAvg(t0,2,2,3)))
        ttl_deoxy_NEU_R_latUp = t( find(dcAvg(ttp_row:tslut_NEU,2,2,3)==max(dcAvg(ttp_row:tslut_NEU,2,2,3)))+ttp_row-1);
        low_deoxy_NEU_R_latUp = max(dcAvg(ttp_row:tslut_NEU,2,2,3))-dcAvg(t0,2,2,3);
    end
end

if abs(max_deoxy_NEU_R_latDown) > abs(min_deoxy_NEU_R_latDown)
   peak_deoxy_NEU_R_latDown = max_deoxy_NEU_R_latDown ;
   ttp_deoxy_NEU_R_latDown = t( dcAvg (:,2,3,3) == max(dcAvg(t0:tslut_NEU,2,3,3))) ;
elseif abs(max_deoxy_NEU_R_latDown) < abs(min_deoxy_NEU_R_latDown)
   peak_deoxy_NEU_R_latDown = min_deoxy_NEU_R_latDown ;
   ttp_deoxy_NEU_R_latDown = t( dcAvg (:,2,3,3) == min(dcAvg(t0:tslut_NEU,2,3,3))) ;
elseif isnan(max_deoxy_NEU_R_latDown)
   peak_deoxy_NEU_R_latDown = NaN;
   ttp_deoxy_NEU_R_latDown = NaN;
elseif abs(max_deoxy_NEU_R_latDown) == abs(min_deoxy_NEU_R_latDown)
   ttp_deoxy_NEU_R_latDown_max = t( dcAvg (:,2,3,3) == max(dcAvg(t0:tslut_NEU,2,3,3))) ;
   ttp_deoxy_NEU_R_latDown_min = t( dcAvg (:,2,3,3) == min(dcAvg(t0:tslut_NEU,2,3,3))) ;
    if ttp_deoxy_NEU_R_latDown_max < ttp_deoxy_NEU_R_latDown_min
       ttp_deoxy_NEU_R_latDown = ttp_deoxy_NEU_R_latDown_max;
       peak_deoxy_NEU_R_latDown = max_deoxy_NEU_R_latDown ; 
    elseif ttp_deoxy_NEU_R_latDown_max > ttp_deoxy_NEU_R_latDown_min
       ttp_deoxy_NEU_R_latDown = ttp_deoxy_NEU_R_latDown_min;
       peak_deoxy_NEU_R_latDown = min_deoxy_NEU_R_latDown ;
    end
end

ttp_row = find(t==ttp_deoxy_NEU_R_latDown);
if isnan(peak_deoxy_NEU_R_latDown)
    ttl_deoxy_NEU_R_latDown = NaN ;
    low_deoxy_NEU_R_latDown = NaN ;
elseif peak_deoxy_NEU_R_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,3,3) < dcAvg(t0,2,3,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,3,3),[ttp_deoxy_NEU_R_latDown 25],[dcAvg(t0,2,3,3) dcAvg(t0,2,3,3)]) ;
        ttl_deoxy_NEU_R_latDown = min(xlin);
        low_deoxy_NEU_R_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,3,3) < dcAvg(t0,2,3,3)))
        ttl_deoxy_NEU_R_latDown = t( find(dcAvg(ttp_row:tslut_NEU,2,3,3)==min(dcAvg(ttp_row:tslut_NEU,2,3,3)))+ttp_row-1);
        low_deoxy_NEU_R_latDown = min(dcAvg(ttp_row:tslut_NEU,2,3,3))-dcAvg(t0,2,3,3);
    end
elseif peak_deoxy_NEU_R_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,3,3) > dcAvg(t0,2,3,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,3,3),[ttp_deoxy_NEU_R_latDown 25],[dcAvg(t0,2,3,3) dcAvg(t0,2,3,3)]) ;
        ttl_deoxy_NEU_R_latDown = min(xlin);
        low_deoxy_NEU_R_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,3,3) > dcAvg(t0,2,3,3)))
        ttl_deoxy_NEU_R_latDown = t( find(dcAvg(ttp_row:tslut_NEU,2,3,3)==max(dcAvg(ttp_row:tslut_NEU,2,3,3)))+ttp_row-1);
        low_deoxy_NEU_R_latDown = max(dcAvg(ttp_row:tslut_NEU,2,3,3))-dcAvg(t0,2,3,3);
    end
end

if abs(max_deoxy_NEU_R_medUp) > abs(min_deoxy_NEU_R_medUp)
   peak_deoxy_NEU_R_medUp = max_deoxy_NEU_R_medUp ;
   ttp_deoxy_NEU_R_medUp = t( dcAvg (:,2,4,3) == max(dcAvg(t0:tslut_NEU,2,4,3))) ;
elseif abs(max_deoxy_NEU_R_medUp) < abs(min_deoxy_NEU_R_medUp)
   peak_deoxy_NEU_R_medUp = min_deoxy_NEU_R_medUp ;
   ttp_deoxy_NEU_R_medUp = t( dcAvg (:,2,4,3) == min(dcAvg(t0:tslut_NEU,2,4,3))) ;
elseif isnan(max_deoxy_NEU_R_medUp)
   ttp_deoxy_NEU_R_medUp = NaN;
   peak_deoxy_NEU_R_medUp = NaN ;
elseif abs(max_deoxy_NEU_R_medUp) == abs(min_deoxy_NEU_R_medUp)
   ttp_deoxy_NEU_R_medUp_max = t( dcAvg (:,2,4,3) == max(dcAvg(t0:tslut_NEU,2,4,3))) ;
   ttp_deoxy_NEU_R_medUp_min = t( dcAvg (:,2,4,3) == min(dcAvg(t0:tslut_NEU,2,4,3))) ;
   if ttp_deoxy_NEU_R_medUp_max < ttp_deoxy_NEU_R_medUp_min
       ttp_deoxy_NEU_R_medUp = ttp_deoxy_NEU_R_medUp_max;
       peak_deoxy_NEU_R_medUp = max_deoxy_NEU_R_medUp ; 
   elseif ttp_deoxy_NEU_R_medUp_max > ttp_deoxy_NEU_R_medUp_min
       ttp_deoxy_NEU_R_medUp = ttp_deoxy_NEU_R_medUp_min;
       peak_deoxy_NEU_R_medUp = min_deoxy_NEU_R_medUp ;
   end
end

ttp_row = find(t==ttp_deoxy_NEU_R_medUp);
if isnan(peak_deoxy_NEU_R_medUp)
    ttl_deoxy_NEU_R_medUp = NaN ;
    low_deoxy_NEU_R_medUp = NaN ;
elseif peak_deoxy_NEU_R_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,4,3) < dcAvg(t0,2,4,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,4,3),[ttp_deoxy_NEU_R_medUp 25],[dcAvg(t0,2,4,3) dcAvg(t0,2,4,3)]) ;
        ttl_deoxy_NEU_R_medUp = min(xlin);
        low_deoxy_NEU_R_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,4,3) < dcAvg(t0,2,4,3)))
        ttl_deoxy_NEU_R_medUp = t( find(dcAvg(ttp_row:tslut_NEU,2,4,3)==min(dcAvg(ttp_row:tslut_NEU,2,4,3)))+ttp_row-1);
        low_deoxy_NEU_R_medUp = min(dcAvg(ttp_row:tslut_NEU,2,4,3))-dcAvg(t0,2,4,3);
    end
elseif peak_deoxy_NEU_R_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,4,3) > dcAvg(t0,2,4,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,4,3),[ttp_deoxy_NEU_R_medUp 25],[dcAvg(t0,2,4,3) dcAvg(t0,2,4,3)]) ;
        ttl_deoxy_NEU_R_medUp = min(xlin);
        low_deoxy_NEU_R_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,4,3) > dcAvg(t0,2,4,3)))
        ttl_deoxy_NEU_R_medUp = t( find(dcAvg(ttp_row:tslut_NEU,2,4,3)==max(dcAvg(ttp_row:tslut_NEU,2,4,3)))+ttp_row-1);
        low_deoxy_NEU_R_medUp = max(dcAvg(ttp_row:tslut_NEU,2,4,3))-dcAvg(t0,2,4,3);
    end
end

if abs(max_deoxy_NEU_R_medDown) > abs(min_deoxy_NEU_R_medDown)
   peak_deoxy_NEU_R_medDown = max_deoxy_NEU_R_medDown ;
   ttp_deoxy_NEU_R_medDown = t( dcAvg (:,2,5,3) == max(dcAvg(t0:tslut_NEU,2,5,3))) ;
elseif abs(max_deoxy_NEU_R_medDown) < abs(min_deoxy_NEU_R_medDown)
   peak_deoxy_NEU_R_medDown = min_deoxy_NEU_R_medDown ;
   ttp_deoxy_NEU_R_medDown = t( dcAvg (:,2,5,3) == min(dcAvg(t0:tslut_NEU,2,5,3))) ;
elseif isnan(max_deoxy_NEU_R_medDown)
   ttp_deoxy_NEU_R_medDown = NaN;
   peak_deoxy_NEU_R_medDown = NaN ;
elseif abs(max_deoxy_NEU_R_medDown) == abs(min_deoxy_NEU_R_medDown)
   ttp_deoxy_NEU_R_medDown_max = t( dcAvg (:,2,5,3) == max(dcAvg(t0:tslut_NEU,2,5,3))) ;
   ttp_deoxy_NEU_R_medDown_min = t( dcAvg (:,2,5,3) == min(dcAvg(t0:tslut_NEU,2,5,3))) ;
   if ttp_deoxy_NEU_R_medDown_max < ttp_deoxy_NEU_R_medDown_min
       ttp_deoxy_NEU_R_medDown = ttp_deoxy_NEU_R_medDown_max;
       peak_deoxy_NEU_R_medDown = max_deoxy_NEU_R_medDown ; 
   elseif ttp_deoxy_NEU_R_medDown_max > ttp_deoxy_NEU_R_medDown_min
       ttp_deoxy_NEU_R_medDown = ttp_deoxy_NEU_R_medDown_min;
       peak_deoxy_NEU_R_medDown = min_deoxy_NEU_R_medDown ;
   end
end

ttp_row = find(t==ttp_deoxy_NEU_R_medDown);
if isnan(peak_deoxy_NEU_R_medDown)
    ttl_deoxy_NEU_R_medDown = NaN ;
    low_deoxy_NEU_R_medDown = NaN ;
elseif peak_deoxy_NEU_R_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,5,3) < dcAvg(t0,2,5,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,5,3),[ttp_deoxy_NEU_R_medDown 25],[dcAvg(t0,2,5,3) dcAvg(t0,2,5,3)]) ;
        ttl_deoxy_NEU_R_medDown = min(xlin);
        low_deoxy_NEU_R_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,5,3) < dcAvg(t0,2,5,3)))
        ttl_deoxy_NEU_R_medDown = t( find(dcAvg(ttp_row:tslut_NEU,2,5,3)==min(dcAvg(ttp_row:tslut_NEU,2,5,3)))+ttp_row-1);
        low_deoxy_NEU_R_medDown = min(dcAvg(ttp_row:tslut_NEU,2,5,3))-dcAvg(t0,2,5,3);
    end
elseif peak_deoxy_NEU_R_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,5,3) > dcAvg(t0,2,5,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,5,3),[ttp_deoxy_NEU_R_medDown 25],[dcAvg(t0,2,5,3) dcAvg(t0,2,5,3)]) ;
        ttl_deoxy_NEU_R_medDown = min(xlin);
        low_deoxy_NEU_R_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,5,3) > dcAvg(t0,2,5,3)))
        ttl_deoxy_NEU_R_medDown = t( find(dcAvg(ttp_row:tslut_NEU,2,5,3)==max(dcAvg(ttp_row:tslut_NEU,2,5,3)))+ttp_row-1);
        low_deoxy_NEU_R_medDown = max(dcAvg(ttp_row:tslut_NEU,2,5,3))-dcAvg(t0,2,5,3);
    end
end

if abs(max_deoxy_NEU_L_medUp) > abs(min_deoxy_NEU_L_medUp)
   peak_deoxy_NEU_L_medUp = max_deoxy_NEU_L_medUp ;
   ttp_deoxy_NEU_L_medUp = t( dcAvg (:,2,8,3) == max(dcAvg(t0:tslut_NEU,2,8,3))) ;
elseif abs(max_deoxy_NEU_L_medUp) < abs(min_deoxy_NEU_L_medUp)
   peak_deoxy_NEU_L_medUp = min_deoxy_NEU_L_medUp ;
   ttp_deoxy_NEU_L_medUp = t( dcAvg (:,2,8,3) == min(dcAvg(t0:tslut_NEU,2,8,3))) ;
elseif isnan(max_deoxy_NEU_L_medUp)
   ttp_deoxy_NEU_L_medUp = NaN;
   peak_deoxy_NEU_L_medUp = NaN ;
elseif abs(max_deoxy_NEU_L_medUp) == abs(min_deoxy_NEU_L_medUp)
   ttp_deoxy_NEU_L_medUp_max = t( dcAvg (:,2,8,3) == max(dcAvg(t0:tslut_NEU,2,8,3))) ;
   ttp_deoxy_NEU_L_medUp_min = t( dcAvg (:,2,8,3) == min(dcAvg(t0:tslut_NEU,2,8,3))) ;
   if ttp_deoxy_NEU_L_medUp_max < ttp_deoxy_NEU_L_medUp_min
       ttp_deoxy_NEU_L_medUp = ttp_deoxy_NEU_L_medUp_max;
       peak_deoxy_NEU_L_medUp = max_deoxy_NEU_L_medUp ; 
   elseif ttp_deoxy_NEU_L_medUp_max > ttp_deoxy_NEU_L_medUp_min
       ttp_deoxy_NEU_L_medUp = ttp_deoxy_NEU_L_medUp_min;
       peak_deoxy_NEU_L_medUp = min_deoxy_NEU_L_medUp ;
   end
end

ttp_row = find(t==ttp_deoxy_NEU_L_medUp);
if isnan(peak_deoxy_NEU_L_medUp)
    ttl_deoxy_NEU_L_medUp = NaN ;
    low_deoxy_NEU_L_medUp = NaN ;
elseif peak_deoxy_NEU_L_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,8,3) < dcAvg(t0,2,8,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,8,3),[ttp_deoxy_NEU_L_medUp 25],[dcAvg(t0,2,8,3) dcAvg(t0,2,8,3)]) ;
        ttl_deoxy_NEU_L_medUp = min(xlin);
        low_deoxy_NEU_L_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,8,3) < dcAvg(t0,2,8,3)))
        ttl_deoxy_NEU_L_medUp = t( find(dcAvg(ttp_row:tslut_NEU,2,8,3)==min(dcAvg(ttp_row:tslut_NEU,2,8,3)))+ttp_row-1);
        low_deoxy_NEU_L_medUp = min(dcAvg(ttp_row:tslut_NEU,2,8,3))-dcAvg(t0,2,8,3);
    end
elseif peak_deoxy_NEU_L_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,8,3) > dcAvg(t0,2,8,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,8,3),[ttp_deoxy_NEU_L_medUp 25],[dcAvg(t0,2,8,3) dcAvg(t0,2,8,3)]) ;
        ttl_deoxy_NEU_L_medUp = min(xlin);
        low_deoxy_NEU_L_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,8,3) > dcAvg(t0,2,8,3)))
        ttl_deoxy_NEU_L_medUp = t( find(dcAvg(ttp_row:tslut_NEU,2,8,3)==max(dcAvg(ttp_row:tslut_NEU,2,8,3)))+ttp_row-1);
        low_deoxy_NEU_L_medUp = max(dcAvg(ttp_row:tslut_NEU,2,8,3))-dcAvg(t0,2,8,3);
    end
end

if abs(max_deoxy_NEU_L_medDown) > abs(min_deoxy_NEU_L_medDown)
   peak_deoxy_NEU_L_medDown = max_deoxy_NEU_L_medDown ;
   ttp_deoxy_NEU_L_medDown = t( dcAvg (:,2,9,3) == max(dcAvg(t0:tslut_NEU,2,9,3))) ;
elseif abs(max_deoxy_NEU_L_medDown) < abs(min_deoxy_NEU_L_medDown)
   peak_deoxy_NEU_L_medDown = min_deoxy_NEU_L_medDown ;
   ttp_deoxy_NEU_L_medDown = t( dcAvg (:,2,9,3) == min(dcAvg(t0:tslut_NEU,2,9,3))) ;
elseif isnan(max_deoxy_NEU_L_medDown)
   ttp_deoxy_NEU_L_medDown = NaN;
   peak_deoxy_NEU_L_medDown = NaN ;
elseif abs(max_deoxy_NEU_L_medDown) == abs(min_deoxy_NEU_L_medDown)
   ttp_deoxy_NEU_L_medDown_max = t( dcAvg (:,2,9,3) == max(dcAvg(t0:tslut_NEU,2,9,3))) ;
   ttp_deoxy_NEU_L_medDown_min = t( dcAvg (:,2,9,3) == min(dcAvg(t0:tslut_NEU,2,9,3))) ;
   if ttp_deoxy_NEU_L_medDown_max < ttp_deoxy_NEU_L_medDown_min
       ttp_deoxy_NEU_L_medDown = ttp_deoxy_NEU_L_medDown_max;
       peak_deoxy_NEU_L_medDown = max_deoxy_NEU_L_medDown ; 
   elseif ttp_deoxy_NEU_L_medDown_max > ttp_deoxy_NEU_L_medDown_min
       ttp_deoxy_NEU_L_medDown = ttp_deoxy_NEU_L_medDown_min;
       peak_deoxy_NEU_L_medDown = min_deoxy_NEU_L_medDown ;
   end
end

ttp_row = find(t==ttp_deoxy_NEU_L_medDown);
if isnan(peak_deoxy_NEU_L_medDown)
    ttl_deoxy_NEU_L_medDown = NaN ;
    low_deoxy_NEU_L_medDown = NaN ;
elseif peak_deoxy_NEU_L_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,9,3) < dcAvg(t0,2,9,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,9,3),[ttp_deoxy_NEU_L_medDown 25],[dcAvg(t0,2,9,3) dcAvg(t0,2,9,3)]) ;
        ttl_deoxy_NEU_L_medDown = min(xlin);
        low_deoxy_NEU_L_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,9,3) < dcAvg(t0,2,9,3)))
        ttl_deoxy_NEU_L_medDown = t( find(dcAvg(ttp_row:tslut_NEU,2,9,3)==min(dcAvg(ttp_row:tslut_NEU,2,9,3)))+ttp_row-1);
        low_deoxy_NEU_L_medDown = min(dcAvg(ttp_row:tslut_NEU,2,9,3))-dcAvg(t0,2,9,3);
    end
elseif peak_deoxy_NEU_L_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,9,3) > dcAvg(t0,2,9,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,9,3),[ttp_deoxy_NEU_L_medDown 25],[dcAvg(t0,2,9,3) dcAvg(t0,2,9,3)]) ;
        ttl_deoxy_NEU_L_medDown = min(xlin);
        low_deoxy_NEU_L_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,9,3) > dcAvg(t0,2,9,3)))
        ttl_deoxy_NEU_L_medDown = t( find(dcAvg(ttp_row:tslut_NEU,2,9,3)==max(dcAvg(ttp_row:tslut_NEU,2,9,3)))+ttp_row-1);
        low_deoxy_NEU_L_medDown = max(dcAvg(ttp_row:tslut_NEU,2,9,3))-dcAvg(t0,2,9,3);
    end
end

if abs(max_deoxy_NEU_L_latUp) > abs(min_deoxy_NEU_L_latUp)
   peak_deoxy_NEU_L_latUp = max_deoxy_NEU_L_latUp ;
   ttp_deoxy_NEU_L_latUp = t( dcAvg (:,2,10,3) == max(dcAvg(t0:tslut_NEU,2,10,3))) ;
elseif abs(max_deoxy_NEU_L_latUp) < abs(min_deoxy_NEU_L_latUp)
   peak_deoxy_NEU_L_latUp = min_deoxy_NEU_L_latUp ;
   ttp_deoxy_NEU_L_latUp = t( dcAvg (:,2,10,3) == min(dcAvg(t0:tslut_NEU,2,10,3))) ;
elseif isnan(max_deoxy_NEU_L_latUp)
   ttp_deoxy_NEU_L_latUp = NaN;
   peak_deoxy_NEU_L_latUp = NaN ;
elseif abs(max_deoxy_NEU_L_latUp) == abs(min_deoxy_NEU_L_latUp)
   ttp_deoxy_NEU_L_latUp_max = t( dcAvg (:,2,10,3) == max(dcAvg(t0:tslut_NEU,2,10,3))) ;
   ttp_deoxy_NEU_L_latUp_min = t( dcAvg (:,2,10,3) == min(dcAvg(t0:tslut_NEU,2,10,3))) ;
   if ttp_deoxy_NEU_L_latUp_max < ttp_deoxy_NEU_L_latUp_min
       ttp_deoxy_NEU_L_latUp = ttp_deoxy_NEU_L_latUp_max;
       peak_deoxy_NEU_L_latUp = max_deoxy_NEU_L_latUp ; 
   elseif ttp_deoxy_NEU_L_latUp_max > ttp_deoxy_NEU_L_latUp_min
       ttp_deoxy_NEU_L_latUp = ttp_deoxy_NEU_L_latUp_min;
       peak_deoxy_NEU_L_latUp = min_deoxy_NEU_L_latUp ;
   end
end
  
ttp_row = find(t==ttp_deoxy_NEU_L_latUp);
if isnan(peak_deoxy_NEU_L_latUp)
    ttl_deoxy_NEU_L_latUp = NaN ;
    low_deoxy_NEU_L_latUp = NaN ;
elseif peak_deoxy_NEU_L_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,10,3) < dcAvg(t0,2,10,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,10,3),[ttp_deoxy_NEU_L_latUp 25],[dcAvg(t0,2,10,3) dcAvg(t0,2,10,3)]) ;
        ttl_deoxy_NEU_L_latUp = min(xlin);
        low_deoxy_NEU_L_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,10,3) < dcAvg(t0,2,10,3)))
        ttl_deoxy_NEU_L_latUp = t( find(dcAvg(ttp_row:tslut_NEU,2,10,3)==min(dcAvg(ttp_row:tslut_NEU,2,10,3)))+ttp_row-1);
        low_deoxy_NEU_L_latUp = min(dcAvg(ttp_row:tslut_NEU,2,10,3))-dcAvg(t0,2,10,3);
    end
elseif peak_deoxy_NEU_L_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,10,3) > dcAvg(t0,2,10,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,10,3),[ttp_deoxy_NEU_L_latUp 25],[dcAvg(t0,2,10,3) dcAvg(t0,2,10,3)]) ;
        ttl_deoxy_NEU_L_latUp = min(xlin);
        low_deoxy_NEU_L_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,10,3) > dcAvg(t0,2,10,3)))
        ttl_deoxy_NEU_L_latUp = t( find(dcAvg(ttp_row:tslut_NEU,2,10,3)==max(dcAvg(ttp_row:tslut_NEU,2,10,3)))+ttp_row-1);
        low_deoxy_NEU_L_latUp = max(dcAvg(ttp_row:tslut_NEU,2,10,3))-dcAvg(t0,2,10,3);
    end
end

if abs(max_deoxy_NEU_L_latDown) > abs(min_deoxy_NEU_L_latDown)
   peak_deoxy_NEU_L_latDown = max_deoxy_NEU_L_latDown ;
   ttp_deoxy_NEU_L_latDown = t( dcAvg (:,2,11,3) == max(dcAvg(t0:tslut_NEU,2,11,3))) ;
elseif abs(max_deoxy_NEU_L_latDown) < abs(min_deoxy_NEU_L_latDown)
   peak_deoxy_NEU_L_latDown = min_deoxy_NEU_L_latDown ;
   ttp_deoxy_NEU_L_latDown = t( dcAvg (:,2,11,3) == min(dcAvg(t0:tslut_NEU,2,11,3))) ;
elseif isnan(max_deoxy_NEU_L_latDown)
   ttp_deoxy_NEU_L_latDown = NaN;
   peak_deoxy_NEU_L_latDown = NaN ;
elseif abs(max_deoxy_NEU_L_latDown) == abs(min_deoxy_NEU_L_latDown)
   ttp_deoxy_NEU_L_latDown_max = t( dcAvg (:,2,11,3) == max(dcAvg(t0:tslut_NEU,2,11,3))) ;
   ttp_deoxy_NEU_L_latDown_min = t( dcAvg (:,2,11,3) == min(dcAvg(t0:tslut_NEU,2,11,3))) ;
   if ttp_deoxy_NEU_L_latDown_max < ttp_deoxy_NEU_L_latDown_min
       ttp_deoxy_NEU_L_latDown = ttp_deoxy_NEU_L_latDown_max;
       peak_deoxy_NEU_L_latDown = max_deoxy_NEU_L_latDown ; 
   elseif ttp_deoxy_NEU_L_latDown_max > ttp_deoxy_NEU_L_latDown_min
       ttp_deoxy_NEU_L_latDown = ttp_deoxy_NEU_L_latDown_min;
       peak_deoxy_NEU_L_latDown = min_deoxy_NEU_L_latDown ;
   end
end

ttp_row = find(t==ttp_deoxy_NEU_L_latDown);
if isnan(peak_deoxy_NEU_L_latDown)
    ttl_deoxy_NEU_L_latDown = NaN ;
    low_deoxy_NEU_L_latDown = NaN ;
elseif peak_deoxy_NEU_L_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,11,3) < dcAvg(t0,2,11,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,11,3),[ttp_deoxy_NEU_L_latDown 25],[dcAvg(t0,2,11,3) dcAvg(t0,2,11,3)]) ;
        ttl_deoxy_NEU_L_latDown = min(xlin);
        low_deoxy_NEU_L_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,11,3) < dcAvg(t0,2,11,3)))
        ttl_deoxy_NEU_L_latDown = t( find(dcAvg(ttp_row:tslut_NEU,2,11,3)==min(dcAvg(ttp_row:tslut_NEU,2,11,3)))+ttp_row-1);
        low_deoxy_NEU_L_latDown = min(dcAvg(ttp_row:tslut_NEU,2,11,3))-dcAvg(t0,2,11,3);
    end
elseif peak_deoxy_NEU_L_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_NEU,2,11,3) > dcAvg(t0,2,11,3)))
        [xlin,~] = intersections(t(ttp_row:tslut_NEU),dcAvg(ttp_row:tslut_NEU,2,11,3),[ttp_deoxy_NEU_L_latDown 25],[dcAvg(t0,2,11,3) dcAvg(t0,2,11,3)]) ;
        ttl_deoxy_NEU_L_latDown = min(xlin);
        low_deoxy_NEU_L_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_NEU,2,11,3) > dcAvg(t0,2,11,3)))
        ttl_deoxy_NEU_L_latDown = t( find(dcAvg(ttp_row:tslut_NEU,2,11,3)==max(dcAvg(ttp_row:tslut_NEU,2,11,3)))+ttp_row-1);
        low_deoxy_NEU_L_latDown = max(dcAvg(ttp_row:tslut_NEU,2,11,3))-dcAvg(t0,2,11,3);
    end
end

%% CONGRUENT OXY
% Finding peak, time-to-peak (TTP), return to baseline (low) and
% time-to-low (TTL)

% Max. oxy during CONGRUENT
max_oxy_CON_R_latUp = max(dcAvg(t0:tslut_INC,1,2,1))-dcAvg(t0,1,2,1); 
max_oxy_CON_R_latDown = max(dcAvg(t0:tslut_INC,1,3,1))-dcAvg(t0,1,3,1);
max_oxy_CON_R_medUp = max(dcAvg(t0:tslut_INC,1,4,1))-dcAvg(t0,1,4,1);
max_oxy_CON_R_medDown = max(dcAvg(t0:tslut_INC,1,5,1))-dcAvg(t0,1,5,1); 
max_oxy_CON_L_medUp = max(dcAvg(t0:tslut_INC,1,8,1))-dcAvg(t0,1,8,1);
max_oxy_CON_L_medDown = max(dcAvg(t0:tslut_INC,1,9,1))-dcAvg(t0,1,9,1);
max_oxy_CON_L_latUp = max(dcAvg(t0:tslut_INC,1,10,1))-dcAvg(t0,1,10,1);
max_oxy_CON_L_latDown = max(dcAvg(t0:tslut_INC,1,11,1))-dcAvg(t0,1,11,1);

% Min. oxy during CONGRUENT
min_oxy_CON_R_latUp = min(dcAvg(t0:tslut_CON,1,2,1))-dcAvg(t0,1,2,1);
min_oxy_CON_R_latDown = min(dcAvg(t0:tslut_CON,1,3,1))-dcAvg(t0,1,3,1);
min_oxy_CON_R_medUp = min(dcAvg(t0:tslut_CON,1,4,1))-dcAvg(t0,1,4,1);
min_oxy_CON_R_medDown = min(dcAvg(t0:tslut_CON,1,5,1))-dcAvg(t0,1,5,1);
min_oxy_CON_L_medUp = min(dcAvg(t0:tslut_CON,1,8,1))-dcAvg(t0,1,8,1);
min_oxy_CON_L_medDown = min(dcAvg(t0:tslut_CON,1,9,1))-dcAvg(t0,1,9,1);
min_oxy_CON_L_latUp = min(dcAvg(t0:tslut_CON,1,10,1))-dcAvg(t0,1,10,1);
min_oxy_CON_L_latDown = min(dcAvg(t0:tslut_CON,1,11,1))-dcAvg(t0,1,11,1);

% Finding peak (min. or max.) and time-to-peak during CONGRUENT 
% If identical, then use the first one
if abs(max_oxy_CON_R_latUp) > abs(min_oxy_CON_R_latUp)
   peak_oxy_CON_R_latUp = max_oxy_CON_R_latUp ;
   ttp_oxy_CON_R_latUp = t( dcAvg (:,1,2,1) == max(dcAvg(t0:tslut_INC,1,2,1))) ;
elseif abs(max_oxy_CON_R_latUp) < abs(min_oxy_CON_R_latUp)
   peak_oxy_CON_R_latUp = min_oxy_CON_R_latUp ;
   ttp_oxy_CON_R_latUp = t( dcAvg (:,1,2,1) == min(dcAvg(t0:tslut_CON,1,2,1))) ;
elseif isnan(max_oxy_CON_R_latUp)
   ttp_oxy_CON_R_latUp = NaN;
   peak_oxy_CON_R_latUp = NaN ;
elseif abs(max_oxy_CON_R_latUp) == abs(min_oxy_CON_R_latUp)
   ttp_oxy_CON_R_latUp_max = t( dcAvg (:,1,2,1) == max(dcAvg(t0:tslut_INC,1,2,1))) ;
   ttp_oxy_CON_R_latUp_min = t( dcAvg (:,1,2,1) == min(dcAvg(t0:tslut_CON,1,2,1))) ;
   if ttp_oxy_CON_R_latUp_max < ttp_oxy_CON_R_latUp_min
       ttp_oxy_CON_R_latUp = ttp_oxy_CON_R_latUp_max;
       peak_oxy_CON_R_latUp = max_oxy_CON_R_latUp ; 
   elseif ttp_oxy_CON_R_latUp_max > ttp_oxy_CON_R_latUp_min
       ttp_oxy_CON_R_latUp = ttp_oxy_CON_R_latUp_min;
       peak_oxy_CON_R_latUp = min_oxy_CON_R_latUp ;
   end
end

ttp_row = find(t==ttp_oxy_CON_R_latUp);
if isnan(peak_oxy_CON_R_latUp)
    ttl_oxy_CON_R_latUp = NaN ;
    low_oxy_CON_R_latUp = NaN ;
elseif peak_oxy_CON_R_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,2,1) < dcAvg(t0,1,2,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,2,1),[ttp_oxy_CON_R_latUp 25],[dcAvg(t0,1,2,1) dcAvg(t0,1,2,1)]) ;
        ttl_oxy_CON_R_latUp = min(xlin);
        low_oxy_CON_R_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,2,1) < dcAvg(t0,1,2,1)))
        ttl_oxy_CON_R_latUp = t( find(dcAvg(ttp_row:tslut_CON,1,2,1)==min(dcAvg(ttp_row:tslut_CON,1,2,1)))+ttp_row-1);
        low_oxy_CON_R_latUp = min(dcAvg(ttp_row:tslut_CON,1,2,1))-dcAvg(t0,1,2,1);
    end
elseif peak_oxy_CON_R_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,2,1) > dcAvg(t0,1,2,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,2,1),[ttp_oxy_CON_R_latUp 25],[dcAvg(t0,1,2,1) dcAvg(t0,1,2,1)]) ;
        ttl_oxy_CON_R_latUp = min(xlin);
        low_oxy_CON_R_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,2,1) > dcAvg(t0,1,2,1)))
        ttl_oxy_CON_R_latUp = t( find(dcAvg(ttp_row:tslut_CON,1,2,1)==max(dcAvg(ttp_row:tslut_CON,1,2,1)))+ttp_row-1);
        low_oxy_CON_R_latUp = max(dcAvg(ttp_row:tslut_CON,1,2,1))-dcAvg(t0,1,2,1);
    end
end

if abs(max_oxy_CON_R_latDown) > abs(min_oxy_CON_R_latDown)
   peak_oxy_CON_R_latDown = max_oxy_CON_R_latDown ;
   ttp_oxy_CON_R_latDown = t( dcAvg (:,1,3,1) == max(dcAvg(t0:tslut_INC,1,3,1))) ;
elseif abs(max_oxy_CON_R_latDown) < abs(min_oxy_CON_R_latDown)
   peak_oxy_CON_R_latDown = min_oxy_CON_R_latDown ;
   ttp_oxy_CON_R_latDown = t( dcAvg (:,1,3,1) == min(dcAvg(t0:tslut_CON,1,3,1))) ;
elseif isnan(max_oxy_CON_R_latDown)
   peak_oxy_CON_R_latDown = NaN;
   ttp_oxy_CON_R_latDown = NaN;
elseif abs(max_oxy_CON_R_latDown) == abs(min_oxy_CON_R_latDown)
   ttp_oxy_CON_R_latDown_max = t( dcAvg (:,1,3,1) == max(dcAvg(t0:tslut_INC,1,3,1))) ;
   ttp_oxy_CON_R_latDown_min = t( dcAvg (:,1,3,1) == min(dcAvg(t0:tslut_CON,1,3,1))) ;
    if ttp_oxy_CON_R_latDown_max < ttp_oxy_CON_R_latDown_min
       ttp_oxy_CON_R_latDown = ttp_oxy_CON_R_latDown_max;
       peak_oxy_CON_R_latDown = max_oxy_CON_R_latDown ; 
    elseif ttp_oxy_CON_R_latDown_max > ttp_oxy_CON_R_latDown_min
       ttp_oxy_CON_R_latDown = ttp_oxy_CON_R_latDown_min;
       peak_oxy_CON_R_latDown = min_oxy_CON_R_latDown ;
    end
end

ttp_row = find(t==ttp_oxy_CON_R_latDown);
if isnan(peak_oxy_CON_R_latDown)
    ttl_oxy_CON_R_latDown = NaN ;
    low_oxy_CON_R_latDown = NaN ;
elseif peak_oxy_CON_R_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,3,1) < dcAvg(t0,1,3,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,3,1),[ttp_oxy_CON_R_latDown 25],[dcAvg(t0,1,3,1) dcAvg(t0,1,3,1)]) ;
        ttl_oxy_CON_R_latDown = min(xlin);
        low_oxy_CON_R_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,3,1) < dcAvg(t0,1,3,1)))
        ttl_oxy_CON_R_latDown = t( find(dcAvg(ttp_row:tslut_CON,1,3,1)==min(dcAvg(ttp_row:tslut_CON,1,3,1)))+ttp_row-1);
        low_oxy_CON_R_latDown = min(dcAvg(ttp_row:tslut_CON,1,3,1))-dcAvg(t0,1,3,1);
    end
elseif peak_oxy_CON_R_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,3,1) > dcAvg(t0,1,3,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,3,1),[ttp_oxy_CON_R_latDown 25],[dcAvg(t0,1,3,1) dcAvg(t0,1,3,1)]) ;
        ttl_oxy_CON_R_latDown = min(xlin);
        low_oxy_CON_R_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,3,1) > dcAvg(t0,1,3,1)))
        ttl_oxy_CON_R_latDown = t( find(dcAvg(ttp_row:tslut_CON,1,3,1)==max(dcAvg(ttp_row:tslut_CON,1,3,1)))+ttp_row-1);
        low_oxy_CON_R_latDown = max(dcAvg(ttp_row:tslut_CON,1,3,1))-dcAvg(t0,1,3,1);
    end
end

if abs(max_oxy_CON_R_medUp) > abs(min_oxy_CON_R_medUp)
   peak_oxy_CON_R_medUp = max_oxy_CON_R_medUp ;
   ttp_oxy_CON_R_medUp = t( dcAvg (:,1,4,1) == max(dcAvg(t0:tslut_INC,1,4,1))) ;
elseif abs(max_oxy_CON_R_medUp) < abs(min_oxy_CON_R_medUp)
   peak_oxy_CON_R_medUp = min_oxy_CON_R_medUp ;
   ttp_oxy_CON_R_medUp = t( dcAvg (:,1,4,1) == min(dcAvg(t0:tslut_CON,1,4,1))) ;
elseif isnan(max_oxy_CON_R_medUp)
   ttp_oxy_CON_R_medUp = NaN;
   peak_oxy_CON_R_medUp = NaN ;
elseif abs(max_oxy_CON_R_medUp) == abs(min_oxy_CON_R_medUp)
   ttp_oxy_CON_R_medUp_max = t( dcAvg (:,1,4,1) == max(dcAvg(t0:tslut_INC,1,4,1))) ;
   ttp_oxy_CON_R_medUp_min = t( dcAvg (:,1,4,1) == min(dcAvg(t0:tslut_CON,1,4,1))) ;
   if ttp_oxy_CON_R_medUp_max < ttp_oxy_CON_R_medUp_min
       ttp_oxy_CON_R_medUp = ttp_oxy_CON_R_medUp_max;
       peak_oxy_CON_R_medUp = max_oxy_CON_R_medUp ; 
   elseif ttp_oxy_CON_R_medUp_max > ttp_oxy_CON_R_medUp_min
       ttp_oxy_CON_R_medUp = ttp_oxy_CON_R_medUp_min;
       peak_oxy_CON_R_medUp = min_oxy_CON_R_medUp ;
   end
end

ttp_row = find(t==ttp_oxy_CON_R_medUp);
if isnan(peak_oxy_CON_R_medUp)
    ttl_oxy_CON_R_medUp = NaN ;
    low_oxy_CON_R_medUp = NaN ;
elseif peak_oxy_CON_R_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,4,1) < dcAvg(t0,1,4,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,4,1),[ttp_oxy_CON_R_medUp 25],[dcAvg(t0,1,4,1) dcAvg(t0,1,4,1)]) ;
        ttl_oxy_CON_R_medUp = min(xlin);
        low_oxy_CON_R_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,4,1) < dcAvg(t0,1,4,1)))
        ttl_oxy_CON_R_medUp = t( find(dcAvg(ttp_row:tslut_CON,1,4,1)==min(dcAvg(ttp_row:tslut_CON,1,4,1)))+ttp_row-1);
        low_oxy_CON_R_medUp = min(dcAvg(ttp_row:tslut_CON,1,4,1))-dcAvg(t0,1,4,1);
    end
elseif peak_oxy_CON_R_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,4,1) > dcAvg(t0,1,4,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,4,1),[ttp_oxy_CON_R_medUp 25],[dcAvg(t0,1,4,1) dcAvg(t0,1,4,1)]) ;
        ttl_oxy_CON_R_medUp = min(xlin);
        low_oxy_CON_R_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,4,1) > dcAvg(t0,1,4,1)))
        ttl_oxy_CON_R_medUp = t( find(dcAvg(ttp_row:tslut_CON,1,4,1)==max(dcAvg(ttp_row:tslut_CON,1,4,1)))+ttp_row-1);
        low_oxy_CON_R_medUp = max(dcAvg(ttp_row:tslut_CON,1,4,1))-dcAvg(t0,1,4,1);
    end
end

if abs(max_oxy_CON_R_medDown) > abs(min_oxy_CON_R_medDown)
   peak_oxy_CON_R_medDown = max_oxy_CON_R_medDown ;
   ttp_oxy_CON_R_medDown = t( dcAvg (:,1,5,1) == max(dcAvg(t0:tslut_INC,1,5,1))) ;
elseif abs(max_oxy_CON_R_medDown) < abs(min_oxy_CON_R_medDown)
   peak_oxy_CON_R_medDown = min_oxy_CON_R_medDown ;
   ttp_oxy_CON_R_medDown = t( dcAvg (:,1,5,1) == min(dcAvg(t0:tslut_CON,1,5,1))) ;
elseif isnan(max_oxy_CON_R_medDown)
   ttp_oxy_CON_R_medDown = NaN;
   peak_oxy_CON_R_medDown = NaN ;
elseif abs(max_oxy_CON_R_medDown) == abs(min_oxy_CON_R_medDown)
   ttp_oxy_CON_R_medDown_max = t( dcAvg (:,1,5,1) == max(dcAvg(t0:tslut_INC,1,5,1))) ;
   ttp_oxy_CON_R_medDown_min = t( dcAvg (:,1,5,1) == min(dcAvg(t0:tslut_CON,1,5,1))) ;
   if ttp_oxy_CON_R_medDown_max < ttp_oxy_CON_R_medDown_min
       ttp_oxy_CON_R_medDown = ttp_oxy_CON_R_medDown_max;
       peak_oxy_CON_R_medDown = max_oxy_CON_R_medDown ; 
   elseif ttp_oxy_CON_R_medDown_max > ttp_oxy_CON_R_medDown_min
       ttp_oxy_CON_R_medDown = ttp_oxy_CON_R_medDown_min;
       peak_oxy_CON_R_medDown = min_oxy_CON_R_medDown ;
   end
end

ttp_row = find(t==ttp_oxy_CON_R_medDown);
if isnan(peak_oxy_CON_R_medDown)
    ttl_oxy_CON_R_medDown = NaN ;
    low_oxy_CON_R_medDown = NaN ;
elseif peak_oxy_CON_R_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,5,1) < dcAvg(t0,1,5,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,5,1),[ttp_oxy_CON_R_medDown 25],[dcAvg(t0,1,5,1) dcAvg(t0,1,5,1)]) ;
        ttl_oxy_CON_R_medDown = min(xlin);
        low_oxy_CON_R_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,5,1) < dcAvg(t0,1,5,1)))
        ttl_oxy_CON_R_medDown = t( find(dcAvg(ttp_row:tslut_CON,1,5,1)==min(dcAvg(ttp_row:tslut_CON,1,5,1)))+ttp_row-1);
        low_oxy_CON_R_medDown = min(dcAvg(ttp_row:tslut_CON,1,5,1))-dcAvg(t0,1,5,1);
    end
elseif peak_oxy_CON_R_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,5,1) > dcAvg(t0,1,5,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,5,1),[ttp_oxy_CON_R_medDown 25],[dcAvg(t0,1,5,1) dcAvg(t0,1,5,1)]) ;
        ttl_oxy_CON_R_medDown = min(xlin);
        low_oxy_CON_R_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,5,1) > dcAvg(t0,1,5,1)))
        ttl_oxy_CON_R_medDown = t( find(dcAvg(ttp_row:tslut_CON,1,5,1)==max(dcAvg(ttp_row:tslut_CON,1,5,1)))+ttp_row-1);
        low_oxy_CON_R_medDown = max(dcAvg(ttp_row:tslut_CON,1,5,1))-dcAvg(t0,1,5,1);
    end
end

if abs(max_oxy_CON_L_medUp) > abs(min_oxy_CON_L_medUp)
   peak_oxy_CON_L_medUp = max_oxy_CON_L_medUp ;
   ttp_oxy_CON_L_medUp = t( dcAvg (:,1,8,1) == max(dcAvg(t0:tslut_INC,1,8,1))) ;
elseif abs(max_oxy_CON_L_medUp) < abs(min_oxy_CON_L_medUp)
   peak_oxy_CON_L_medUp = min_oxy_CON_L_medUp ;
   ttp_oxy_CON_L_medUp = t( dcAvg (:,1,8,1) == min(dcAvg(t0:tslut_CON,1,8,1))) ;
elseif isnan(max_oxy_CON_L_medUp)
   ttp_oxy_CON_L_medUp = NaN;
   peak_oxy_CON_L_medUp = NaN ;
elseif abs(max_oxy_CON_L_medUp) == abs(min_oxy_CON_L_medUp)
   ttp_oxy_CON_L_medUp_max = t( dcAvg (:,1,8,1) == max(dcAvg(t0:tslut_INC,1,8,1))) ;
   ttp_oxy_CON_L_medUp_min = t( dcAvg (:,1,8,1) == min(dcAvg(t0:tslut_CON,1,8,1))) ;
   if ttp_oxy_CON_L_medUp_max < ttp_oxy_CON_L_medUp_min
       ttp_oxy_CON_L_medUp = ttp_oxy_CON_L_medUp_max;
       peak_oxy_CON_L_medUp = max_oxy_CON_L_medUp ; 
   elseif ttp_oxy_CON_L_medUp_max > ttp_oxy_CON_L_medUp_min
       ttp_oxy_CON_L_medUp = ttp_oxy_CON_L_medUp_min;
       peak_oxy_CON_L_medUp = min_oxy_CON_L_medUp ;
   end
end

ttp_row = find(t==ttp_oxy_CON_L_medUp);
if isnan(peak_oxy_CON_L_medUp)
    ttl_oxy_CON_L_medUp = NaN ;
    low_oxy_CON_L_medUp = NaN ;
elseif peak_oxy_CON_L_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,8,1) < dcAvg(t0,1,8,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,8,1),[ttp_oxy_CON_L_medUp 25],[dcAvg(t0,1,8,1) dcAvg(t0,1,8,1)]) ;
        ttl_oxy_CON_L_medUp = min(xlin);
        low_oxy_CON_L_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,8,1) < dcAvg(t0,1,8,1)))
        ttl_oxy_CON_L_medUp = t( find(dcAvg(ttp_row:tslut_CON,1,8,1)==min(dcAvg(ttp_row:tslut_CON,1,8,1)))+ttp_row-1);
        low_oxy_CON_L_medUp = min(dcAvg(ttp_row:tslut_CON,1,8,1))-dcAvg(t0,1,8,1);
    end
elseif peak_oxy_CON_L_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,8,1) > dcAvg(t0,1,8,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,8,1),[ttp_oxy_CON_L_medUp 25],[dcAvg(t0,1,8,1) dcAvg(t0,1,8,1)]) ;
        ttl_oxy_CON_L_medUp = min(xlin);
        low_oxy_CON_L_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,8,1) > dcAvg(t0,1,8,1)))
        ttl_oxy_CON_L_medUp = t( find(dcAvg(ttp_row:tslut_CON,1,8,1)==max(dcAvg(ttp_row:tslut_CON,1,8,1)))+ttp_row-1);
        low_oxy_CON_L_medUp = max(dcAvg(ttp_row:tslut_CON,1,8,1))-dcAvg(t0,1,8,1);
    end
end

if abs(max_oxy_CON_L_medDown) > abs(min_oxy_CON_L_medDown)
   peak_oxy_CON_L_medDown = max_oxy_CON_L_medDown ;
   ttp_oxy_CON_L_medDown = t( dcAvg (:,1,9,1) == max(dcAvg(t0:tslut_INC,1,9,1))) ;
elseif abs(max_oxy_CON_L_medDown) < abs(min_oxy_CON_L_medDown)
   peak_oxy_CON_L_medDown = min_oxy_CON_L_medDown ;
   ttp_oxy_CON_L_medDown = t( dcAvg (:,1,9,1) == min(dcAvg(t0:tslut_CON,1,9,1))) ;
elseif isnan(max_oxy_CON_L_medDown)
   ttp_oxy_CON_L_medDown = NaN;
   peak_oxy_CON_L_medDown = NaN ;
elseif abs(max_oxy_CON_L_medDown) == abs(min_oxy_CON_L_medDown)
   ttp_oxy_CON_L_medDown_max = t( dcAvg (:,1,9,1) == max(dcAvg(t0:tslut_INC,1,9,1))) ;
   ttp_oxy_CON_L_medDown_min = t( dcAvg (:,1,9,1) == min(dcAvg(t0:tslut_CON,1,9,1))) ;
   if ttp_oxy_CON_L_medDown_max < ttp_oxy_CON_L_medDown_min
       ttp_oxy_CON_L_medDown = ttp_oxy_CON_L_medDown_max;
       peak_oxy_CON_L_medDown = max_oxy_CON_L_medDown ; 
   elseif ttp_oxy_CON_L_medDown_max > ttp_oxy_CON_L_medDown_min
       ttp_oxy_CON_L_medDown = ttp_oxy_CON_L_medDown_min;
       peak_oxy_CON_L_medDown = min_oxy_CON_L_medDown ;
   end
end

ttp_row = find(t==ttp_oxy_CON_L_medDown);
if isnan(peak_oxy_CON_L_medDown)
    ttl_oxy_CON_L_medDown = NaN ;
    low_oxy_CON_L_medDown = NaN ;
elseif peak_oxy_CON_L_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,9,1) < dcAvg(t0,1,9,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,9,1),[ttp_oxy_CON_L_medDown 25],[dcAvg(t0,1,9,1) dcAvg(t0,1,9,1)]) ;
        ttl_oxy_CON_L_medDown = min(xlin);
        low_oxy_CON_L_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,9,1) < dcAvg(t0,1,9,1)))
        ttl_oxy_CON_L_medDown = t( find(dcAvg(ttp_row:tslut_CON,1,9,1)==min(dcAvg(ttp_row:tslut_CON,1,9,1)))+ttp_row-1);
        low_oxy_CON_L_medDown = min(dcAvg(ttp_row:tslut_CON,1,9,1))-dcAvg(t0,1,9,1);
    end
elseif peak_oxy_CON_L_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,9,1) > dcAvg(t0,1,9,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,9,1),[ttp_oxy_CON_L_medDown 25],[dcAvg(t0,1,9,1) dcAvg(t0,1,9,1)]) ;
        ttl_oxy_CON_L_medDown = min(xlin);
        low_oxy_CON_L_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,9,1) > dcAvg(t0,1,9,1)))
        ttl_oxy_CON_L_medDown = t( find(dcAvg(ttp_row:tslut_CON,1,9,1)==max(dcAvg(ttp_row:tslut_CON,1,9,1)))+ttp_row-1);
        low_oxy_CON_L_medDown = max(dcAvg(ttp_row:tslut_CON,1,9,1))-dcAvg(t0,1,9,1);
    end
end

if abs(max_oxy_CON_L_latUp) > abs(min_oxy_CON_L_latUp)
   peak_oxy_CON_L_latUp = max_oxy_CON_L_latUp ;
   ttp_oxy_CON_L_latUp = t( dcAvg (:,1,10,1) == max(dcAvg(t0:tslut_INC,1,10,1))) ;
elseif abs(max_oxy_CON_L_latUp) < abs(min_oxy_CON_L_latUp)
   peak_oxy_CON_L_latUp = min_oxy_CON_L_latUp ;
   ttp_oxy_CON_L_latUp = t( dcAvg (:,1,10,1) == min(dcAvg(t0:tslut_CON,1,10,1))) ;
elseif isnan(max_oxy_CON_L_latUp)
   ttp_oxy_CON_L_latUp = NaN;
   peak_oxy_CON_L_latUp = NaN ;
elseif abs(max_oxy_CON_L_latUp) == abs(min_oxy_CON_L_latUp)
   ttp_oxy_CON_L_latUp_max = t( dcAvg (:,1,10,1) == max(dcAvg(t0:tslut_INC,1,10,1))) ;
   ttp_oxy_CON_L_latUp_min = t( dcAvg (:,1,10,1) == min(dcAvg(t0:tslut_CON,1,10,1))) ;
   if ttp_oxy_CON_L_latUp_max < ttp_oxy_CON_L_latUp_min
       ttp_oxy_CON_L_latUp = ttp_oxy_CON_L_latUp_max;
       peak_oxy_CON_L_latUp = max_oxy_CON_L_latUp ; 
   elseif ttp_oxy_CON_L_latUp_max > ttp_oxy_CON_L_latUp_min
       ttp_oxy_CON_L_latUp = ttp_oxy_CON_L_latUp_min;
       peak_oxy_CON_L_latUp = min_oxy_CON_L_latUp ;
   end
end
  
ttp_row = find(t==ttp_oxy_CON_L_latUp);
if isnan(peak_oxy_CON_L_latUp)
    ttl_oxy_CON_L_latUp = NaN ;
    low_oxy_CON_L_latUp = NaN ;
elseif peak_oxy_CON_L_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,10,1) < dcAvg(t0,1,10,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,10,1),[ttp_oxy_CON_L_latUp 25],[dcAvg(t0,1,10,1) dcAvg(t0,1,10,1)]) ;
        ttl_oxy_CON_L_latUp = min(xlin);
        low_oxy_CON_L_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,10,1) < dcAvg(t0,1,10,1)))
        ttl_oxy_CON_L_latUp = t( find(dcAvg(ttp_row:tslut_CON,1,10,1)==min(dcAvg(ttp_row:tslut_CON,1,10,1)))+ttp_row-1);
        low_oxy_CON_L_latUp = min(dcAvg(ttp_row:tslut_CON,1,10,1))-dcAvg(t0,1,10,1);
    end
elseif peak_oxy_CON_L_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,10,1) > dcAvg(t0,1,10,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,10,1),[ttp_oxy_CON_L_latUp 25],[dcAvg(t0,1,10,1) dcAvg(t0,1,10,1)]) ;
        ttl_oxy_CON_L_latUp = min(xlin);
        low_oxy_CON_L_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,10,1) > dcAvg(t0,1,10,1)))
        ttl_oxy_CON_L_latUp = t( find(dcAvg(ttp_row:tslut_CON,1,10,1)==max(dcAvg(ttp_row:tslut_CON,1,10,1)))+ttp_row-1);
        low_oxy_CON_L_latUp = max(dcAvg(ttp_row:tslut_CON,1,10,1))-dcAvg(t0,1,10,1);
    end
end

if abs(max_oxy_CON_L_latDown) > abs(min_oxy_CON_L_latDown)
   peak_oxy_CON_L_latDown = max_oxy_CON_L_latDown ;
   ttp_oxy_CON_L_latDown = t( dcAvg (:,1,11,1) == max(dcAvg(t0:tslut_INC,1,11,1))) ;
elseif abs(max_oxy_CON_L_latDown) < abs(min_oxy_CON_L_latDown)
   peak_oxy_CON_L_latDown = min_oxy_CON_L_latDown ;
   ttp_oxy_CON_L_latDown = t( dcAvg (:,1,11,1) == min(dcAvg(t0:tslut_CON,1,11,1))) ;
elseif isnan(max_oxy_CON_L_latDown)
   ttp_oxy_CON_L_latDown = NaN;
   peak_oxy_CON_L_latDown = NaN ;
elseif abs(max_oxy_CON_L_latDown) == abs(min_oxy_CON_L_latDown)
   ttp_oxy_CON_L_latDown_max = t( dcAvg (:,1,11,1) == max(dcAvg(t0:tslut_INC,1,11,1))) ;
   ttp_oxy_CON_L_latDown_min = t( dcAvg (:,1,11,1) == min(dcAvg(t0:tslut_CON,1,11,1))) ;
   if ttp_oxy_CON_L_latDown_max < ttp_oxy_CON_L_latDown_min
       ttp_oxy_CON_L_latDown = ttp_oxy_CON_L_latDown_max;
       peak_oxy_CON_L_latDown = max_oxy_CON_L_latDown ; 
   elseif ttp_oxy_CON_L_latDown_max > ttp_oxy_CON_L_latDown_min
       ttp_oxy_CON_L_latDown = ttp_oxy_CON_L_latDown_min;
       peak_oxy_CON_L_latDown = min_oxy_CON_L_latDown ;
   end
end

ttp_row = find(t==ttp_oxy_CON_L_latDown);
if isnan(peak_oxy_CON_L_latDown)
    ttl_oxy_CON_L_latDown = NaN ;
    low_oxy_CON_L_latDown = NaN ;
elseif peak_oxy_CON_L_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,11,1) < dcAvg(t0,1,11,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,11,1),[ttp_oxy_CON_L_latDown 25],[dcAvg(t0,1,11,1) dcAvg(t0,1,11,1)]) ;
        ttl_oxy_CON_L_latDown = min(xlin);
        low_oxy_CON_L_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,11,1) < dcAvg(t0,1,11,1)))
        ttl_oxy_CON_L_latDown = t( find(dcAvg(ttp_row:tslut_CON,1,11,1)==min(dcAvg(ttp_row:tslut_CON,1,11,1)))+ttp_row-1);
        low_oxy_CON_L_latDown = min(dcAvg(ttp_row:tslut_CON,1,11,1))-dcAvg(t0,1,11,1);
    end
elseif peak_oxy_CON_L_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,1,11,1) > dcAvg(t0,1,11,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,1,11,1),[ttp_oxy_CON_L_latDown 25],[dcAvg(t0,1,11,1) dcAvg(t0,1,11,1)]) ;
        ttl_oxy_CON_L_latDown = min(xlin);
        low_oxy_CON_L_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,1,11,1) > dcAvg(t0,1,11,1)))
        ttl_oxy_CON_L_latDown = t( find(dcAvg(ttp_row:tslut_CON,1,11,1)==max(dcAvg(ttp_row:tslut_CON,1,11,1)))+ttp_row-1);
        low_oxy_CON_L_latDown = max(dcAvg(ttp_row:tslut_CON,1,11,1))-dcAvg(t0,1,11,1);
    end
end

%% CONGRUENT DEOXY
% Finding peak, time-to-peak (TTP), return to baseline (low) and
% time-to-low (TTL)

% Max. deoxy during CONGRUENT
max_deoxy_CON_R_latUp = max(dcAvg(t0:tslut_INC,2,2,1))-dcAvg(t0,2,2,1); 
max_deoxy_CON_R_latDown = max(dcAvg(t0:tslut_INC,2,3,1))-dcAvg(t0,2,3,1);
max_deoxy_CON_R_medUp = max(dcAvg(t0:tslut_INC,2,4,1))-dcAvg(t0,2,4,1);
max_deoxy_CON_R_medDown = max(dcAvg(t0:tslut_INC,2,5,1))-dcAvg(t0,2,5,1); 
max_deoxy_CON_L_medUp = max(dcAvg(t0:tslut_INC,2,8,1))-dcAvg(t0,2,8,1);
max_deoxy_CON_L_medDown = max(dcAvg(t0:tslut_INC,2,9,1))-dcAvg(t0,2,9,1);
max_deoxy_CON_L_latUp = max(dcAvg(t0:tslut_INC,2,10,1))-dcAvg(t0,2,10,1);
max_deoxy_CON_L_latDown = max(dcAvg(t0:tslut_INC,2,11,1))-dcAvg(t0,2,11,1);

% Min. deoxy during CONGRUENT
min_deoxy_CON_R_latUp = min(dcAvg(t0:tslut_CON,2,2,1))-dcAvg(t0,2,2,1);
min_deoxy_CON_R_latDown = min(dcAvg(t0:tslut_CON,2,3,1))-dcAvg(t0,2,3,1);
min_deoxy_CON_R_medUp = min(dcAvg(t0:tslut_CON,2,4,1))-dcAvg(t0,2,4,1);
min_deoxy_CON_R_medDown = min(dcAvg(t0:tslut_CON,2,5,1))-dcAvg(t0,2,5,1);
min_deoxy_CON_L_medUp = min(dcAvg(t0:tslut_CON,2,8,1))-dcAvg(t0,2,8,1);
min_deoxy_CON_L_medDown = min(dcAvg(t0:tslut_CON,2,9,1))-dcAvg(t0,2,9,1);
min_deoxy_CON_L_latUp = min(dcAvg(t0:tslut_CON,2,10,1))-dcAvg(t0,2,10,1);
min_deoxy_CON_L_latDown = min(dcAvg(t0:tslut_CON,2,11,1))-dcAvg(t0,2,11,1);

% Finding peak (min. or max.) and time-to-peak during CONGRUENT 
% If identical, then use the first one
if abs(max_deoxy_CON_R_latUp) > abs(min_deoxy_CON_R_latUp)
   peak_deoxy_CON_R_latUp = max_deoxy_CON_R_latUp ;
   ttp_deoxy_CON_R_latUp = t( dcAvg (:,2,2,1) == max(dcAvg(t0:tslut_INC,2,2,1))) ;
elseif abs(max_deoxy_CON_R_latUp) < abs(min_deoxy_CON_R_latUp)
   peak_deoxy_CON_R_latUp = min_deoxy_CON_R_latUp ;
   ttp_deoxy_CON_R_latUp = t( dcAvg (:,2,2,1) == min(dcAvg(t0:tslut_CON,2,2,1))) ;
elseif isnan(max_deoxy_CON_R_latUp)
   ttp_deoxy_CON_R_latUp = NaN;
   peak_deoxy_CON_R_latUp = NaN ;
elseif abs(max_deoxy_CON_R_latUp) == abs(min_deoxy_CON_R_latUp)
   ttp_deoxy_CON_R_latUp_max = t( dcAvg (:,2,2,1) == max(dcAvg(t0:tslut_INC,2,2,1))) ;
   ttp_deoxy_CON_R_latUp_min = t( dcAvg (:,2,2,1) == min(dcAvg(t0:tslut_CON,2,2,1))) ;
   if ttp_deoxy_CON_R_latUp_max < ttp_deoxy_CON_R_latUp_min
       ttp_deoxy_CON_R_latUp = ttp_deoxy_CON_R_latUp_max;
       peak_deoxy_CON_R_latUp = max_deoxy_CON_R_latUp ; 
   elseif ttp_deoxy_CON_R_latUp_max > ttp_deoxy_CON_R_latUp_min
       ttp_deoxy_CON_R_latUp = ttp_deoxy_CON_R_latUp_min;
       peak_deoxy_CON_R_latUp = min_deoxy_CON_R_latUp ;
   end
end

ttp_row = find(t==ttp_deoxy_CON_R_latUp);
if isnan(peak_deoxy_CON_R_latUp)
    ttl_deoxy_CON_R_latUp = NaN ;
    low_deoxy_CON_R_latUp = NaN ;
elseif peak_deoxy_CON_R_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,2,1) < dcAvg(t0,2,2,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,2,1),[ttp_deoxy_CON_R_latUp 25],[dcAvg(t0,2,2,1) dcAvg(t0,2,2,1)]) ;
        ttl_deoxy_CON_R_latUp = min(xlin);
        low_deoxy_CON_R_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,2,1) < dcAvg(t0,2,2,1)))
        ttl_deoxy_CON_R_latUp = t( find(dcAvg(ttp_row:tslut_CON,2,2,1)==min(dcAvg(ttp_row:tslut_CON,2,2,1)))+ttp_row-1);
        low_deoxy_CON_R_latUp = min(dcAvg(ttp_row:tslut_CON,2,2,1))-dcAvg(t0,2,2,1);
    end
elseif peak_deoxy_CON_R_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,2,1) > dcAvg(t0,2,2,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,2,1),[ttp_deoxy_CON_R_latUp 25],[dcAvg(t0,2,2,1) dcAvg(t0,2,2,1)]) ;
        ttl_deoxy_CON_R_latUp = min(xlin);
        low_deoxy_CON_R_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,2,1) > dcAvg(t0,2,2,1)))
        ttl_deoxy_CON_R_latUp = t( find(dcAvg(ttp_row:tslut_CON,2,2,1)==max(dcAvg(ttp_row:tslut_CON,2,2,1)))+ttp_row-1);
        low_deoxy_CON_R_latUp = max(dcAvg(ttp_row:tslut_CON,2,2,1))-dcAvg(t0,2,2,1);
    end
end

if abs(max_deoxy_CON_R_latDown) > abs(min_deoxy_CON_R_latDown)
   peak_deoxy_CON_R_latDown = max_deoxy_CON_R_latDown ;
   ttp_deoxy_CON_R_latDown = t( dcAvg (:,2,3,1) == max(dcAvg(t0:tslut_INC,2,3,1))) ;
elseif abs(max_deoxy_CON_R_latDown) < abs(min_deoxy_CON_R_latDown)
   peak_deoxy_CON_R_latDown = min_deoxy_CON_R_latDown ;
   ttp_deoxy_CON_R_latDown = t( dcAvg (:,2,3,1) == min(dcAvg(t0:tslut_CON,2,3,1))) ;
elseif isnan(max_deoxy_CON_R_latDown)
   peak_deoxy_CON_R_latDown = NaN; 
   ttp_deoxy_CON_R_latDown = NaN;
elseif abs(max_deoxy_CON_R_latDown) == abs(min_deoxy_CON_R_latDown)
   ttp_deoxy_CON_R_latDown_max = t( dcAvg (:,2,3,1) == max(dcAvg(t0:tslut_INC,2,3,1))) ;
   ttp_deoxy_CON_R_latDown_min = t( dcAvg (:,2,3,1) == min(dcAvg(t0:tslut_CON,2,3,1))) ;
    if ttp_deoxy_CON_R_latDown_max < ttp_deoxy_CON_R_latDown_min
       ttp_deoxy_CON_R_latDown = ttp_deoxy_CON_R_latDown_max;
       peak_deoxy_CON_R_latDown = max_deoxy_CON_R_latDown ; 
    elseif ttp_deoxy_CON_R_latDown_max > ttp_deoxy_CON_R_latDown_min
       ttp_deoxy_CON_R_latDown = ttp_deoxy_CON_R_latDown_min;
       peak_deoxy_CON_R_latDown = min_deoxy_CON_R_latDown ;
    end
end

ttp_row = find(t==ttp_deoxy_CON_R_latDown);
if isnan(peak_deoxy_CON_R_latDown)
    ttl_deoxy_CON_R_latDown = NaN ;
    low_deoxy_CON_R_latDown = NaN ;
elseif peak_deoxy_CON_R_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,3,1) < dcAvg(t0,2,3,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,3,1),[ttp_deoxy_CON_R_latDown 25],[dcAvg(t0,2,3,1) dcAvg(t0,2,3,1)]) ;
        ttl_deoxy_CON_R_latDown = min(xlin);
        low_deoxy_CON_R_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,3,1) < dcAvg(t0,2,3,1)))
        ttl_deoxy_CON_R_latDown = t( find(dcAvg(ttp_row:tslut_CON,2,3,1)==min(dcAvg(ttp_row:tslut_CON,2,3,1)))+ttp_row-1);
        low_deoxy_CON_R_latDown = min(dcAvg(ttp_row:tslut_CON,2,3,1))-dcAvg(t0,2,3,1);
    end
elseif peak_deoxy_CON_R_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,3,1) > dcAvg(t0,2,3,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,3,1),[ttp_deoxy_CON_R_latDown 25],[dcAvg(t0,2,3,1) dcAvg(t0,2,3,1)]) ;
        ttl_deoxy_CON_R_latDown = min(xlin);
        low_deoxy_CON_R_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,3,1) > dcAvg(t0,2,3,1)))
        ttl_deoxy_CON_R_latDown = t( find(dcAvg(ttp_row:tslut_CON,2,3,1)==max(dcAvg(ttp_row:tslut_CON,2,3,1)))+ttp_row-1);
        low_deoxy_CON_R_latDown = max(dcAvg(ttp_row:tslut_CON,2,3,1))-dcAvg(t0,2,3,1);
    end
end

if abs(max_deoxy_CON_R_medUp) > abs(min_deoxy_CON_R_medUp)
   peak_deoxy_CON_R_medUp = max_deoxy_CON_R_medUp ;
   ttp_deoxy_CON_R_medUp = t( dcAvg (:,2,4,1) == max(dcAvg(t0:tslut_INC,2,4,1))) ;
elseif abs(max_deoxy_CON_R_medUp) < abs(min_deoxy_CON_R_medUp)
   peak_deoxy_CON_R_medUp = min_deoxy_CON_R_medUp ;
   ttp_deoxy_CON_R_medUp = t( dcAvg (:,2,4,1) == min(dcAvg(t0:tslut_CON,2,4,1))) ;
elseif isnan(max_deoxy_CON_R_medUp)
   ttp_deoxy_CON_R_medUp = NaN;
   peak_deoxy_CON_R_medUp = NaN ;
elseif abs(max_deoxy_CON_R_medUp) == abs(min_deoxy_CON_R_medUp)
   ttp_deoxy_CON_R_medUp_max = t( dcAvg (:,2,4,1) == max(dcAvg(t0:tslut_INC,2,4,1))) ;
   ttp_deoxy_CON_R_medUp_min = t( dcAvg (:,2,4,1) == min(dcAvg(t0:tslut_CON,2,4,1))) ;
   if ttp_deoxy_CON_R_medUp_max < ttp_deoxy_CON_R_medUp_min
       ttp_deoxy_CON_R_medUp = ttp_deoxy_CON_R_medUp_max;
       peak_deoxy_CON_R_medUp = max_deoxy_CON_R_medUp ; 
   elseif ttp_deoxy_CON_R_medUp_max > ttp_deoxy_CON_R_medUp_min
       ttp_deoxy_CON_R_medUp = ttp_deoxy_CON_R_medUp_min;
       peak_deoxy_CON_R_medUp = min_deoxy_CON_R_medUp ;
   end
end

ttp_row = find(t==ttp_deoxy_CON_R_medUp);
if isnan(peak_deoxy_CON_R_medUp)
    ttl_deoxy_CON_R_medUp = NaN ;
    low_deoxy_CON_R_medUp = NaN ;
elseif peak_deoxy_CON_R_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,4,1) < dcAvg(t0,2,4,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,4,1),[ttp_deoxy_CON_R_medUp 25],[dcAvg(t0,2,4,1) dcAvg(t0,2,4,1)]) ;
        ttl_deoxy_CON_R_medUp = min(xlin);
        low_deoxy_CON_R_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,4,1) < dcAvg(t0,2,4,1)))
        ttl_deoxy_CON_R_medUp = t( find(dcAvg(ttp_row:tslut_CON,2,4,1)==min(dcAvg(ttp_row:tslut_CON,2,4,1)))+ttp_row-1);
        low_deoxy_CON_R_medUp = min(dcAvg(ttp_row:tslut_CON,2,4,1))-dcAvg(t0,2,4,1);
    end
elseif peak_deoxy_CON_R_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,4,1) > dcAvg(t0,2,4,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,4,1),[ttp_deoxy_CON_R_medUp 25],[dcAvg(t0,2,4,1) dcAvg(t0,2,4,1)]) ;
        ttl_deoxy_CON_R_medUp = min(xlin);
        low_deoxy_CON_R_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,4,1) > dcAvg(t0,2,4,1)))
        ttl_deoxy_CON_R_medUp = t( find(dcAvg(ttp_row:tslut_CON,2,4,1)==max(dcAvg(ttp_row:tslut_CON,2,4,1)))+ttp_row-1);
        low_deoxy_CON_R_medUp = max(dcAvg(ttp_row:tslut_CON,2,4,1))-dcAvg(t0,2,4,1);
    end
end

if abs(max_deoxy_CON_R_medDown) > abs(min_deoxy_CON_R_medDown)
   peak_deoxy_CON_R_medDown = max_deoxy_CON_R_medDown ;
   ttp_deoxy_CON_R_medDown = t( dcAvg (:,2,5,1) == max(dcAvg(t0:tslut_INC,2,5,1))) ;
elseif abs(max_deoxy_CON_R_medDown) < abs(min_deoxy_CON_R_medDown)
   peak_deoxy_CON_R_medDown = min_deoxy_CON_R_medDown ;
   ttp_deoxy_CON_R_medDown = t( dcAvg (:,2,5,1) == min(dcAvg(t0:tslut_CON,2,5,1))) ;
elseif isnan(max_deoxy_CON_R_medDown)
   ttp_deoxy_CON_R_medDown = NaN;
   peak_deoxy_CON_R_medDown = NaN ;
elseif abs(max_deoxy_CON_R_medDown) == abs(min_deoxy_CON_R_medDown)
   ttp_deoxy_CON_R_medDown_max = t( dcAvg (:,2,5,1) == max(dcAvg(t0:tslut_INC,2,5,1))) ;
   ttp_deoxy_CON_R_medDown_min = t( dcAvg (:,2,5,1) == min(dcAvg(t0:tslut_CON,2,5,1))) ;
   if ttp_deoxy_CON_R_medDown_max < ttp_deoxy_CON_R_medDown_min
       ttp_deoxy_CON_R_medDown = ttp_deoxy_CON_R_medDown_max;
       peak_deoxy_CON_R_medDown = max_deoxy_CON_R_medDown ; 
   elseif ttp_deoxy_CON_R_medDown_max > ttp_deoxy_CON_R_medDown_min
       ttp_deoxy_CON_R_medDown = ttp_deoxy_CON_R_medDown_min;
       peak_deoxy_CON_R_medDown = min_deoxy_CON_R_medDown ;
   end
end

ttp_row = find(t==ttp_deoxy_CON_R_medDown);
if isnan(peak_deoxy_CON_R_medDown)
    ttl_deoxy_CON_R_medDown = NaN ;
    low_deoxy_CON_R_medDown = NaN ;
elseif peak_deoxy_CON_R_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,5,1) < dcAvg(t0,2,5,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,5,1),[ttp_deoxy_CON_R_medDown 25],[dcAvg(t0,2,5,1) dcAvg(t0,2,5,1)]) ;
        ttl_deoxy_CON_R_medDown = min(xlin);
        low_deoxy_CON_R_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,5,1) < dcAvg(t0,2,5,1)))
        ttl_deoxy_CON_R_medDown = t( find(dcAvg(ttp_row:tslut_CON,2,5,1)==min(dcAvg(ttp_row:tslut_CON,2,5,1)))+ttp_row-1);
        low_deoxy_CON_R_medDown = min(dcAvg(ttp_row:tslut_CON,2,5,1))-dcAvg(t0,2,5,1);
    end
elseif peak_deoxy_CON_R_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,5,1) > dcAvg(t0,2,5,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,5,1),[ttp_deoxy_CON_R_medDown 25],[dcAvg(t0,2,5,1) dcAvg(t0,2,5,1)]) ;
        ttl_deoxy_CON_R_medDown = min(xlin);
        low_deoxy_CON_R_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,5,1) > dcAvg(t0,2,5,1)))
        ttl_deoxy_CON_R_medDown = t( find(dcAvg(ttp_row:tslut_CON,2,5,1)==max(dcAvg(ttp_row:tslut_CON,2,5,1)))+ttp_row-1);
        low_deoxy_CON_R_medDown = max(dcAvg(ttp_row:tslut_CON,2,5,1))-dcAvg(t0,2,5,1);
    end
end

if abs(max_deoxy_CON_L_medUp) > abs(min_deoxy_CON_L_medUp)
   peak_deoxy_CON_L_medUp = max_deoxy_CON_L_medUp ;
   ttp_deoxy_CON_L_medUp = t( dcAvg (:,2,8,1) == max(dcAvg(t0:tslut_INC,2,8,1))) ;
elseif abs(max_deoxy_CON_L_medUp) < abs(min_deoxy_CON_L_medUp)
   peak_deoxy_CON_L_medUp = min_deoxy_CON_L_medUp ;
   ttp_deoxy_CON_L_medUp = t( dcAvg (:,2,8,1) == min(dcAvg(t0:tslut_CON,2,8,1))) ;
elseif isnan(max_deoxy_CON_L_medUp)
   ttp_deoxy_CON_L_medUp = NaN;
   peak_deoxy_CON_L_medUp = NaN ;
elseif abs(max_deoxy_CON_L_medUp) == abs(min_deoxy_CON_L_medUp)
   ttp_deoxy_CON_L_medUp_max = t( dcAvg (:,2,8,1) == max(dcAvg(t0:tslut_INC,2,8,1))) ;
   ttp_deoxy_CON_L_medUp_min = t( dcAvg (:,2,8,1) == min(dcAvg(t0:tslut_CON,2,8,1))) ;
   if ttp_deoxy_CON_L_medUp_max < ttp_deoxy_CON_L_medUp_min
       ttp_deoxy_CON_L_medUp = ttp_deoxy_CON_L_medUp_max;
       peak_deoxy_CON_L_medUp = max_deoxy_CON_L_medUp ; 
   elseif ttp_deoxy_CON_L_medUp_max > ttp_deoxy_CON_L_medUp_min
       ttp_deoxy_CON_L_medUp = ttp_deoxy_CON_L_medUp_min;
       peak_deoxy_CON_L_medUp = min_deoxy_CON_L_medUp ;
   end
end

ttp_row = find(t==ttp_deoxy_CON_L_medUp);
if isnan(peak_deoxy_CON_L_medUp)
    ttl_deoxy_CON_L_medUp = NaN ;
    low_deoxy_CON_L_medUp = NaN ;
elseif peak_deoxy_CON_L_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,8,1) < dcAvg(t0,2,8,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,8,1),[ttp_deoxy_CON_L_medUp 25],[dcAvg(t0,2,8,1) dcAvg(t0,2,8,1)]) ;
        ttl_deoxy_CON_L_medUp = min(xlin);
        low_deoxy_CON_L_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,8,1) < dcAvg(t0,2,8,1)))
        ttl_deoxy_CON_L_medUp = t( find(dcAvg(ttp_row:tslut_CON,2,8,1)==min(dcAvg(ttp_row:tslut_CON,2,8,1)))+ttp_row-1);
        low_deoxy_CON_L_medUp = min(dcAvg(ttp_row:tslut_CON,2,8,1))-dcAvg(t0,2,8,1);
    end
elseif peak_deoxy_CON_L_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,8,1) > dcAvg(t0,2,8,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,8,1),[ttp_deoxy_CON_L_medUp 25],[dcAvg(t0,2,8,1) dcAvg(t0,2,8,1)]) ;
        ttl_deoxy_CON_L_medUp = min(xlin);
        low_deoxy_CON_L_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,8,1) > dcAvg(t0,2,8,1)))
        ttl_deoxy_CON_L_medUp = t( find(dcAvg(ttp_row:tslut_CON,2,8,1)==max(dcAvg(ttp_row:tslut_CON,2,8,1)))+ttp_row-1);
        low_deoxy_CON_L_medUp = max(dcAvg(ttp_row:tslut_CON,2,8,1))-dcAvg(t0,2,8,1);
    end
end

if abs(max_deoxy_CON_L_medDown) > abs(min_deoxy_CON_L_medDown)
   peak_deoxy_CON_L_medDown = max_deoxy_CON_L_medDown ;
   ttp_deoxy_CON_L_medDown = t( dcAvg (:,2,9,1) == max(dcAvg(t0:tslut_INC,2,9,1))) ;
elseif abs(max_deoxy_CON_L_medDown) < abs(min_deoxy_CON_L_medDown)
   peak_deoxy_CON_L_medDown = min_deoxy_CON_L_medDown ;
   ttp_deoxy_CON_L_medDown = t( dcAvg (:,2,9,1) == min(dcAvg(t0:tslut_CON,2,9,1))) ;
elseif isnan(max_deoxy_CON_L_medDown)
   ttp_deoxy_CON_L_medDown = NaN;
   peak_deoxy_CON_L_medDown = NaN ;
elseif abs(max_deoxy_CON_L_medDown) == abs(min_deoxy_CON_L_medDown)
   ttp_deoxy_CON_L_medDown_max = t( dcAvg (:,2,9,1) == max(dcAvg(t0:tslut_INC,2,9,1))) ;
   ttp_deoxy_CON_L_medDown_min = t( dcAvg (:,2,9,1) == min(dcAvg(t0:tslut_CON,2,9,1))) ;
   if ttp_deoxy_CON_L_medDown_max < ttp_deoxy_CON_L_medDown_min
       ttp_deoxy_CON_L_medDown = ttp_deoxy_CON_L_medDown_max;
       peak_deoxy_CON_L_medDown = max_deoxy_CON_L_medDown ; 
   elseif ttp_deoxy_CON_L_medDown_max > ttp_deoxy_CON_L_medDown_min
       ttp_deoxy_CON_L_medDown = ttp_deoxy_CON_L_medDown_min;
       peak_deoxy_CON_L_medDown = min_deoxy_CON_L_medDown ;
   end
end

ttp_row = find(t==ttp_deoxy_CON_L_medDown);
if isnan(peak_deoxy_CON_L_medDown)
    ttl_deoxy_CON_L_medDown = NaN ;
    low_deoxy_CON_L_medDown = NaN ;
elseif peak_deoxy_CON_L_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,9,1) < dcAvg(t0,2,9,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,9,1),[ttp_deoxy_CON_L_medDown 25],[dcAvg(t0,2,9,1) dcAvg(t0,2,9,1)]) ;
        ttl_deoxy_CON_L_medDown = min(xlin);
        low_deoxy_CON_L_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,9,1) < dcAvg(t0,2,9,1)))
        ttl_deoxy_CON_L_medDown = t( find(dcAvg(ttp_row:tslut_CON,2,9,1)==min(dcAvg(ttp_row:tslut_CON,2,9,1)))+ttp_row-1);
        low_deoxy_CON_L_medDown = min(dcAvg(ttp_row:tslut_CON,2,9,1))-dcAvg(t0,2,9,1);
    end
elseif peak_deoxy_CON_L_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,9,1) > dcAvg(t0,2,9,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,9,1),[ttp_deoxy_CON_L_medDown 25],[dcAvg(t0,2,9,1) dcAvg(t0,2,9,1)]) ;
        ttl_deoxy_CON_L_medDown = min(xlin);
        low_deoxy_CON_L_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,9,1) > dcAvg(t0,2,9,1)))
        ttl_deoxy_CON_L_medDown = t( find(dcAvg(ttp_row:tslut_CON,2,9,1)==max(dcAvg(ttp_row:tslut_CON,2,9,1)))+ttp_row-1);
        low_deoxy_CON_L_medDown = max(dcAvg(ttp_row:tslut_CON,2,9,1))-dcAvg(t0,2,9,1);
    end
end

if abs(max_deoxy_CON_L_latUp) > abs(min_deoxy_CON_L_latUp)
   peak_deoxy_CON_L_latUp = max_deoxy_CON_L_latUp ;
   ttp_deoxy_CON_L_latUp = t( dcAvg (:,2,10,1) == max(dcAvg(t0:tslut_INC,2,10,1))) ;
elseif abs(max_deoxy_CON_L_latUp) < abs(min_deoxy_CON_L_latUp)
   peak_deoxy_CON_L_latUp = min_deoxy_CON_L_latUp ;
   ttp_deoxy_CON_L_latUp = t( dcAvg (:,2,10,1) == min(dcAvg(t0:tslut_CON,2,10,1))) ;
elseif isnan(max_deoxy_CON_L_latUp)
   ttp_deoxy_CON_L_latUp = NaN;
   peak_deoxy_CON_L_latUp = NaN ;
elseif abs(max_deoxy_CON_L_latUp) == abs(min_deoxy_CON_L_latUp)
   ttp_deoxy_CON_L_latUp_max = t( dcAvg (:,2,10,1) == max(dcAvg(t0:tslut_INC,2,10,1))) ;
   ttp_deoxy_CON_L_latUp_min = t( dcAvg (:,2,10,1) == min(dcAvg(t0:tslut_CON,2,10,1))) ;
   if ttp_deoxy_CON_L_latUp_max < ttp_deoxy_CON_L_latUp_min
       ttp_deoxy_CON_L_latUp = ttp_deoxy_CON_L_latUp_max;
       peak_deoxy_CON_L_latUp = max_deoxy_CON_L_latUp ; 
   elseif ttp_deoxy_CON_L_latUp_max > ttp_deoxy_CON_L_latUp_min
       ttp_deoxy_CON_L_latUp = ttp_deoxy_CON_L_latUp_min;
       peak_deoxy_CON_L_latUp = min_deoxy_CON_L_latUp ;
   end
end
  
ttp_row = find(t==ttp_deoxy_CON_L_latUp);
if isnan(peak_deoxy_CON_L_latUp)
    ttl_deoxy_CON_L_latUp = NaN ;
    low_deoxy_CON_L_latUp = NaN ;
elseif peak_deoxy_CON_L_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,10,1) < dcAvg(t0,2,10,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,10,1),[ttp_deoxy_CON_L_latUp 25],[dcAvg(t0,2,10,1) dcAvg(t0,2,10,1)]) ;
        ttl_deoxy_CON_L_latUp = min(xlin);
        low_deoxy_CON_L_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,10,1) < dcAvg(t0,2,10,1)))
        ttl_deoxy_CON_L_latUp = t( find(dcAvg(ttp_row:tslut_CON,2,10,1)==min(dcAvg(ttp_row:tslut_CON,2,10,1)))+ttp_row-1);
        low_deoxy_CON_L_latUp = min(dcAvg(ttp_row:tslut_CON,2,10,1))-dcAvg(t0,2,10,1);
    end
elseif peak_deoxy_CON_L_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,10,1) > dcAvg(t0,2,10,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,10,1),[ttp_deoxy_CON_L_latUp 25],[dcAvg(t0,2,10,1) dcAvg(t0,2,10,1)]) ;
        ttl_deoxy_CON_L_latUp = min(xlin);
        low_deoxy_CON_L_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,10,1) > dcAvg(t0,2,10,1)))
        ttl_deoxy_CON_L_latUp = t( find(dcAvg(ttp_row:tslut_CON,2,10,1)==max(dcAvg(ttp_row:tslut_CON,2,10,1)))+ttp_row-1);
        low_deoxy_CON_L_latUp = max(dcAvg(ttp_row:tslut_CON,2,10,1))-dcAvg(t0,2,10,1);
    end
end

if abs(max_deoxy_CON_L_latDown) > abs(min_deoxy_CON_L_latDown)
   peak_deoxy_CON_L_latDown = max_deoxy_CON_L_latDown ;
   ttp_deoxy_CON_L_latDown = t( dcAvg (:,2,11,1) == max(dcAvg(t0:tslut_INC,2,11,1))) ;
elseif abs(max_deoxy_CON_L_latDown) < abs(min_deoxy_CON_L_latDown)
   peak_deoxy_CON_L_latDown = min_deoxy_CON_L_latDown ;
   ttp_deoxy_CON_L_latDown = t( dcAvg (:,2,11,1) == min(dcAvg(t0:tslut_CON,2,11,1))) ;
elseif isnan(max_deoxy_CON_L_latDown)
   ttp_deoxy_CON_L_latDown = NaN;
   peak_deoxy_CON_L_latDown = NaN ;
elseif abs(max_deoxy_CON_L_latDown) == abs(min_deoxy_CON_L_latDown)
   ttp_deoxy_CON_L_latDown_max = t( dcAvg (:,2,11,1) == max(dcAvg(t0:tslut_INC,2,11,1))) ;
   ttp_deoxy_CON_L_latDown_min = t( dcAvg (:,2,11,1) == min(dcAvg(t0:tslut_CON,2,11,1))) ;
   if ttp_deoxy_CON_L_latDown_max < ttp_deoxy_CON_L_latDown_min
       ttp_deoxy_CON_L_latDown = ttp_deoxy_CON_L_latDown_max;
       peak_deoxy_CON_L_latDown = max_deoxy_CON_L_latDown ; 
   elseif ttp_deoxy_CON_L_latDown_max > ttp_deoxy_CON_L_latDown_min
       ttp_deoxy_CON_L_latDown = ttp_deoxy_CON_L_latDown_min;
       peak_deoxy_CON_L_latDown = min_deoxy_CON_L_latDown ;
   end
end

ttp_row = find(t==ttp_deoxy_CON_L_latDown);
if isnan(peak_deoxy_CON_L_latDown)
    ttl_deoxy_CON_L_latDown = NaN ;
    low_deoxy_CON_L_latDown = NaN ;
elseif peak_deoxy_CON_L_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,11,1) < dcAvg(t0,2,11,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,11,1),[ttp_deoxy_CON_L_latDown 25],[dcAvg(t0,2,11,1) dcAvg(t0,2,11,1)]) ;
        ttl_deoxy_CON_L_latDown = min(xlin);
        low_deoxy_CON_L_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,11,1) < dcAvg(t0,2,11,1)))
        ttl_deoxy_CON_L_latDown = t( find(dcAvg(ttp_row:tslut_CON,2,11,1)==min(dcAvg(ttp_row:tslut_CON,2,11,1)))+ttp_row-1);
        low_deoxy_CON_L_latDown = min(dcAvg(ttp_row:tslut_CON,2,11,1))-dcAvg(t0,2,11,1);
    end
elseif peak_deoxy_CON_L_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_CON,2,11,1) > dcAvg(t0,2,11,1)))
        [xlin,~] = intersections(t(ttp_row:tslut_CON),dcAvg(ttp_row:tslut_CON,2,11,1),[ttp_deoxy_CON_L_latDown 25],[dcAvg(t0,2,11,1) dcAvg(t0,2,11,1)]) ;
        ttl_deoxy_CON_L_latDown = min(xlin);
        low_deoxy_CON_L_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_CON,2,11,1) > dcAvg(t0,2,11,1)))
        ttl_deoxy_CON_L_latDown = t( find(dcAvg(ttp_row:tslut_CON,2,11,1)==max(dcAvg(ttp_row:tslut_CON,2,11,1)))+ttp_row-1);
        low_deoxy_CON_L_latDown = max(dcAvg(ttp_row:tslut_CON,2,11,1))-dcAvg(t0,2,11,1);
    end
end

%% INCONGRUENT OXY
% Finding peak, time-to-peak (TTP), return to baseline (low) and
% time-to-low (TTL)

% Max. oxy during INCONGRUENT
max_oxy_INC_R_latUp = max(dcAvg(t0:tslut_INC,1,2,2))-dcAvg(t0,1,2,2); 
max_oxy_INC_R_latDown = max(dcAvg(t0:tslut_INC,1,3,2))-dcAvg(t0,1,3,2);
max_oxy_INC_R_medUp = max(dcAvg(t0:tslut_INC,1,4,2))-dcAvg(t0,1,4,2);
max_oxy_INC_R_medDown = max(dcAvg(t0:tslut_INC,1,5,2))-dcAvg(t0,1,3,2); 
max_oxy_INC_L_medUp = max(dcAvg(t0:tslut_INC,1,8,2))-dcAvg(t0,1,8,2);
max_oxy_INC_L_medDown = max(dcAvg(t0:tslut_INC,1,9,2))-dcAvg(t0,1,9,2);
max_oxy_INC_L_latUp = max(dcAvg(t0:tslut_INC,1,10,2))-dcAvg(t0,1,10,2);
max_oxy_INC_L_latDown = max(dcAvg(t0:tslut_INC,1,11,2))-dcAvg(t0,1,11,2);

% MIn. oxy during INCONGRUENT
min_oxy_INC_R_latUp = min(dcAvg(t0:tslut_INC,1,2,2))-dcAvg(t0,1,2,2);
min_oxy_INC_R_latDown = min(dcAvg(t0:tslut_INC,1,3,2))-dcAvg(t0,1,3,2);
min_oxy_INC_R_medUp = min(dcAvg(t0:tslut_INC,1,4,2))-dcAvg(t0,1,4,2);
min_oxy_INC_R_medDown = min(dcAvg(t0:tslut_INC,1,5,2))-dcAvg(t0,1,5,2);
min_oxy_INC_L_medUp = min(dcAvg(t0:tslut_INC,1,8,2))-dcAvg(t0,1,8,2);
min_oxy_INC_L_medDown = min(dcAvg(t0:tslut_INC,1,9,2))-dcAvg(t0,1,9,2);
min_oxy_INC_L_latUp = min(dcAvg(t0:tslut_INC,1,10,2))-dcAvg(t0,1,10,2);
min_oxy_INC_L_latDown = min(dcAvg(t0:tslut_INC,1,11,2))-dcAvg(t0,1,11,2);

% Finding peak (min. or max.) and time-to-peak during INCONGRUENT 
% If identical, then use the first one
if abs(max_oxy_INC_R_latUp) > abs(min_oxy_INC_R_latUp)
   peak_oxy_INC_R_latUp = max_oxy_INC_R_latUp ;
   ttp_oxy_INC_R_latUp = t( dcAvg (:,1,2,2) == max(dcAvg(t0:tslut_INC,1,2,2))) ;
elseif abs(max_oxy_INC_R_latUp) < abs(min_oxy_INC_R_latUp)
   peak_oxy_INC_R_latUp = min_oxy_INC_R_latUp ;
   ttp_oxy_INC_R_latUp = t( dcAvg (:,1,2,2) == min(dcAvg(t0:tslut_INC,1,2,2))) ;
elseif isnan(max_oxy_INC_R_latUp)
   ttp_oxy_INC_R_latUp = NaN;
   peak_oxy_INC_R_latUp = NaN ;
elseif abs(max_oxy_INC_R_latUp) == abs(min_oxy_INC_R_latUp)
   ttp_oxy_INC_R_latUp_max = t( dcAvg (:,1,2,2) == max(dcAvg(t0:tslut_INC,1,2,2))) ;
   ttp_oxy_INC_R_latUp_min = t( dcAvg (:,1,2,2) == min(dcAvg(t0:tslut_INC,1,2,2))) ;
   if ttp_oxy_INC_R_latUp_max < ttp_oxy_INC_R_latUp_min
       ttp_oxy_INC_R_latUp = ttp_oxy_INC_R_latUp_max;
       peak_oxy_INC_R_latUp = max_oxy_INC_R_latUp ; 
   elseif ttp_oxy_INC_R_latUp_max > ttp_oxy_INC_R_latUp_min
       ttp_oxy_INC_R_latUp = ttp_oxy_INC_R_latUp_min;
       peak_oxy_INC_R_latUp = min_oxy_INC_R_latUp ;
   end
end
  
ttp_row = find(t==ttp_oxy_INC_R_latUp);
if isnan(peak_oxy_INC_R_latUp)
    ttl_oxy_INC_R_latUp = NaN ;
    low_oxy_INC_R_latUp = NaN ;
elseif peak_oxy_INC_R_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,2,2) < dcAvg(t0,1,2,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,2,2),[ttp_oxy_INC_R_latUp 25],[dcAvg(t0,1,2,2) dcAvg(t0,1,2,2)]) ;
        ttl_oxy_INC_R_latUp = min(xlin);
        low_oxy_INC_R_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,2,2) < dcAvg(t0,1,2,2)))
        ttl_oxy_INC_R_latUp = t( find(dcAvg(ttp_row:tslut_INC,1,2,2)==min(dcAvg(ttp_row:tslut_INC,1,2,2)))+ttp_row-1);
        low_oxy_INC_R_latUp = min(dcAvg(ttp_row:tslut_INC,1,2,2))-dcAvg(t0,1,2,2);
    end
elseif peak_oxy_INC_R_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,2,2) > dcAvg(t0,1,2,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,2,2),[ttp_oxy_INC_R_latUp 25],[dcAvg(t0,1,2,2) dcAvg(t0,1,2,2)]) ;
        ttl_oxy_INC_R_latUp = min(xlin);
        low_oxy_INC_R_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,2,2) > dcAvg(t0,1,2,2)))
        ttl_oxy_INC_R_latUp = t( find(dcAvg(ttp_row:tslut_INC,1,2,2)==max(dcAvg(ttp_row:tslut_INC,1,2,2)))+ttp_row-1);
        low_oxy_INC_R_latUp = max(dcAvg(ttp_row:tslut_INC,1,2,2))-dcAvg(t0,1,2,2);
    end
end

if abs(max_oxy_INC_R_latDown) > abs(min_oxy_INC_R_latDown)
   peak_oxy_INC_R_latDown = max_oxy_INC_R_latDown ;
   ttp_oxy_INC_R_latDown = t( dcAvg (:,1,3,2) == max(dcAvg(t0:tslut_INC,1,3,2))) ;
elseif abs(max_oxy_INC_R_latDown) < abs(min_oxy_INC_R_latDown)
   peak_oxy_INC_R_latDown = min_oxy_INC_R_latDown ;
   ttp_oxy_INC_R_latDown = t( dcAvg (:,1,3,2) == min(dcAvg(t0:tslut_INC,1,3,2))) ;
elseif isnan(max_oxy_INC_R_latDown)
   peak_oxy_INC_R_latDown = NaN;
   ttp_oxy_INC_R_latDown = NaN;
elseif abs(max_oxy_INC_R_latDown) == abs(min_oxy_INC_R_latDown)
   ttp_oxy_INC_R_latDown_max = t( dcAvg (:,1,3,2) == max(dcAvg(t0:tslut_INC,1,3,2))) ;
   ttp_oxy_INC_R_latDown_min = t( dcAvg (:,1,3,2) == min(dcAvg(t0:tslut_INC,1,3,2))) ;
    if ttp_oxy_INC_R_latDown_max < ttp_oxy_INC_R_latDown_min
       ttp_oxy_INC_R_latDown = ttp_oxy_INC_R_latDown_max;
       peak_oxy_INC_R_latDown = max_oxy_INC_R_latDown ; 
    elseif ttp_oxy_INC_R_latDown_max > ttp_oxy_INC_R_latDown_min
       ttp_oxy_INC_R_latDown = ttp_oxy_INC_R_latDown_min;
       peak_oxy_INC_R_latDown = min_oxy_INC_R_latDown ;
    end
end

ttp_row = find(t==ttp_oxy_INC_R_latDown);
if isnan(peak_oxy_INC_R_latDown)
    ttl_oxy_INC_R_latDown = NaN ;
    low_oxy_INC_R_latDown = NaN ;
elseif peak_oxy_INC_R_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,3,2) < dcAvg(t0,1,3,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,3,2),[ttp_oxy_INC_R_latDown 25],[dcAvg(t0,1,3,2) dcAvg(t0,1,3,2)]) ;
        ttl_oxy_INC_R_latDown = min(xlin);
        low_oxy_INC_R_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,3,2) < dcAvg(t0,1,3,2)))
        ttl_oxy_INC_R_latDown = t( find(dcAvg(ttp_row:tslut_INC,1,3,2)==min(dcAvg(ttp_row:tslut_INC,1,3,2)))+ttp_row-1);
        low_oxy_INC_R_latDown = min(dcAvg(ttp_row:tslut_INC,1,3,2))-dcAvg(t0,1,3,2);
    end
elseif peak_oxy_INC_R_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,3,2) > dcAvg(t0,1,3,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,3,2),[ttp_oxy_INC_R_latDown 25],[dcAvg(t0,1,3,2) dcAvg(t0,1,3,2)]) ;
        ttl_oxy_INC_R_latDown = min(xlin);
        low_oxy_INC_R_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,3,2) > dcAvg(t0,1,3,2)))
        ttl_oxy_INC_R_latDown = t( find(dcAvg(ttp_row:tslut_INC,1,3,2)==max(dcAvg(ttp_row:tslut_INC,1,3,2)))+ttp_row-1);
        low_oxy_INC_R_latDown = max(dcAvg(ttp_row:tslut_INC,1,3,2))-dcAvg(t0,1,3,2);
    end
end

if abs(max_oxy_INC_R_medUp) > abs(min_oxy_INC_R_medUp)
   peak_oxy_INC_R_medUp = max_oxy_INC_R_medUp ;
   ttp_oxy_INC_R_medUp = t( dcAvg (:,1,4,2) == max(dcAvg(t0:tslut_INC,1,4,2))) ;
elseif abs(max_oxy_INC_R_medUp) < abs(min_oxy_INC_R_medUp)
   peak_oxy_INC_R_medUp = min_oxy_INC_R_medUp ;
   ttp_oxy_INC_R_medUp = t( dcAvg (:,1,4,2) == min(dcAvg(t0:tslut_INC,1,4,2))) ;
elseif isnan(max_oxy_INC_R_medUp)
   ttp_oxy_INC_R_medUp = NaN;
   peak_oxy_INC_R_medUp = NaN ;
elseif abs(max_oxy_INC_R_medUp) == abs(min_oxy_INC_R_medUp)
   ttp_oxy_INC_R_medUp_max = t( dcAvg (:,1,4,2) == max(dcAvg(t0:tslut_INC,1,4,2))) ;
   ttp_oxy_INC_R_medUp_min = t( dcAvg (:,1,4,2) == min(dcAvg(t0:tslut_INC,1,4,2))) ;
   if ttp_oxy_INC_R_medUp_max < ttp_oxy_INC_R_medUp_min
       ttp_oxy_INC_R_medUp = ttp_oxy_INC_R_medUp_max;
       peak_oxy_INC_R_medUp = max_oxy_INC_R_medUp ; 
   elseif ttp_oxy_INC_R_medUp_max > ttp_oxy_INC_R_medUp_min
       ttp_oxy_INC_R_medUp = ttp_oxy_INC_R_medUp_min;
       peak_oxy_INC_R_medUp = min_oxy_INC_R_medUp ;
   end
end

ttp_row = find(t==ttp_oxy_INC_R_medUp);
if isnan(peak_oxy_INC_R_medUp)
    ttl_oxy_INC_R_medUp = NaN ;
    low_oxy_INC_R_medUp = NaN ;
elseif peak_oxy_INC_R_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,4,2) < dcAvg(t0,1,4,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,4,2),[ttp_oxy_INC_R_medUp 25],[dcAvg(t0,1,4,2) dcAvg(t0,1,4,2)]) ;
        ttl_oxy_INC_R_medUp = min(xlin);
        low_oxy_INC_R_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,4,2) < dcAvg(t0,1,4,2)))
        ttl_oxy_INC_R_medUp = t( find(dcAvg(ttp_row:tslut_INC,1,4,2)==min(dcAvg(ttp_row:tslut_INC,1,4,2)))+ttp_row-1);
        low_oxy_INC_R_medUp = min(dcAvg(ttp_row:tslut_INC,1,4,2))-dcAvg(t0,1,4,2);
    end
elseif peak_oxy_INC_R_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,4,2) > dcAvg(t0,1,4,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,4,2),[ttp_oxy_INC_R_medUp 25],[dcAvg(t0,1,4,2) dcAvg(t0,1,4,2)]) ;
        ttl_oxy_INC_R_medUp = min(xlin);
        low_oxy_INC_R_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,4,2) > dcAvg(t0,1,4,2)))
        ttl_oxy_INC_R_medUp = t( find(dcAvg(ttp_row:tslut_INC,1,4,2)==max(dcAvg(ttp_row:tslut_INC,1,4,2)))+ttp_row-1);
        low_oxy_INC_R_medUp = max(dcAvg(ttp_row:tslut_INC,1,4,2))-dcAvg(t0,1,4,2);
    end
end

if abs(max_oxy_INC_R_medDown) > abs(min_oxy_INC_R_medDown)
   peak_oxy_INC_R_medDown = max_oxy_INC_R_medDown ;
   ttp_oxy_INC_R_medDown = t( dcAvg (:,1,5,2) == max(dcAvg(t0:tslut_INC,1,5,2))) ;
elseif abs(max_oxy_INC_R_medDown) < abs(min_oxy_INC_R_medDown)
   peak_oxy_INC_R_medDown = min_oxy_INC_R_medDown ;
   ttp_oxy_INC_R_medDown = t( dcAvg (:,1,5,2) == min(dcAvg(t0:tslut_INC,1,5,2))) ;
elseif isnan(max_oxy_INC_R_medDown)
   ttp_oxy_INC_R_medDown = NaN;
   peak_oxy_INC_R_medDown = NaN ;
elseif abs(max_oxy_INC_R_medDown) == abs(min_oxy_INC_R_medDown)
   ttp_oxy_INC_R_medDown_max = t( dcAvg (:,1,5,2) == max(dcAvg(t0:tslut_INC,1,5,2))) ;
   ttp_oxy_INC_R_medDown_min = t( dcAvg (:,1,5,2) == min(dcAvg(t0:tslut_INC,1,5,2))) ;
   if ttp_oxy_INC_R_medDown_max < ttp_oxy_INC_R_medDown_min
       ttp_oxy_INC_R_medDown = ttp_oxy_INC_R_medDown_max;
       peak_oxy_INC_R_medDown = max_oxy_INC_R_medDown ; 
   elseif ttp_oxy_INC_R_medDown_max > ttp_oxy_INC_R_medDown_min
       ttp_oxy_INC_R_medDown = ttp_oxy_INC_R_medDown_min;
       peak_oxy_INC_R_medDown = min_oxy_INC_R_medDown ;
   end
end

ttp_row = find(t==ttp_oxy_INC_R_medDown);
if isnan(peak_oxy_INC_R_medDown)
    ttl_oxy_INC_R_medDown = NaN ;
    low_oxy_INC_R_medDown = NaN ;
elseif peak_oxy_INC_R_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,5,2) < dcAvg(t0,1,5,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,5,2),[ttp_oxy_INC_R_medDown 25],[dcAvg(t0,1,5,2) dcAvg(t0,1,5,2)]) ;
        ttl_oxy_INC_R_medDown = min(xlin);
        low_oxy_INC_R_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,5,2) < dcAvg(t0,1,5,2)))
        ttl_oxy_INC_R_medDown = t( find(dcAvg(ttp_row:tslut_INC,1,5,2)==min(dcAvg(ttp_row:tslut_INC,1,5,2)))+ttp_row-1);
        low_oxy_INC_R_medDown = min(dcAvg(ttp_row:tslut_INC,1,5,2))-dcAvg(t0,1,5,2);
    end
elseif peak_oxy_INC_R_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,5,2) > dcAvg(t0,1,5,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,5,2),[ttp_oxy_INC_R_medDown 25],[dcAvg(t0,1,5,2) dcAvg(t0,1,5,2)]) ;
        ttl_oxy_INC_R_medDown = min(xlin);
        low_oxy_INC_R_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,5,2) > dcAvg(t0,1,5,2)))
        ttl_oxy_INC_R_medDown = t( find(dcAvg(ttp_row:tslut_INC,1,5,2)==max(dcAvg(ttp_row:tslut_INC,1,5,2)))+ttp_row-1);
        low_oxy_INC_R_medDown = max(dcAvg(ttp_row:tslut_INC,1,5,2))-dcAvg(t0,1,5,2);
    end
end

if abs(max_oxy_INC_L_medUp) > abs(min_oxy_INC_L_medUp)
   peak_oxy_INC_L_medUp = max_oxy_INC_L_medUp ;
   ttp_oxy_INC_L_medUp = t( dcAvg (:,1,8,2) == max(dcAvg(t0:tslut_INC,1,8,2))) ;
elseif abs(max_oxy_INC_L_medUp) < abs(min_oxy_INC_L_medUp)
   peak_oxy_INC_L_medUp = min_oxy_INC_L_medUp ;
   ttp_oxy_INC_L_medUp = t( dcAvg (:,1,8,2) == min(dcAvg(t0:tslut_INC,1,8,2))) ;
elseif isnan(max_oxy_INC_L_medUp)
   ttp_oxy_INC_L_medUp = NaN;
   peak_oxy_INC_L_medUp = NaN ;
elseif abs(max_oxy_INC_L_medUp) == abs(min_oxy_INC_L_medUp)
   ttp_oxy_INC_L_medUp_max = t( dcAvg (:,1,8,2) == max(dcAvg(t0:tslut_INC,1,8,2))) ;
   ttp_oxy_INC_L_medUp_min = t( dcAvg (:,1,8,2) == min(dcAvg(t0:tslut_INC,1,8,2))) ;
   if ttp_oxy_INC_L_medUp_max < ttp_oxy_INC_L_medUp_min
       ttp_oxy_INC_L_medUp = ttp_oxy_INC_L_medUp_max;
       peak_oxy_INC_L_medUp = max_oxy_INC_L_medUp ; 
   elseif ttp_oxy_INC_L_medUp_max > ttp_oxy_INC_L_medUp_min
       ttp_oxy_INC_L_medUp = ttp_oxy_INC_L_medUp_min;
       peak_oxy_INC_L_medUp = min_oxy_INC_L_medUp ;
   end
end

ttp_row = find(t==ttp_oxy_INC_L_medUp);
if isnan(peak_oxy_INC_L_medUp)
    ttl_oxy_INC_L_medUp = NaN ;
    low_oxy_INC_L_medUp = NaN ;
elseif peak_oxy_INC_L_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,8,2) < dcAvg(t0,1,8,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,8,2),[ttp_oxy_INC_L_medUp 25],[dcAvg(t0,1,8,2) dcAvg(t0,1,8,2)]) ;
        ttl_oxy_INC_L_medUp = min(xlin);
        low_oxy_INC_L_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,8,2) < dcAvg(t0,1,8,2)))
        ttl_oxy_INC_L_medUp = t( find(dcAvg(ttp_row:tslut_INC,1,8,2)==min(dcAvg(ttp_row:tslut_INC,1,8,2)))+ttp_row-1);
        low_oxy_INC_L_medUp = min(dcAvg(ttp_row:tslut_INC,1,8,2))-dcAvg(t0,1,8,2);
    end
elseif peak_oxy_INC_L_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,8,2) > dcAvg(t0,1,8,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,8,2),[ttp_oxy_INC_L_medUp 25],[dcAvg(t0,1,8,2) dcAvg(t0,1,8,2)]) ;
        ttl_oxy_INC_L_medUp = min(xlin);
        low_oxy_INC_L_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,8,2) > dcAvg(t0,1,8,2)))
        ttl_oxy_INC_L_medUp = t( find(dcAvg(ttp_row:tslut_INC,1,8,2)==max(dcAvg(ttp_row:tslut_INC,1,8,2)))+ttp_row-1);
        low_oxy_INC_L_medUp = max(dcAvg(ttp_row:tslut_INC,1,8,2))-dcAvg(t0,1,8,2);
    end
end

if abs(max_oxy_INC_L_medDown) > abs(min_oxy_INC_L_medDown)
   peak_oxy_INC_L_medDown = max_oxy_INC_L_medDown ;
   ttp_oxy_INC_L_medDown = t( dcAvg (:,1,9,2) == max(dcAvg(t0:tslut_INC,1,9,2))) ;
elseif abs(max_oxy_INC_L_medDown) < abs(min_oxy_INC_L_medDown)
   peak_oxy_INC_L_medDown = min_oxy_INC_L_medDown ;
   ttp_oxy_INC_L_medDown = t( dcAvg (:,1,9,2) == min(dcAvg(t0:tslut_INC,1,9,2))) ;
elseif isnan(max_oxy_INC_L_medDown)
   ttp_oxy_INC_L_medDown = NaN;
   peak_oxy_INC_L_medDown = NaN ;
elseif abs(max_oxy_INC_L_medDown) == abs(min_oxy_INC_L_medDown)
   ttp_oxy_INC_L_medDown_max = t( dcAvg (:,1,9,2) == max(dcAvg(t0:tslut_INC,1,9,2))) ;
   ttp_oxy_INC_L_medDown_min = t( dcAvg (:,1,9,2) == min(dcAvg(t0:tslut_INC,1,9,2))) ;
   if ttp_oxy_INC_L_medDown_max < ttp_oxy_INC_L_medDown_min
       ttp_oxy_INC_L_medDown = ttp_oxy_INC_L_medDown_max;
       peak_oxy_INC_L_medDown = max_oxy_INC_L_medDown ; 
   elseif ttp_oxy_INC_L_medDown_max > ttp_oxy_INC_L_medDown_min
       ttp_oxy_INC_L_medDown = ttp_oxy_INC_L_medDown_min;
       peak_oxy_INC_L_medDown = min_oxy_INC_L_medDown ;
   end
end

ttp_row = find(t==ttp_oxy_INC_L_medDown);
if isnan(peak_oxy_INC_L_medDown)
    ttl_oxy_INC_L_medDown = NaN ;
    low_oxy_INC_L_medDown = NaN ;
elseif peak_oxy_INC_L_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,9,2) < dcAvg(t0,1,9,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,9,2),[ttp_oxy_INC_L_medDown 25],[dcAvg(t0,1,9,2) dcAvg(t0,1,9,2)]) ;
        ttl_oxy_INC_L_medDown = min(xlin);
        low_oxy_INC_L_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,9,2) < dcAvg(t0,1,9,2)))
        ttl_oxy_INC_L_medDown = t( find(dcAvg(ttp_row:tslut_INC,1,9,2)==min(dcAvg(ttp_row:tslut_INC,1,9,2)))+ttp_row-1);
        low_oxy_INC_L_medDown = min(dcAvg(ttp_row:tslut_INC,1,9,2))-dcAvg(t0,1,9,2);
    end
elseif peak_oxy_INC_L_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,9,2) > dcAvg(t0,1,9,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,9,2),[ttp_oxy_INC_L_medDown 25],[dcAvg(t0,1,9,2) dcAvg(t0,1,9,2)]) ;
        ttl_oxy_INC_L_medDown = min(xlin);
        low_oxy_INC_L_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,9,2) > dcAvg(t0,1,9,2)))
        ttl_oxy_INC_L_medDown = t( find(dcAvg(ttp_row:tslut_INC,1,9,2)==max(dcAvg(ttp_row:tslut_INC,1,9,2)))+ttp_row-1);
        low_oxy_INC_L_medDown = max(dcAvg(ttp_row:tslut_INC,1,9,2))-dcAvg(t0,1,9,2);
    end
end

if abs(max_oxy_INC_L_latUp) > abs(min_oxy_INC_L_latUp)
   peak_oxy_INC_L_latUp = max_oxy_INC_L_latUp ;
   ttp_oxy_INC_L_latUp = t( dcAvg (:,1,10,2) == max(dcAvg(t0:tslut_INC,1,10,2))) ;
elseif abs(max_oxy_INC_L_latUp) < abs(min_oxy_INC_L_latUp)
   peak_oxy_INC_L_latUp = min_oxy_INC_L_latUp ;
   ttp_oxy_INC_L_latUp = t( dcAvg (:,1,10,2) == min(dcAvg(t0:tslut_INC,1,10,2))) ;
elseif isnan(max_oxy_INC_L_latUp)
   ttp_oxy_INC_L_latUp = NaN;
   peak_oxy_INC_L_latUp = NaN ;
elseif abs(max_oxy_INC_L_latUp) == abs(min_oxy_INC_L_latUp)
   ttp_oxy_INC_L_latUp_max = t( dcAvg (:,1,10,2) == max(dcAvg(t0:tslut_INC,1,10,2))) ;
   ttp_oxy_INC_L_latUp_min = t( dcAvg (:,1,10,2) == min(dcAvg(t0:tslut_INC,1,10,2))) ;
   if ttp_oxy_INC_L_latUp_max < ttp_oxy_INC_L_latUp_min
       ttp_oxy_INC_L_latUp = ttp_oxy_INC_L_latUp_max;
       peak_oxy_INC_L_latUp = max_oxy_INC_L_latUp ; 
   elseif ttp_oxy_INC_L_latUp_max > ttp_oxy_INC_L_latUp_min
       ttp_oxy_INC_L_latUp = ttp_oxy_INC_L_latUp_min;
       peak_oxy_INC_L_latUp = min_oxy_INC_L_latUp ;
   end
end
  
ttp_row = find(t==ttp_oxy_INC_L_latUp);
if isnan(peak_oxy_INC_L_latUp)
    ttl_oxy_INC_L_latUp = NaN ;
    low_oxy_INC_L_latUp = NaN ;
elseif peak_oxy_INC_L_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,10,2) < dcAvg(t0,1,10,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,10,2),[ttp_oxy_INC_L_latUp 25],[dcAvg(t0,1,10,2) dcAvg(t0,1,10,2)]) ;
        ttl_oxy_INC_L_latUp = min(xlin);
        low_oxy_INC_L_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,10,2) < dcAvg(t0,1,10,2)))
        ttl_oxy_INC_L_latUp = t( find(dcAvg(ttp_row:tslut_INC,1,10,2)==min(dcAvg(ttp_row:tslut_INC,1,10,2)))+ttp_row-1);
        low_oxy_INC_L_latUp = min(dcAvg(ttp_row:tslut_INC,1,10,2))-dcAvg(t0,1,10,2);
    end
elseif peak_oxy_INC_L_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,10,2) > dcAvg(t0,1,10,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,10,2),[ttp_oxy_INC_L_latUp 25],[dcAvg(t0,1,10,2) dcAvg(t0,1,10,2)]) ;
        ttl_oxy_INC_L_latUp = min(xlin);
        low_oxy_INC_L_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,10,2) > dcAvg(t0,1,10,2)))
        ttl_oxy_INC_L_latUp = t( find(dcAvg(ttp_row:tslut_INC,1,10,2)==max(dcAvg(ttp_row:tslut_INC,1,10,2)))+ttp_row-1);
        low_oxy_INC_L_latUp = max(dcAvg(ttp_row:tslut_INC,1,10,2))-dcAvg(t0,1,10,2);
    end
end

if abs(max_oxy_INC_L_latDown) > abs(min_oxy_INC_L_latDown)
   peak_oxy_INC_L_latDown = max_oxy_INC_L_latDown ;
   ttp_oxy_INC_L_latDown = t( dcAvg (:,1,11,2) == max(dcAvg(t0:tslut_INC,1,11,2))) ;
elseif abs(max_oxy_INC_L_latDown) < abs(min_oxy_INC_L_latDown)
   peak_oxy_INC_L_latDown = min_oxy_INC_L_latDown ;
   ttp_oxy_INC_L_latDown = t( dcAvg (:,1,11,2) == min(dcAvg(t0:tslut_INC,1,11,2))) ;
elseif isnan(max_oxy_INC_L_latDown)
   ttp_oxy_INC_L_latDown = NaN;
   peak_oxy_INC_L_latDown = NaN ;
elseif abs(max_oxy_INC_L_latDown) == abs(min_oxy_INC_L_latDown)
   ttp_oxy_INC_L_latDown_max = t( dcAvg (:,1,11,2) == max(dcAvg(t0:tslut_INC,1,11,2))) ;
   ttp_oxy_INC_L_latDown_min = t( dcAvg (:,1,11,2) == min(dcAvg(t0:tslut_INC,1,11,2))) ;
   if ttp_oxy_INC_L_latDown_max < ttp_oxy_INC_L_latDown_min
       ttp_oxy_INC_L_latDown = ttp_oxy_INC_L_latDown_max;
       peak_oxy_INC_L_latDown = max_oxy_INC_L_latDown ; 
   elseif ttp_oxy_INC_L_latDown_max > ttp_oxy_INC_L_latDown_min
       ttp_oxy_INC_L_latDown = ttp_oxy_INC_L_latDown_min;
       peak_oxy_INC_L_latDown = min_oxy_INC_L_latDown ;
   end
end

ttp_row = find(t==ttp_oxy_INC_L_latDown);
if isnan(peak_oxy_INC_L_latDown)
    ttl_oxy_INC_L_latDown = NaN ;
    low_oxy_INC_L_latDown = NaN ;
elseif peak_oxy_INC_L_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,11,2) < dcAvg(t0,1,11,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,11,2),[ttp_oxy_INC_L_latDown 25],[dcAvg(t0,1,11,2) dcAvg(t0,1,11,2)]) ;
        ttl_oxy_INC_L_latDown = min(xlin);
        low_oxy_INC_L_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,11,2) < dcAvg(t0,1,11,2)))
        ttl_oxy_INC_L_latDown = t( find(dcAvg(ttp_row:tslut_INC,1,11,2)==min(dcAvg(ttp_row:tslut_INC,1,11,2)))+ttp_row-1);
        low_oxy_INC_L_latDown = min(dcAvg(ttp_row:tslut_INC,1,11,2))-dcAvg(t0,1,11,2);
    end
elseif peak_oxy_INC_L_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,1,11,2) > dcAvg(t0,1,11,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,1,11,2),[ttp_oxy_INC_L_latDown 25],[dcAvg(t0,1,11,2) dcAvg(t0,1,11,2)]) ;
        ttl_oxy_INC_L_latDown = min(xlin);
        low_oxy_INC_L_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,1,11,2) > dcAvg(t0,1,11,2)))
        ttl_oxy_INC_L_latDown = t( find(dcAvg(ttp_row:tslut_INC,1,11,2)==max(dcAvg(ttp_row:tslut_INC,1,11,2)))+ttp_row-1);
        low_oxy_INC_L_latDown = max(dcAvg(ttp_row:tslut_INC,1,11,2))-dcAvg(t0,1,11,2);
    end
end

%% INCONGRUENT DEOXY
% Finding peak, time-to-peak (TTP), return to baseline (low) and
% time-to-low (TTL)

% Max. deoxy during INCONGRUENT
max_deoxy_INC_R_latUp = max(dcAvg(t0:tslut_INC,2,2,2))-dcAvg(t0,2,2,2); 
max_deoxy_INC_R_latDown = max(dcAvg(t0:tslut_INC,2,3,2))-dcAvg(t0,2,3,2);
max_deoxy_INC_R_medUp = max(dcAvg(t0:tslut_INC,2,4,2))-dcAvg(t0,2,4,2);
max_deoxy_INC_R_medDown = max(dcAvg(t0:tslut_INC,2,5,2))-dcAvg(t0,2,3,2); 
max_deoxy_INC_L_medUp = max(dcAvg(t0:tslut_INC,2,8,2))-dcAvg(t0,2,8,2);
max_deoxy_INC_L_medDown = max(dcAvg(t0:tslut_INC,2,9,2))-dcAvg(t0,2,9,2);
max_deoxy_INC_L_latUp = max(dcAvg(t0:tslut_INC,2,10,2))-dcAvg(t0,2,10,2);
max_deoxy_INC_L_latDown = max(dcAvg(t0:tslut_INC,2,11,2))-dcAvg(t0,2,11,2);

% Min. deoxy during INCONGRUENT
min_deoxy_INC_R_latUp = min(dcAvg(t0:tslut_INC,2,2,2))-dcAvg(t0,2,2,2);
min_deoxy_INC_R_latDown = min(dcAvg(t0:tslut_INC,2,3,2))-dcAvg(t0,2,3,2);
min_deoxy_INC_R_medUp = min(dcAvg(t0:tslut_INC,2,4,2))-dcAvg(t0,2,4,2);
min_deoxy_INC_R_medDown = min(dcAvg(t0:tslut_INC,2,5,2))-dcAvg(t0,2,5,2);
min_deoxy_INC_L_medUp = min(dcAvg(t0:tslut_INC,2,8,2))-dcAvg(t0,2,8,2);
min_deoxy_INC_L_medDown = min(dcAvg(t0:tslut_INC,2,9,2))-dcAvg(t0,2,9,2);
min_deoxy_INC_L_latUp = min(dcAvg(t0:tslut_INC,2,10,2))-dcAvg(t0,2,10,2);
min_deoxy_INC_L_latDown = min(dcAvg(t0:tslut_INC,2,11,2))-dcAvg(t0,2,11,2);

% Finding peak (min. or max.) and time-to-peak during INCONGRUENT 
% If identical, then use the first one
if abs(max_deoxy_INC_R_latUp) > abs(min_deoxy_INC_R_latUp)
   peak_deoxy_INC_R_latUp = max_deoxy_INC_R_latUp ;
   ttp_deoxy_INC_R_latUp = t( dcAvg (:,2,2,2) == max(dcAvg(t0:tslut_INC,2,2,2))) ;
elseif abs(max_deoxy_INC_R_latUp) < abs(min_deoxy_INC_R_latUp)
   peak_deoxy_INC_R_latUp = min_deoxy_INC_R_latUp ;
   ttp_deoxy_INC_R_latUp = t( dcAvg (:,2,2,2) == min(dcAvg(t0:tslut_INC,2,2,2))) ;
elseif isnan(max_deoxy_INC_R_latUp)
   ttp_deoxy_INC_R_latUp = NaN;
   peak_deoxy_INC_R_latUp = NaN ;
elseif abs(max_deoxy_INC_R_latUp) == abs(min_deoxy_INC_R_latUp)
   ttp_deoxy_INC_R_latUp_max = t( dcAvg (:,2,2,2) == max(dcAvg(t0:tslut_INC,2,2,2))) ;
   ttp_deoxy_INC_R_latUp_min = t( dcAvg (:,2,2,2) == min(dcAvg(t0:tslut_INC,2,2,2))) ;
   if ttp_deoxy_INC_R_latUp_max < ttp_deoxy_INC_R_latUp_min
       ttp_deoxy_INC_R_latUp = ttp_deoxy_INC_R_latUp_max;
       peak_deoxy_INC_R_latUp = max_deoxy_INC_R_latUp ; 
   elseif ttp_deoxy_INC_R_latUp_max > ttp_deoxy_INC_R_latUp_min
       ttp_deoxy_INC_R_latUp = ttp_deoxy_INC_R_latUp_min;
       peak_deoxy_INC_R_latUp = min_deoxy_INC_R_latUp ;
   end
end
  
ttp_row = find(t==ttp_deoxy_INC_R_latUp);
if isnan(peak_deoxy_INC_R_latUp)
    ttl_deoxy_INC_R_latUp = NaN ;
    low_deoxy_INC_R_latUp = NaN ;
elseif peak_deoxy_INC_R_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,2,2) < dcAvg(t0,2,2,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,2,2),[ttp_deoxy_INC_R_latUp 25],[dcAvg(t0,2,2,2) dcAvg(t0,2,2,2)]) ;
        ttl_deoxy_INC_R_latUp = min(xlin);
        low_deoxy_INC_R_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,2,2) < dcAvg(t0,2,2,2)))
        ttl_deoxy_INC_R_latUp = t( find(dcAvg(ttp_row:tslut_INC,2,2,2)==min(dcAvg(ttp_row:tslut_INC,2,2,2)))+ttp_row-1);
        low_deoxy_INC_R_latUp = min(dcAvg(ttp_row:tslut_INC,2,2,2))-dcAvg(t0,2,2,2);
    end
elseif peak_deoxy_INC_R_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,2,2) > dcAvg(t0,2,2,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,2,2),[ttp_deoxy_INC_R_latUp 25],[dcAvg(t0,2,2,2) dcAvg(t0,2,2,2)]) ;
        ttl_deoxy_INC_R_latUp = min(xlin);
        low_deoxy_INC_R_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,2,2) > dcAvg(t0,2,2,2)))
        ttl_deoxy_INC_R_latUp = t( find(dcAvg(ttp_row:tslut_INC,2,2,2)==max(dcAvg(ttp_row:tslut_INC,2,2,2)))+ttp_row-1);
        low_deoxy_INC_R_latUp = max(dcAvg(ttp_row:tslut_INC,2,2,2))-dcAvg(t0,2,2,2);
    end
end

if abs(max_deoxy_INC_R_latDown) > abs(min_deoxy_INC_R_latDown)
   peak_deoxy_INC_R_latDown = max_deoxy_INC_R_latDown ;
   ttp_deoxy_INC_R_latDown = t( dcAvg (:,2,3,2) == max(dcAvg(t0:tslut_INC,2,3,2))) ;
elseif abs(max_deoxy_INC_R_latDown) < abs(min_deoxy_INC_R_latDown)
   peak_deoxy_INC_R_latDown = min_deoxy_INC_R_latDown ;
   ttp_deoxy_INC_R_latDown = t( dcAvg (:,2,3,2) == min(dcAvg(t0:tslut_INC,2,3,2))) ;
elseif isnan(max_deoxy_INC_R_latDown)
   peak_deoxy_INC_R_latDown = NaN;
   ttp_deoxy_INC_R_latDown = NaN;
elseif abs(max_deoxy_INC_R_latDown) == abs(min_deoxy_INC_R_latDown)
   ttp_deoxy_INC_R_latDown_max = t( dcAvg (:,2,3,2) == max(dcAvg(t0:tslut_INC,2,3,2))) ;
   ttp_deoxy_INC_R_latDown_min = t( dcAvg (:,2,3,2) == min(dcAvg(t0:tslut_INC,2,3,2))) ;
    if ttp_deoxy_INC_R_latDown_max < ttp_deoxy_INC_R_latDown_min
       ttp_deoxy_INC_R_latDown = ttp_deoxy_INC_R_latDown_max;
       peak_deoxy_INC_R_latDown = max_deoxy_INC_R_latDown ; 
    elseif ttp_deoxy_INC_R_latDown_max > ttp_deoxy_INC_R_latDown_min
       ttp_deoxy_INC_R_latDown = ttp_deoxy_INC_R_latDown_min;
       peak_deoxy_INC_R_latDown = min_deoxy_INC_R_latDown ;
    end
end

ttp_row = find(t==ttp_deoxy_INC_R_latDown);
if isnan(peak_deoxy_INC_R_latDown)
    ttl_deoxy_INC_R_latDown = NaN ;
    low_deoxy_INC_R_latDown = NaN ;
elseif peak_deoxy_INC_R_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,3,2) < dcAvg(t0,2,3,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,3,2),[ttp_deoxy_INC_R_latDown 25],[dcAvg(t0,2,3,2) dcAvg(t0,2,3,2)]) ;
        ttl_deoxy_INC_R_latDown = min(xlin);
        low_deoxy_INC_R_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,3,2) < dcAvg(t0,2,3,2)))
        ttl_deoxy_INC_R_latDown = t( find(dcAvg(ttp_row:tslut_INC,2,3,2)==min(dcAvg(ttp_row:tslut_INC,2,3,2)))+ttp_row-1);
        low_deoxy_INC_R_latDown = min(dcAvg(ttp_row:tslut_INC,2,3,2))-dcAvg(t0,2,3,2);
    end
elseif peak_deoxy_INC_R_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,3,2) > dcAvg(t0,2,3,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,3,2),[ttp_deoxy_INC_R_latDown 25],[dcAvg(t0,2,3,2) dcAvg(t0,2,3,2)]) ;
        ttl_deoxy_INC_R_latDown = min(xlin);
        low_deoxy_INC_R_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,3,2) > dcAvg(t0,2,3,2)))
        ttl_deoxy_INC_R_latDown = t( find(dcAvg(ttp_row:tslut_INC,2,3,2)==max(dcAvg(ttp_row:tslut_INC,2,3,2)))+ttp_row-1);
        low_deoxy_INC_R_latDown = max(dcAvg(ttp_row:tslut_INC,2,3,2))-dcAvg(t0,2,3,2);
    end
end

if abs(max_deoxy_INC_R_medUp) > abs(min_deoxy_INC_R_medUp)
   peak_deoxy_INC_R_medUp = max_deoxy_INC_R_medUp ;
   ttp_deoxy_INC_R_medUp = t( dcAvg (:,2,4,2) == max(dcAvg(t0:tslut_INC,2,4,2))) ;
elseif abs(max_deoxy_INC_R_medUp) < abs(min_deoxy_INC_R_medUp)
   peak_deoxy_INC_R_medUp = min_deoxy_INC_R_medUp ;
   ttp_deoxy_INC_R_medUp = t( dcAvg (:,2,4,2) == min(dcAvg(t0:tslut_INC,2,4,2))) ;
elseif isnan(max_deoxy_INC_R_medUp)
   ttp_deoxy_INC_R_medUp = NaN;
   peak_deoxy_INC_R_medUp = NaN ;
elseif abs(max_deoxy_INC_R_medUp) == abs(min_deoxy_INC_R_medUp)
   ttp_deoxy_INC_R_medUp_max = t( dcAvg (:,2,4,2) == max(dcAvg(t0:tslut_INC,2,4,2))) ;
   ttp_deoxy_INC_R_medUp_min = t( dcAvg (:,2,4,2) == min(dcAvg(t0:tslut_INC,2,4,2))) ;
   if ttp_deoxy_INC_R_medUp_max < ttp_deoxy_INC_R_medUp_min
       ttp_deoxy_INC_R_medUp = ttp_deoxy_INC_R_medUp_max;
       peak_deoxy_INC_R_medUp = max_deoxy_INC_R_medUp ; 
   elseif ttp_deoxy_INC_R_medUp_max > ttp_deoxy_INC_R_medUp_min
       ttp_deoxy_INC_R_medUp = ttp_deoxy_INC_R_medUp_min;
       peak_deoxy_INC_R_medUp = min_deoxy_INC_R_medUp ;
   end
end

ttp_row = find(t==ttp_deoxy_INC_R_medUp);
if isnan(peak_deoxy_INC_R_medUp)
    ttl_deoxy_INC_R_medUp = NaN ;
    low_deoxy_INC_R_medUp = NaN ;
elseif peak_deoxy_INC_R_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,4,2) < dcAvg(t0,2,4,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,4,2),[ttp_deoxy_INC_R_medUp 25],[dcAvg(t0,2,4,2) dcAvg(t0,2,4,2)]) ;
        ttl_deoxy_INC_R_medUp = min(xlin);
        low_deoxy_INC_R_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,4,2) < dcAvg(t0,2,4,2)))
        ttl_deoxy_INC_R_medUp = t( find(dcAvg(ttp_row:tslut_INC,2,4,2)==min(dcAvg(ttp_row:tslut_INC,2,4,2)))+ttp_row-1);
        low_deoxy_INC_R_medUp = min(dcAvg(ttp_row:tslut_INC,2,4,2))-dcAvg(t0,2,4,2);
    end
elseif peak_deoxy_INC_R_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,4,2) > dcAvg(t0,2,4,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,4,2),[ttp_deoxy_INC_R_medUp 25],[dcAvg(t0,2,4,2) dcAvg(t0,2,4,2)]) ;
        ttl_deoxy_INC_R_medUp = min(xlin);
        low_deoxy_INC_R_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,4,2) > dcAvg(t0,2,4,2)))
        ttl_deoxy_INC_R_medUp = t( find(dcAvg(ttp_row:tslut_INC,2,4,2)==max(dcAvg(ttp_row:tslut_INC,2,4,2)))+ttp_row-1);
        low_deoxy_INC_R_medUp = max(dcAvg(ttp_row:tslut_INC,2,4,2))-dcAvg(t0,2,4,2);
    end
end

if abs(max_deoxy_INC_R_medDown) > abs(min_deoxy_INC_R_medDown)
   peak_deoxy_INC_R_medDown = max_deoxy_INC_R_medDown ;
   ttp_deoxy_INC_R_medDown = t( dcAvg (:,2,5,2) == max(dcAvg(t0:tslut_INC,2,5,2))) ;
elseif abs(max_deoxy_INC_R_medDown) < abs(min_deoxy_INC_R_medDown)
   peak_deoxy_INC_R_medDown = min_deoxy_INC_R_medDown ;
   ttp_deoxy_INC_R_medDown = t( dcAvg (:,2,5,2) == min(dcAvg(t0:tslut_INC,2,5,2))) ;
elseif isnan(max_deoxy_INC_R_medDown)
   ttp_deoxy_INC_R_medDown = NaN;
   peak_deoxy_INC_R_medDown = NaN ;
elseif abs(max_deoxy_INC_R_medDown) == abs(min_deoxy_INC_R_medDown)
   ttp_deoxy_INC_R_medDown_max = t( dcAvg (:,2,5,2) == max(dcAvg(t0:tslut_INC,2,5,2))) ;
   ttp_deoxy_INC_R_medDown_min = t( dcAvg (:,2,5,2) == min(dcAvg(t0:tslut_INC,2,5,2))) ;
   if ttp_deoxy_INC_R_medDown_max < ttp_deoxy_INC_R_medDown_min
       ttp_deoxy_INC_R_medDown = ttp_deoxy_INC_R_medDown_max;
       peak_deoxy_INC_R_medDown = max_deoxy_INC_R_medDown ; 
   elseif ttp_deoxy_INC_R_medDown_max > ttp_deoxy_INC_R_medDown_min
       ttp_deoxy_INC_R_medDown = ttp_deoxy_INC_R_medDown_min;
       peak_deoxy_INC_R_medDown = min_deoxy_INC_R_medDown ;
   end
end

ttp_row = find(t==ttp_deoxy_INC_R_medDown);
if isnan(peak_deoxy_INC_R_medDown)
    ttl_deoxy_INC_R_medDown = NaN ;
    low_deoxy_INC_R_medDown = NaN ;
elseif peak_deoxy_INC_R_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,5,2) < dcAvg(t0,2,5,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,5,2),[ttp_deoxy_INC_R_medDown 25],[dcAvg(t0,2,5,2) dcAvg(t0,2,5,2)]) ;
        ttl_deoxy_INC_R_medDown = min(xlin);
        low_deoxy_INC_R_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,5,2) < dcAvg(t0,2,5,2)))
        ttl_deoxy_INC_R_medDown = t( find(dcAvg(ttp_row:tslut_INC,2,5,2)==min(dcAvg(ttp_row:tslut_INC,2,5,2)))+ttp_row-1);
        low_deoxy_INC_R_medDown = min(dcAvg(ttp_row:tslut_INC,2,5,2))-dcAvg(t0,2,5,2);
    end
elseif peak_deoxy_INC_R_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,5,2) > dcAvg(t0,2,5,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,5,2),[ttp_deoxy_INC_R_medDown 25],[dcAvg(t0,2,5,2) dcAvg(t0,2,5,2)]) ;
        ttl_deoxy_INC_R_medDown = min(xlin);
        low_deoxy_INC_R_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,5,2) > dcAvg(t0,2,5,2)))
        ttl_deoxy_INC_R_medDown = t( find(dcAvg(ttp_row:tslut_INC,2,5,2)==max(dcAvg(ttp_row:tslut_INC,2,5,2)))+ttp_row-1);
        low_deoxy_INC_R_medDown = max(dcAvg(ttp_row:tslut_INC,2,5,2))-dcAvg(t0,2,5,2);
    end
end

if abs(max_deoxy_INC_L_medUp) > abs(min_deoxy_INC_L_medUp)
   peak_deoxy_INC_L_medUp = max_deoxy_INC_L_medUp ;
   ttp_deoxy_INC_L_medUp = t( dcAvg (:,2,8,2) == max(dcAvg(t0:tslut_INC,2,8,2))) ;
elseif abs(max_deoxy_INC_L_medUp) < abs(min_deoxy_INC_L_medUp)
   peak_deoxy_INC_L_medUp = min_deoxy_INC_L_medUp ;
   ttp_deoxy_INC_L_medUp = t( dcAvg (:,2,8,2) == min(dcAvg(t0:tslut_INC,2,8,2))) ;
elseif isnan(max_deoxy_INC_L_medUp)
   ttp_deoxy_INC_L_medUp = NaN;
   peak_deoxy_INC_L_medUp = NaN ;
elseif abs(max_deoxy_INC_L_medUp) == abs(min_deoxy_INC_L_medUp)
   ttp_deoxy_INC_L_medUp_max = t( dcAvg (:,2,8,2) == max(dcAvg(t0:tslut_INC,2,8,2))) ;
   ttp_deoxy_INC_L_medUp_min = t( dcAvg (:,2,8,2) == min(dcAvg(t0:tslut_INC,2,8,2))) ;
   if ttp_deoxy_INC_L_medUp_max < ttp_deoxy_INC_L_medUp_min
       ttp_deoxy_INC_L_medUp = ttp_deoxy_INC_L_medUp_max;
       peak_deoxy_INC_L_medUp = max_deoxy_INC_L_medUp ; 
   elseif ttp_deoxy_INC_L_medUp_max > ttp_deoxy_INC_L_medUp_min
       ttp_deoxy_INC_L_medUp = ttp_deoxy_INC_L_medUp_min;
       peak_deoxy_INC_L_medUp = min_deoxy_INC_L_medUp ;
   end
end

ttp_row = find(t==ttp_deoxy_INC_L_medUp);
if isnan(peak_deoxy_INC_L_medUp)
    ttl_deoxy_INC_L_medUp = NaN ;
    low_deoxy_INC_L_medUp = NaN ;
elseif peak_deoxy_INC_L_medUp > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,8,2) < dcAvg(t0,2,8,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,8,2),[ttp_deoxy_INC_L_medUp 25],[dcAvg(t0,2,8,2) dcAvg(t0,2,8,2)]) ;
        ttl_deoxy_INC_L_medUp = min(xlin);
        low_deoxy_INC_L_medUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,8,2) < dcAvg(t0,2,8,2)))
        ttl_deoxy_INC_L_medUp = t( find(dcAvg(ttp_row:tslut_INC,2,8,2)==min(dcAvg(ttp_row:tslut_INC,2,8,2)))+ttp_row-1);
        low_deoxy_INC_L_medUp = min(dcAvg(ttp_row:tslut_INC,2,8,2))-dcAvg(t0,2,8,2);
    end
elseif peak_deoxy_INC_L_medUp < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,8,2) > dcAvg(t0,2,8,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,8,2),[ttp_deoxy_INC_L_medUp 25],[dcAvg(t0,2,8,2) dcAvg(t0,2,8,2)]) ;
        ttl_deoxy_INC_L_medUp = min(xlin);
        low_deoxy_INC_L_medUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,8,2) > dcAvg(t0,2,8,2)))
        ttl_deoxy_INC_L_medUp = t( find(dcAvg(ttp_row:tslut_INC,2,8,2)==max(dcAvg(ttp_row:tslut_INC,2,8,2)))+ttp_row-1);
        low_deoxy_INC_L_medUp = max(dcAvg(ttp_row:tslut_INC,2,8,2))-dcAvg(t0,2,8,2);
    end
end

if abs(max_deoxy_INC_L_medDown) > abs(min_deoxy_INC_L_medDown)
   peak_deoxy_INC_L_medDown = max_deoxy_INC_L_medDown ;
   ttp_deoxy_INC_L_medDown = t( dcAvg (:,2,9,2) == max(dcAvg(t0:tslut_INC,2,9,2))) ;
elseif abs(max_deoxy_INC_L_medDown) < abs(min_deoxy_INC_L_medDown)
   peak_deoxy_INC_L_medDown = min_deoxy_INC_L_medDown ;
   ttp_deoxy_INC_L_medDown = t( dcAvg (:,2,9,2) == min(dcAvg(t0:tslut_INC,2,9,2))) ;
elseif isnan(max_deoxy_INC_L_medDown)
   ttp_deoxy_INC_L_medDown = NaN;
   peak_deoxy_INC_L_medDown = NaN ;
elseif abs(max_deoxy_INC_L_medDown) == abs(min_deoxy_INC_L_medDown)
   ttp_deoxy_INC_L_medDown_max = t( dcAvg (:,2,9,2) == max(dcAvg(t0:tslut_INC,2,9,2))) ;
   ttp_deoxy_INC_L_medDown_min = t( dcAvg (:,2,9,2) == min(dcAvg(t0:tslut_INC,2,9,2))) ;
   if ttp_deoxy_INC_L_medDown_max < ttp_deoxy_INC_L_medDown_min
       ttp_deoxy_INC_L_medDown = ttp_deoxy_INC_L_medDown_max;
       peak_deoxy_INC_L_medDown = max_deoxy_INC_L_medDown ; 
   elseif ttp_deoxy_INC_L_medDown_max > ttp_deoxy_INC_L_medDown_min
       ttp_deoxy_INC_L_medDown = ttp_deoxy_INC_L_medDown_min;
       peak_deoxy_INC_L_medDown = min_deoxy_INC_L_medDown ;
   end
end

ttp_row = find(t==ttp_deoxy_INC_L_medDown);
if isnan(peak_deoxy_INC_L_medDown)
    ttl_deoxy_INC_L_medDown = NaN ;
    low_deoxy_INC_L_medDown = NaN ;
elseif peak_deoxy_INC_L_medDown > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,9,2) < dcAvg(t0,2,9,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,9,2),[ttp_deoxy_INC_L_medDown 25],[dcAvg(t0,2,9,2) dcAvg(t0,2,9,2)]) ;
        ttl_deoxy_INC_L_medDown = min(xlin);
        low_deoxy_INC_L_medDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,9,2) < dcAvg(t0,2,9,2)))
        ttl_deoxy_INC_L_medDown = t( find(dcAvg(ttp_row:tslut_INC,2,9,2)==min(dcAvg(ttp_row:tslut_INC,2,9,2)))+ttp_row-1);
        low_deoxy_INC_L_medDown = min(dcAvg(ttp_row:tslut_INC,2,9,2))-dcAvg(t0,2,9,2);
    end
elseif peak_deoxy_INC_L_medDown < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,9,2) > dcAvg(t0,2,9,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,9,2),[ttp_deoxy_INC_L_medDown 25],[dcAvg(t0,2,9,2) dcAvg(t0,2,9,2)]) ;
        ttl_deoxy_INC_L_medDown = min(xlin);
        low_deoxy_INC_L_medDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,9,2) > dcAvg(t0,2,9,2)))
        ttl_deoxy_INC_L_medDown = t( find(dcAvg(ttp_row:tslut_INC,2,9,2)==max(dcAvg(ttp_row:tslut_INC,2,9,2)))+ttp_row-1);
        low_deoxy_INC_L_medDown = max(dcAvg(ttp_row:tslut_INC,2,9,2))-dcAvg(t0,2,9,2);
    end
end

if abs(max_deoxy_INC_L_latUp) > abs(min_deoxy_INC_L_latUp)
   peak_deoxy_INC_L_latUp = max_deoxy_INC_L_latUp ;
   ttp_deoxy_INC_L_latUp = t( dcAvg (:,2,10,2) == max(dcAvg(t0:tslut_INC,2,10,2))) ;
elseif abs(max_deoxy_INC_L_latUp) < abs(min_deoxy_INC_L_latUp)
   peak_deoxy_INC_L_latUp = min_deoxy_INC_L_latUp ;
   ttp_deoxy_INC_L_latUp = t( dcAvg (:,2,10,2) == min(dcAvg(t0:tslut_INC,2,10,2))) ;
elseif isnan(max_deoxy_INC_L_latUp)
   ttp_deoxy_INC_L_latUp = NaN;
   peak_deoxy_INC_L_latUp = NaN ;
elseif abs(max_deoxy_INC_L_latUp) == abs(min_deoxy_INC_L_latUp)
   ttp_deoxy_INC_L_latUp_max = t( dcAvg (:,2,10,2) == max(dcAvg(t0:tslut_INC,2,10,2))) ;
   ttp_deoxy_INC_L_latUp_min = t( dcAvg (:,2,10,2) == min(dcAvg(t0:tslut_INC,2,10,2))) ;
   if ttp_deoxy_INC_L_latUp_max < ttp_deoxy_INC_L_latUp_min
       ttp_deoxy_INC_L_latUp = ttp_deoxy_INC_L_latUp_max;
       peak_deoxy_INC_L_latUp = max_deoxy_INC_L_latUp ; 
   elseif ttp_deoxy_INC_L_latUp_max > ttp_deoxy_INC_L_latUp_min
       ttp_deoxy_INC_L_latUp = ttp_deoxy_INC_L_latUp_min;
       peak_deoxy_INC_L_latUp = min_deoxy_INC_L_latUp ;
   end
end

ttp_row = find(t==ttp_deoxy_INC_L_latUp);
if isnan(peak_deoxy_INC_L_latUp)
    ttl_deoxy_INC_L_latUp = NaN ;
    low_deoxy_INC_L_latUp = NaN ;
elseif peak_deoxy_INC_L_latUp > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,10,2) < dcAvg(t0,2,10,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,10,2),[ttp_deoxy_INC_L_latUp 25],[dcAvg(t0,2,10,2) dcAvg(t0,2,10,2)]) ;
        ttl_deoxy_INC_L_latUp = min(xlin);
        low_deoxy_INC_L_latUp = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,10,2) < dcAvg(t0,2,10,2)))
        ttl_deoxy_INC_L_latUp = t( find(dcAvg(ttp_row:tslut_INC,2,10,2)==min(dcAvg(ttp_row:tslut_INC,2,10,2)))+ttp_row-1);
        low_deoxy_INC_L_latUp = min(dcAvg(ttp_row:tslut_INC,2,10,2))-dcAvg(t0,2,10,2);
    end
elseif peak_deoxy_INC_L_latUp < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,10,2) > dcAvg(t0,2,10,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,10,2),[ttp_deoxy_INC_L_latUp 25],[dcAvg(t0,2,10,2) dcAvg(t0,2,10,2)]) ;
        ttl_deoxy_INC_L_latUp = min(xlin);
        low_deoxy_INC_L_latUp = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,10,2) > dcAvg(t0,2,10,2)))
        ttl_deoxy_INC_L_latUp = t( find(dcAvg(ttp_row:tslut_INC,2,10,2)==max(dcAvg(ttp_row:tslut_INC,2,10,2)))+ttp_row-1);
        low_deoxy_INC_L_latUp = max(dcAvg(ttp_row:tslut_INC,2,10,2))-dcAvg(t0,2,10,2);
    end
end
  
if abs(max_deoxy_INC_L_latDown) > abs(min_deoxy_INC_L_latDown)
   peak_deoxy_INC_L_latDown = max_deoxy_INC_L_latDown ;
   ttp_deoxy_INC_L_latDown = t( dcAvg (:,2,11,2) == max(dcAvg(t0:tslut_INC,2,11,2))) ;
elseif abs(max_deoxy_INC_L_latDown) < abs(min_deoxy_INC_L_latDown)
   peak_deoxy_INC_L_latDown = min_deoxy_INC_L_latDown ;
   ttp_deoxy_INC_L_latDown = t( dcAvg (:,2,11,2) == min(dcAvg(t0:tslut_INC,2,11,2))) ;
elseif isnan(max_deoxy_INC_L_latDown)
   ttp_deoxy_INC_L_latDown = NaN;
   peak_deoxy_INC_L_latDown = NaN ;
elseif abs(max_deoxy_INC_L_latDown) == abs(min_deoxy_INC_L_latDown)
   ttp_deoxy_INC_L_latDown_max = t( dcAvg (:,2,11,2) == max(dcAvg(t0:tslut_INC,2,11,2))) ;
   ttp_deoxy_INC_L_latDown_min = t( dcAvg (:,2,11,2) == min(dcAvg(t0:tslut_INC,2,11,2))) ;
   if ttp_deoxy_INC_L_latDown_max < ttp_deoxy_INC_L_latDown_min
       ttp_deoxy_INC_L_latDown = ttp_deoxy_INC_L_latDown_max;
       peak_deoxy_INC_L_latDown = max_deoxy_INC_L_latDown ; 
   elseif ttp_deoxy_INC_L_latDown_max > ttp_deoxy_INC_L_latDown_min
       ttp_deoxy_INC_L_latDown = ttp_deoxy_INC_L_latDown_min;
       peak_deoxy_INC_L_latDown = min_deoxy_INC_L_latDown ;
   end
end

ttp_row = find(t==ttp_deoxy_INC_L_latDown);
if isnan(peak_deoxy_INC_L_latDown)
    ttl_deoxy_INC_L_latDown = NaN ;
    low_deoxy_INC_L_latDown = NaN ;
elseif peak_deoxy_INC_L_latDown > 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,11,2) < dcAvg(t0,2,11,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,11,2),[ttp_deoxy_INC_L_latDown 25],[dcAvg(t0,2,11,2) dcAvg(t0,2,11,2)]) ;
        ttl_deoxy_INC_L_latDown = min(xlin);
        low_deoxy_INC_L_latDown = 0;
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,11,2) < dcAvg(t0,2,11,2)))
        ttl_deoxy_INC_L_latDown = t( find(dcAvg(ttp_row:tslut_INC,2,11,2)==min(dcAvg(ttp_row:tslut_INC,2,11,2)))+ttp_row-1);
        low_deoxy_INC_L_latDown = min(dcAvg(ttp_row:tslut_INC,2,11,2))-dcAvg(t0,2,11,2);
    end
elseif peak_deoxy_INC_L_latDown < 0
    if logical(any(dcAvg(ttp_row:tslut_INC,2,11,2) > dcAvg(t0,2,11,2)))
        [xlin,~] = intersections(t(ttp_row:tslut_INC),dcAvg(ttp_row:tslut_INC,2,11,2),[ttp_deoxy_INC_L_latDown 25],[dcAvg(t0,2,11,2) dcAvg(t0,2,11,2)]) ;
        ttl_deoxy_INC_L_latDown = min(xlin);
        low_deoxy_INC_L_latDown = 0;    
    elseif ~logical(any(dcAvg(ttp_row:tslut_INC,2,11,2) > dcAvg(t0,2,11,2)))
        ttl_deoxy_INC_L_latDown = t( find(dcAvg(ttp_row:tslut_INC,2,11,2)==max(dcAvg(ttp_row:tslut_INC,2,11,2)))+ttp_row-1);
        low_deoxy_INC_L_latDown = max(dcAvg(ttp_row:tslut_INC,2,11,2))-dcAvg(t0,2,11,2);
    end
end

%% Average coefficient to peak

Slope1_oxy_NEU_R_latUp = peak_oxy_NEU_R_latUp / ttp_oxy_NEU_R_latUp;
Slope1_oxy_NEU_R_latDown = peak_oxy_NEU_R_latDown / ttp_oxy_NEU_R_latDown;
Slope1_oxy_NEU_R_medUp = peak_oxy_NEU_R_medUp / ttp_oxy_NEU_R_medUp;
Slope1_oxy_NEU_R_medDown = peak_oxy_NEU_R_medDown / ttp_oxy_NEU_R_medDown;
Slope1_oxy_NEU_L_latUp = peak_oxy_NEU_L_latUp / ttp_oxy_NEU_L_latUp;
Slope1_oxy_NEU_L_latDown = peak_oxy_NEU_L_latDown / ttp_oxy_NEU_L_latDown;
Slope1_oxy_NEU_L_medUp = peak_oxy_NEU_L_medUp / ttp_oxy_NEU_L_medUp;
Slope1_oxy_NEU_L_medDown = peak_oxy_NEU_L_medDown / ttp_oxy_NEU_L_medDown;

Slope1_deoxy_NEU_R_latUp = peak_deoxy_NEU_R_latUp / ttp_deoxy_NEU_R_latUp;
Slope1_deoxy_NEU_R_latDown = peak_deoxy_NEU_R_latDown / ttp_deoxy_NEU_R_latDown;
Slope1_deoxy_NEU_R_medUp = peak_deoxy_NEU_R_medUp / ttp_deoxy_NEU_R_medUp;
Slope1_deoxy_NEU_R_medDown = peak_deoxy_NEU_R_medDown / ttp_deoxy_NEU_R_medDown;
Slope1_deoxy_NEU_L_latUp = peak_deoxy_NEU_L_latUp / ttp_deoxy_NEU_L_latUp;
Slope1_deoxy_NEU_L_latDown = peak_deoxy_NEU_L_latDown / ttp_deoxy_NEU_L_latDown;
Slope1_deoxy_NEU_L_medUp = peak_deoxy_NEU_L_medUp / ttp_deoxy_NEU_L_medUp;
Slope1_deoxy_NEU_L_medDown = peak_deoxy_NEU_L_medDown / ttp_deoxy_NEU_L_medDown;

Slope1_oxy_CON_R_latUp = peak_oxy_CON_R_latUp / ttp_oxy_CON_R_latUp;
Slope1_oxy_CON_R_latDown = peak_oxy_CON_R_latDown / ttp_oxy_CON_R_latDown;
Slope1_oxy_CON_R_medUp = peak_oxy_CON_R_medUp / ttp_oxy_CON_R_medUp;
Slope1_oxy_CON_R_medDown = peak_oxy_CON_R_medDown / ttp_oxy_CON_R_medDown;
Slope1_oxy_CON_L_latUp = peak_oxy_CON_L_latUp / ttp_oxy_CON_L_latUp;
Slope1_oxy_CON_L_latDown = peak_oxy_CON_L_latDown / ttp_oxy_CON_L_latDown;
Slope1_oxy_CON_L_medUp = peak_oxy_CON_L_medUp / ttp_oxy_CON_L_medUp;
Slope1_oxy_CON_L_medDown = peak_oxy_CON_L_medDown / ttp_oxy_CON_L_medDown;

Slope1_deoxy_CON_R_latUp = peak_deoxy_CON_R_latUp / ttp_deoxy_CON_R_latUp;
Slope1_deoxy_CON_R_latDown = peak_deoxy_CON_R_latDown / ttp_deoxy_CON_R_latDown;
Slope1_deoxy_CON_R_medUp = peak_deoxy_CON_R_medUp / ttp_deoxy_CON_R_medUp;
Slope1_deoxy_CON_R_medDown = peak_deoxy_CON_R_medDown / ttp_deoxy_CON_R_medDown;
Slope1_deoxy_CON_L_latUp = peak_deoxy_CON_L_latUp / ttp_deoxy_CON_L_latUp;
Slope1_deoxy_CON_L_latDown = peak_deoxy_CON_L_latDown / ttp_deoxy_CON_L_latDown;
Slope1_deoxy_CON_L_medUp = peak_deoxy_CON_L_medUp / ttp_deoxy_CON_L_medUp;
Slope1_deoxy_CON_L_medDown = peak_deoxy_CON_L_medDown / ttp_deoxy_CON_L_medDown;

Slope1_oxy_INC_R_latUp = peak_oxy_INC_R_latUp / ttp_oxy_INC_R_latUp;
Slope1_oxy_INC_R_latDown = peak_oxy_INC_R_latDown / ttp_oxy_INC_R_latDown;
Slope1_oxy_INC_R_medUp = peak_oxy_INC_R_medUp / ttp_oxy_INC_R_medUp;
Slope1_oxy_INC_R_medDown = peak_oxy_INC_R_medDown / ttp_oxy_INC_R_medDown;
Slope1_oxy_INC_L_latUp = peak_oxy_INC_L_latUp / ttp_oxy_INC_L_latUp;
Slope1_oxy_INC_L_latDown = peak_oxy_INC_L_latDown / ttp_oxy_INC_L_latDown;
Slope1_oxy_INC_L_medUp = peak_oxy_INC_L_medUp / ttp_oxy_INC_L_medUp;
Slope1_oxy_INC_L_medDown = peak_oxy_INC_L_medDown / ttp_oxy_INC_L_medDown;

Slope1_deoxy_INC_R_latUp = peak_deoxy_INC_R_latUp / ttp_deoxy_INC_R_latUp;
Slope1_deoxy_INC_R_latDown = peak_deoxy_INC_R_latDown / ttp_deoxy_INC_R_latDown;
Slope1_deoxy_INC_R_medUp = peak_deoxy_INC_R_medUp / ttp_deoxy_INC_R_medUp;
Slope1_deoxy_INC_R_medDown = peak_deoxy_INC_R_medDown / ttp_deoxy_INC_R_medDown;
Slope1_deoxy_INC_L_latUp = peak_deoxy_INC_L_latUp / ttp_deoxy_INC_L_latUp;
Slope1_deoxy_INC_L_latDown = peak_deoxy_INC_L_latDown / ttp_deoxy_INC_L_latDown;
Slope1_deoxy_INC_L_medUp = peak_deoxy_INC_L_medUp / ttp_deoxy_INC_L_medUp;
Slope1_deoxy_INC_L_medDown = peak_deoxy_INC_L_medDown / ttp_deoxy_INC_L_medDown;

%% Average slope from peak to baseline return

Slope2_oxy_NEU_R_latUp = (low_oxy_NEU_R_latUp - peak_oxy_NEU_R_latUp)/ttl_oxy_NEU_R_latUp;
Slope2_oxy_NEU_R_latDown = (low_oxy_NEU_R_latDown - peak_oxy_NEU_R_latDown)/ttl_oxy_NEU_R_latDown;
Slope2_oxy_NEU_R_medUp = (low_oxy_NEU_R_medUp - peak_oxy_NEU_R_medUp)/ttl_oxy_NEU_R_medUp;
Slope2_oxy_NEU_R_medDown = (low_oxy_NEU_R_medDown - peak_oxy_NEU_R_medDown)/ttl_oxy_NEU_R_medDown;
Slope2_oxy_NEU_L_latUp = (low_oxy_NEU_L_latUp - peak_oxy_NEU_L_latUp)/ttl_oxy_NEU_L_latUp;
Slope2_oxy_NEU_L_latDown = (low_oxy_NEU_L_latDown - peak_oxy_NEU_L_latDown)/ttl_oxy_NEU_L_latDown;
Slope2_oxy_NEU_L_medUp = (low_oxy_NEU_L_medUp  - peak_oxy_NEU_L_medUp)/ttl_oxy_NEU_L_medUp;
Slope2_oxy_NEU_L_medDown = (low_oxy_NEU_L_medDown  - peak_oxy_NEU_L_medDown)/ttl_oxy_NEU_L_medDown;

Slope2_deoxy_NEU_R_latUp = (low_deoxy_NEU_R_latUp - peak_deoxy_NEU_R_latUp)/ttl_deoxy_NEU_R_latUp;
Slope2_deoxy_NEU_R_latDown = (low_deoxy_NEU_R_latDown - peak_deoxy_NEU_R_latDown)/ttl_deoxy_NEU_R_latDown;
Slope2_deoxy_NEU_R_medUp = (low_deoxy_NEU_R_medUp - peak_deoxy_NEU_R_medUp)/ttl_deoxy_NEU_R_medUp;
Slope2_deoxy_NEU_R_medDown = (low_deoxy_NEU_R_medDown  - peak_deoxy_NEU_R_medDown)/ttl_deoxy_NEU_R_medDown;
Slope2_deoxy_NEU_L_latUp = (low_deoxy_NEU_L_latUp - peak_deoxy_NEU_L_latUp)/ttl_deoxy_NEU_L_latUp;
Slope2_deoxy_NEU_L_latDown = (low_deoxy_NEU_L_latDown - peak_deoxy_NEU_L_latDown)/ttl_deoxy_NEU_L_latDown;
Slope2_deoxy_NEU_L_medUp = (low_deoxy_NEU_L_medUp - peak_deoxy_NEU_L_medUp)/ttl_deoxy_NEU_L_medUp;
Slope2_deoxy_NEU_L_medDown = (low_deoxy_NEU_L_medDown - peak_deoxy_NEU_L_medDown)/ttl_deoxy_NEU_L_medDown;

Slope2_oxy_CON_R_latUp = (low_oxy_CON_R_latUp - peak_oxy_CON_R_latUp)/ttl_oxy_CON_R_latUp;
Slope2_oxy_CON_R_latDown = (low_oxy_CON_R_latDown - peak_oxy_CON_R_latDown)/ttl_oxy_CON_R_latDown;
Slope2_oxy_CON_R_medUp = (low_oxy_CON_R_medUp - peak_oxy_CON_R_medUp)/ttl_oxy_CON_R_medUp;
Slope2_oxy_CON_R_medDown = (low_oxy_CON_R_medDown - peak_oxy_CON_R_medDown)/ttl_oxy_CON_R_medDown;
Slope2_oxy_CON_L_latUp = (low_oxy_CON_L_latUp - peak_oxy_CON_L_latUp)/ttl_oxy_CON_L_latUp;
Slope2_oxy_CON_L_latDown = (low_oxy_CON_L_latDown - peak_oxy_CON_L_latDown)/ttl_oxy_CON_L_latDown;
Slope2_oxy_CON_L_medUp = (low_oxy_CON_L_medUp - peak_oxy_CON_L_medUp)/ttl_oxy_CON_L_medUp;
Slope2_oxy_CON_L_medDown = (low_oxy_CON_L_medDown - peak_oxy_CON_L_medDown)/ttl_oxy_CON_L_medDown;

Slope2_deoxy_CON_R_latUp = (low_deoxy_CON_R_latUp - peak_deoxy_CON_R_latUp)/ttl_deoxy_CON_R_latUp;
Slope2_deoxy_CON_R_latDown = (low_deoxy_CON_R_latDown - peak_deoxy_CON_R_latDown)/ttl_deoxy_CON_R_latDown;
Slope2_deoxy_CON_R_medUp = (low_deoxy_CON_R_medUp - peak_deoxy_CON_R_medUp)/ttl_deoxy_CON_R_medUp;
Slope2_deoxy_CON_R_medDown = (low_deoxy_CON_R_medDown - peak_deoxy_CON_R_medDown)/ttl_deoxy_CON_R_medDown;
Slope2_deoxy_CON_L_latUp = (low_deoxy_CON_L_latUp - peak_deoxy_CON_L_latUp)/ttl_deoxy_CON_L_latUp;
Slope2_deoxy_CON_L_latDown = (low_deoxy_CON_L_latDown - peak_deoxy_CON_L_latDown)/ttl_deoxy_CON_L_latDown;
Slope2_deoxy_CON_L_medUp = (low_deoxy_CON_L_medUp  - peak_deoxy_CON_L_medUp)/ttl_deoxy_CON_L_medUp;
Slope2_deoxy_CON_L_medDown = (low_deoxy_CON_L_medDown - peak_deoxy_CON_L_medDown)/ttl_deoxy_CON_L_medDown;

Slope2_oxy_INC_R_latUp = (low_oxy_INC_R_latUp  - peak_oxy_INC_R_latUp)/ttl_oxy_INC_R_latUp;
Slope2_oxy_INC_R_latDown = (low_oxy_INC_R_latDown - peak_oxy_INC_R_latDown)/ttl_oxy_INC_R_latDown;
Slope2_oxy_INC_R_medUp = (low_oxy_INC_R_medUp - peak_oxy_INC_R_medUp)/ttl_oxy_INC_R_medUp;
Slope2_oxy_INC_R_medDown = (low_oxy_INC_R_medDown - peak_oxy_INC_R_medDown)/ttl_oxy_INC_R_medDown;
Slope2_oxy_INC_L_latUp = (low_oxy_INC_L_latUp - peak_oxy_INC_L_latUp)/ttl_oxy_INC_L_latUp;
Slope2_oxy_INC_L_latDown = (low_oxy_INC_L_latDown - peak_oxy_INC_L_latDown)/ttl_oxy_INC_L_latDown;
Slope2_oxy_INC_L_medUp = (low_oxy_INC_L_medUp - peak_oxy_INC_L_medUp)/ttl_oxy_INC_L_medUp;
Slope2_oxy_INC_L_medDown = (low_oxy_INC_L_medDown - peak_oxy_INC_L_medDown)/ttl_oxy_INC_L_medDown;

Slope2_deoxy_INC_R_latUp = (low_deoxy_INC_R_latUp - peak_deoxy_INC_R_latUp)/ttl_deoxy_INC_R_latUp;
Slope2_deoxy_INC_R_latDown = (low_deoxy_INC_R_latDown - peak_deoxy_INC_R_latDown)/ttl_deoxy_INC_R_latDown;
Slope2_deoxy_INC_R_medUp = (low_deoxy_INC_R_medUp  - peak_deoxy_INC_R_medUp)/ttl_deoxy_INC_R_medUp;
Slope2_deoxy_INC_R_medDown = (low_deoxy_INC_R_medDown - peak_deoxy_INC_R_medDown)/ttl_deoxy_INC_R_medDown;
Slope2_deoxy_INC_L_latUp = (low_deoxy_INC_L_latUp - peak_deoxy_INC_L_latUp)/ttl_deoxy_INC_L_latUp;
Slope2_deoxy_INC_L_latDown = (low_deoxy_INC_L_latDown - peak_deoxy_INC_L_latDown)/ttl_deoxy_INC_L_latDown;
Slope2_deoxy_INC_L_medUp = (low_deoxy_INC_L_medUp - peak_deoxy_INC_L_medUp)/ttl_deoxy_INC_L_medUp;
Slope2_deoxy_INC_L_medDown = (low_deoxy_INC_L_medDown - peak_deoxy_INC_L_medDown)/ttl_deoxy_INC_L_medDown;

%% Activational pattern OXY NEU
qt = tinv(0.975,4) ;
qts = qt/sqrt(5) ;

% If 0 < abs(peak) - SD*qts corresponding to 0 < abs(peak) - qt*SEM
% I.e. if 0 is outside 95% CI


if abs(peak_oxy_NEU_R_latUp) < qts*dcSD(t==ttp_oxy_NEU_R_latUp,1,2,3)
   act_oxy_NEU_R_latUp=0;
elseif abs(peak_oxy_NEU_R_latUp) > qts*dcSD(t==ttp_oxy_NEU_R_latUp,1,2,3)
   if peak_oxy_NEU_R_latUp > 0
      act_oxy_NEU_R_latUp=1;
   elseif peak_oxy_NEU_R_latUp < 0
      act_oxy_NEU_R_latUp=-1;
   end
elseif isnan(peak_oxy_NEU_R_latUp)
    act_oxy_NEU_R_latUp = NaN;
end

if abs(peak_oxy_NEU_R_latDown) < qts*dcSD(t==ttp_oxy_NEU_R_latDown,1,3,3)
   act_oxy_NEU_R_latDown=0;
elseif abs(peak_oxy_NEU_R_latDown) > qts*dcSD(t==ttp_oxy_NEU_R_latDown,1,3,3)
   if peak_oxy_NEU_R_latDown > 0
      act_oxy_NEU_R_latDown=1;
   elseif peak_oxy_NEU_R_latDown < 0
      act_oxy_NEU_R_latDown=-1;
   end
elseif isnan(peak_oxy_NEU_R_latDown)
    act_oxy_NEU_R_latDown = NaN;
end

if abs(peak_oxy_NEU_R_medUp) < qts*dcSD(t==ttp_oxy_NEU_R_medUp,1,4,3)
   act_oxy_NEU_R_medUp=0;
elseif abs(peak_oxy_NEU_R_medUp) > qts*dcSD(t==ttp_oxy_NEU_R_medUp,1,4,3)
   if peak_oxy_NEU_R_medUp > 0
      act_oxy_NEU_R_medUp=1;
   elseif peak_oxy_NEU_R_medUp < 0
      act_oxy_NEU_R_medUp=-1;
   end
elseif isnan(peak_oxy_NEU_R_medUp)
    act_oxy_NEU_R_medUp = NaN;
end

if abs(peak_oxy_NEU_R_medDown) < qts*dcSD(t==ttp_oxy_NEU_R_medDown,1,5,3)
   act_oxy_NEU_R_medDown=0;
elseif abs(peak_oxy_NEU_R_medDown) > qts*dcSD(t==ttp_oxy_NEU_R_medDown,1,5,3)
   if peak_oxy_NEU_R_medDown > 0
      act_oxy_NEU_R_medDown=1;
   elseif peak_oxy_NEU_R_medDown < 0
      act_oxy_NEU_R_medDown=-1;
   end
elseif isnan(peak_oxy_NEU_R_medDown)
    act_oxy_NEU_R_medDown = NaN;
end

if abs(peak_oxy_NEU_L_medUp) < qts*dcSD(t==ttp_oxy_NEU_L_medUp,1,8,3)
   act_oxy_NEU_L_medUp=0;
elseif abs(peak_oxy_NEU_L_medUp) > qts*dcSD(t==ttp_oxy_NEU_L_medUp,1,8,3)
   if peak_oxy_NEU_L_medUp > 0
      act_oxy_NEU_L_medUp=1;
   elseif peak_oxy_NEU_L_medUp < 0
      act_oxy_NEU_L_medUp=-1;
   end
elseif isnan(peak_oxy_NEU_L_medUp)
    act_oxy_NEU_L_medUp = NaN;
end

if abs(peak_oxy_NEU_L_medDown) < qts*dcSD(t==ttp_oxy_NEU_L_medDown,1,9,3)
   act_oxy_NEU_L_medDown=0;
elseif abs(peak_oxy_NEU_L_medDown) > qts*dcSD(t==ttp_oxy_NEU_L_medDown,1,9,3)
   if peak_oxy_NEU_L_medDown > 0
      act_oxy_NEU_L_medDown=1;
   elseif peak_oxy_NEU_L_medDown < 0
      act_oxy_NEU_L_medDown=-1;
   end
elseif isnan(peak_oxy_NEU_L_medDown)
    act_oxy_NEU_L_medDown = NaN;
end

if abs(peak_oxy_NEU_L_latUp) < qts*dcSD(t==ttp_oxy_NEU_L_latUp,1,10,3)
   act_oxy_NEU_L_latUp=0;
elseif abs(peak_oxy_NEU_L_latUp) > qts*dcSD(t==ttp_oxy_NEU_L_latUp,1,10,3)
   if peak_oxy_NEU_L_latUp > 0
      act_oxy_NEU_L_latUp=1;
   elseif peak_oxy_NEU_L_latUp < 0
      act_oxy_NEU_L_latUp=-1;
   end
elseif isnan(peak_oxy_NEU_L_latUp)
    act_oxy_NEU_L_latUp = NaN;
end

if abs(peak_oxy_NEU_L_latDown) < qts*dcSD(t==ttp_oxy_NEU_L_latDown,1,11,3)
   act_oxy_NEU_L_latDown=0;
elseif abs(peak_oxy_NEU_L_latDown) > qts*dcSD(t==ttp_oxy_NEU_L_latDown,1,11,3)
   if peak_oxy_NEU_L_latDown > 0
      act_oxy_NEU_L_latDown=1;
   elseif peak_oxy_NEU_L_latDown < 0
      act_oxy_NEU_L_latDown=-1;
   end
elseif isnan(peak_oxy_NEU_L_latDown)
    act_oxy_NEU_L_latDown = NaN;
end

%% Activational pattern DEOXY NEU

if abs(peak_deoxy_NEU_R_latUp) < qts*dcSD(t==ttp_deoxy_NEU_R_latUp,2,2,3)
   act_deoxy_NEU_R_latUp=0;
elseif abs(peak_deoxy_NEU_R_latUp) > qts*dcSD(t==ttp_deoxy_NEU_R_latUp,2,2,3)
   if peak_deoxy_NEU_R_latUp > 0
      act_deoxy_NEU_R_latUp=1;
   elseif peak_deoxy_NEU_R_latUp < 0
      act_deoxy_NEU_R_latUp=-1;
   end
elseif isnan(peak_deoxy_NEU_R_latUp)
    act_deoxy_NEU_R_latUp = NaN;
end

if abs(peak_deoxy_NEU_R_latDown) < qts*dcSD(t==ttp_deoxy_NEU_R_latDown,2,3,3)
   act_deoxy_NEU_R_latDown=0;
elseif abs(peak_deoxy_NEU_R_latDown) > qts*dcSD(t==ttp_deoxy_NEU_R_latDown,2,3,3)
   if peak_deoxy_NEU_R_latDown > 0
      act_deoxy_NEU_R_latDown=1;
   elseif peak_deoxy_NEU_R_latDown < 0
      act_deoxy_NEU_R_latDown=-1;
   end
elseif isnan(peak_deoxy_NEU_R_latDown)
    act_deoxy_NEU_R_latDown = NaN;
end

if abs(peak_deoxy_NEU_R_medUp) < qts*dcSD(t==ttp_deoxy_NEU_R_medUp,2,4,3)
   act_deoxy_NEU_R_medUp=0;
elseif abs(peak_deoxy_NEU_R_medUp) > qts*dcSD(t==ttp_deoxy_NEU_R_medUp,2,4,3)
   if peak_deoxy_NEU_R_medUp > 0
      act_deoxy_NEU_R_medUp=1;
   elseif peak_deoxy_NEU_R_medUp < 0
      act_deoxy_NEU_R_medUp=-1;
   end
elseif isnan(peak_deoxy_NEU_R_medUp)
    act_deoxy_NEU_R_medUp = NaN;
end

if abs(peak_deoxy_NEU_R_medDown) < qts*dcSD(t==ttp_deoxy_NEU_R_medDown,2,5,3)
   act_deoxy_NEU_R_medDown=0;
elseif abs(peak_deoxy_NEU_R_medDown) > qts*dcSD(t==ttp_deoxy_NEU_R_medDown,2,5,3)
   if peak_deoxy_NEU_R_medDown > 0
      act_deoxy_NEU_R_medDown=1;
   elseif peak_deoxy_NEU_R_medDown < 0
      act_deoxy_NEU_R_medDown=-1;
   end
elseif isnan(peak_deoxy_NEU_R_medDown)
    act_deoxy_NEU_R_medDown = NaN;
end

if abs(peak_deoxy_NEU_L_medUp) < qts*dcSD(t==ttp_deoxy_NEU_L_medUp,2,8,3)
   act_deoxy_NEU_L_medUp=0;
elseif abs(peak_deoxy_NEU_L_medUp) > qts*dcSD(t==ttp_deoxy_NEU_L_medUp,2,8,3)
   if peak_deoxy_NEU_L_medUp > 0
      act_deoxy_NEU_L_medUp=1;
   elseif peak_deoxy_NEU_L_medUp < 0
      act_deoxy_NEU_L_medUp=-1;
   end
elseif isnan(peak_deoxy_NEU_L_medUp)
    act_deoxy_NEU_L_medUp = NaN;
end

if abs(peak_deoxy_NEU_L_medDown) < qts*dcSD(t==ttp_deoxy_NEU_L_medDown,2,9,3)
   act_deoxy_NEU_L_medDown=0;
elseif abs(peak_deoxy_NEU_L_medDown) > qts*dcSD(t==ttp_deoxy_NEU_L_medDown,2,9,3)
   if peak_deoxy_NEU_L_medDown > 0
      act_deoxy_NEU_L_medDown=1;
   elseif peak_deoxy_NEU_L_medDown < 0
      act_deoxy_NEU_L_medDown=-1;
   end
elseif isnan(peak_deoxy_NEU_L_medDown)
    act_deoxy_NEU_L_medDown = NaN;
end

if abs(peak_deoxy_NEU_L_latUp) < qts*dcSD(t==ttp_deoxy_NEU_L_latUp,2,10,3)
   act_deoxy_NEU_L_latUp=0;
elseif abs(peak_deoxy_NEU_L_latUp) > qts*dcSD(t==ttp_deoxy_NEU_L_latUp,2,10,3)
   if peak_deoxy_NEU_L_latUp > 0
      act_deoxy_NEU_L_latUp=1;
   elseif peak_deoxy_NEU_L_latUp < 0
      act_deoxy_NEU_L_latUp=-1;
   end
elseif isnan(peak_deoxy_NEU_L_latUp)
    act_deoxy_NEU_L_latUp = NaN;
end

if abs(peak_deoxy_NEU_L_latDown) < qts*dcSD(t==ttp_deoxy_NEU_L_latDown,2,11,3)
   act_deoxy_NEU_L_latDown=0;
elseif abs(peak_deoxy_NEU_L_latDown) > qts*dcSD(t==ttp_deoxy_NEU_L_latDown,2,11,3)
   if peak_deoxy_NEU_L_latDown > 0
      act_deoxy_NEU_L_latDown=1;
   elseif peak_deoxy_NEU_L_latDown < 0
      act_deoxy_NEU_L_latDown=-1;
   end
elseif isnan(peak_deoxy_NEU_L_latDown)
    act_deoxy_NEU_L_latDown = NaN;
end

%% Activational pattern OXY CON

if abs(peak_oxy_CON_R_latUp) < qts*dcSD(t==ttp_oxy_CON_R_latUp,1,2,1)
   act_oxy_CON_R_latUp=0;
elseif abs(peak_oxy_CON_R_latUp) > qts*dcSD(t==ttp_oxy_CON_R_latUp,1,2,1)
   if peak_oxy_CON_R_latUp > 0
      act_oxy_CON_R_latUp=1;
   elseif peak_oxy_CON_R_latUp < 0
      act_oxy_CON_R_latUp=-1;
   end
elseif isnan(peak_oxy_CON_R_latUp)
    act_oxy_CON_R_latUp = NaN;
end

if abs(peak_oxy_CON_R_latDown) < qts*dcSD(t==ttp_oxy_CON_R_latDown,1,3,1)
   act_oxy_CON_R_latDown=0;
elseif abs(peak_oxy_CON_R_latDown) > qts*dcSD(t==ttp_oxy_CON_R_latDown,1,3,1)
   if peak_oxy_CON_R_latDown > 0
      act_oxy_CON_R_latDown=1;
   elseif peak_oxy_CON_R_latDown < 0
      act_oxy_CON_R_latDown=-1;
   end
elseif isnan(peak_oxy_CON_R_latDown)
    act_oxy_CON_R_latDown = NaN;
end

if abs(peak_oxy_CON_R_medUp) < qts*dcSD(t==ttp_oxy_CON_R_medUp,1,4,1)
   act_oxy_CON_R_medUp=0;
elseif abs(peak_oxy_CON_R_medUp) > qts*dcSD(t==ttp_oxy_CON_R_medUp,1,4,1)
   if peak_oxy_CON_R_medUp > 0
      act_oxy_CON_R_medUp=1;
   elseif peak_oxy_CON_R_medUp < 0
      act_oxy_CON_R_medUp=-1;
   end
elseif isnan(peak_oxy_CON_R_medUp)
    act_oxy_CON_R_medUp = NaN;
end

if abs(peak_oxy_CON_R_medDown) < qts*dcSD(t==ttp_oxy_CON_R_medDown,1,5,1)
   act_oxy_CON_R_medDown=0;
elseif abs(peak_oxy_CON_R_medDown) > qts*dcSD(t==ttp_oxy_CON_R_medDown,1,5,1)
   if peak_oxy_CON_R_medDown > 0
      act_oxy_CON_R_medDown=1;
   elseif peak_oxy_CON_R_medDown < 0
      act_oxy_CON_R_medDown=-1;
   end
elseif isnan(peak_oxy_CON_R_medDown)
    act_oxy_CON_R_medDown = NaN;
end

if abs(peak_oxy_CON_L_medUp) < qts*dcSD(t==ttp_oxy_CON_L_medUp,1,8,1)
   act_oxy_CON_L_medUp=0;
elseif abs(peak_oxy_CON_L_medUp) > qts*dcSD(t==ttp_oxy_CON_L_medUp,1,8,1)
   if peak_oxy_CON_L_medUp > 0
      act_oxy_CON_L_medUp=1;
   elseif peak_oxy_CON_L_medUp < 0
      act_oxy_CON_L_medUp=-1;
   end
elseif isnan(peak_oxy_CON_L_medUp)
    act_oxy_CON_L_medUp = NaN;
end

if abs(peak_oxy_CON_L_medDown) < qts*dcSD(t==ttp_oxy_CON_L_medDown,1,9,1)
   act_oxy_CON_L_medDown=0;
elseif abs(peak_oxy_CON_L_medDown) > qts*dcSD(t==ttp_oxy_CON_L_medDown,1,9,1)
   if peak_oxy_CON_L_medDown > 0
      act_oxy_CON_L_medDown=1;
   elseif peak_oxy_CON_L_medDown < 0
      act_oxy_CON_L_medDown=-1;
   end
elseif isnan(peak_oxy_CON_L_medDown)
    act_oxy_CON_L_medDown = NaN;
end

if abs(peak_oxy_CON_L_latUp) < qts*dcSD(t==ttp_oxy_CON_L_latUp,1,10,1)
   act_oxy_CON_L_latUp=0;
elseif abs(peak_oxy_CON_L_latUp) > qts*dcSD(t==ttp_oxy_CON_L_latUp,1,10,1)
   if peak_oxy_CON_L_latUp > 0
      act_oxy_CON_L_latUp=1;
   elseif peak_oxy_CON_L_latUp < 0
      act_oxy_CON_L_latUp=-1;
   end
elseif isnan(peak_oxy_CON_L_latUp)
    act_oxy_CON_L_latUp = NaN;
end

if abs(peak_oxy_CON_L_latDown) < qts*dcSD(t==ttp_oxy_CON_L_latDown,1,11,1)
   act_oxy_CON_L_latDown=0;
elseif abs(peak_oxy_CON_L_latDown) > qts*dcSD(t==ttp_oxy_CON_L_latDown,1,11,1)
   if peak_oxy_CON_L_latDown > 0
      act_oxy_CON_L_latDown=1;
   elseif peak_oxy_CON_L_latDown < 0
      act_oxy_CON_L_latDown=-1;
   end
elseif isnan(peak_oxy_CON_L_latDown)
    act_oxy_CON_L_latDown = NaN;
end

%% Activational pattern DEOXY CON

if abs(peak_deoxy_CON_R_latUp) < qts*dcSD(t==ttp_deoxy_CON_R_latUp,2,2,1)
   act_deoxy_CON_R_latUp=0;
elseif abs(peak_deoxy_CON_R_latUp) > qts*dcSD(t==ttp_deoxy_CON_R_latUp,2,2,1)
   if peak_deoxy_CON_R_latUp > 0
      act_deoxy_CON_R_latUp=1;
   elseif peak_deoxy_CON_R_latUp < 0
      act_deoxy_CON_R_latUp=-1;
   end
elseif isnan(peak_deoxy_CON_R_latUp)
    act_deoxy_CON_R_latUp = NaN;
end

if abs(peak_deoxy_CON_R_latDown) < qts*dcSD(t==ttp_deoxy_CON_R_latDown,2,3,1)
   act_deoxy_CON_R_latDown=0;
elseif abs(peak_deoxy_CON_R_latDown) > qts*dcSD(t==ttp_deoxy_CON_R_latDown,2,3,1)
   if peak_deoxy_CON_R_latDown > 0
      act_deoxy_CON_R_latDown=1;
   elseif peak_deoxy_CON_R_latDown < 0
      act_deoxy_CON_R_latDown=-1;
   end
elseif isnan(peak_deoxy_CON_R_latDown)
    act_deoxy_CON_R_latDown = NaN;
end

if abs(peak_deoxy_CON_R_medUp) < qts*dcSD(t==ttp_deoxy_CON_R_medUp,2,4,1)
   act_deoxy_CON_R_medUp=0;
elseif abs(peak_deoxy_CON_R_medUp) > qts*dcSD(t==ttp_deoxy_CON_R_medUp,2,4,1)
   if peak_deoxy_CON_R_medUp > 0
      act_deoxy_CON_R_medUp=1;
   elseif peak_deoxy_CON_R_medUp < 0
      act_deoxy_CON_R_medUp=-1;
   end
elseif isnan(peak_deoxy_CON_R_medUp)
    act_deoxy_CON_R_medUp = NaN;
end

if abs(peak_deoxy_CON_R_medDown) < qts*dcSD(t==ttp_deoxy_CON_R_medDown,2,5,1)
   act_deoxy_CON_R_medDown=0;
elseif abs(peak_deoxy_CON_R_medDown) > qts*dcSD(t==ttp_deoxy_CON_R_medDown,2,5,1)
   if peak_deoxy_CON_R_medDown > 0
      act_deoxy_CON_R_medDown=1;
   elseif peak_deoxy_CON_R_medDown < 0
      act_deoxy_CON_R_medDown=-1;
   end
elseif isnan(peak_deoxy_CON_R_medDown)
    act_deoxy_CON_R_medDown = NaN;
end

if abs(peak_deoxy_CON_L_medUp) < qts*dcSD(t==ttp_deoxy_CON_L_medUp,2,8,1)
   act_deoxy_CON_L_medUp=0;
elseif abs(peak_deoxy_CON_L_medUp) > qts*dcSD(t==ttp_deoxy_CON_L_medUp,2,8,1)
   if peak_deoxy_CON_L_medUp > 0
      act_deoxy_CON_L_medUp=1;
   elseif peak_deoxy_CON_L_medUp < 0
      act_deoxy_CON_L_medUp=-1;
   end
elseif isnan(peak_deoxy_CON_L_medUp)
    act_deoxy_CON_L_medUp = NaN;
end

if abs(peak_deoxy_CON_L_medDown) < qts*dcSD(t==ttp_deoxy_CON_L_medDown,2,9,1)
   act_deoxy_CON_L_medDown=0;
elseif abs(peak_deoxy_CON_L_medDown) > qts*dcSD(t==ttp_deoxy_CON_L_medDown,2,9,1)
   if peak_deoxy_CON_L_medDown > 0
      act_deoxy_CON_L_medDown=1;
   elseif peak_deoxy_CON_L_medDown < 0
      act_deoxy_CON_L_medDown=-1;
   end
elseif isnan(peak_deoxy_CON_L_medDown)
    act_deoxy_CON_L_medDown = NaN;
end

if abs(peak_deoxy_CON_L_latUp) < qts*dcSD(t==ttp_deoxy_CON_L_latUp,2,10,1)
   act_deoxy_CON_L_latUp=0;
elseif abs(peak_deoxy_CON_L_latUp) > qts*dcSD(t==ttp_deoxy_CON_L_latUp,2,10,1)
   if peak_deoxy_CON_L_latUp > 0
      act_deoxy_CON_L_latUp=1;
   elseif peak_deoxy_CON_L_latUp < 0
      act_deoxy_CON_L_latUp=-1;
   end
elseif isnan(peak_deoxy_CON_L_latUp)
    act_deoxy_CON_L_latUp = NaN;
end

if abs(peak_deoxy_CON_L_latDown) < qts*dcSD(t==ttp_deoxy_CON_L_latDown,2,11,1)
   act_deoxy_CON_L_latDown=0;
elseif abs(peak_deoxy_CON_L_latDown) > qts*dcSD(t==ttp_deoxy_CON_L_latDown,2,11,1)
   if peak_deoxy_CON_L_latDown > 0
      act_deoxy_CON_L_latDown=1;
   elseif peak_deoxy_CON_L_latDown < 0
      act_deoxy_CON_L_latDown=-1;
   end
elseif isnan(peak_deoxy_CON_L_latDown)
    act_deoxy_CON_L_latDown = NaN;
end

%% Activational pattern OXY INC

if abs(peak_oxy_INC_R_latUp) < qts*dcSD(t==ttp_oxy_INC_R_latUp,1,2,2)
   act_oxy_INC_R_latUp=0;
elseif abs(peak_oxy_INC_R_latUp) > qts*dcSD(t==ttp_oxy_INC_R_latUp,1,2,2)
   if peak_oxy_INC_R_latUp > 0
      act_oxy_INC_R_latUp=1;
   elseif peak_oxy_INC_R_latUp < 0
      act_oxy_INC_R_latUp=-1;
   end
elseif isnan(peak_oxy_INC_R_latUp)
    act_oxy_INC_R_latUp = NaN;
end

if abs(peak_oxy_INC_R_latDown) < qts*dcSD(t==ttp_oxy_INC_R_latDown,1,3,2)
   act_oxy_INC_R_latDown=0;
elseif abs(peak_oxy_INC_R_latDown) > qts*dcSD(t==ttp_oxy_INC_R_latDown,1,3,2)
   if peak_oxy_INC_R_latDown > 0
      act_oxy_INC_R_latDown=1;
   elseif peak_oxy_INC_R_latDown < 0
      act_oxy_INC_R_latDown=-1;
   end
elseif isnan(peak_oxy_INC_R_latDown)
    act_oxy_INC_R_latDown = NaN;
end

if abs(peak_oxy_INC_R_medUp) < qts*dcSD(t==ttp_oxy_INC_R_medUp,1,4,2)
   act_oxy_INC_R_medUp=0;
elseif abs(peak_oxy_INC_R_medUp) > qts*dcSD(t==ttp_oxy_INC_R_medUp,1,4,2)
   if peak_oxy_INC_R_medUp > 0
      act_oxy_INC_R_medUp=1;
   elseif peak_oxy_INC_R_medUp < 0
      act_oxy_INC_R_medUp=-1;
   end
elseif isnan(peak_oxy_INC_R_medUp)
    act_oxy_INC_R_medUp = NaN;
end

if abs(peak_oxy_INC_R_medDown) < qts*dcSD(t==ttp_oxy_INC_R_medDown,1,5,2)
   act_oxy_INC_R_medDown=0;
elseif abs(peak_oxy_INC_R_medDown) > qts*dcSD(t==ttp_oxy_INC_R_medDown,1,5,2)
   if peak_oxy_INC_R_medDown > 0
      act_oxy_INC_R_medDown=1;
   elseif peak_oxy_INC_R_medDown < 0
      act_oxy_INC_R_medDown=-1;
   end
elseif isnan(peak_oxy_INC_R_medDown)
    act_oxy_INC_R_medDown = NaN;
end

if abs(peak_oxy_INC_L_medUp) < qts*dcSD(t==ttp_oxy_INC_L_medUp,1,8,2)
   act_oxy_INC_L_medUp=0;
elseif abs(peak_oxy_INC_L_medUp) > qts*dcSD(t==ttp_oxy_INC_L_medUp,1,8,2)
   if peak_oxy_INC_L_medUp > 0
      act_oxy_INC_L_medUp=1;
   elseif peak_oxy_INC_L_medUp < 0
      act_oxy_INC_L_medUp=-1;
   end
elseif isnan(peak_oxy_INC_L_medUp)
    act_oxy_INC_L_medUp = NaN;
end

if abs(peak_oxy_INC_L_medDown) < qts*dcSD(t==ttp_oxy_INC_L_medDown,1,9,2)
   act_oxy_INC_L_medDown=0;
elseif abs(peak_oxy_INC_L_medDown) > qts*dcSD(t==ttp_oxy_INC_L_medDown,1,9,2)
   if peak_oxy_INC_L_medDown > 0
      act_oxy_INC_L_medDown=1;
   elseif peak_oxy_INC_L_medDown < 0
      act_oxy_INC_L_medDown=-1;
   end
elseif isnan(peak_oxy_INC_L_medDown)
    act_oxy_INC_L_medDown = NaN;
end

if abs(peak_oxy_INC_L_latUp) < qts*dcSD(t==ttp_oxy_INC_L_latUp,1,10,2)
   act_oxy_INC_L_latUp=0;
elseif abs(peak_oxy_INC_L_latUp) > qts*dcSD(t==ttp_oxy_INC_L_latUp,1,10,2)
   if peak_oxy_INC_L_latUp > 0
      act_oxy_INC_L_latUp=1;
   elseif peak_oxy_INC_L_latUp < 0
      act_oxy_INC_L_latUp=-1;
   end
elseif isnan(peak_oxy_INC_L_latUp)
    act_oxy_INC_L_latUp = NaN;
end

if abs(peak_oxy_INC_L_latDown) < qts*dcSD(t==ttp_oxy_INC_L_latDown,1,11,2)
   act_oxy_INC_L_latDown=0;
elseif abs(peak_oxy_INC_L_latDown) > qts*dcSD(t==ttp_oxy_INC_L_latDown,1,11,2)
   if peak_oxy_INC_L_latDown > 0
      act_oxy_INC_L_latDown=1;
   elseif peak_oxy_INC_L_latDown < 0
      act_oxy_INC_L_latDown=-1;
   end
elseif isnan(peak_oxy_INC_L_latDown)
    act_oxy_INC_L_latDown = NaN;
end

%% Activational pattern DEOXY INC

if abs(peak_deoxy_INC_R_latUp) < qts*dcSD(t==ttp_deoxy_INC_R_latUp,2,2,2)
   act_deoxy_INC_R_latUp=0;
elseif abs(peak_deoxy_INC_R_latUp) > qts*dcSD(t==ttp_deoxy_INC_R_latUp,2,2,2)
   if peak_deoxy_INC_R_latUp > 0
      act_deoxy_INC_R_latUp=1;
   elseif peak_deoxy_INC_R_latUp < 0
      act_deoxy_INC_R_latUp=-1;
   end
elseif isnan(peak_deoxy_INC_R_latUp)
    act_deoxy_INC_R_latUp = NaN;
end

if abs(peak_deoxy_INC_R_latDown) < qts*dcSD(t==ttp_deoxy_INC_R_latDown,2,3,2)
   act_deoxy_INC_R_latDown=0;
elseif abs(peak_deoxy_INC_R_latDown) > qts*dcSD(t==ttp_deoxy_INC_R_latDown,2,3,2)
   if peak_deoxy_INC_R_latDown > 0
      act_deoxy_INC_R_latDown=1;
   elseif peak_deoxy_INC_R_latDown < 0
      act_deoxy_INC_R_latDown=-1;
   end
elseif isnan(peak_deoxy_INC_R_latDown)
    act_deoxy_INC_R_latDown = NaN;
end

if abs(peak_deoxy_INC_R_medUp) < qts*dcSD(t==ttp_deoxy_INC_R_medUp,2,4,2)
   act_deoxy_INC_R_medUp=0;
elseif abs(peak_deoxy_INC_R_medUp) > qts*dcSD(t==ttp_deoxy_INC_R_medUp,2,4,2)
   if peak_deoxy_INC_R_medUp > 0
      act_deoxy_INC_R_medUp=1;
   elseif peak_deoxy_INC_R_medUp < 0
      act_deoxy_INC_R_medUp=-1;
   end
elseif isnan(peak_deoxy_INC_R_medUp)
    act_deoxy_INC_R_medUp = NaN;
end

if abs(peak_deoxy_INC_R_medDown) < qts*dcSD(t==ttp_deoxy_INC_R_medDown,2,5,2)
   act_deoxy_INC_R_medDown=0;
elseif abs(peak_deoxy_INC_R_medDown) > qts*dcSD(t==ttp_deoxy_INC_R_medDown,2,5,2)
   if peak_deoxy_INC_R_medDown > 0
      act_deoxy_INC_R_medDown=1;
   elseif peak_deoxy_INC_R_medDown < 0
      act_deoxy_INC_R_medDown=-1;
   end
elseif isnan(peak_deoxy_INC_R_medDown)
    act_deoxy_INC_R_medDown = NaN;
end

if abs(peak_deoxy_INC_L_medUp) < qts*dcSD(t==ttp_deoxy_INC_L_medUp,2,8,2)
   act_deoxy_INC_L_medUp=0;
elseif abs(peak_deoxy_INC_L_medUp) > qts*dcSD(t==ttp_deoxy_INC_L_medUp,2,8,2)
   if peak_deoxy_INC_L_medUp > 0
      act_deoxy_INC_L_medUp=1;
   elseif peak_deoxy_INC_L_medUp < 0
      act_deoxy_INC_L_medUp=-1;
   end
elseif isnan(peak_deoxy_INC_L_medUp)
    act_deoxy_INC_L_medUp = NaN;
end

if abs(peak_deoxy_INC_L_medDown) < qts*dcSD(t==ttp_deoxy_INC_L_medDown,2,9,2)
   act_deoxy_INC_L_medDown=0;
elseif abs(peak_deoxy_INC_L_medDown) > qts*dcSD(t==ttp_deoxy_INC_L_medDown,2,9,2)
   if peak_deoxy_INC_L_medDown > 0
      act_deoxy_INC_L_medDown=1;
   elseif peak_deoxy_INC_L_medDown < 0
      act_deoxy_INC_L_medDown=-1;
   end
elseif isnan(peak_deoxy_INC_L_medDown)
    act_deoxy_INC_L_medDown = NaN;
end

if abs(peak_deoxy_INC_L_latUp) < qts*dcSD(t==ttp_deoxy_INC_L_latUp,2,10,2)
   act_deoxy_INC_L_latUp=0;
elseif abs(peak_deoxy_INC_L_latUp) > qts*dcSD(t==ttp_deoxy_INC_L_latUp,2,10,2)
   if peak_deoxy_INC_L_latUp > 0
      act_deoxy_INC_L_latUp=1;
   elseif peak_deoxy_INC_L_latUp < 0
      act_deoxy_INC_L_latUp=-1;
   end
elseif isnan(peak_deoxy_INC_L_latUp)
    act_deoxy_INC_L_latUp = NaN;
end

if abs(peak_deoxy_INC_L_latDown) < qts*dcSD(t==ttp_deoxy_INC_L_latDown,2,11,2)
   act_deoxy_INC_L_latDown=0;
elseif abs(peak_deoxy_INC_L_latDown) > qts*dcSD(t==ttp_deoxy_INC_L_latDown,2,11,2)
   if peak_deoxy_INC_L_latDown > 0
      act_deoxy_INC_L_latDown=1;
   elseif peak_deoxy_INC_L_latDown < 0
      act_deoxy_INC_L_latDown=-1;
   end
elseif isnan(peak_deoxy_INC_L_latDown)
    act_deoxy_INC_L_latDown = NaN;
end

%% Load temporary tables and concenate

load('Peak_TTP_Slope1_Control.mat')
load('Downslope_Time2Minimum_Control.mat')
load('Act_Control.mat')

Peak_TTP_Slope1_concenate = table(peak_oxy_NEU_R_latUp,     ttp_oxy_NEU_R_latUp,     Slope1_oxy_NEU_R_latUp,     ...
        peak_deoxy_NEU_R_latUp,   ttp_deoxy_NEU_R_latUp,   Slope1_deoxy_NEU_R_latUp,   ...
        peak_oxy_NEU_R_latDown,   ttp_oxy_NEU_R_latDown,   Slope1_oxy_NEU_R_latDown,   ...
        peak_deoxy_NEU_R_latDown, ttp_deoxy_NEU_R_latDown, Slope1_deoxy_NEU_R_latDown, ...
        peak_oxy_NEU_R_medUp,     ttp_oxy_NEU_R_medUp,     Slope1_oxy_NEU_R_medUp,     ...
        peak_deoxy_NEU_R_medUp,   ttp_deoxy_NEU_R_medUp,   Slope1_deoxy_NEU_R_medUp,   ...
        peak_oxy_NEU_R_medDown,   ttp_oxy_NEU_R_medDown,   Slope1_oxy_NEU_R_medDown,   ...
        peak_deoxy_NEU_R_medDown, ttp_deoxy_NEU_R_medDown, Slope1_deoxy_NEU_R_medDown, ...
        peak_oxy_NEU_L_medUp,     ttp_oxy_NEU_L_medUp,     Slope1_oxy_NEU_L_medUp,     ...
        peak_deoxy_NEU_L_medUp,   ttp_deoxy_NEU_L_medUp,   Slope1_deoxy_NEU_L_medUp,   ...
        peak_oxy_NEU_L_medDown,   ttp_oxy_NEU_L_medDown,   Slope1_oxy_NEU_L_medDown,   ...
        peak_deoxy_NEU_L_medDown, ttp_deoxy_NEU_L_medDown, Slope1_deoxy_NEU_L_medDown, ...
        peak_oxy_NEU_L_latUp,     ttp_oxy_NEU_L_latUp,     Slope1_oxy_NEU_L_latUp,     ...
        peak_deoxy_NEU_L_latUp,   ttp_deoxy_NEU_L_latUp,   Slope1_deoxy_NEU_L_latUp,   ...
        peak_oxy_NEU_L_latDown,   ttp_oxy_NEU_L_latDown,   Slope1_oxy_NEU_L_latDown,   ...
        peak_deoxy_NEU_L_latDown, ttp_deoxy_NEU_L_latDown, Slope1_deoxy_NEU_L_latDown, ...
                                                                                       ...
        peak_oxy_CON_R_latUp,     ttp_oxy_CON_R_latUp,     Slope1_oxy_CON_R_latUp,     ...
        peak_deoxy_CON_R_latUp,   ttp_deoxy_CON_R_latUp,   Slope1_deoxy_CON_R_latUp,   ...
        peak_oxy_CON_R_latDown,   ttp_oxy_CON_R_latDown,   Slope1_oxy_CON_R_latDown,   ...
        peak_deoxy_CON_R_latDown, ttp_deoxy_CON_R_latDown, Slope1_deoxy_CON_R_latDown, ...
        peak_oxy_CON_R_medUp,     ttp_oxy_CON_R_medUp,     Slope1_oxy_CON_R_medUp,     ...
        peak_deoxy_CON_R_medUp,   ttp_deoxy_CON_R_medUp,   Slope1_deoxy_CON_R_medUp,   ...
        peak_oxy_CON_R_medDown,   ttp_oxy_CON_R_medDown,   Slope1_oxy_CON_R_medDown,   ...
        peak_deoxy_CON_R_medDown, ttp_deoxy_CON_R_medDown, Slope1_deoxy_CON_R_medDown, ...
        peak_oxy_CON_L_medUp,     ttp_oxy_CON_L_medUp,     Slope1_oxy_CON_L_medUp,     ...
        peak_deoxy_CON_L_medUp,   ttp_deoxy_CON_L_medUp,   Slope1_deoxy_CON_L_medUp,   ...
        peak_oxy_CON_L_medDown,   ttp_oxy_CON_L_medDown,   Slope1_oxy_CON_L_medDown,   ...
        peak_deoxy_CON_L_medDown, ttp_deoxy_CON_L_medDown, Slope1_deoxy_CON_L_medDown, ...
        peak_oxy_CON_L_latUp,     ttp_oxy_CON_L_latUp,     Slope1_oxy_CON_L_latUp,     ...
        peak_deoxy_CON_L_latUp,   ttp_deoxy_CON_L_latUp,   Slope1_deoxy_CON_L_latUp,   ...
        peak_oxy_CON_L_latDown,   ttp_oxy_CON_L_latDown,   Slope1_oxy_CON_L_latDown,   ...
        peak_deoxy_CON_L_latDown, ttp_deoxy_CON_L_latDown, Slope1_deoxy_CON_L_latDown, ...
                                                                                       ...
        peak_oxy_INC_R_latUp,     ttp_oxy_INC_R_latUp,     Slope1_oxy_INC_R_latUp,     ...
        peak_deoxy_INC_R_latUp,   ttp_deoxy_INC_R_latUp,   Slope1_deoxy_INC_R_latUp,   ...
        peak_oxy_INC_R_latDown,   ttp_oxy_INC_R_latDown,   Slope1_oxy_INC_R_latDown,   ...
        peak_deoxy_INC_R_latDown, ttp_deoxy_INC_R_latDown, Slope1_deoxy_INC_R_latDown, ...
        peak_oxy_INC_R_medUp,     ttp_oxy_INC_R_medUp,     Slope1_oxy_INC_R_medUp,     ...
        peak_deoxy_INC_R_medUp,   ttp_deoxy_INC_R_medUp,   Slope1_deoxy_INC_R_medUp,   ...
        peak_oxy_INC_R_medDown,   ttp_oxy_INC_R_medDown,   Slope1_oxy_INC_R_medDown,   ...
        peak_deoxy_INC_R_medDown, ttp_deoxy_INC_R_medDown, Slope1_deoxy_INC_R_medDown, ...
        peak_oxy_INC_L_medUp,     ttp_oxy_INC_L_medUp,     Slope1_oxy_INC_L_medUp,     ...
        peak_deoxy_INC_L_medUp,   ttp_deoxy_INC_L_medUp,   Slope1_deoxy_INC_L_medUp,   ...
        peak_oxy_INC_L_medDown,   ttp_oxy_INC_L_medDown,   Slope1_oxy_INC_L_medDown,   ...
        peak_deoxy_INC_L_medDown, ttp_deoxy_INC_L_medDown, Slope1_deoxy_INC_L_medDown, ...
        peak_oxy_INC_L_latUp,     ttp_oxy_INC_L_latUp,     Slope1_oxy_INC_L_latUp,     ...
        peak_deoxy_INC_L_latUp,   ttp_deoxy_INC_L_latUp,   Slope1_deoxy_INC_L_latUp,   ...
        peak_oxy_INC_L_latDown,   ttp_oxy_INC_L_latDown,   Slope1_oxy_INC_L_latDown,   ...
        peak_deoxy_INC_L_latDown, ttp_deoxy_INC_L_latDown, Slope1_deoxy_INC_L_latDown, ...
    'VariableNames',                                            ...
        {'Peak_oxy_NEU_R_latUp',     'TTP_oxy_NEU_R_latUp',    'Slope1_oxy_NEU_R_latUp',      ...
         'Peak_deoxy_NEU_R_latUp',   'TTP_deoxy_NEU_R_latUp',  'Slope1_deoxy_NEU_R_latUp',    ...
         'Peak_oxy_NEU_R_latDown',   'TTP_oxy_NEU_R_latDown',  'Slope1_oxy_NEU_R_latDown',    ...
         'Peak_deoxy_NEU_R_latDown', 'TTP_deoxy_NEU_R_latDown','Slope1_deoxy_NEU_R_latDown',  ...
         'Peak_oxy_NEU_R_medUp',     'TTP_oxy_NEU_R_medUp',    'Slope1_oxy_NEU_R_medUp',      ...
         'Peak_deoxy_NEU_R_medUp',   'TTP_deoxy_NEU_R_medUp',  'Slope1_deoxy_NEU_R_medUp',    ...
         'Peak_oxy_NEU_R_medDown',   'TTP_oxy_NEU_R_medDown',  'Slope1_oxy_NEU_R_medDown',    ...
         'Peak_deoxy_NEU_R_medDown', 'TTP_deoxy_NEU_R_medDown','Slope1_deoxy_NEU_R_medDown',  ...
         'Peak_oxy_NEU_L_medUp',     'TTP_oxy_NEU_L_medUp',    'Slope1_oxy_NEU_L_medUp',      ...
         'Peak_deoxy_NEU_L_medUp',   'TTP_deoxy_NEU_L_medUp',  'Slope1_deoxy_NEU_L_medUp',    ...
         'Peak_oxy_NEU_L_medDown',   'TTP_oxy_NEU_L_medDown',  'Slope1_oxy_NEU_L_medDown',    ...
         'Peak_deoxy_NEU_L_medDown', 'TTP_deoxy_NEU_L_medDown','Slope1_deoxy_NEU_L_medDown',  ...
         'Peak_oxy_NEU_L_latUp',     'TTP_oxy_NEU_L_latUp',    'Slope1_oxy_NEU_L_latUp',      ...
         'Peak_deoxy_NEU_L_latUp',   'TTP_deoxy_NEU_L_latUp',  'Slope1_deoxy_NEU_L_latUp',    ...
         'Peak_oxy_NEU_L_latDown',   'TTP_oxy_NEU_L_latDown',  'Slope1_oxy_NEU_L_latDown',    ...
         'Peak_deoxy_NEU_L_latDown', 'TTP_deoxy_NEU_L_latDown','Slope1_deoxy_NEU_L_latDown',  ...
                                                                                           ...
         'Peak_oxy_CON_R_latUp',     'TTP_oxy_CON_R_latUp',    'Slope1_oxy_CON_R_latUp',      ...
         'Peak_deoxy_CON_R_latUp',   'TTP_deoxy_CON_R_latUp',  'Slope1_deoxy_CON_R_latUp',    ...
         'Peak_oxy_CON_R_latDown',   'TTP_oxy_CON_R_latDown',  'Slope1_oxy_CON_R_latDown',    ...
         'Peak_deoxy_CON_R_latDown', 'TTP_deoxy_CON_R_latDown','Slope1_deoxy_CON_R_latDown',  ...
         'Peak_oxy_CON_R_medUp',     'TTP_oxy_CON_R_medUp',    'Slope1_oxy_CON_R_medUp',      ...
         'Peak_deoxy_CON_R_medUp',   'TTP_deoxy_CON_R_medUp',  'Slope1_deoxy_CON_R_medUp',    ...
         'Peak_oxy_CON_R_medDown',   'TTP_oxy_CON_R_medDown',  'Slope1_oxy_CON_R_medDown',    ...
         'Peak_deoxy_CON_R_medDown', 'TTP_deoxy_CON_R_medDown','Slope1_deoxy_CON_R_medDown',  ...
         'Peak_oxy_CON_L_medUp',     'TTP_oxy_CON_L_medUp',    'Slope1_oxy_CON_L_medUp',      ...
         'Peak_deoxy_CON_L_medUp',   'TTP_deoxy_CON_L_medUp',  'Slope1_deoxy_CON_L_medUp',    ...
         'Peak_oxy_CON_L_medDown',   'TTP_oxy_CON_L_medDown',  'Slope1_oxy_CON_L_medDown',    ...
         'Peak_deoxy_CON_L_medDown', 'TTP_deoxy_CON_L_medDown','Slope1_deoxy_CON_L_medDown',  ...
         'Peak_oxy_CON_L_latUp',     'TTP_oxy_CON_L_latUp',    'Slope1_oxy_CON_L_latUp',      ...
         'Peak_deoxy_CON_L_latUp',   'TTP_deoxy_CON_L_latUp',  'Slope1_deoxy_CON_L_latUp',    ...
         'Peak_oxy_CON_L_latDown',   'TTP_oxy_CON_L_latDown',  'Slope1_oxy_CON_L_latDown',    ...
         'Peak_deoxy_CON_L_latDown', 'TTP_deoxy_CON_L_latDown','Slope1_deoxy_CON_L_latDown',  ...
                                                                                           ...
         'Peak_oxy_INC_R_latUp',     'TTP_oxy_INC_R_latUp',    'Slope1_oxy_INC_R_latUp',      ...
         'Peak_deoxy_INC_R_latUp',   'TTP_deoxy_INC_R_latUp',  'Slope1_deoxy_INC_R_latUp',    ...
         'Peak_oxy_INC_R_latDown',   'TTP_oxy_INC_R_latDown',  'Slope1_oxy_INC_R_latDown',    ...
         'Peak_deoxy_INC_R_latDown', 'TTP_deoxy_INC_R_latDown','Slope1_deoxy_INC_R_latDown',  ...
         'Peak_oxy_INC_R_medUp',     'TTP_oxy_INC_R_medUp',    'Slope1_oxy_INC_R_medUp',      ...
         'Peak_deoxy_INC_R_medUp',   'TTP_deoxy_INC_R_medUp',  'Slope1_deoxy_INC_R_medUp',    ...
         'Peak_oxy_INC_R_medDown',   'TTP_oxy_INC_R_medDown',  'Slope1_oxy_INC_R_medDown',    ...
         'Peak_deoxy_INC_R_medDown', 'TTP_deoxy_INC_R_medDown','Slope1_deoxy_INC_R_medDown',  ...
         'Peak_oxy_INC_L_medUp',     'TTP_oxy_INC_L_medUp',    'Slope1_oxy_INC_L_medUp',      ...
         'Peak_deoxy_INC_L_medUp',   'TTP_deoxy_INC_L_medUp',  'Slope1_deoxy_INC_L_medUp',    ...
         'Peak_oxy_INC_L_medDown',   'TTP_oxy_INC_L_medDown',  'Slope1_oxy_INC_L_medDown',    ...
         'Peak_deoxy_INC_L_medDown', 'TTP_deoxy_INC_L_medDown','Slope1_deoxy_INC_L_medDown',  ...
         'Peak_oxy_INC_L_latUp',     'TTP_oxy_INC_L_latUp',    'Slope1_oxy_INC_L_latUp',      ...
         'Peak_deoxy_INC_L_latUp',   'TTP_deoxy_INC_L_latUp',  'Slope1_deoxy_INC_L_latUp',    ...
         'Peak_oxy_INC_L_latDown',   'TTP_oxy_INC_L_latDown',  'Slope1_oxy_INC_L_latDown',    ...
         'Peak_deoxy_INC_L_latDown', 'TTP_deoxy_INC_L_latDown','Slope1_deoxy_INC_L_latDown' },...
    'RowNames',{subj}) ;

Peak_TTP_Slope1 = [Peak_TTP_Slope1 ; Peak_TTP_Slope1_concenate] ;


Downslope_Time2Minimum_concenate = table(ttl_oxy_NEU_R_latUp,     Slope2_oxy_NEU_R_latUp,     ...
        ttl_deoxy_NEU_R_latUp,   Slope2_deoxy_NEU_R_latUp,   ...
        ttl_oxy_NEU_R_latDown,   Slope2_oxy_NEU_R_latDown,   ...
        ttl_deoxy_NEU_R_latDown, Slope2_deoxy_NEU_R_latDown, ...
        ttl_oxy_NEU_R_medUp,     Slope2_oxy_NEU_R_medUp,     ...
        ttl_deoxy_NEU_R_medUp,   Slope2_deoxy_NEU_R_medUp,   ...
        ttl_oxy_NEU_R_medDown,   Slope2_oxy_NEU_R_medDown,   ...
        ttl_deoxy_NEU_R_medDown, Slope2_deoxy_NEU_R_medDown, ...
        ttl_oxy_NEU_L_medUp,     Slope2_oxy_NEU_L_medUp,     ...
        ttl_deoxy_NEU_L_medUp,   Slope2_deoxy_NEU_L_medUp,   ...
        ttl_oxy_NEU_L_medDown,   Slope2_oxy_NEU_L_medDown,   ...
        ttl_deoxy_NEU_L_medDown, Slope2_deoxy_NEU_L_medDown, ...
        ttl_oxy_NEU_L_latUp,     Slope2_oxy_NEU_L_latUp,     ...
        ttl_deoxy_NEU_L_latUp,   Slope2_deoxy_NEU_L_latUp,   ...
        ttl_oxy_NEU_L_latDown,   Slope2_oxy_NEU_L_latDown,   ...
        ttl_deoxy_NEU_L_latDown, Slope2_deoxy_NEU_L_latDown, ...
                                                                                       ...
        ttl_oxy_CON_R_latUp,     Slope2_oxy_CON_R_latUp,     ...
        ttl_deoxy_CON_R_latUp,   Slope2_deoxy_CON_R_latUp,   ...
        ttl_oxy_CON_R_latDown,   Slope2_oxy_CON_R_latDown,   ...
        ttl_deoxy_CON_R_latDown, Slope2_deoxy_CON_R_latDown, ...
        ttl_oxy_CON_R_medUp,     Slope2_oxy_CON_R_medUp,     ...
        ttl_deoxy_CON_R_medUp,   Slope2_deoxy_CON_R_medUp,   ...
        ttl_oxy_CON_R_medDown,   Slope2_oxy_CON_R_medDown,   ...
        ttl_deoxy_CON_R_medDown, Slope2_deoxy_CON_R_medDown, ...
        ttl_oxy_CON_L_medUp,     Slope2_oxy_CON_L_medUp,     ...
        ttl_deoxy_CON_L_medUp,   Slope2_deoxy_CON_L_medUp,   ...
        ttl_oxy_CON_L_medDown,   Slope2_oxy_CON_L_medDown,   ...
        ttl_deoxy_CON_L_medDown, Slope2_deoxy_CON_L_medDown, ...
        ttl_oxy_CON_L_latUp,     Slope2_oxy_CON_L_latUp,     ...
        ttl_deoxy_CON_L_latUp,   Slope2_deoxy_CON_L_latUp,   ...
        ttl_oxy_CON_L_latDown,   Slope2_oxy_CON_L_latDown,   ...
        ttl_deoxy_CON_L_latDown, Slope2_deoxy_CON_L_latDown, ...
                                                                                       ...
        ttl_oxy_INC_R_latUp,     Slope2_oxy_INC_R_latUp,     ...
        ttl_deoxy_INC_R_latUp,   Slope2_deoxy_INC_R_latUp,   ...
        ttl_oxy_INC_R_latDown,   Slope2_oxy_INC_R_latDown,   ...
        ttl_deoxy_INC_R_latDown, Slope2_deoxy_INC_R_latDown, ...
        ttl_oxy_INC_R_medUp,     Slope2_oxy_INC_R_medUp,     ...
        ttl_deoxy_INC_R_medUp,   Slope2_deoxy_INC_R_medUp,   ...
        ttl_oxy_INC_R_medDown,   Slope2_oxy_INC_R_medDown,   ...
        ttl_deoxy_INC_R_medDown, Slope2_deoxy_INC_R_medDown, ...
        ttl_oxy_INC_L_medUp,     Slope2_oxy_INC_L_medUp,     ...
        ttl_deoxy_INC_L_medUp,   Slope2_deoxy_INC_L_medUp,   ...
        ttl_oxy_INC_L_medDown,   Slope2_oxy_INC_L_medDown,   ...
        ttl_deoxy_INC_L_medDown, Slope2_deoxy_INC_L_medDown, ...
        ttl_oxy_INC_L_latUp,     Slope2_oxy_INC_L_latUp,     ...
        ttl_deoxy_INC_L_latUp,   Slope2_deoxy_INC_L_latUp,   ...
        ttl_oxy_INC_L_latDown,   Slope2_oxy_INC_L_latDown,   ...
        ttl_deoxy_INC_L_latDown, Slope2_deoxy_INC_L_latDown, ...
    'VariableNames',                                            ...
        {'T2M_oxy_NEU_R_latUp',    'Slope2_oxy_NEU_R_latUp',      ...
         'T2M_deoxy_NEU_R_latUp',  'Slope2_deoxy_NEU_R_latUp',    ...
         'T2M_oxy_NEU_R_latDown',  'Slope2_oxy_NEU_R_latDown',    ...
         'T2M_deoxy_NEU_R_latDown','Slope2_deoxy_NEU_R_latDown',  ...
         'T2M_oxy_NEU_R_medUp',    'Slope2_oxy_NEU_R_medUp',      ...
         'T2M_deoxy_NEU_R_medUp',  'Slope2_deoxy_NEU_R_medUp',    ...
         'T2M_oxy_NEU_R_medDown',  'Slope2_oxy_NEU_R_medDown',    ...
         'T2M_deoxy_NEU_R_medDown','Slope2_deoxy_NEU_R_medDown',  ...
         'T2M_oxy_NEU_L_medUp',    'Slope2_oxy_NEU_L_medUp',      ...
         'T2M_deoxy_NEU_L_medUp',  'Slope2_deoxy_NEU_L_medUp',    ...
         'T2M_oxy_NEU_L_medDown',  'Slope2_oxy_NEU_L_medDown',    ...
         'T2M_deoxy_NEU_L_medDown','Slope2_deoxy_NEU_L_medDown',  ...
         'T2M_oxy_NEU_L_latUp',    'Slope2_oxy_NEU_L_latUp',      ...
         'T2M_deoxy_NEU_L_latUp',  'Slope2_deoxy_NEU_L_latUp',    ...
         'T2M_oxy_NEU_L_latDown',  'Slope2_oxy_NEU_L_latDown',    ...
         'T2M_deoxy_NEU_L_latDown','Slope2_deoxy_NEU_L_latDown',  ...
                                                                                           ...
         'T2M_oxy_CON_R_latUp',    'Slope2_oxy_CON_R_latUp',      ...
         'T2M_deoxy_CON_R_latUp',  'Slope2_deoxy_CON_R_latUp',    ...
         'T2M_oxy_CON_R_latDown',  'Slope2_oxy_CON_R_latDown',    ...
         'T2M_deoxy_CON_R_latDown','Slope2_deoxy_CON_R_latDown',  ...
         'T2M_oxy_CON_R_medUp',    'Slope2_oxy_CON_R_medUp',      ...
         'T2M_deoxy_CON_R_medUp',  'Slope2_deoxy_CON_R_medUp',    ...
         'T2M_oxy_CON_R_medDown',  'Slope2_oxy_CON_R_medDown',    ...
         'T2M_deoxy_CON_R_medDown','Slope2_deoxy_CON_R_medDown',  ...
         'T2M_oxy_CON_L_medUp',    'Slope2_oxy_CON_L_medUp',      ...
         'T2M_deoxy_CON_L_medUp',  'Slope2_deoxy_CON_L_medUp',    ...
         'T2M_oxy_CON_L_medDown',  'Slope2_oxy_CON_L_medDown',    ...
         'T2M_deoxy_CON_L_medDown','Slope2_deoxy_CON_L_medDown',  ...
         'T2M_oxy_CON_L_latUp',    'Slope2_oxy_CON_L_latUp',      ...
         'T2M_deoxy_CON_L_latUp',  'Slope2_deoxy_CON_L_latUp',    ...
         'T2M_oxy_CON_L_latDown',  'Slope2_oxy_CON_L_latDown',    ...
         'T2M_deoxy_CON_L_latDown','Slope2_deoxy_CON_L_latDown',  ...
                                                                                           ...
         'T2M_oxy_INC_R_latUp',    'Slope2_oxy_INC_R_latUp',      ...
         'T2M_deoxy_INC_R_latUp',  'Slope2_deoxy_INC_R_latUp',    ...
         'T2M_oxy_INC_R_latDown',  'Slope2_oxy_INC_R_latDown',    ...
         'T2M_deoxy_INC_R_latDown','Slope2_deoxy_INC_R_latDown',  ...
         'T2M_oxy_INC_R_medUp',    'Slope2_oxy_INC_R_medUp',      ...
         'T2M_deoxy_INC_R_medUp',  'Slope2_deoxy_INC_R_medUp',    ...
         'T2M_oxy_INC_R_medDown',  'Slope2_oxy_INC_R_medDown',    ...
         'T2M_deoxy_INC_R_medDown','Slope2_deoxy_INC_R_medDown',  ...
         'T2M_oxy_INC_L_medUp',    'Slope2_oxy_INC_L_medUp',      ...
         'T2M_deoxy_INC_L_medUp',  'Slope2_deoxy_INC_L_medUp',    ...
         'T2M_oxy_INC_L_medDown',  'Slope2_oxy_INC_L_medDown',    ...
         'T2M_deoxy_INC_L_medDown','Slope2_deoxy_INC_L_medDown',  ...
         'T2M_oxy_INC_L_latUp',    'Slope2_oxy_INC_L_latUp',      ...
         'T2M_deoxy_INC_L_latUp',  'Slope2_deoxy_INC_L_latUp',    ...
         'T2M_oxy_INC_L_latDown',  'Slope2_oxy_INC_L_latDown',    ...
         'T2M_deoxy_INC_L_latDown','Slope2_deoxy_INC_L_latDown' },...
    'RowNames',{subj}) ;

Downslope_Time2Minimum = [Downslope_Time2Minimum ; Downslope_Time2Minimum_concenate] ;


Act_concenate = table(act_oxy_NEU_R_latUp, act_deoxy_NEU_R_latUp,    ...
        act_oxy_NEU_R_latDown,    act_deoxy_NEU_R_latDown,           ...
        act_oxy_NEU_R_medUp,      act_deoxy_NEU_R_medUp,             ...
        act_oxy_NEU_R_medDown,    act_deoxy_NEU_R_medDown,           ...
        act_oxy_NEU_L_medUp,      act_deoxy_NEU_L_medUp,             ...
        act_oxy_NEU_L_medDown,    act_deoxy_NEU_L_medDown,           ...
        act_oxy_NEU_L_latUp,      act_deoxy_NEU_L_latUp,             ...
        act_oxy_NEU_L_latDown,    act_deoxy_NEU_L_latDown,           ...
                                                                     ...
        act_oxy_CON_R_latUp,      act_deoxy_CON_R_latUp,             ...
        act_oxy_CON_R_latDown,    act_deoxy_CON_R_latDown,           ...
        act_oxy_CON_R_medUp,      act_deoxy_CON_R_medUp,             ...
        act_oxy_CON_R_medDown,    act_deoxy_CON_R_medDown,           ...
        act_oxy_CON_L_medUp,      act_deoxy_CON_L_medUp,             ...
        act_oxy_CON_L_medDown,    act_deoxy_CON_L_medDown,           ...
        act_oxy_CON_L_latUp,      act_deoxy_CON_L_latUp,             ...
        act_oxy_CON_L_latDown,    act_deoxy_CON_L_latDown,           ...
                                                                     ...
        act_oxy_INC_R_latUp,      act_deoxy_INC_R_latUp,             ...
        act_oxy_INC_R_latDown,    act_deoxy_INC_R_latDown,           ...
        act_oxy_INC_R_medUp,      act_deoxy_INC_R_medUp,             ...
        act_oxy_INC_R_medDown,    act_deoxy_INC_R_medDown,           ...
        act_oxy_INC_L_medUp,      act_deoxy_INC_L_medUp,             ...
        act_oxy_INC_L_medDown,    act_deoxy_INC_L_medDown,           ...
        act_oxy_INC_L_latUp,      act_deoxy_INC_L_latUp,             ...
        act_oxy_INC_L_latDown,    act_deoxy_INC_L_latDown,           ...
    'VariableNames',                                                 ...
        {'Act_oxy_NEU_R_latUp'      'Act_deoxy_NEU_R_latUp',          ...
         'Act_oxy_NEU_R_latDown',  'Act_deoxy_NEU_R_latDown',        ...
         'Act_oxy_NEU_R_medUp',    'Act_deoxy_NEU_R_medUp',          ...
         'Act_oxy_NEU_R_medDown',  'Act_deoxy_NEU_R_medDown',        ...
         'Act_oxy_NEU_L_medUp',    'Act_deoxy_NEU_L_medUp',          ...
         'Act_oxy_NEU_L_medDown',  'Act_deoxy_NEU_L_medDown',        ...
         'Act_oxy_NEU_L_latUp',    'Act_deoxy_NEU_L_latUp',          ...
         'Act_oxy_NEU_L_latDown',  'Act_deoxy_NEU_L_latDown',        ...
                                                                     ...
         'Act_oxy_CON_R_latUp',    'Act_deoxy_CON_R_latUp',          ...
         'Act_oxy_CON_R_latDown',  'Act_deoxy_CON_R_latDown',        ...
         'Act_oxy_CON_R_medUp',    'Act_deoxy_CON_R_medUp',          ...
         'Act_oxy_CON_R_medDown',  'Act_deoxy_CON_R_medDown',        ...
         'Act_oxy_CON_L_medUp',    'Act_deoxy_CON_L_medUp',          ...
         'Act_oxy_CON_L_medDown',  'Act_deoxy_CON_L_medDown',        ...
         'Act_oxy_CON_L_latUp',    'Act_deoxy_CON_L_latUp',          ...
         'Act_oxy_CON_L_latDown',  'Act_deoxy_CON_L_latDown',        ...
                                                                     ...
         'Act_oxy_INC_R_latUp',    'Act_deoxy_INC_R_latUp',          ...
         'Act_oxy_INC_R_latDown',  'Act_deoxy_INC_R_latDown',        ...
         'Act_oxy_INC_R_medUp',    'Act_deoxy_INC_R_medUp',          ...
         'Act_oxy_INC_R_medDown',  'Act_deoxy_INC_R_medDown',        ...
         'Act_oxy_INC_L_medUp',    'Act_deoxy_INC_L_medUp',          ...
         'Act_oxy_INC_L_medDown',  'Act_deoxy_INC_L_medDown',        ...
         'Act_oxy_INC_L_latUp',    'Act_deoxy_INC_L_latUp',          ...
         'Act_oxy_INC_L_latDown',  'Act_deoxy_INC_L_latDown',      },...
    'RowNames',{subj}) ;

Act = [Act ; Act_concenate] ;

% For the first run
% Peak_TTP_Slope1 = Peak_TTP_Slope1_concenate ;
% Downslope_Time2Minimum = Downslope_Time2Minimum_concenate ;
% Act = Act_concenate ; 

%% Save output

save('Peak_TTP_Slope1_Control.mat', '-mat', 'Peak_TTP_Slope1')
save('Downslope_Time2Minimum_Control.mat', '-mat', 'Downslope_Time2Minimum') 
save('Act_Control.mat', '-mat', 'Act')


%% Delete other variables than group (substantial loading time)
clearvars ('-except', 'group')

