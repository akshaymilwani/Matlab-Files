
%% Q1
clear all
clc

C1 = xlsread('HW8data', 'P1-Data')
b = [C1(1,2) C1(1,3)]
e = [C1(2:46,2) C1(2:46,3)]
XY = [C1(:,2) C1(:,3)]
sh = vec2struct('b',1,'e',2:size(XY,1));

q = C1(2:46,5)./2000
s = C1(2:46,4)
sh = vec2struct(sh,'q',q,'s',s)
tr = struct('b',1,'e',1);
tr = vec2struct(tr,'Kwt',25,'Kcu',2750)
sdisp(sh)

% Use Delaunay triangulation as road network
tri = delaunay(XY(:,1),XY(:,2));  % Delaunay triangulation
IJ = tri2list(tri);  % Convert to arc list
d = diag(dists(XY(IJ(:,1),:),XY(abs(IJ(:,2)),:),'mi'));
d = d * 1.2;  % Convert great circle to estimated road distances
d = round(d);
IJD = [IJ d];
pplot(IJD,XY,'m-')
pplot(IJD,num2cellstr(d),XY)

% Get road network
expansionAroundXY = 0.12;
[XY2,IJD,isXY,isIJD] = subgraph(usrdnode('XY'),...
   isinrect(usrdnode('XY'),boundrect(XY,expansionAroundXY)),...
   usrdlink('IJD'));

% Label type of road
s10 = usrdlink(isIJD);
isI = s10.Type == 'I';         % Interstate highways
isIR = isI & s10.Urban == ' '; % Rural Interstate highways
isIU = isI & ~isIR;          % Urban Interstate highways
isR = s10.Urban == ' ' & ~isI; % Rural non-Interstate roads
isU = ~isI & ~isR;           % Urban non-Interstate roads

% Plot roads
makemap(XY2,0.12)  % 3% expansion
h = [];  % Keep handle to each plot for legend
h = [h pplot(IJD(isR,:),XY2,'r-','DisplayName','Rural Roads')];
h = [h pplot(IJD(isU,:),XY2,'k-','DisplayName','Urban Roads')];
h = [h pplot(IJD(isI,:),XY2,'c-','DisplayName','Interstate Roads')];

% Add connector roads from cities to road network
[IJD11,IJD12,IJD22] = addconnector(XY,XY2,IJD);
h = [h pplot(IJD12,[XY; XY2],'b-','DisplayName','Connector Roads')];
h = [h pplot(XY,'g.','DisplayName','Cities')];

% Convert road distances to travel times (needs to be after ADDCONNECTOR)
v.IR = 70;  % Rural Interstate highways average speed (mph)
v.IU = 50;  % Urban Interstate highways average speed (mph)
v.R = 45;   % Rural non-Interstate roads average speed (mph)
v.U = 20;   % Urban non-Interstate roads average speed (mph)
v.C = 15;   % Facility to road connector average speed (mph)

IJT = IJD;
IJT(isIR,3) = IJD(isIR,3)/v.IR;
IJT(isIU,3) = IJD(isIU,3)/v.IU;
IJT(isR,3) = IJD(isR,3)/v.R;
IJT(isU,3) = IJD(isU,3)/v.U;

IJT22 = IJD22;                % road to road
IJT22(:,3) = IJT(:,3);
IJT12 = IJD12;                % facility to road
IJT12(:,3) = IJD12(:,3)/v.C;  % (IJD11 facility to facility arcs ignored)

% Shortest time routes
n = size(XY,1);
tic
[T,P] = dijk(list2adj([IJT12; IJT22]),1:n,1:n);
toc

D = dists(XY,XY,'mi')*1.2;
Tmax = 10;
T = T + 5/60;   % 5 min positioning time at each stop
tU = q * 3/60; % 3 min/ton unloading time (ignore loading time at depot)
temin = 7;  
temax = 18; 
sh = vec2struct('b',1,'e',2:size(XY,1));
sh = vec2struct(sh,'q',q,'s',s);
sh = vec2struct(sh,'tU',tU,'temin',temin,'temax',temax)
tr = struct('b',1,'e',1,'Kwt',25,'Kcu',2750, 'tbmin',7, 'tbmax',18, 'temin',7, 'temax',18)
sdisp(sh)

% Construct & improve routes:
rTDh0 = @(rte) rteTC(rte,sh,T,tr);
rTDh = @(rte) myrteTT(rte,rTDh0,Tmax);
ph = @(rte) plotshmt(sh,XY,rte,tr);
IJS = pairwisesavings(rTDh,sh);

r6 = twoopt(savings(rTDh,sh,IJS,ph),rTDh,ph);
pplot(XY(1,:),'ks')
plotshmt(sh,XY,r6,tr)

% Display route output structure
[TC6,Xflg6,out6] = rTDh0(r6);
for i = 1:length(out6), sdisp(out6(i),false,i), end


%% Q2
clear all
XY = xlsread('HW8data', 'P2-Locations');

R = xlsread('HW8data', 'P2-Requests');
Pickup = R(:,1);
Dropoff = R(:,2);
Arrive = R(:,3);
Depart = R(:,4);

PickupXY = [XY(Pickup,1), XY(Pickup,2)];
DropoffXY = [XY(Dropoff,1), XY(Dropoff,2)];

XY0 = [XY(1,1) XY(1,2)]
XY = [XY0; PickupXY; DropoffXY]

D = dists(XY,XY,'mi')*1.2;

%For onward journey

sh = vec2struct('b',Pickup,'e',Dropoff)
tr = struct('b',1,'e',1,'maxCust',8);
sdisp(sh)

T = D/50;       % 50 mph average travel speed
tU = 3/60; % 3 min unloading time
tL =3/60;  % 3 min loading time   
sh = vec2struct(sh,'tU',tU,'tL',tL,'tbmin',Arrive-2,'temax',Arrive);
sdisp(sh)

% Construct & improve routes:
maxCust = 8;
rTDh0 = @(rte) rteTC(rte,sh,T,tr);
rTDh = @(rte) myrteCust(rte,rTDh0,maxCust);
ph = @(rte) plotshmt(sh,XY,rte,tr);
IJS = pairwisesavings(rTDh,sh);
r6 = twoopt(savings(rTDh,sh,IJS,ph),rTDh,ph);
pplot(XY(1,:),'ks')
plotshmt(sh,XY,r6,tr)

% Display route output structure
[TC6,Xflg6,out6] = rTDh0(r6);
for i = 1:length(out6), sdisp(out6(i),false,i), end

%For return journey

sh = vec2struct('b',Dropoff,'e',Pickup)
tr = struct('b',1,'e',1,'maxCust',8);
sdisp(sh)

T = D/50;       % 50 mph average travel speed
tU = 3/60; % 3 min unloading time
tL =3/60;  % 3 min loading time   
sh = vec2struct(sh,'tU',tU,'tL',tL,'tbmin',Depart,'temax',Depart+2);
sdisp(sh)

% Construct & improve routes:
maxCust = 8;
rTDh0 = @(rte) rteTC(rte,sh,T,tr);
rTDh = @(rte) myrteCust(rte,rTDh0,maxCust);
ph = @(rte) plotshmt(sh,XY,rte,tr);
IJS = pairwisesavings(rTDh,sh);
r6 = twoopt(savings(rTDh,sh,IJS,ph),rTDh,ph);
pplot(XY(1,:),'ks')
plotshmt(sh,XY,r6,tr)

% Display route output structure
[TC6,Xflg6,out6] = rTDh0(r6);
for i = 1:length(out6), sdisp(out6(i),false,i), end


%% Question 3
clear, close all
load shmtNC30
tr = struct('r',2,'Kwt',25,'Kcu',2750);
sh = vec2struct('b',b([1 3 26 5]),'e',e([1 3 26 5]));
sh = vec2struct(sh,'d',diag(D([sh.b],[sh.e])));
sdisp(sh,1)
r = 2
rTCh = @(rte) rteTC(rte,sh,D*r);
n = length(sh)

% Independent transport charge
c0 = [sh.d]*r

% Min incremental charge for all possible routes
R = perms(1:n)
R = sortrows(R,1:n)
C = zeros(size(R));
for i = 1:size(C,1)
   for j = 1:size(C,2)
      Rj = perms(R(i,1:j));  % Try all permutations to get optimal
      TC(j) = Inf;
      for k = 1:size(Rj,1)
         [~,TCj] = insertimprove(Rj(k,:),rTCh,sh);
         if TCj < TC(j), TC(j) = TCj; end
      end
   end
   C(i,:) = TC;
   TC = diff([0 TC]);
   C(i,:) = TC(invperm(R(i,:)));
end
mdisp(C,sum(R.*repmat(10.^[n-1:-1:0],size(R,1),1),2))

% Equal charge allocation
TCc = min(sum(C,2))
c_equal = repmat(TCc/n,1,n)
pct_reduct = round(100*(1 - c_equal./c0))

% Equal savings allocation
Sn = sum(c0) - TCc
c_eq_sav = c0 - Sn/n
pct_reduct = round(100*(1 - c_eq_sav./c0))

% Exact Shapley allocation
c_Shap_exact = mean(C,1)
pct_reduct = round(100*(1 - c_Shap_exact./c0))

% Pairwise approximate Shapely allocation
[~,S2] = pairwisesavings(rTCh,sh)
c_Shap_approx = c0 - (Sn/n + sum(S2)/(n-1) - sum(sum(S2))/(n*(n-1)))
pct_reduct = round(100*(1 - c_Shap_approx./c0))

%% Comparison
vdisp('c0,c_equal,c_eq_sav,c_Shap_exact,c_Shap_approx',true,true)






