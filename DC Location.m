
%% Q1
sh.d = 532;
uwt = 70;
ucu = 20;
ud = 25
q = ud*(uwt/2000)
sh.s = uwt/ucu
tr.Kwt = 25; tr.Kcu = 2750;
qmax = maxpayld(sh,tr)
qLTL = q - qmax*floor(q/qmax)
q = qLTL; class = 250; disc = 0; MC = 95.23;
qB = [0 0.25 0.5 1 2.5 5 10 15 20 Inf];
i = find(qB(1:end-1) <= q & q < qB(2:end))
ODi = 124.23
ci = ODi*20*q
ODiplus1 = 101.83
qBi = qB(i+1)
ciplus1 = ODiplus1*20*qBi
c_tar = (1 - disc)*max(MC,min(ci,ciplus1))

qTL = qmax*floor(q/qmax)
ppiLTL = 144.3; % Jan 2018 (P)
cLTL = transcharge(qLTL,sh,[],ppiLTL)
rLTL = rateLTL(qLTL,sh.s,sh.d,ppiLTL)
cLTL = rLTL * sh.d * qLTL

Diff = c_tar - cLTL


%% Q2
clear all
sh.d = 625;
sh.s = 12;
tr.Kwt = 25; tr.Kcu = 2750;
ppiTL = 123.4; % Jan 2018 (P)
ppiLTL = 141.4; % Jan 2018 (P)
tr.r = 2 * (ppiTL/102.7)
sh.f = 75
sh = vec2struct(sh,'a',1,'v',12000,'h',.3)
[TLC,qTL] = minTLC(sh,tr)
[TLC,TC,IC] = totlogcost(qTL,transcharge(qTL,sh,tr),sh)

% Considering one shipment per week
q1 = 75/52;
n1 = sh.f/q1;
t1 = (1/n1)*52;
[c,isLTL,cTL,cLTL] = transcharge(q1,sh,tr,ppiLTL);
[TLC1,TC1,IC1] = totlogcost(q1,c,sh);
vdisp('TLC,TC,IC');
vdisp('TLC1,TC1,IC1');

Diff = TLC1-TLC    %Increase in total logistic cost


%% Q3 

clear all
% Considering independent P2P truckloads for A and B
% For component A
sh.d = 500;
uwt = 24; ucu = 12;
sh.s = uwt/ucu
tr.Kwt = 25; tr.Kcu = 2750;
tr.r = 2
qAmax = maxpayld(sh,tr)

sh.f = 24*1000/2000
nA = sh.f/qAmax

tA = qAmax/sh.f, tA = 1/nA
daysA = tA * 365.25

rFTLA = tr.r/qAmax
TC_FTLA = sh.f*rFTLA*sh.d

a = (2000/24)*30
h = (.80/4)+0.04+0.06
sh = vec2struct(sh,'a',0.5,'v',2500,'h',.30)
IC_FTLA = sh.a*sh.v*sh.h * qAmax

TLC_FTLA = TC_FTLA + IC_FTLA

% For component B

sh.d = 500;
uwt = 45; ucu = 4;
sh.s = uwt/ucu
tr.Kwt = 25; tr.Kcu = 2750;
tr.r = 2
qBmax = maxpayld(sh,tr)

sh.f = 45*1200/2000
nB = sh.f/qBmax

tB = qBmax/sh.f, tB = 1/nB
daysB = tB * 365.25

rFTLB = tr.r/qBmax
TC_FTLB = sh.f*rFTLB*sh.d

a = (2000/45)*45
h = (.80/6)+0.04+0.06
sh = vec2struct(sh,'a',.5,'v',2000,'h',.2333)

IC_FTLB = sh.a*sh.v*sh.h * qBmax
TLC_FTLB = TC_FTLB + IC_FTLB

TLC_FTL = TLC_FTLA + TLC_FTLB

% Considering combined P2P truckloads for A and B

sh = vec2struct('d', [500 500],'f',[12 27],'s',[2 11.5], 'a', [0.5 0.5],'v', [2500 2000], 'h', [0.3000 0.2333], 'qmax', [qAmax qBmax], 'TLC', [TLC_FTLA TLC_FTLB])

ash = aggshmt(sh);
sdisp([sh ash],[],'sh12,agg3')
[TLC_AB,qAB,isLTL] = minTLC(ash,tr)

Diff = TLC_AB - TLC_FTL  %Decrease in total logistic cost with combined shipment.

%% Q4
clear all
city2lonlat = @(city,st) ...
       uscity('XY',mand(city,uscity('Name'),st,uscity('ST'))) %Defining city2lonlat function handle
    DCcity = 'Winston-Salem'; %DC location
    DC = city2lonlat(DCcity,'NC');
    Ccity = {'Statesville','Wilmington'}; %Customer's location
    Cst = {'NC','NC'};
    CXY = city2lonlat(Ccity,Cst);
    Scity = {'Asheville','Raleigh'}; %Plant's location
    Sst = {'NC','NC'};
    SXY = city2lonlat(Scity,Sst);

    %Creating inbound and outbound shipments
    aS = 0;                % alpha for batch production at suppliers
    aC = 0.5;              % alpha for constant consumption at customers
    aO = aS+0.5;           % inbound alpha for no coordination   
    aD = 0 + aC;           % outbound alpha for no coordination
    ppiTL = 125;           % Given
    ppiLTL = 140;          % Given
    tr = struct('r',2*(ppiTL/102.7),'Kwt',25,'Kcu',2750);
    sS = [150/15 250/10]; % density (lb/ft^3) 
    fS = [4000*150/2000 1200*250/2000];        % annual demand (ton/yr)
    pct = [40 60]/100;   % Plant demand percentage
    uval = [800 450];      % unit value ($)
    uwt = [150 250];         %Unit Weight (lb)
    v = uval./(uwt/2000);   % product cost ($/ton)
    shS = vec2struct(...
       'f',fS,'s',sS,'a',aO,'v',v,'h',0.3,'d',dists(SXY,DC,'mi')*1.2);
    shC = vec2struct(...
       aggshmt(shS),'f',sum(fS)*pct,'a',aD,'d',dists(DC,CXY,'mi')*1.2);
    sh = [shS shC];
    qmax = maxpayld(sh,tr);
    sdisp(sh)

    %Determining total Logistics Cost
    [TLC,q2,isLTL2] = minTLC(sh,tr,ppiLTL);
    n2 = [sh.f]./q2;
    t2 = 365.25./n2;
    vdisp('TLC,q2,qmax,n2,t2,isLTL2',true,true)
    disp('Total Logistics Cost is :') 
    disp (sum(TLC))

%% Q5 
clear all
T = readtable('HW6data.xlsx');
S = table2struct(T);
s = uszip5;

tr.Kwt = 25; tr.Kcu = 2750;
tr.r = 2;
ucu = [S.cu];
uwt = [S.wt];
ud = [S.ud];
fS = [ud .* uwt]./2000;
zip = [S.zip];
uval = [S.uc];

Z = uszip5('XY',...
   mand(zip,uszip5('Code5')));

city2lonlat = @(city,st) ...
   uscity('XY',mand(city,uscity('Name'),st,uscity('ST')));

CXY = city2lonlat('Petersburg','VA');
a = 0.5; 
shS = vec2struct('f',fS,'s',uwt./ucu,'d',dists(CXY,Z),'a',a,'v',uval./(uwt/2000),'h',.3);

%Assuming min TLC
[TC0,q0] = minTLC(shS,tr);
TC0tot = sum(TC0)

% Finding qmin for shipment interval less than 1 month
qmin = (1/12).*[shS.f];

TC1 = totlogcost(qmin,transcharge(qmin,shS,tr),shS);

TC1tot = sum(TC1)


%% Q6
clear all
    city2lonlat = @(city,st) ...
       uscity('XY',mand(city,uscity('Name'),st,uscity('ST'))) %Defining city2lonlat function handle
    Ccity = {'Atlanta','Nashville','Boston','Denver','Chicago','Columbia'}; %Location of Plants
    Cst = {'GA','TN','MA','CO','IL','SC'};
    CXY = city2lonlat(Ccity,Cst);
    Scity = {'San Diego','Seattle','Portland'}; %Location of Suppliers
    Sst = {'CA','WA','ME'};
    SXY = city2lonlat(Scity,Sst);

    % Different Products Supplied to 6 plants
    ppiTL = 131.4; % Jan 2018 (P)
    ppiLTL = 179.4; % Jan 2018 (P)
    tr = struct('r',2*(ppiTL/102.7),'Kwt',25,'Kcu',2750);
    sS = [60/3 20/5 84/12];    % density (lb/ft^3) 
    fS = [32800*60/2000 4800*20/2000 13400*84/2000];         % demand (ton/yr)
    aS = 0;                    % alpha for batch production at suppliers
    aC = 0.5;                  % alpha for constant consumption at customers
    pct = [21 15 24 12 18 10]/100;   % customer demand percentage
    uval = [4500 700 3500];      % unit value ($)
    uwt = [60 20 84];         %Unit Weight (lb)
    v = uval./(uwt/2000);   % product cost ($/ton)

    %Creating Shipments
    F = fS(:)*pct
    S = repmat(sS(:),1,size(CXY,1));
    V = repmat(v(:),1,size(CXY,1));
    D = dists(SXY,CXY,'mi')*1.2;
    sh = vec2struct('f',F(:),'s',S(:),'a',aS+aC,'v',V(:),'h',.3,'d',D(:));
    sdisp(sh)

    % Determine optimal shipment sizes and cost
    [TLC1,q1,isLTL1] = minTLC(sh,tr,ppiLTL);
    qmax1 = maxpayld(sh,tr);
    n1 = [sh.f]./q1;
    t1 = 365.25./n1;
    tmax = 1/12;
    nmin = 1/tmax;
    ndisp = max(n1,nmin)
    vdisp('TLC1,q1,qmax1,ndisp,t1,isLTL1',true,true)
   
    disp (sum(TLC1)) %Total Logistics Cost in direct shipment

    DCcity = 'Kansas City';
    DC = city2lonlat(DCcity,'KS');

    % Create shipments inbound and outbound to the DC
    aO = aS+0.5;           % inbound alpha for no coordination   
    aD = 0 + aC;           % outbound alpha for no coordination
    shS = vec2struct(...
       'f',fS,'s',sS,'a',aO,'v',v,'h',0.3,'d',dists(SXY,DC,'mi')*1.2);
    shC = vec2struct(...
       aggshmt(shS),'f',sum(fS)*pct,'a',aD,'d',dists(DC,CXY,'mi')*1.2);
    sh = [shS shC];
    qmax = maxpayld(sh,tr);
    sdisp(sh)

    % Determine initial DC Location using FTL approximation
    qmax = maxpayld(sh,tr)
    w = [sh.f]./qmax
    DCftl = minisumloc([SXY; CXY],w,'mi')
    [~,~,~,dstr] = lonlat2city(DCftl); disp(dstr);

    % Optimal DC Location with No Cooridination
    % Determine optimal shipment size as subroutine in location procedure
    aO = aS+0.5;           % inbound alpha for no coordination   
    aD = 0 + aC;           % outbound alpha for no coordination
    shS = vec2struct(...
       'f',fS,'s',sS,'a',aO,'v',v,'h',0.3,'d',dists(SXY,DC,'mi')*1.2);
    shC = vec2struct(...
       aggshmt(shS),'f',sum(fS)*pct,'a',aD,'d',dists(DC,CXY,'mi')*1.2);
    [shS.a] = deal(0.5); sh = [shS shC];
    shd_h = @(xy) vec2struct(sh,'d',[dists(SXY,[xy(1) xy(2)],'mi')'*1.2 ...
       dists([xy(1) xy(2)],CXY,'mi')*1.2]);
    TLC5h = @(xy) sum(minTLC(shd_h(xy),tr,ppiLTL));
    [DC5,TLC5] = fminsearch(TLC5h,DCftl)
    [~,~,~,dstr] = lonlat2city(DC5); disp(dstr)

    [TLC5,q5,isLTL5] = minTLC(shd_h(DC5),tr,ppiLTL);
    n5 = [sh.f]./q5;
    t5 = 365.25./n5;
    vdisp('TLC5,q5,qmax,n5,t5,isLTL5',true,true)

    % Optimal DC Location with Perfect Cross-Docking
    % Determine optimal common shipment interval along with location (3-D)
    [shS.a] = deal(0);   % set inbound alpha to 0
    sh = [shS shC];
    shd_h = @(xy) vec2struct(sh,'d',[dists(SXY,[xy(1) xy(2)],'mi')' ...
       dists([xy(1) xy(2)],CXY,'mi')]*1.2);
    TLC60h = @(xyt) totlogcost([sh.f]*xyt(3),...
       transcharge([sh.f]*xyt(3),shd_h(xyt(1:2)),tr,ppiLTL),shd_h(xyt(1:2)));
    TLC6h = @(xyt) iff([sh.f]*xyt(3) <= qmax, TLC60h(xyt), Inf);
    tx0 = min(qmax./[sh.f]);
    [xyt6,TLC6] = fminsearch(@(xyt) sum(TLC6h(xyt)),[DCftl tx0])
    [~,~,~,dstr] = lonlat2city(xyt6(1:2)); disp(dstr)

    t6 = xyt6(3);
    q6 = [sh.f]*t6;
    [~,isLTL6] = transcharge(q6,shd_h(xyt6(1:2)),tr,ppiLTL);
    TLC6 = TLC6h(xyt6);
    n6 = [sh.f]./q6;
    t6 = 365.25./n6;
    vdisp('TLC6,q6,qmax,n6,t6,isLTL6',true,true)

    %Final Results
   
    disp (sum(TLC1))  % Total Logistics Cost in direct shipment

    disp (sum(TLC5)) % Total Logistics Cost in uncordinated shipment with DC near Hendersonville, TN

    disp (sum(TLC6)) % Total Logistics Cost in perfect cross-docking with DC near Bowling Green, KY

    disp(dstr) %The best location for DC.

    disp(sum(TLC1)-sum(TLC6)) %Best savings are obtained during perfect cross docking scenario with location at Bowling Green, KY.

