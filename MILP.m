
%% Q1a

clear all
s1 = uscity;
DC = uszip5('XY',mand([20548 26149 36317],uszip5('Code5')))
Cust = uszip5('XY',mand([30669 38339 30732 23830 23154],uszip5('Code5')))
dem = [40 55 35 70 25]

k = rand(1,1)
C = k.*dists(DC,Cust,'mi')

s = [Inf Inf Inf]
[F,TC] = trans(C,s,dem)

Tot1 = sum(sum(F.*dists(DC,Cust,'mi'),2))

%% Q1b

s1 = uscity;
DC = uszip5('XY',mand([20548 26149 36317],uszip5('Code5')))
Cust = uszip5('XY',mand([30669 38339 30732 23830 23154],uszip5('Code5')))
dem = [40 55 35 70 25]

k = rand(1,1)
C = k.*dists(DC,Cust,'mi')

s = [80 80 80]
[F,TC] = trans(C,s,dem)

Tot2 = sum(sum(F.*dists(DC,Cust,'mi'),2))
Diff1 = Tot2 - Tot1  % Change in ton-miles

%% Q1c 

s1 = uscity;
DC = uszip5('XY',mand([20548 26149 36317],uszip5('Code5')))
Cust = uszip5('XY',mand([30669 38339 30732 23830 23154],uszip5('Code5')))
dem = [40 55 35 70 25]

k = rand(1,1)
C = k.*dists(DC,Cust,'mi')

s = [80 80 80]

IJC = adj2list(lev2adj(C));

A(1:15)=35;
A1 = A';

IJC = [IJC, A1];

s = [s -dem];

lp = mcnf2lp(IJC,s);
[x,TC,XFlg,out] = lplog(lp{:})
[f,TC,nf] = lp2mcnf(x,IJC,s)

IJF = [IJC(:,[1 2]) f]
AF = list2adj(IJF)
F = adj2lev(AF,size(C))

Tot3 = sum(sum(F.*dists(DC,Cust,'mi'),2))
Diff2 = Tot3 - Tot2   % Change in ton-miles


%% Q1d 

s1 = uscity;
Plant = uszip5('XY',mand([28124 27115 37421 27513],uszip5('Code5')))
DC = uszip5('XY',mand([20548 26149 36317],uszip5('Code5')))
Cust = uszip5('XY',mand([30669 38339 30732 23830 23154],uszip5('Code5')))

s = [60 60 60 60 0 0 0 -dem] 
C23=C;
C12 =k.*dists(Plant,DC,'mi')
W = lev2adj(C12,C23)
IJC = adj2list(W)

lp = mcnf2lp(IJC,s);
[x,TC,XFlg,out] = lplog(lp{:})
[f,TC,nf] = lp2mcnf(x,IJC,s)

IJF = [IJC(:,[1 2]) f]
AF = list2adj(IJF)

[F12,F23] = adj2lev(AF,[4 3 5]) 
% F12 represents the product supplied from each plant to the DCs and F23 represents the product supplied from each DC to the customers 

W = lev2adj(C12,C23);
mdisp(W)


%% Q2 Shortest Road Travel Time

clear all
clc

Ral = uscity('XY',mand({'Raleigh'},uscity('Name'),{'NC'},uscity('ST')));
[Name,XY1]=uscity10k('Name','XY',(dists(Ral,uscity10k('XY'),'mi')<200)' ...
   & uscity10k('Pop')>30000);

% XY1 = [-80.2527 37.1354] 
XY1(33,1) = [-80.2527]
XY1(33,2) = [37.1354]
Name(33) = {'Virginia Tech'}

% Get road network
expansionAroundXY = 0.1;
[XY2,IJD,isXY,isIJD] = subgraph(usrdnode('XY'),...
   isinrect(usrdnode('XY'),boundrect(XY1,expansionAroundXY)),...
   usrdlink('IJD'));

% Label type of road
s = usrdlink(isIJD);
isI = s.Type == 'I';         % Interstate highways
isIR = isI & s.Urban == ' '; % Rural Interstate highways
isIU = isI & ~isIR;          % Urban Interstate highways
isR = s.Urban == ' ' & ~isI; % Rural non-Interstate roads
isU = ~isI & ~isR;           % Urban non-Interstate roads

% Plot roads
makemap(XY2,0.03)  % 3% expansion
h = [];  % Keep handle to each plot for legend
h = [h pplot(IJD(isR,:),XY2,'r-','DisplayName','Rural Roads')];
h = [h pplot(IJD(isU,:),XY2,'k-','DisplayName','Urban Roads')];
h = [h pplot(IJD(isI,:),XY2,'c-','DisplayName','Interstate Roads')];

% Add connector roads from cities to road network
[IJD11,IJD12,IJD22] = addconnector(XY1,XY2,IJD);
h = [h pplot(IJD12,[XY1; XY2],'b-','DisplayName','Connector Roads')];
h = [h pplot(XY1,'g.','DisplayName','Cities')];

% Convert road distances to travel times (needs to be after ADDCONNECTOR)
v.IR = 75;  % Rural Interstate highways average speed (mph)
v.IU = 65;  % Urban Interstate highways average speed (mph)
v.R = 50;   % Rural non-Interstate roads average speed (mph)
v.U = 25;   % Urban non-Interstate roads average speed (mph)
v.C = 20;   % Facility to road connector average speed (mph)

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
n = size(XY1,1);
tic
[T,P] = dijk(list2adj([IJT12; IJT22]),1:n,1:n);
toc

% Distance of shortest time route
W = list2adj([IJD12; IJD22]);
D = zeros(n);  
for i = 1:n
   for j = 1:n
      D(i,j) = locTC(pred2path(P,i,j),W);
   end
end

% Find shortest path from Raleigh (node 21) to Virginia (node 33)
idx1 = 21; idx2 = 33; 

[t,p] = dijk(list2adj([IJT12; IJT22]),idx1,idx2);
h = [h ...
   pplot({p},[XY1;XY2],'y-','LineWidth',2,'DisplayName','Shortest Path')];
pplot(XY1([idx1 idx2],:),Name([idx1 idx2]))
title(sprintf(...
   'From %s to %s: Distance %.2f mi, Time = %d hr %d min',...
   Name{idx1},Name{idx2},D(idx1,idx2),floor(t),round(60*(t-floor(t)))));
legend(h),shg


%% Q3 2-product, 3-stage, 26-week Input data

clear all
clc

K = [60 50;               % capacity for product g in stage m (ton)
     55 45;
     50 35];
T = 26, rng(196),
D = round([gamrnd(6,4,T,1) gamrnd(4,3,T,1)]);

Cp = [12  20;           % production cost of product g at stage m ($/ton)
      75 130;
      35  60];
h = 0.3/12
Ci = cumsum(Cp,1)*h       % inventory cost of product g for stage m ($/ton)

Cs = [400 600;         % stage-m product-g setup cost ($)
       90  110;
       50  60];
yinit = [0  0;           % initial product g inventory at stage m (ton)
         0  0;
         0 11];
yfinal = zeros(3,2);     % final product g inventory at stage m (ton)
k0 = [1 0;               % initial setup at stage m for product g
      1 0;
      1 0];
M = size(K,1);           % number of production stages = 3         
G = size(K,2);           % number of products produced = 2
T = size(D,1);

Cp = reshape(repmat(Cp,[T 1 1]),M,T,G)     % create M x T x G array (3-D)
Ci = reshape(repmat(Ci,[T+1 1 1]),M,T+1,G) % create M x (T+1) x G array
Ci(:,1,:) = 0   % intital inventory cost already accounted for last period
Cs = reshape(repmat(Cs,[T 1 1]),M,T,G)     % create M x T x G array
mp = Milp('PPlan');
mp.addobj('min',Cp,Ci,Cs,zeros(M,T,G))     % zeros(M,T,G) dummy array for k
for g = 1:G
   for t = 1:T
      for m = 1:M-1
         mp.addcstr({[1 -1],{[m m+1],t,g}},{[1 -1],{m,[t t+1],g}},0,0,'=',0)
      end
      mp.addcstr({M,t,g},{[1 -1],{M,[t t+1],g}},0,0,'=',D(t,g))
      for m = 1:M
         mp.addcstr({m,t,g},0,0,'<=',{K(m,g),{m,t,g}})
      end
   end
   for m = 1:M
      mp.addcstr(0,0,{-1,{m,1,g}},{m,1,g},'<=',k0(m,g))
      for t = 2:T
         mp.addcstr(0,0,{-1,{m,t,g}},{[1 -1],{m,[t t-1],g}},'<=',0)
      end
   end
end
for m = 1:M, for t = 1:T, mp.addcstr(0,0,0,{m,t,':'},'=',1), end, end
mp.addlb(0,horzcat(reshape(yinit,M,1,G),zeros(M,T-1,G),reshape(yfinal,M,1,G)),0,0)
mp.addub(Inf,horzcat(reshape(yinit,M,1,G),inf(M,T-1,G),reshape(yfinal,M,1,G)),1,1)
mp.addctype('C','C','B','B')
spy(mp.Model.A),shg

ilp = mp.milp2ilp;
[x,TC,exitflag,output] = intlinprog(ilp{:});
TC,output
x = mp.namesolution(x)

Fp = x.Cp;
Fi = x.Ci;
Fs = x.Cs;
Fk = x.arg4;
for g = 1:G
   mdisp(D(:,g)',[],[],['D' num2str(g)])           
   mdisp(Fp(:,:,g),[],[],['Fp' num2str(g)])
   mdisp(Fi(:,:,g),[],[],['Fi' num2str(g)])
   mdisp(Fs(:,:,g),[],[],['Fs' num2str(g)])
   mdisp(Fk(:,:,g),[],[],['Fk' num2str(g)])
end

%% Q4 
clear all
clc

% Input data
K = [26 15;
     26 15;
     26 15];               % capacity for product g in stage m (ton)
T = 10;
rng(286); 
D = round([gamrnd(6,4,T,1) gamrnd(4,3,T,1)])   % demand for product g in period t (ton)

Cp = [120 200;
        20 20;
        32 32];% production/transportation cost of product g at stage m ($/ton)

  
h = 0.15/12
Ci = cumsum(Cp,1)*h      % inventory cost of product g for stage m ($/ton)

Cs = [0 0;         % stage-m product-g setup cost ($)
      0 0;
      0 0];
yinit = [0  0;           % initial product g inventory at stage m (ton)
         0  0;
         0 5];
yfinal = [0 0;          % final product g inventory at stage m (ton)
          0 0;
          0 5];
k0 = [1 0;               % initial setup at stage m for product g
      1 0;
      1 0];
M = size(K,1);           % number of production stages = 3
T = size(D,1);           % number of periods of production = 10
G = size(K,2);           % number of products produced = 2

% Create MILP model
Cp = reshape(repmat(Cp,[T 1 1]),M,T,G)     % create M x T x G array (3-D)
Ci = reshape(repmat(Ci,[T+1 1 1]),M,T+1,G) % create M x (T+1) x G array
Ci(:,1,:) = 0   % intital inventory cost already accounted for last period
Cs = reshape(repmat(Cs,[T 1 1]),M,T,G)     % create M x T x G array
mp = Milp('PPlan');
mp.addobj('min',Cp,Ci,Cs,zeros(M,T,G))     % zeros(M,T,G) dummy array for k
for g = 1:G
   for t = 1:T
      for m = 1:M-1
         mp.addcstr({[1 -1],{[m m+1],t,g}},{[1 -1],{m,[t t+1],g}},0,0,'=',0)
      end
      mp.addcstr({M,t,g},{[1 -1],{M,[t t+1],g}},0,0,'=',D(t,g))
      for m = 1:M
         mp.addcstr({m,t,g},0,0,'<=',{K(m,g),{m,t,g}})
      end
   end
   for m = 1:M
      mp.addcstr(0,0,{-1,{m,1,g}},{m,1,g},'<=',k0(m,g))
      for t = 2:T
         mp.addcstr(0,0,{-1,{m,t,g}},{[1 -1],{m,[t t-1],g}},'<=',0)
      end
   end
end
%for m = 1:M, for t = 1:T, mp.addcstr(0,0,0,{m,t,':'},'=',1), end, end
mp.addlb(0,horzcat(reshape(yinit,M,1,G),zeros(M,T-1,G),reshape(yfinal,M,1,G)),0,0)
mp.addub(Inf,horzcat(reshape(yinit,M,1,G),inf(M,T-1,G),reshape(yfinal,M,1,G)),1,1)
mp.addctype('C','C','B','B')

% Display constraint matrix
spy(mp.Model.A),shg

% Solve using INTLINPROG (k,z are binary variables)
ilp = mp.milp2ilp;
[x,TC,exitflag,output] = intlinprog(ilp{:});
TC,output
x = mp.namesolution(x)

Fp = x.Cp;
Fi = x.Ci;
Fs = x.Cs;
Fk = x.arg4;
for g = 1:G
   mdisp(D(:,g)',[],[],['D' num2str(g)])           
   mdisp(Fp(:,:,g),[],[],['Fp' num2str(g)])
   mdisp(Fi(:,:,g),[],[],['Fi' num2str(g)])
   mdisp(Fs(:,:,g),[],[],['Fs' num2str(g)])
   mdisp(Fk(:,:,g),[],[],['Fk' num2str(g)])
end

