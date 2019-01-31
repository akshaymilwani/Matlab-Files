
%% Q1b: Solve using UFLADD in Matlog 
k = [30 30 30 30] % k is same for each NF
C = [0 92 30 46; 92 0 40 94; 30 40 0 18; 46 94 18 0]
[y,TC,X] = ufladd(k,C);
y,TC,mdisp(X)

%% Q1c: Solve using MILP
mp = Milp('UFL');
mp.addobj('min',k,C)
[n m] = size(C);
 for j = 1:m
   mp.addcstr(0,{':',j},'=',1)
 end
 for i = 1:n
   mp.addcstr({m,{i}},'>=',{i,':'})  
 end
 mp.addub(Inf,1)
 mp.addctype('B','C')
[x,TC,nevals,XFlg] = milplog(mp); TC;nevals;XFlg;
 x = mp.namesolution(x), xC = x.C
 TC = k*x.k' + sum(sum(C.*xC))

%% Q2
clear all
s=uszip3;

T = readtable('PopcoData.xlsx')
S = table2struct(T)
DC = [[S.lon]' [S.lat]'];
fDC = [S.demand]';
TPC = [S.prod_cost]';
TDC = [S.dist_cost]';

figure
plot(fDC,TPC,'r.'),shg

x = fDC;
y = TPC;
yest = @(x,p) p(1) + p(2)*x;
fh = @(p) sum((y - yest(x,p)).^2);
ab = fminsearch(fh,[0 1])
k = ab(1), cp = ab(2)
hold on, fplot(@(x) yest(x,ab),[0 max(x)],'k-'), shg, hold off

[P,q,a] = uszip3('XY','Pop','LandArea', uszip3('Pop') > 0);
dafh = @(XY1,a1,XY2,a2) max(1.2*dists(XY1,XY2,'mi'),...
   0.675*max(sqrt(a1(:)),sqrt(a2(:)')));

idx = zeros(1,length(q));
for i = 1:length(q)
   [di,idxi] = min(dafh(P(i,:),a(i),DC,0));   % Actual DCs have 0 area
   dmax = 200;   % Max plant-to-customer distance for round-trip travel
   if di <= dmax
      idx(i) = idxi;
   end
end

P(idx==0,:) = [];
q(idx==0) = [];
a(idx==0) = [];
idx(idx==0) = [];

idx0 = setdiff(1:length(fDC),unique(idx));

D = dafh(DC,0,P,a);

f = zeros(length(q),1);
for i = 1:length(fDC)
   Mi{i} = find(argmin(D) == i);
   f(Mi{i}) = fDC(i)*(q(Mi{i})/sum(q(Mi{i})));
end

F = sparse(argmin(D),1:length(f),f);
r = sum(TDC)/sum(sum(F.*D))  
% per capita demand
pcd=s.Pop*(sum(fDC)/sum(q)); 
P=s.XY;
a=s.LandArea;
D = dafh([P; DC],0,P,a);   
C = r*(pcd(:)'.*D);

[y,TC,X] = ufl(k,C);
nNF = length(y) %(Total number of facilities required across continental US)
TDC_new = full(sum(sum(C.*X)))
TC_new = length(y)*k + TDC_new

NFsites = [P; DC];
makemap(NFsites);

pplot(NFsites(y,:),'g.')

%% 3a
clear all
x = uszip5(mor({'NC'},uszip5('ST')) & uszip5('Pop') > 10000)
P = x.XY;
a = x.LandArea;
dafh = @(XY1,a1,XY2,a2) max(1.2*dists(XY1,XY2,'mi'),...
   0.675*max(sqrt(a1(:)),sqrt(a2(:)')));
Da = dafh(P,a,P,a);  
f = x.Pop';          
r = 1/10000;         
C = r*Da.*f;         
[y,TC] = pmedian(10,C)

%% 3b
k = 0;                
mp = Milp('UFL')
mp.Model
[n m] = size(C)
kn = iff(isscalar(k),repmat(k,1,n),k(:)');  
mp.addobj('min',kn,C) 
n + n*m
for j = 1:m
   mp.addcstr(0,{':',j},'=',1)   
end
for i = 1:n
   mp.addcstr({m,{i}},'>=',{i,':'})  
end
mp.addub(1,1)
mp.addctype('B','C')         
mp.Model

mp.Model.name = 'pMedian'
mp.addcstr(1,0,'=',10)        
mp.Model

clear model params
model = mp.milp2gb;
params.outputflag = 1;
result = gurobi(model,params);
x = result.x;
TC = result.objval;

x = mp.namesolution(x)
TC
TCcalc = sum(kn.*x.kn) + full(sum(sum(C.*x.C)))
idxNF = find(x.kn)
nNF = sum(x.kn)


%% 4
clear all
X = xlsread('HW5data.xlsx','Clinics');
DC = X(:,[3 2]);
pplot(DC,'r.');
D = dists(DC,DC,'mi') * 1.2;
A = false(size(D));
A(D < 10) = true;
c = ones(1,size(A,2));
mp = Milp('Set Covering')
mp.addobj('min',c)
mp.addcstr(A,'>=',1)
mp.addctype('B')
spy(A)
clear model params
model = mp.milp2gb;
params.outputflag = 1;
result = gurobi(model,params);
x = result.x;
nNF = sum(x)