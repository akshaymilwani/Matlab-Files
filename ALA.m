

%% Q2a 
city2lonlat = @(city,st) uscity('XY',mand(city,uscity('Name'),st,uscity('ST')));
X1 = city2lonlat({'Gainesville','Warren'},{'GA','MI'});
X2 = uszip5('XY',mand([10020 17112 27606 32606 48234 56123],uszip5('Code5')));

P = [X1;X2];
w = [2.232 0.484 1.2 4.8 4.2 1.8 3.6 3];

[xy,TC] = minisumloc(P,w,'mi')

s1 = uscity50k
s1.Name{argmin(dists(xy,s1.XY,'mi'))}  % Closest city

%% Q2b

W = [2.232 0.484 0 0 0 0 0 0; 0 0 1.2 4.8 4.2 0 0 0; 0 0 0 0 0 1.8 3.6 3], V = [0 3.4 2.8; 0 0 0; 0 0 0]
[X,TC1] = minisumloc(P,W,'mi',V)

s1 = uscity50k
s1.Name{argmin(dists([-80.1163   40.4321],s1.XY,'mi'))}  % Closest city to widget factory
s1.Name{argmin(dists([-77.0105   40.2024],s1.XY,'mi'))}  % Closest city to DC1
s1.Name{argmin(dists([-83.0395   42.4312],s1.XY,'mi'))}  % Closest city to DC2

%% Q2c 
TC
TC1

d = TC - TC1  % Difference in cost
p = d/TC*100  % Percentage change (decrease) in cost 

%% Q3b
P = [50 150 190 220 270 295 420]'
m = size(P,1)
w = ones(1,7)
X = [60 125 130]'
n = size(X,1)
p = 1
[X,TC] = ala(X,w,P,p) %Optimal locations are at X = 50,190,295 where TC = 220

clear all
clc

%% Q4a
s = uscity;
s1 = uscity10k;
ST = uscity10k('ST');

city2lonlat = @(city,st) uscity10k('XY',mand(city,uscity10k('Name'),st,uscity10k('ST')));
Q = city2lonlat({'Richmond'},{'VA'})

s2 = uscity10k(strcmp('NC',ST) | strcmp('SC',ST) | strcmp('VA',ST) & s1.XY(:,2) <= Q(:,2))

city2lonlat1 = @(city,st) uscity10k('XY',mand(city,uscity10k('Name'),st,uscity10k('ST')));
x1 = city2lonlat1({'Rocky Mount'}, {'NC'})

W = (s2.Pop)';
P = s2.XY;
TD = W*dists(x1,P,'mi')'

TC = 6700000
k = TC/TD

W1 = k*W;

[x,TC1] = minisumloc(P,W1,'mi')

c1 = TC - TC1    % Maximum expected reduction in annual outbound transportation costs if the DC could be re-located to any other location 

%% Q4b 

P = s2.XY;
W = (s2.Pop)';
k = TC/TD
W1 = k*W;
n = 2
p = 'mi'
[x2,TC2] = ala(randX(P,n),W1,P,p)
c2 = TC - TC2   %  Maximum expected reduction in annual outbound transportation costs if two DCs could be located anywhere and the existing DC would be closed 


%% Q4c

m = size(P,1);
TCh = @(W,X)sum(sum(W(2,:).*dists(X,P,p)))...
        + sum(sum(W(1,:).*dists(x1,P,p)))
n = 1; 
alloc_h = @(X) full(sparse(argmin(dists(X,P,p)),1:m,W,2,m))
loc_h = @(W,X0) fminsearch(@(X) TCh(W,X),X0);
rng(5638)
X0 = randX(P,n);
X2 = randX(P,n);

Xnew = [x1,X2]
TC3 = Inf;
done = false;
while ~done
   Xnew = [x1; X2]
   W2i = alloc_h(Xnew);     % allocate
   X2i = loc_h(W2i,X2);   % locate
   TC3i = TCh(W2i,X2i);
   if TC3i < TC3
      TC3 = TC3i; X2 = X2i; W2 = W2i;
   else
      done = true;
   end
end
X2
TC3 = k*TC3
c3 = TC3-TC2  % Increase in cost of (c) as compared to (b) is 2.1553e+05. 

%% Q5 
clear XY Area Pop ST
z = uszip3;
xycog = @(XY,w)[w(:)'*XY(:,1) w(:)'*XY(:,2)]/sum(w);
ST = unique(z.ST);  

for i = 1:length(ST)
    zi = uszip3(mor(ST(i),z.ST));
    XY(i,:) = xycog(zi.XY,zi.Pop);
    Area(i) = sum(zi.LandArea);
    Pop(i) = sum(zi.Pop);
end
save usstate ST XY Area Pop
makemap(XY)
pplot(XY,'r.')
