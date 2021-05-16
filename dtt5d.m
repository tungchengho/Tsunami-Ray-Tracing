% This script calculates tsunami travel time from a point to a 5x5 grid.
% 
% Input:
% bathy: 5x5 matrix of bathymetry, positive: depth
% xgd: 5x5 matrix of x-coordinate
% ygd: 5x5 matrix of y-coordinate
% arrt_in: input grid of arrival time.

function Tart5 = dtt5d(bathy, xgd, ygd)
% dtt5.m
% g0 = 9.81;
% g0 -> gv(ygd(cy));
% g0 -> gv(ygd(cj));

Tart5 = inf*ones(5,5);  % travel time
Tart5(3,3) = 0;

% xgd->x_dist
% ygd->y_dist
x_dist = zeros(5,5);
for jj = 1:5
    x_dist(jj,:) = xgd*2*pi*6371.008/360.*cosd(ygd(jj))*1000;
end
y_dist = (ygd*2*pi*6371.008/360)'*1000;

% Center point
cy = 3; cx = 3;

% Slowness
slo = zeros(5,5);
for jj = 1:5
    slo(jj,:) = 1./sqrt(gv(ygd(jj)).*bathy(jj,:));
end

%% For center 9 elements
%  o o o o o
%  o • • • o
%  o • @ • o
%  o • • • o
%  o o o o o
% 

Tart5(2:4,2:4) = dtt3s( slo(2:4,2:4), xgd(2:4), ygd(2:4));

%% For L step elements in 4 quadrant
% @: center point, •: L-corner grid point
% 8 points, each quadrant has 2 points: (p1y,p1x) and (p2y,p2x)
% R1_half, R2_half: half distance between @ and •
% dep_mid1, dep_mid2: depth at the middle between @ and •
%
%  o • o • o
%  • o o o •
%  o o @ o o
%  • o o o •
%  o • o • o
% 

for ofx = [1 -1]
    for ofy = [1 -1]
        % calculation
        p1y = cy + 2*ofy; p1x = cx +   ofx;
        p2y = cy +   ofy; p2x = cx + 2*ofx;
        R1_half = 0.5*deg2m( xgd(cx), ygd(cy), xgd(p1x), ygd(p1y));
        R2_half = 0.5*deg2m( xgd(cx), ygd(cy), xgd(p2x), ygd(p2y));
        v_p1_half = sqrt( gv(ygd(cy+ofy))*0.5*( bathy(cy+ofy,cx+ofx) + bathy(cy+ofy,cx    ) ) );
        v_p2_half = sqrt( gv(ygd(cy))*0.5*( bathy(cy+ofy,cx+ofx) + bathy(cy    ,cx+ofx) ) );
        Tart5(p1y,p1x) = 0.5*R1_half*slo(cy,cx) + R1_half./v_p1_half + 0.5*R1_half*slo(p1y,p1x);
        Tart5(p2y,p2x) = 0.5*R2_half*slo(cy,cx) + R2_half./v_p2_half + 0.5*R2_half*slo(p2y,p2x);
    end
end

%% For center 9 elements of 4 diagonal elements
% 
%  o o o o o
%  o o o o o
%  • • • o o
%  • @ • o o
%  • • • o o
% 
% cj = 2; ci = 2;
subi = 1:3; subj = 1:3; 
Tart5(subj,subi) = dtt3s( slo(subj,subi), xgd(subi), ygd(subj), Tart5(subj,subi));

% 
%  o o o o o
%  o o o o o
%  o o • • •
%  o o • @ •
%  o o • • •
% 
% cj = 2; ci = 4;
subi = 3:5; subj = 1:3; 
Tart5(subj,subi) = dtt3s( slo(subj,subi), xgd(subi), ygd(subj), Tart5(subj,subi));

% 
%  o o • • •
%  o o • @ •
%  o o • • •
%  o o o o o
%  o o o o o
% 
% cj = 4; ci = 4;
subi = 3:5; subj = 3:5; 
Tart5(subj,subi) = dtt3s( slo(subj,subi), xgd(subi), ygd(subj), Tart5(subj,subi));

% 
%  • • • o o
%  • @ • o o
%  • • • o o
%  o o o o o
%  o o o o o
% 
% cj = 4; ci = 2;
subi = 1:3; subj = 3:5; 
Tart5(subj,subi) = dtt3s( slo(subj,subi), xgd(subi), ygd(subj), Tart5(subj,subi));


%% For center 9 elements of 4 connecting elements
% 
%  o o o o o
%  o o o o o
%  o • • • o
%  o • @ • o
%  • • • • •
% 
% ## cj = 2; ci = 3; ##

% ===== 9 elements =====
subj = 1:3; subi = 2:4; 
Tart5(subj,subi) = dtt3s( slo(subj,subi), xgd(subi), ygd(subj), Tart5(subj,subi));

% ===== L-corner =====
% ----- Down-side ----- ofy = -1
cj = 2; ci = 3;
p1y = cj - 1; p1x = ci + 2;
p2y = cj - 1; p2x = ci - 2;
R1_half = 0.5*deg2m( xgd(ci), ygd(cj), xgd(p1x), ygd(p1y));
R2_half = 0.5*deg2m( xgd(ci), ygd(cj), xgd(p2x), ygd(p2y));
v_p1_half = sqrt( gv(ygd(cj))*0.5*( bathy(cj-1,ci+1) + bathy(cj,ci+1) ) );
v_p2_half = sqrt( gv(ygd(cj))*0.5*( bathy(cj-1,ci-1) + bathy(cj,ci-1) ) );
Tart5(p1y,p1x) = min([Tart5(p1y,p1x), Tart5(cj,ci) + 0.5*R1_half*slo(cj,ci) + R1_half./v_p1_half + 0.5*R1_half*slo(p1y,p1x)]);
Tart5(p2y,p2x) = min([Tart5(p2y,p2x), Tart5(cj,ci) + 0.5*R2_half*slo(cj,ci) + R2_half./v_p2_half + 0.5*R2_half*slo(p2y,p2x)]);

%
%  • • • • •
%  o • @ • o
%  o • • • o
%  o o o o o
%  o o o o o
% 
% ## cj = 4; ci = 3; ##

% ===== 9 elements =====
subj = 3:5; subi = 2:4; 
Tart5(subj,subi) = dtt3s( slo(subj,subi), xgd(subi), ygd(subj), Tart5(subj,subi));

% ===== L-corner =====
% ----- Up-side ----- ofy = 1;
cj = 4; ci = 3;
p1y = cj + 1; p1x = ci + 2;
p2y = cj + 1; p2x = ci - 2;
R1_half = 0.5*deg2m( xgd(ci), ygd(cj), xgd(p1x), ygd(p1y));
R2_half = 0.5*deg2m( xgd(ci), ygd(cj), xgd(p2x), ygd(p2y));
v_p1_half = sqrt( gv(ygd(cj))*0.5*( bathy(cj+1,ci+1) + bathy(cj,ci+1) ) );
v_p2_half = sqrt( gv(ygd(cj))*0.5*( bathy(cj+1,ci-1) + bathy(cj,ci-1) ) );
Tart5(p1y,p1x) = min([Tart5(p1y,p1x), Tart5(cj,ci) + 0.5*R1_half*slo(cj,ci) + R1_half./v_p1_half + 0.5*R1_half*slo(p1y,p1x)]);
Tart5(p2y,p2x) = min([Tart5(p2y,p2x), Tart5(cj,ci) + 0.5*R2_half*slo(cj,ci) + R2_half./v_p2_half + 0.5*R2_half*slo(p2y,p2x)]);

%
%  • o o o o
%  • • • o o
%  • @ • o o
%  • • • o o
%  • o o o o
% 
% ## cj = 3; ci = 2; ##

% ===== 9 elements =====
subj = 2:4; subi = 1:3; 
Tart5(subj,subi) = dtt3s( slo(subj,subi), xgd(subi), ygd(subj), Tart5(subj,subi));

% ===== L-corner =====
% ----- Left-side ----- ofx = -1
cj = 3; ci = 2;
p1y = cj + 2; p1x = ci - 1;
p2y = cj - 2; p2x = ci - 1;
R1_half = 0.5*deg2m( xgd(ci), ygd(cj), xgd(p1x), ygd(p1y));
R2_half = 0.5*deg2m( xgd(ci), ygd(cj), xgd(p2x), ygd(p2y));
v_p1_half = sqrt( gv(ygd(cj))*0.5*( bathy(cj+1,ci-1) + bathy(cj+1,ci) ) );
v_p2_half = sqrt( gv(ygd(cj))*0.5*( bathy(cj-1,ci-1) + bathy(cj-1,ci) ) );
Tart5(p1y,p1x) = min([Tart5(p1y,p1x), Tart5(cj,ci) + 0.5*R1_half*slo(cj,ci) + R1_half./v_p1_half + 0.5*R1_half*slo(p1y,p1x)]);
Tart5(p2y,p2x) = min([Tart5(p2y,p2x), Tart5(cj,ci) + 0.5*R2_half*slo(cj,ci) + R2_half./v_p2_half + 0.5*R2_half*slo(p2y,p2x)]);

%
%  o o o o •
%  o o • • •
%  o o • @ •
%  o o • • •
%  o o o o •
%
% ## cj = 3; ci = 4; ##

% ===== 9 elements =====
subj = 2:4; subi = 3:5; 
Tart5(subj,subi) = dtt3s( slo(subj,subi), xgd(subi), ygd(subj), Tart5(subj,subi));

% ===== L-corner =====
% ----- Right-side ----- ofx = 1;
cj = 3; ci = 4;
p1y = cj + 2; p1x = ci + 1;
p2y = cj - 2; p2x = ci + 1;
R1_half = 0.5*deg2m( xgd(ci), ygd(cj), xgd(p1x), ygd(p1y));
R2_half = 0.5*deg2m( xgd(ci), ygd(cj), xgd(p2x), ygd(p2y));
v_p1_half = sqrt( gv(ygd(cj))*0.5*( bathy(cj+1,ci+1) + bathy(cj+1,ci) ) );
v_p2_half = sqrt( gv(ygd(cj))*0.5*( bathy(cj-1,ci+1) + bathy(cj-1,ci) ) );
Tart5(p1y,p1x) = min([Tart5(p1y,p1x), Tart5(cj,ci) + 0.5*R1_half*slo(cj,ci) + R1_half./v_p1_half + 0.5*R1_half*slo(p1y,p1x)]);
Tart5(p2y,p2x) = min([Tart5(p2y,p2x), Tart5(cj,ci) + 0.5*R2_half*slo(cj,ci) + R2_half./v_p2_half + 0.5*R2_half*slo(p2y,p2x)]);

end     % end of function dtt5