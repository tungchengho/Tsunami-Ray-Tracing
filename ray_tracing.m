% This script traces the ray path 
% Input:
%     xgd: domain x coordinate
%     ygd: domain y coordinate
%     dep: water depth (>0: topography. <0:bathymetry)
%     T_A: tsunami travel time map from source
%     T_B: tsunami travel time map from receiver
%     case_name: filename of output files
% 
% Output
%     path_loc: [lon, lat]. ray path in Longitude and Latitude
%     path_dist: distance of each segment


function ray_tracing(xgd, ygd, dep, T_A, T_B, case_name)

%% Initialization
% Convert travel time from hour to second
T_A = T_A*3600;
T_B = T_B*3600;
% Find location (index) of Point 1
[epiy, epix] = find( T_A == min(T_A(:)) );
% Find location (index) of Point 2
[stn_y, stn_x] = find( T_B == min(T_B(:)) );

fprintf('    TB(A) = %10.5f, TA(B) = %10.5f\n', T_B(epiy(1),epix(1))/60, T_A(stn_y,stn_x)/60)

%% Calculating
fprintf('Now calculating the ray path\n')

% Trace the path from Point 2 to Point 1
% Get first point (Point 2)
cnt_x = stn_x; cnt_y = stn_y;
path_ind = [cnt_x, cnt_y];
path_dist = [0, 0, 0];
itr = 1;
T = T_A + T_B;    % T >= 0, 
rd = 2;
mk_in = false(5,5);
mk_in(2:4,2:4) = true;
mk_out = ~mk_in;
mk_in(3,3) = false;
while min((cnt_x - epix).^2) + min((cnt_y - epiy).^2) > 0 && itr <= size(dep, 1) + size(dep, 2)
    itr = itr + 1;
    itvx = cnt_x - rd : cnt_x + rd; itvy = cnt_y - rd : cnt_y + rd;
    Ts = T(itvy, itvx);  % subdomain of T

    % ets: subdomain of et
    T_1s = T_A(itvy, itvx);
    % etsc: time offsets relative to center point (idxx,idyy) of ets. >0: close to station
    T_1s0 = T_1s - T_1s(rd+1, rd+1);
    
%% Calculate travel time from center of current point
    % EsT: estimated travel time
    EsT = dtt5d( -dep(itvy,itvx) , xgd(itvx) , ygd(itvy) );
    D_arrt = EsT + T_1s0;   % D_arrt: difference of arrival time
    dtr = abs(D_arrt./EsT);   % dtr=0, wave direction
%     fprintf('%d ', 2*min(dtr(:))>0.1)
    Ts(~(dtr<=0.2)) = nan;

    % Find the subscripts of least dt
    [~,ind] = min(abs(Ts(:)));
    [yy,xx] = ind2sub(size(Ts),ind);

    % Replace current idxx and idyy with new one
    cnt_x = cnt_x-rd-1+xx;
    cnt_y = cnt_y-rd-1+yy;

    % Save location and distance
    path_ind(itr,:) = [cnt_x, cnt_y];
    path_dist(itr-1,:) = [ -dep(cnt_y,cnt_x), deg2m( xgd(path_ind(itr,1)), ygd(path_ind(itr,2)), ...
    xgd(path_ind(itr-1,1)), ygd(path_ind(itr-1,2))), T_A(cnt_y,cnt_x)/60];

end % while out of specific area

path_loc = [xgd(path_ind(:,1))', ygd(path_ind(:,2))'];

% Output path loc and distance
dlmwrite([ case_name '_nodes.txt'], path_loc, 'precision', '%10.3f', 'delimiter', ' ')
dlmwrite([ case_name '_dist.txt'], path_dist, 'precision', '%10.3f', 'delimiter', ' ')


% end
