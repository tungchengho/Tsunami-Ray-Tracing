% This script calculates tsunami travel time for a 3x3 grid.

function Tart3 = dtt3s(slo, xgd, ygd, arrt_in)
if nargin < 4
    Tart3 = inf*ones(3,3);
    Tart3(2,2) = 0;
else
    Tart3 = arrt_in;
end

for ii = 1:3
    jj = 1;
    Tart3(jj,ii) = min([ Tart3(2,2) + deg2m( xgd(2), ygd(2), xgd(ii), ygd(jj) ) ...
                         * 0.5* (slo(2,2)+slo(jj,ii)) , Tart3(jj,ii) ]);
                     
    jj = 3;
    Tart3(jj,ii) = min([ Tart3(2,2) + deg2m( xgd(2), ygd(2), xgd(ii), ygd(jj) ) ...
                         * 0.5* (slo(2,2)+slo(jj,ii)) , Tart3(jj,ii) ]);
end

jj = 2; ii = 1;
Tart3(jj,ii) = min([ Tart3(2,2) + deg2m( xgd(2), ygd(2), xgd(ii), ygd(jj) ) ...
                         * 0.5* (slo(2,2)+slo(jj,ii)) , Tart3(jj,ii) ]);
                     
jj = 2; ii = 3;
Tart3(jj,ii) = min([ Tart3(2,2) + deg2m( xgd(2), ygd(2), xgd(ii), ygd(jj) ) ...
                         * 0.5* (slo(2,2)+slo(jj,ii)) , Tart3(jj,ii) ]);
end
