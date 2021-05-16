function g0 = gv(lat)

g0 = 9.780327 .* (1 + 0.0052792 .* sind(lat) .* sind(lat) ...
    + 0.0000232 .* sind(lat) .* sind(lat) .* sind(lat) .* sind(lat) );

end