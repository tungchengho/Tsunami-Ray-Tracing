function dd=deg2m(x1,y1,x2,y2)

dr=acosd(cosd(90-y1).*cosd(90-y2) + sind(90-y1).*sind(90-y2).*cosd(x2-x1));
dd=dr.*pi.*6371.008./180*1000;


