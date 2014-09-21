function [xyCoordinates,arcLengths] = discretizePoly(poly,desiredLength)
%Returns the coordinates of the discretized polynomial. Discretization is
%500/desiredLength steps. xyCoordinates are the x and y coordinates of the discretized
%polynomial. arcLengths contains:
%   Column 1: Distance between xyCoordinates(index+1) and
%   xyCoordinates(index) at arcLengths (index,1)
%        ex: arcLengths(21,1) is the distance between xyCoordinates(21) and
%        xyCoordinates(22)
%
%   Column 2: Cumulative distance, or arc length, at xyCoordinate(index+1)
%   is given by arcLength(index,2)
%        ex: arcLengths(21,2) is the arc length from xyCoordinates(0) to
%        xyCoordinates(22)

    polynomial = poly;
    xVal = zeros(10000,1);
    yVal = zeros(10000,1);
    distance = zeros(9999,1);
    arcLength = zeros(9999,1);

    xVal(1) = 0;
    yVal(1) = polyval(polynomial,xVal(1));
    xVal(2) = 2*desiredLength/500;
    yVal(2) = polyval(polynomial,xVal(2));

    distance(1) = pdist([xVal(2),yVal(2);xVal(1),yVal(1)]);
    arcLength(1) = distance(1);

    i = 3;
    while arcLength(i-2) < desiredLength
        xVal(i) = i*desiredLength/500;
        yVal(i) = polyval(polynomial,xVal(i));
        distance(i-1) = pdist([xVal(i),yVal(i);xVal(i-1),yVal(i-1)]);
        arcLength(i-1) = arcLength(i-2) + distance(i-1);
        i = i+1;
    end
    xVal = xVal(1:i-1);
    yVal = yVal(1:i-1);
    distance = distance(1:i-2);
    arcLength = arcLength(1:i-2);
    xyCoordinates = [xVal,yVal];
    arcLengths = [distance,arcLength];
    
    
    %WHY DOES THIS RETURN x ~~ desired length? this makes no sense
end
