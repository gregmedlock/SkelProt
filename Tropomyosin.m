classdef Tropomyosin<handle
    %Creates an object for an individual protein from ProteinSummary.
    %
    %Defines each tropomyosin object by starting with an image input and
    %a starting coordinate at either the end or beginning of the protein.
    %The constructor builds out from the starting coordinate by scanning
    %neighboring pixels in the binary image and saving each coordinate.
    %Analysis Coordinates are created so the original Skeleton Coordinates
    %can be preserved if modifications need to be made to finish analysis.
    
    properties (SetAccess = public)
        Img
        Scale
        Start
        End
        SkelCoordinates
        Length
        AnalysisCoordinates
        Polynomial
    end
    
    methods
        
        %Constructs a Tropomyosin object from an image using either the
        %first coordinate on either end of the protein.
        function obj = Tropomyosin(I, coordinate1, scale)
            obj.Scale = scale;
            obj.Img = I;
            obj.Start = coordinate1;
            point1 = obj.Start;
            obj.SkelCoordinates(1,:) = point1;
            point2 = [0, 0];
            
            %n and m describe the spatial arrangement of the 8 points
            %around the current point in (n, m) = (x, y) format
            n = [-1 -1 -1 0 0 1 1 1];
            m = [-1 0 1 -1 1 -1 0 1];
            %adjacents holds a 1 in the indices of n and m where a 
            %neighboring pixel is occupied in the n and m matrices
            adjacents = [0 0 0 0 0 0 0 0];
            i = 1;
            
            %Loops through each of the 8 neighboring pixels
            while i <= 8
                newPoint = NextPoint(point1, n(i),m(i));
                if I(newPoint(1),newPoint(2)) == 1 && isequal(newPoint,point2) == 0
                    skelSize = size(obj.SkelCoordinates);
                    skelRows = skelSize(1);
                    for j = 1:skelRows
                        skelPoint = obj.SkelCoordinates(j,1:2);
                        if newPoint(1) == skelPoint(1) && newPoint(2) == skelPoint(2)
                            same = 1;
                            break
                        else
                            same = 0;
                        end
                    end
                    
                    if same == 0

                        adjacents(i) = 1;
                        if nnz(adjacents) > 2
                            adjacents = [0 0 0 0 0 0 0 0];
                            error('Image is either branched or has width > 1')
                        end
                    end
                end
                
                %When all neighbors have been scanned, makes sure there are
                %at most 2 neighbors, then picks the closest neighbor as
                %the next pixel.
                if i == 8
                    if nnz(adjacents) > 2
                        adjacents = [0 0 0 0 0 0 0 0];
                        error('Image is either branched or has width > 1')
                    elseif nnz(adjacents) == 2 || nnz(adjacents) == 1
                        location = find(adjacents);
                        distances = zeros(1,length(location));
                        newPoint = zeros(length(location),2);
                        for q = 1:length(location)
                            index = location(q);
                            point = NextPoint(point1, n(index),m(index));
                            newPoint(q,1:2) = point;
                            distance = pdist([point1;newPoint(q,:)]);
                            distances(q) = distance;
                        end
                        [minimum,ind] = min(distances);
                        obj.SkelCoordinates = cat(1,obj.SkelCoordinates,newPoint(ind,:));
                        point2 = point1;
                        point1 = newPoint(ind,:);
                        i = 0;
                    end                   
                    adjacents = [0 0 0 0 0 0 0 0];
                end
                i = i + 1;   
            end
            obj.End = obj.SkelCoordinates(end,:);
            obj.Length = length(obj.SkelCoordinates);
            obj.AnalysisCoordinates = obj.SkelCoordinates;
            %make analysis coordinates anchored at the origin and
            %differentiable.
            obj.anchorCoordinates();
            obj.makeDifferentiable();
            obj.Polynomial = createPolyfit(obj);
        end
        
        %% get methods
        
        %returns the arc length of the Tropomyosin object
        function arcLength = getArcLength(obj)
            arcLength = 0;
            for i = 1:length(obj.SkelCoordinates) - 1
                point1 = obj.SkelCoordinates(i,1:2);
                point2 = obj.SkelCoordinates(i+1,1:2);
                arcLength = arcLength + obj.Scale*pdist([point1;point2]);
            end
        end
        
        %returns the end coordinate of the tropomyosin object
        function endPoint = getEnd(obj)
            endPoint = obj.End;
        end
        
        %returns the end-to-end length of the Tropomyosin object
        function endToEnd = getEndToEnd(obj)
            endToEnd = obj.Scale*pdist([obj.Start;obj.End]);
        end
        
        %return a binary image containing the tropomyosin object
        %coordinates
        function image = getImage(obj)
            canvas = zeros(size(obj.Img));
            for i = 1:length(obj.SkelCoordinates)
                canvas(obj.SkelCoordinates(i,1),obj.SkelCoordinates(i,2)) = 1;
            end
            image = canvas;
        end
        
        %returns the start coordinate of the tropomyosin object
        function startPoint = getStart(obj)
            startPoint = obj.Start;
        end
        
        %% Persistence Length Calculation
        
        %Standard step size = 1nm
        %in development, need to pretreat coordinates
        function arcStep = discretizedArcLengths(obj,stepSize,anchoredDims)
            polynomialFit = polyfit(anchoredDims(:,1),anchoredDims(:,2),5);
            xVal(1) = 0;
            yVal(1) = polyval(polynomial,xVal(1));
            xVal(2) = .01;
            yVal(2) = polyval(polynomial,xVal(2));
            
            distance(1) = pdist([xVal(2),yVal(2);xVal(1),yVal(1)]);
            arcLength(1) = distance(1);
            xRange = round(min(anchoredDims(:,1))):round(max(anchoredDims(:,1)));
            for i = 3:10000
                xVal(i) = i*.01;
                yVal(i) = polyval(polynomial,xVal(i));
                distance(i-1) = pdist([xVal(i),yVal(i);xVal(i-1),yVal(i-1)]);
                arcLength(i-1) = arcLength(i-2) + distance(i-1);
            end
        end
        
        %Calculate persistence length using the tangent correlation
%         function Lp = persistenceLengthTangent(obj, varargin)
%             avgCos = getAvgDeviation(obj);
%             if nargin == 2
%                 startLength = varargin{1};
%                 endLength = varargin{2};
%             else
%                 startLength = 0;
%                 endLength = 
%             end
%             arcLength = avgCos(:,1);
%             angleCos = avgCos(:,2);
%             logCos = log(angleCos);
%             lineFit = polyfit(arcLength,logCos,1);
%             slope = lineFit(1);
%             Lp = 1/(-2*slope);
%         end
        
            %Calculate persistence length using the tangent correlation
        function Lp = persistenceLengthTangent(obj)
            avgCos = getAvgDeviation(obj);
            arcLength = avgCos(:,1);
            angleCos = avgCos(:,2);
            logCos = log(angleCos);
            lineFit = polyfit(arcLength,logCos,1);
            slope = lineFit(1);
            Lp = 1/(-2*slope);
        end
        
        %Calculate persistence length using the second moment of tangent
        %angles
        function Lp = persistenceLengthMoment(obj)
            avgCos = getAvgDeviation(obj);
            arcLength = avgCos(:,1);
            angleCos = avgCos(:,2);
            Lp = polyfit(arcLength,acosd(angleCos).^2,1);
            Lp = 1/Lp(1);
        end
        

        function avgCos = getAvgDeviation(obj)
            valuesOfL = obj.Length;
            avgCos = zeros(valuesOfL*2,2);
            vectorGap = 1; %nm
            i = 1;
            while vectorGap < obj.getArcLength - 3
                angles = angleDeviation(obj,vectorGap);
                avgCos(i,1) = vectorGap;
                avgCos(i,2) = mean(cosd(angles(:,2)));
                vectorGap = vectorGap + 1; %add 1nm to vector spacing
                i = i+1;
            end
            %Need to get rid of any values at end that were not assigned
            avgCos = avgCos(1:i-1,1:2);
        end
        
        %Continuos deviation from polyfit of AnalysisCoordinates
        %
        %Returns column vector [s,angle] where s is the arc length of the
        %initial tangent vector for which angle is the angle of deviation with the
        %second tangent vector at arc length vectorGap.
        function angle = angleDeviation(obj,vectorGap)
            xVal1 = 0;
            s1 = 0;
            polynomialFit = polyfit(obj.AnalysisCoordinates(:,1),obj.AnalysisCoordinates(:,2),5);
            derivative = polyder(polynomialFit);
            tangent1 = polyval(derivative,xVal1);
            vector1 = [1,tangent1];
            magnitude1 = sqrt(vector1(1)^2 + vector1(2)^2);
            angle = zeros(3*length(obj.SkelCoordinates)-1,2);
            i = 0;
            s2 = vectorGap;
            desiredLength = obj.getArcLength();
            [xyCoordinates,arcLengths,xValAtDesiredLength] = constrainPolyLength(obj,desiredLength);
            xVal2 = getContourLengthXValue(obj,xyCoordinates,arcLengths,s2);
            while s2 < obj.getArcLength()
                i = i+1;
                tangent2 = polyval(derivative,xVal2);
                vector2 = [1,tangent2];
                magnitude2 = sqrt(vector2(1)^2 + vector2(2)^2);
                angle(i,2) = acosd(dot(vector1,vector2)/(magnitude1*magnitude2));
                angle(i,1) = s2;
                %advance vector 1, recalculate tangent & magnitude
                s2 = s2 + 1; % advance by 1nm
                s1 = s1 + 1; % advance by 1nm
                xVal2 = getContourLengthXValue(obj,xyCoordinates,arcLengths,s2);
                xVal1 = getContourLengthXValue(obj,xyCoordinates,arcLengths,s1);
                tangent1 = polyval(derivative,xVal1);
                vector1 = [1,tangent1];
                magnitude1 = sqrt(vector1(1)^2 + vector1(2)^2); 
            end
            angle = angle(1:i,1:2);
        end
        
        %Returns the x value at which a desired contour length
        %occurs in a set of skeleton coordinates. xyCoordinates and
        %arcLengths can be obtained with a call to constrainPolyLength
        function xVal = getContourLengthXValue(obj,xyCoordinates,arcLengths,desiredLength)
            i = 1;
            while desiredLength > arcLengths(i,2)
               i = i + 1;
               if i > length(arcLengths)
                   i = i - 1;
                   break
               end
            end
            xVal = xyCoordinates(i,1);
        end
        
        
        function closerX = getCloserValue(obj,fun,funHandle,derivative,derivative2Func,xVal,value,desiredValue)
            %derivativeofDerivative = fun(1:end-1).* (4:-1:1);
            %THIS DOESNT WORKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK       
            evaluatedDerivative = derivative2Func(xVal);
            %if sign(desiredValue) == sign(value)
            differenceConstant = .05;
            %else
            %    differenceConstant = 1
            %end
            
            %derivativeofDerivativeVal = polyval(derivativeOfDerivative,xVal)
            if desiredValue - value < 5 && ...
                    desiredValue - value > -5 || desiredValue - value > 1000
                closerX = xVal;
            elseif desiredValue - value > 0 && ...
                    evaluatedDerivative > 0%Contour length is too small, derivative is positive
                newX = xVal +differenceConstant*abs(xVal);
                length = integral(funHandle,0,newX);
                'value too small, derivative is positive';
                closerX = getCloserValue(obj,fun,funHandle,derivative,derivative2Func,newX,length,desiredValue);
            elseif desiredValue - value < 0 && ...
                    evaluatedDerivative > 0%Contour length is too large,derivative is positive
                newX = xVal -differenceConstant*abs(xVal);
                length = integral(funHandle,0,newX);
                'value too large, derivative is positive';
                closerX = getCloserValue(obj,fun,funHandle,derivative,derivative2Func,newX,length,desiredValue);
            elseif desiredValue - value > 0 && ...
                    evaluatedDerivative < 0%Contour length is too small, derivative is negative
                newX = xVal -differenceConstant*abs(xVal);
                length = integral(funHandle,0,newX);
                'value too small, derivative is negative';
                closerX = getCloserValue(obj,fun,funHandle,derivative,derivative2Func,newX,length,desiredValue);
            else %Contour length is too large,derivative is negative
                %Condition is desiredValue - value < 0 && ...
                %    polyval(derivative,xVal) < 0
                newX = xVal +differenceConstant*abs(xVal);
                length = integral(funHandle,0,newX);
                'value too large, derivative is negative';
                closerX = getCloserValue(obj,fun,funHandle,derivative,derivative2Func,newX,length,desiredValue);
                
            end
        end
        
        %% Helper & Test Functions
        
        %Anchor coordinates to origin to make polyfit more accurate
        %Need to adjust other functions to accommodate AnalysisCoordinates
        %now instead of returning new coordinates with anchorcoordinates
        function obj = anchorCoordinates(obj)
            for i = 1:length(obj.SkelCoordinates())
                coordinateX = obj.SkelCoordinates(i,1) - obj.SkelCoordinates(1,1);
                obj.AnalysisCoordinates(i,1) = coordinateX;
                coordinateY = obj.SkelCoordinates(i,2) - obj.SkelCoordinates(1,2);
                obj.AnalysisCoordinates(i,2) = coordinateY;
            end
        end
        
        %Rotates AnalysisCoordinates so that it is differentiable.
        %This function ensures that AnalysisCoordinate x values never
        %decrease along y (i.e. the skeleton passes the line test)
        function obj = makeDifferentiable(obj)
            rotationMatrix = [0 1;-1 0]; %rotates a coordinate pi/4 counterclockwise
            bool = obj.checkX();
            while bool == false
                for i = 1:length(obj.AnalysisCoordinates)
                    obj.AnalysisCoordinates(i,1:2) = obj.AnalysisCoordinates(i,1:2)*rotationMatrix;    
                end
                bool = obj.checkX();
                if nnz(obj.AnalysisCoordinates - (obj.SkelCoordinates - obj.SkelCoordinates(1,1:2))) == 0
                    bool = true;
                end
            end
        end
        
        %Returns a 5th order polynomial fit of the skeleton coordinates
        function poly = createPolyfit(obj)
            poly = polyfit(obj.AnalysisCoordinates(:,1),obj.AnalysisCoordinates(:,2),5);
        end
        
        %Returns x value for the polynomial when constrained to the length
        %of actual Tm in nm.
        function [xyCoordinates,arcLengths,xValAtDesiredLength] = constrainPolyLength(obj,desiredLength)
            [xyCoordinates,arcLengths] = discretizePoly(obj.Polynomial,desiredLength);
            i = 1;
            while arcLengths(i,2) < desiredLength
                i = i+1;
                if i > length(arcLengths)
                    i = i-1;
                    break
                end
            end
            xValAtDesiredLength = xyCoordinates(i,1);
            xyCoordinates = xyCoordinates(1:i,1:2);
        end
        
        
        %Returns True if x values never decrease thoughout a skeleton
        function bool = checkX(obj)
           bool = true;
           numPoints = length(obj.AnalysisCoordinates);
           point1 = obj.AnalysisCoordinates(1);
           for point = 2:numPoints
               point2 = obj.AnalysisCoordinates(point);
               if point2(1) < point1(1)
                   bool = false;
               end
               point1 = point2;
           end
        end
        
        %returns the arc length up to the index of coordinate s, for use in
        %determining persistence length as a function of s
        function arcLength = getArcLengthToPoint(obj,sValue)
           arcLength = 0;
           for i = 1:sValue - 1
                point1 = obj.SkelCoordinates(i,1:2);
                point2 = obj.SkelCoordinates(i+1,1:2);
                arcLength = arcLength + obj.Scale*pdist([point1;point2]);
            end
        end
        
        %returns 1 if the tropoyosin object contains point
        function bool = containsPoint(obj, point)
            bool = false;
            for i = 1:length(obj.SkelCoordinates)
                if obj.SkelCoordinates(i,1) == point(1,1) && obj.SkelCoordinates(i,2) == point(1,2)
                    bool = true;
                    [obj.SkelCoordinates(i,1),point(1,1);obj.SkelCoordinates(i,2),point(1,2)];
                    break
                end
            end
        end
        
        %function endToEndLp = getEndToEndLp(obj,scale)
        %    contourLength = getArcLength(obj,scale);
        %    endToEndLength = getEndToEnd(obj,scale);
        %    
        %    endToEndLp =  
        %end
        

        %% Discrete Persistence Length (DEPRECATED, USE CONTINUOUS)
        function persistenceLength = getPersistenceLength(obj,pointGap)
            angleDeviations = getAllAnglesContinuous(obj,pointGap);
            arcLength = angleDeviations(:,1);
            angleCos = angleDeviations(:,2);
            logCos = log(angleCos);
            lineFit = polyfit(arcLength,logCos,1);
            slope = lineFit(1);
            persistenceLength = -2*slope;
        end
        
        %returns a vector of angles of deviation for <cos(theta)>
        %in: pointGap: Distance in Pixels between the points that make up
        %each vector
        %in: vectorGap: Distance in pixels between the 2nd point in vector
        %1 and the first point in vector 2.
        function cosDeviation = angleDeviationDiscrete(obj,pointGap,vectorGap)
            point1 = obj.Start;
            point2 = obj.SkelCoordinates(2,1:2);
            vector1 = [point1-point1;point2-point1]; %anchor to origin
            mag1 = sqrt(vector1(2,1)^2 + vector1(2,2)^2);
            angle = zeros(length(obj.SkelCoordinates-1),1);
            for i = 1:(pointGap+1):length(obj.SkelCoordinates) - (vectorGap+2)
                point3 = obj.SkelCoordinates(i+vectorGap,1:2);
                point4 = obj.SkelCoordinates(i+vectorGap+pointGap,1:2);
                vector2 = [point3-point3;point4-point3];
                mag2 = sqrt(vector2(2,1)^2 + vector2(2,2)^2);
                %angle given by dot product
                angle(i) = acosd(dot(vector1(2,1:2),vector2(2,1:2))/(mag1*mag2));
                %advance vector 1 by two pixels, recalc magnitude
                point1 = obj.SkelCoordinates(i+2,1:2);
                point2 = obj.SkelCoordinates(i+3,1:2);
                vector1 = [point1-point1;point2-point1];
                mag1 = sqrt(vector1(2,1)^2 + vector1(2,2)^2); 
            end
            cosDeviation = angle;
        end
        
        %Need to calculate cos(theta) as a function of vectorGap
        %NEED TO HANDLE AND NORMALIZE BY CONTOUR LENGTH
        function avgCos = getAllAngles(obj,pointGap)
            %vectorGap is a minimum of 1 and must have room at end for 2
            %points and any pointGap space required
            valuesOfL = length(obj.SkelCoordinates)-(pointGap+2);
            avgCos = zeros(valuesOfL,2);
            for i = 1:valuesOfL
                angles = angleDeviationDiscrete(obj,pointGap,i);
                vectorGapLength = obj.Scale*i;
                avgCos(i,1) = mean(cos(angles));
                avgCos(i,2) = vectorGapLength;
            end
        end
        
    end
    
end

