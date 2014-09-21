classdef ProteinSummary < handle
    %A ProteinSummary object contains information for a single
    %Rotary-Shadowed EM image and method for analyzing the Protein objects,
    %in this case Tropomyosin objects, contained.

    
    properties
        Proteins
        Image
        Scale %in nm/pixel
        Dims
    end
    
    methods
        %constructor takes image scale, image dimensions, and protein as
        %arguments, respectively. Constructed with an empty Proteins field
        %unless specified.
        function obj = ProteinSummary(varargin)
            if nargin == 4
                obj.Proteins = [varargin{4}];
                obj.Scale = varargin{2};
                obj.Dims = varargin{3};
                obj.Image = varargin{1};
            else
                obj.Proteins = [];
                obj.Scale = varargin{2};
                obj.Dims = varargin{3};
                obj.Image = varargin{1};
            end
        end

        %Adds a Tropomyosin object to the Proteins field
        %
        %input: Tropomyosin object
        %output: None
        function obj = addProtein(obj,prot)
            if isempty(obj.Proteins) == 1 
               obj.Proteins = prot;
            else
               obj.Proteins = [obj.Proteins prot];
            end
        end
        
        %Returns the average end to end length of all Tropomyosin objects
        %
        %input: None
        %output: Double
        function avgEndToEnd = getAvgEndToEnd(obj)
            numProts = length(obj.Proteins);
            lengthTotal = 0;
            for i = 1:numProts
                lengthTotal = lengthTotal + obj.Proteins(i).getEndToEnd(obj.Scale);
            end
            avgEndToEnd = lengthTotal/numProts;
        end
        
        %Returns the average arc length of all Tropomyosin objects
        %
        %input: None
        %output: Double
        function avgArcLength = getAvgArcLength(obj)
            numProts = length(obj.Proteins);
            lengthTotal = 0;
            for i = 1:numProts
                lengthTotal = lengthTotal + obj.Proteins(i).getArcLength(obj.Scale);
            end
            avgArcLength = lengthTotal/numProts;
        end
        
        %Returns a matrix containing arc lengths in odd columns and 
        %cos(AngleDeviations) for that arc length in even columns for 
        %all proteins in ProteinSummary.
        %
        %input: None
        %output: Matrix[[Double], [Double]]
        function avgDeviations = getAvgDeviations(obj)
           avgDeviations = zeros(100,100);
           for protein = 2:2:length(obj.Proteins)*2
               deviation = obj.Proteins(protein/2).getAvgDeviation();
               arcLengths = deviation(:,1);
               cosAngleDeviations = deviation(:,2);
               vectorLength = length(arcLengths);
               avgDeviations(1:vectorLength,protein-1) = arcLengths;
               avgDeviations(1:vectorLength,protein) = cosAngleDeviations;
           end
        end
        
        %Displays all Tropomyosin Objects
        %
        %input: None
        %output: None (figure displayed)
        function show(obj)
           display = zeros(obj.Dims);
           for i = 1:length(obj.Proteins)
              display = display + obj.Proteins(i).getImage();
           end
           imshow(display)
        end
        
        %Returns the binary map of Tropomyosin objects
        %
        %input: None
        %output: Matrix
        function display = map(obj)
           display = zeros(obj.Dims);
           for i = 1:length(obj.Proteins)
              display = display + obj.Proteins(i).getImage();
           end
        end
        
        %Displays the skeleton for all Tropomyosin objects within
        %ProteinSummary over the original image.
        %
        %input: None
        %output: None (Displays Figure)
        function obj = visualProteinShow(obj)
           i = 1;
           figure(1)
           subplotRow = round(length(obj.Proteins)/5);
           hold on
           while i<length(obj.Proteins)
               prot = obj.Proteins(i);
               proteinImage = prot.getImage();
               [xCoordinates,yCoordinates] = ind2sub(size(proteinImage),find(proteinImage));
               minX = min(xCoordinates);
               minY = min(yCoordinates);
               maxX = max(xCoordinates);
               maxY = max(yCoordinates);
               baseImage = obj.Image;
               imageBorder = 10;
               resizedBase = baseImage(minX-imageBorder:maxX+imageBorder,minY-imageBorder:maxY+imageBorder);
               resizedProtein = proteinImage(minX-imageBorder:maxX+imageBorder,minY-imageBorder:maxY+imageBorder);
               subplot(subplotRow,5,i)
               imshow(highlight(resizedBase,resizedProtein));
               i = i+1;
           end
        end
        
        %To be implemented. Will display all proteins in ProteinSummary,
        %allowing each to be checked/unchecked. Unchecked Proteins will be
        %removed, allowing manual exclusion prior to analysis.
        %
        %IMPLEMENT UITABLE TO CHECK ONE PROTEIN AT A TIME. NEED TO FIGURE
        %OUT HOW TO ADD IMAGE TO UITABLE.
        function obj = visualProteinCheck(obj)
            f = figure('Position', [100 100 500 500]);
            table = uitable('Parent', f, 'Position', [25 25 400 400]);
            i = 1;
            figure(1)
            hold on
            while i<length(obj.Proteins)
                prot = obj.Proteins(i);
                proteinImage = prot.getImage();
                [xCoordinates,yCoordinates] = ind2sub(size(proteinImage),find(proteinImage));
                minX = min(xCoordinates);
                minY = min(yCoordinates);
                maxX = max(xCoordinates);
                maxY = max(yCoordinates);
                baseImage = obj.Image;
                imageBorder = 10;
                resizedBase = baseImage(minX-imageBorder:maxX+imageBorder,minY-imageBorder:maxY+imageBorder);
                resizedProtein = proteinImage(minX-imageBorder:maxX+imageBorder,minY-imageBorder:maxY+imageBorder);
                set(t, 'Data', [imshow(highlight(resizedBase,resizedProtein)), false]);
                imshow(highlight(resizedBase,resizedProtein));
                i = i+1;
            end
        end
        
        %For remove methods, add varargin for specific lengths
        
        %Removes proteins with an end to end length of 35nm or less
        %
        %input: None
        %output: None
        function obj = removeSmallProteins(obj)
           i = 1;
           while i<length(obj.Proteins)
               if obj.Proteins(i).getEndToEnd() < 35
                   obj.Proteins(i) = [];
               else
                   i = i+1;
               end
           end
        end
        
        %Removes proteins with an end to end length of 55nm or more
        %
        %input: None
        %output: None
        function obj = removeLargeProteins(obj)
           i = 1;
           while i<length(obj.Proteins)
               if obj.Proteins(i).getEndToEnd() > 45
                   obj.Proteins(i) = [];
               else
                   i = i+1;
               end
           end
        end
        
        % Currently, ProteinSummary will allow two proteins to be created
        % for one unique protein; that is, there will be one protein object
        % with endpoints [end1,end2] and one with [end2,end1], which are
        % identical. removeOverlaps removes such copies.
        %
        %input: None
        %output: None
        function obj = removeOverlaps(obj)
            i = 1;
            len = length(obj.Proteins);
            while i <= len
                j = i + 1;
                prot = obj.Proteins(i);
                while j <= len
                    
                    k = 1;
                    while k <= length(obj.Proteins(j).SkelCoordinates)
                        if prot.containsPoint(obj.Proteins(j).SkelCoordinates(k,1:2))
                            [i,j];
                            if isequal(obj.Proteins(i).SkelCoordinates, flipud(obj.Proteins(j).SkelCoordinates)) == 1
                                obj.Proteins(i) = [];
                                len = len - 1;
                            else
                                obj.Proteins(i) = [];
                                obj.Proteins(j - 1) = [];
                                len = len - 2;
                            end
                            k = 0;
                            j = i+1;
                            if j > length(obj.Proteins) 
                                break
                            end
                            
                            prot = obj.Proteins(i);
                        end
                        k = k+1;
                    end
                    j = j+1;
                end
                i = i+1;
            end
            
        end
        
        %Method no longer needed because of adjacents during Tropomyosin
        %Constructor. Deprecated.
        function obj = removeBranches(obj)
            prots = obj.map;
            branchPoints = bwmorph(prots, 'branchpoints');
            nonZeroCoordinates = find(branchPoints);
            [x,y] = ind2sub(size(branchPoints),nonZeroCoordinates);
            branchPointCoordinates = [x , y];
            i = 1;
            len = length(obj.Proteins);
            while i <= len
                j = 1;
                prot = obj.Proteins(i);
                while j <= length(x)
                    
                    if prot.containsPoint(branchPointCoordinates(j,1:2))
                        length(obj.Proteins)
                        len = len - 1;
                        obj.Proteins(i) = [];
                        j = 1;
                        if i > length(obj.Proteins)
                            break
                        end
                        prot = obj.Proteins(i);
                    end
                    j = j+1;
                end
                i = i+1;
            end
            
        end
        
        %Returns counts from histogram of Image in ProteinSummary.
        %
        %input: None
        %output: [[Double], [Double]]
        function binCounts = returnHistogram(obj)
            binCounts = histc(obj.Image,0:256);
            
        end
        
    end
    
end

