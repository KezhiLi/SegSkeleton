function skeleton =  segWorm_ske(img, frame, isNormalized, verbose, varargin)
%SEGWORM Segment the worm in an image and organize the information in a
%   structure.
%
%   WORM = SEGWORM(IMG, FRAME, ISNORMALIZED, VERBOSE)
%
%   WORM = SEGWORM(IMG, FRAME, ISNORMALIZED, VERBOSE, NUMERODE, NUMDILATE)
%
%   WORM = SEGWORM(IMG, FRAME, ISNORMALIZED, VERBOSE, NUMERODE, NUMDILATE,
%                  SAMPLES, ISINTERP)
%
%   Inputs:
%       img          - the image to segment
%       frame        - the frame number (if the image comes from video)
%       isNormalized - is the image already normalized (i.e., all pixel
%                      values are between 0 to 1, inclusive)?
%       verbose      - verbose mode shows the results in a figure
%       numErode     - the number of time to erode the binary worm image;
%                      if empty or 0, the binary image is left unchanged
%       numDilate    - the number of time to dilate the binary worm image;
%                      if empty or 0, the binary image is left unchanged
%       samples      - the number of samples to use in verbose mode;
%                      if empty, all the worm is used.
%       isInterp     - when downsampling, should we interpolate the missing
%                      data or copy it from the original worm;
%                      if empty, we interpolate the missing data.
%
%   Output:
%       worm - the worm information organized in a structure
%              This structure contains 8 sub-structures,
%              6 sub-sub-structures, and 4 sub-sub-sub-structures:
%
%              * Video *
%              video = {frame}
%
%              * Contour *
%              contour = {pixels, touchI, inI, outI, angles, headI, tailI}
%
%              * Skeleton *
%              skeleton = {pixels, touchI, inI, outI, inOutI, angles,
%                          length, chainCodeLengths, widths}
%
%              Note: positive skeleton angles bulge towards the side
%              clockwise from the worm's head (unless the worm is flipped).
%
%              * Head *
%              head = {bounds, pixels, area,
%                      cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%              head.bounds{contour.left (indices for [start end]),
%                          contour.right (indices for [start end]),
%                          skeleton indices for [start end]}
%
%              * Tail *
%              tail = {bounds, pixels, area,
%                      cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%              tail.bounds{contour.left (indices for [start end]),
%                          contour.right (indices for [start end]),
%                          skeleton indices for [start end]}
%
%              * Left Side (Counter Clockwise from the Head) *
%              left = {bounds, pixels, area,
%                      cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%              left.bounds{contour (indices for [start end]),
%                          skeleton (indices for [start end])}
%
%              * Right Side (Clockwise from the Head) *
%              right = {bounds, pixels, area,
%                       cdf (at [2.5% 25% 50% 75% 97.5%]), stdev}
%              right.bounds{contour (indices for [start end]),
%                           skeleton (indices for [start end])}
%
%              * Orientation *
%              orientation = {head, vulva}
%              orientation.head = {isFlipped,
%                                  confidence.head, confidence.tail}
%              orientation.vulva = {isClockwiseFromHead,
%                                  confidence.vulva, confidence.nonVulva}
%
%       errNum - the error number if segmentation failed
%                (see also WORMFRAMEANNOTATION)
%       errMsg - the error message if segmentation failed
%
%   See also WORM2STRUCT, NORMWORMS
%
%
% � Medical Research Council 2012
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% Are we eroding and/or dilating the worm?
numErode = [];
numDilate = [];
if ~isempty(varargin)
    numErode = varargin{1};
end
if length(varargin) > 1
    numDilate = varargin{2};
end

% Convert the image to grayscale.
worm = [];
vWorm = []; % the verbose-mode worm
if (size(img,3) == 3)
    img = rgb2gray(img);
end

% Store the original then binarize the image.
oImg = img;
img = otsuImg(img, isNormalized);
%img = ~img;

% Erode and dilate the binary image.
if ~isempty(numErode) && numErode > 0
    img = imerode(img, strel('disk', numErode));
end
if ~isempty(numDilate) && numDilate > 0
    img = imdilate(img, strel('disk', numDilate));
end

% Find the worm.
errNum = [];
errMsg = [];
cc = bwconncomp(img);
wormPixels = [];
if ~isempty(cc.PixelIdxList)
    maxCCIdx = 0;
    maxCCSize = 0;
    for i = 1:length(cc.PixelIdxList)
        ccSize = length(cc.PixelIdxList{i});
        if ccSize > maxCCSize
            maxCCSize = ccSize;
            maxCCIdx = i;
        end
    end
    wormPixels = cc.PixelIdxList{maxCCIdx};
end

% No worm found.
if isempty(wormPixels)
    errNum = 101;
    errMsg = 'No worm was found.';
    
    % Show the failure.
    if verbose
        warning('segWorm:NoWormFound', ['Frame %d: ' errMsg], frame);

        % Open a big figure.
        figure('OuterPosition', [50 50 1280 960]);
        set(gcf, 'Color', [1 .5 .5]);
        
        % Show the original image.
        imshow(oImg);
        title('Original Image');
    end
    return;
end

% Find a point on the contour.
[y x] = ind2sub(size(img), min(wormPixels));

% Trace the contour clockwise.
contour = bwClockTrace(img, [x y], true);

% The contour touches a boundary.
if min(contour(:,1)) == 1 || min(contour(:,2)) == 1 || ...
        max(contour(:,1)) == size(img, 1) || max(contour(:,2)) == size(img, 2)
    errNum = 102;
    errMsg = 'The worm contour touches the image boundary.';
    
    % Show the failure.
    if verbose
        warning('segWorm:ContourTouchesBoundary', ...
            ['Frame %d: ' errMsg], frame);
        
        % Open a big figure.
        figure('OuterPosition', [50 50 1280 960]);
        set(gcf, 'Color', [1 .5 .5]);
        
        % Show the original image.
        subplot(1,2,1), imshow(oImg);
        title('Original Image');
        
        % Show the thresholded image.
        subplot(1,2,2), imshow(img);
        title('Thresholded Image');
    end
    return;
end

% The worm is roughly divided into 24 segments of musculature (i.e., hinges
% that represent degrees of freedom) on each side. Therefore, 48 segments
% around a 2-D contour.
% Note: "In C. elegans the 95 rhomboid-shaped body wall muscle cells are
% arranged as staggered pairs in four longitudinal bundles located in four
% quadrants. Three of these bundles (DL, DR, VR) contain 24 cells each,
% whereas VL bundle contains 23 cells." - www.wormatlas.org
sWormSegs = 24;
cWormSegs = 2 * sWormSegs;

% The contour is too small.
if size(contour, 1) < cWormSegs
    errNum = 103;
    errMsg = 'The worm contour is too small.';
    
    % Show the failure.
    if verbose
        warning('segWorm:ContourTooSmall', ['Frame %d: ' errMsg], frame);
        
        % Open a big figure.
        figure('OuterPosition', [50 50 1280 960]);
        set(gcf, 'Color', [1 .5 .5]);
        
        % Show the original image.
        subplot(1,2,1), imshow(oImg);
        title('Original Image');
        
        % Show the thresholded image.
        subplot(1,2,2), imshow(img);
        title('Thresholded Image');
    end
    return;
end

% Clean up the worm's contour.
if verbose
    roughContour = contour;
end
contour = cleanWorm(contour, size(contour, 1) / cWormSegs);

% The contour is too small.
if size(contour, 1) < cWormSegs
    errNum = 103;
    errMsg = 'The worm contour is too small.';
    
    % Show the failure.
    if verbose
        warning('segWorm:ContourTooSmall', ['Frame %d: ' errMsg], frame);
        
        % Open a big figure.
        figure('OuterPosition', [50 50 1280 960]);
        set(gcf, 'Color', [1 .5 .5]);
        
        % Show the original image.
        subplot(1,2,1), imshow(oImg);
        title('Original Image');
        
        % Show the thresholded image.
        subplot(1,2,2), imshow(img);
        title('Thresholded Image');
    end
    return;
end

% Compute the contour's local high/low-frequency curvature.
% Note: worm body muscles are arranged and innervated as staggered pairs.
% Therefore, 2 segments have one theoretical degree of freedom (i.e. one
% approximation of a hinge). In the head, muscles are innervated
% individually. Therefore, we sample the worm head's curvature at twice the
% frequency of its body.
% Note 2: we ignore Nyquist sampling theorem (sampling at twice the
% frequency) since the worm's cuticle constrains its mobility and practical
% degrees of freedom.
cCCLengths = circComputeChainCodeLengths(contour);
wormSegLength = (cCCLengths(1) + cCCLengths(end)) / cWormSegs;
hfAngleEdgeLength = wormSegLength;
hfCAngles = circCurvature(contour, hfAngleEdgeLength, cCCLengths);
lfAngleEdgeLength = 2 * hfAngleEdgeLength;
lfCAngles = circCurvature(contour, lfAngleEdgeLength, cCCLengths);

% Blur the contour's local high-frequency curvature.
% Note: on a small scale, noise causes contour imperfections that shift an
% angle from its correct location. Therefore, blurring angles by averaging
% them with their neighbors can localize them better.
wormSegSize = size(contour, 1) / cWormSegs;
hfAngleEdgeSize = wormSegSize;
hfBlurSize = ceil(hfAngleEdgeSize / 2);
hfBlurWin(1:hfBlurSize) = 1 / hfBlurSize;
mhfCAngles = circConv(hfCAngles, hfBlurWin);

% Compute the contour's local high/low-frequency curvature maxima.
[mhfCMaxP mhfCMaxI] = maxPeaksCircDist(mhfCAngles, hfAngleEdgeLength, ...
    cCCLengths);
[lfCMaxP lfCMaxI] = maxPeaksCircDist(lfCAngles, lfAngleEdgeLength, ...
    cCCLengths);

% Are there too many possible head/tail points?
lfHT = lfCMaxP > 90;
lfHTSize = sum(lfHT);
if lfHTSize > 2
    errNum = 104;
    errMsg = ['The worm has 3 or more low-frequency sampled convexities' ...
        'sharper than 90 degrees (possible head/tail points).'];
    
    % Organize the available worm information.
    if verbose
        warning('segWorm:TooManyEnds', ['Frame %d: ' errMsg], frame);
        vWorm = worm2struct(frame, contour, [], [], [], lfCAngles, [], ...
            [], cCCLengths, [], [], [], [], [], [], [], [], [], [], [], ...
            [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], ...
            [], [], [], [], [], [], [], [], [], 0, [], [], 0, [], []);
    else
        return;
    end
end

% Are the head and tail on the outer contour?
mhfHT = mhfCMaxP > 60;
mhfHTSize = sum(mhfHT);
if mhfHTSize < 2
    errNum = 105;
    errMsg = ['The worm contour has less than 2 high-frequency sampled '...
        'convexities sharper than 60 degrees (the head and tail). ' ...
        'Therefore, the worm is coiled or obscured and cannot be segmented.'];
    
    % Organize the available worm information.
    if verbose
        warning('segWorm:TooFewEnds', ['Frame %d: ' errMsg], frame);
        vWorm = worm2struct(frame, contour, [], [], [], lfCAngles, [], ...
            [], cCCLengths, [], [], [], [], [], [], [], [], [], [], [], ...
            [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], ...
            [], [], [], [], [], [], [], [], [], 0, [], [], 0, [], []);
    else
        return;
    end

% The head and tail are on the outer contour.
else
    
    % The low-frequency sampling identified the head and tail. 
    if lfHTSize > 1
        
        % Find the head and tail convexities in the low-frequency sampling.
        % Note: the tail should have a sharper angle.
        lfHTI = lfCMaxI(lfHT);
        lfHTP = lfCMaxP(lfHT);
        if lfHTP(1) <= lfHTP(2)
            headI = lfHTI(1);
            tailI = lfHTI(2);
        else
            headI = lfHTI(2);
            tailI = lfHTI(1);
        end
        
        % Localize the head by finding its nearest, sharpest (but blurred),
        % high-frequency convexity.
        mhfHTI = mhfCMaxI(mhfHT);
        dhfHeadI = abs(cCCLengths(headI) - cCCLengths(mhfHTI));
        dhfHeadI = min(dhfHeadI, cCCLengths(end) - dhfHeadI);
        [~, hfHeadI] = min(dhfHeadI);
        headI = mhfHTI(hfHeadI);
        
        % Localize the tail by finding its nearest, sharpest (but blurred),
        % high-frequency convexity.
        dhfTailI = abs(cCCLengths(tailI) - cCCLengths(mhfHTI));
        dhfTailI = min(dhfTailI, cCCLengths(end) - dhfTailI);
        [~, hfTailI] = min(dhfTailI);
        tailI = mhfHTI(hfTailI);
        
    % The high-frequency sampling identifies the head and tail. 
    elseif mhfHTSize < 3
        
        % Find the head and tail convexities in the high-frequency sampling.
        % Note: the tail should have a sharper angle.
        mhfHTI = mhfCMaxI(mhfHT);
        mhfHTP = mhfCMaxP(mhfHT);
        if mhfHTP(1) <= mhfHTP(2)
            headI = mhfHTI(1);
            tailI = mhfHTI(2);
        else
            headI = mhfHTI(2);
            tailI = mhfHTI(1);
        end
        
    % The high-frequency sampling identifies several, potential heads/tails. 
    else
        
        % Initialize our head and tail choicse.
        mhfHTI = mhfCMaxI(mhfHT);
        mhfHTI1 = mhfHTI(1);
        mhfHTI2 = mhfHTI(2);
        
        % How far apart are the head and tail?
        dmhfHTI12 = abs(cCCLengths(mhfHTI(1)) - cCCLengths(mhfHTI(2)));
        dmhfHTI12 = min(dmhfHTI12, cCCLengths(end) - dmhfHTI12);
        
        % Search for the 2 sharp convexities that are furthest apart.
        for i = 1:(mhfHTSize - 1)
            for j = (i + 1):mhfHTSize
                
                % How far apart are these 2 convexities?
                dmhfHTIij = abs(cCCLengths(mhfHTI(i)) - ...
                    cCCLengths(mhfHTI(j)));
                dmhfHTIij = min(dmhfHTIij, cCCLengths(end) - dmhfHTIij);
                
                % These 2 convexities are better head and tail choices.
                if dmhfHTIij > dmhfHTI12
                    mhfHTI1 = mhfHTI(i);
                    mhfHTI2 = mhfHTI(j);
                    dmhfHTI12 = dmhfHTIij;
                end
            end
        end
        
        % Which convexity is the head and which is the tail?
        % Note: the tail should have a sharper angle.
        if mhfCAngles(mhfHTI1) < mhfCAngles(mhfHTI2)
            headI = mhfHTI1;
            tailI = mhfHTI2;
        else
            headI = mhfHTI2;
            tailI = mhfHTI1;
        end            
    end
    
    % Find the length of each side.
    if headI > tailI
        size1 = cCCLengths(headI) - cCCLengths(tailI);
        size2 = cCCLengths(end) - cCCLengths(headI) + cCCLengths(tailI);
    else
        size1 = cCCLengths(tailI) - cCCLengths(headI);
        size2 = cCCLengths(end) - cCCLengths(tailI) + cCCLengths(headI);
    end
    
    % Are the sides within 50% of each others size?
    % Note: if a worm's length from head to tail is at least twice larger
    % on one side (relative to the other), than the worm must be touching
    % itself.
    if min(size1, size2)/ max(size1, size2) <= .5
        errNum = 106;
        errMsg = ['The worm length, from head to tail, is more than ' ...
            'twice as large on one side than it is on the other. ' ...
            'Therefore, the worm is coiled or obscured and cannot be segmented.'];
        
        % Organize the available worm information.
        if verbose
            warning('segWorm:DoubleLengthSide', ['Frame %d: ' errMsg], frame);
            vWorm = worm2struct(frame, contour, [], [], [], lfCAngles, ...
                headI, tailI, cCCLengths, [], [], [], [], [], [], [], ...
                [], [], [], [], [], [], [], [], [], [], [], [], [], [], ...
                [], [], [], [], [], [], [], [], [], [], [], [], [], [], ...
                0, [], [], 0, [], []);
        else
            return;
        end
    end

%if 0
    % Orient the contour and angles at the maximum curvature (the head or tail).
    if headI > 1
        contour = [contour(headI:end,:); contour(1:(headI - 1),:)];
        cCCLengths = [cCCLengths(headI:end) - cCCLengths(headI - 1); ...
            cCCLengths(1:(headI - 1)) + ...
            (cCCLengths(end) - cCCLengths(headI - 1))];
        %hfCAngles = [hfCAngles(headI:end); hfCAngles(1:(headI - 1))];
        lfCAngles = [lfCAngles(headI:end); lfCAngles(1:(headI - 1))];
        lfCMaxI = lfCMaxI - headI + 1;
        wrap = lfCMaxI < 1;
        lfCMaxI(wrap) = lfCMaxI(wrap) + length(lfCAngles);
        tailI = tailI - headI + 1;
        headI = 1;
        if tailI < 1
            tailI = tailI + size(contour, 1);
        end
    end
%end

    % Compute the contour's local low-frequency curvature minima.
    [lfCMinP lfCMinI] = minPeaksCircDist(lfCAngles, lfAngleEdgeLength, ...
        cCCLengths);
    
    % Compute the worm's skeleton.
    [skeleton cWidths] = linearSkeleton(headI, tailI, lfCMinP, lfCMinI, ...
        lfCMaxP, lfCMaxI, contour, wormSegLength, cCCLengths);
    
    % Measure the skeleton's chain code length.
    sCCLengths = computeChainCodeLengths(skeleton);
    sLength = sCCLengths(end);
    
    % Compute the worm's head and tail (at this point, we cannot
    % distinguish between the two). The worm's head and tail occupy,
    % approximately, 4 muscle segments each, on the skeleton and either
    % side of the contour.
    % Note: "The first two muscle cells in the two ventral and two dorsal
    % rows [of the head] are smaller than their lateral counterparts,
    % giving a stagger to the packing of the two rows of cells in a
    % quadrant. The first four muscles in each quadrant are innervated
    % exclusively by motoneurons in the nerve ring. The second block of
    % four muscles is dually innervated, receiving synaptic input from
    % motoneurons in the nerve ring and the anterior ventral cord. The rest
    % of the muscles in the body are exclusively innervated by NMJs in the
    % dorsal and ventral cords." - The Structure of the Nervous System of
    % the Nematode C. elegans, on www.wormatlas.org
    htSSegLength = sCCLengths(end) * (4 / sWormSegs);
    [head hlcBounds hrcBounds hsBounds] = ...
        worm2poly(1, chainCodeLength2Index(htSSegLength, sCCLengths), ...
        skeleton, headI, tailI, contour, false, sCCLengths, cCCLengths);
    [tail tlcBounds trcBounds tsBounds] = ...
        worm2poly(size(skeleton, 1), ...
        chainCodeLength2Index(sCCLengths(end) - htSSegLength, sCCLengths), ...
        skeleton, headI, tailI, contour, false, sCCLengths, cCCLengths);

    % Compute the contour's local low-frequency curvature minima.
    [lfCMinP, lfCMinI] = minPeaksCircDist(lfCAngles, lfAngleEdgeLength, ...
        cCCLengths);

    % Is the worm coiled?
    % If there are no large concavities, the worm is not coiled.
    lfCBendI = lfCMinI(lfCMinP < -30);
    
end

end
