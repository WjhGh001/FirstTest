%参见链接(里面有相应翻译)：https://ww2.mathworks.cn/matlabcentral/answers/
%156301-why-i-got-error-when-executing-below-code-what-can-i-do

% SLIC Simple Linear Iterative Clustering SuperPixels
%
% Implementation of Achanta, Shaji, Smith, Lucchi, Fua and Susstrunk's
% SLIC Superpixels
%
% Usage:   [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw)
%
% Arguments:  im - Image to be segmented.
%              k - Number of desired superpixels. Note that this is nominal
%                  the actual number of superpixels generated will generally
%                  be a bit larger, espiecially if parameter m is small.
%                  所需 superpixels 的数量。请注意, 这是名义上生成的数量，实际
%                  数量一般会更大一点, 尤其是当参数 m 小的时候。
%              m - Weighting factor between colour and spatial
%                  differences. Values from about 5 to 40 are useful.  Use a
%                  large value to enforce superpixels with more regular and
%                  smoother shapes. Try a value of 10 to start with.
%                  颜色和空间差异之间的加权因子。在CIELAB空间，范围在[1,40].
%                  从大约5到40是有用的。 使用大值来强制 superpixels 的形状更加
%                  整齐和平滑(使用较大的值可以使超级像素具有更规则和更平滑的形状)。
%                  请尝试使用值10开始。
%       seRadius - Regions morphologically smaller than this are merged with
%                  adjacent regions. Try a value of 1 or 1.5.  Use 0 to
%                  disable.
%                  在形态上小于此的区域与相邻区域合并。尝试值1或1.5。使用0禁用。
%         colopt - String 'mean' or 'median' indicating how the cluster
%                  colour centre should be computed. Defaults to 'mean'。
%                  字符串 "平均值" 或 "中间值", 指示如何计算簇颜色中心。默认为 "平均值"
%             mw - Optional median filtering window size.  Image compression
%                  can result in noticeable artifacts in the a*b* components
%                  of the image.  Median filtering can reduce this. mw can be
%                  a single value in which case the same median filtering is
%                  applied to each L* a* and b* components.  Alternatively it
%                  can be a 2-vector where mw(1) specifies the median
%                  filtering window to be applied to L* and mw(2) is the
%                  median filtering window to be applied to a* and b*.
%                  可选的中值过滤窗口大小。图像压缩可能会在a * b *组件中导致明显的伪像
%                  mw可以是单个值，在这种情况下，将相同的中值滤波应用于每个L * a *和b *分量。
%                  或者，它可以是2维向量，其中mw（1）指定要应用于L *的中值滤波窗口，
%                  而mw（2）是要应用于a *和b *的中值滤波窗口。

% Returns:     l - Labeled image of superpixels. Labels range from 1 to k.
%                  标记的超像素图像。标签范围从1到k。
%             Am - Adjacency matrix of segments.  Am(i, j) indicates whether
%                  segments labeled i and j are connected/adjacent
%                  段的邻接矩阵。 Am（i，j）指示标记为i和j的段是否已连接/相邻
%             Sp - Superpixel attribute structure array with fields:
%                   L  - Mean L* value
%                   a  - Mean a* value
%                   b  - Mean b* value
%                   r  - Mean row value
%                   c  - Mean column value
%                   stdL  - Standard deviation of L* 标准偏差 L *
%                   stda  - Standard deviation of a* 标准偏差 a *
%                   stdb  - Standard deviation of b* 标准偏差 b *
%                   N - Number of pixels
%                   edges - List of edge numbers that bound each
%                           superpixel. This field is allocated, but not set,
%                           by SLIC. Use SPEDGES(spedges) for this.
%              d - Distance image giving the distance each pixel is from its
%                  associated superpixel centre.
%
% It is suggested that use of this function is followed by SPDBSCAN to perform a
% DBSCAN clustering of superpixels.  This results in a simple and fast
% segmentation of an image.
% 
% Minor variations from the original algorithm as defined in Achanta et al's
% paper:
% 
% - SuperPixel centres are initialised on a hexagonal grid rather than a square
%   one. This results in a segmentation that will be nominally 6-connected
%   which hopefully facilitates any subsequent post-processing that seeks to
%   merge superpixels.
% - Initial cluster positions are not shifted to point of lowest gradient
%   within a 3x3 neighbourhood because this will be rendered irrelevant the
%   first time cluster centres are updated.
% 
% Reference: R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua and
% S. Susstrunk. "SLIC Superpixels Compared to State-of-the-Art Superpixel
% Methods"  PAMI. Vol 34 No 11.  November 2012. pp 2274-2281.
% 
% See also: SPDBSCAN, MCLEANUPREGIONS, REGIONADJACENCY, DRAWREGIONBOUNDARIES, RGB2LAB
% 
% Copyright (c) 2013 Peter Kovesi
% Centre for Exploration Targeting
% School of Earth and Environment
% The University of Western Australia
% peter.kovesi at uwa edu au
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
% 
% The Software is provided "as is", without warranty of any kind.
% 
% Feb  2013
% July 2013 Super pixel attributes returned as a structure array
% 
% Note that most of the computation time is not in the clustering, but rather
% in the region cleanup process.


function [l, Am, Sp, d] = slic(im, k, m, seRadius, colopt, mw, nItr, eim, We)
%% 
    if ~exist('colopt','var') || isempty(colopt), colopt = 'mean'; end  %设置默认值，~exist('colopt','var')是检测colopt变量是否存在
    if ~exist('mw','var')     || isempty(mw),         mw = 0;      end
    if ~exist('nItr','var')   || isempty(nItr),     nItr = 10;     end
    
    if exist('eim', 'var'), USEDIST = 0; else, USEDIST = 1; end
        
    MEANCENTRE = 1;
    MEDIANCENTRE = 2;
    
    if strcmp(colopt, 'mean')%比较函数
        centre = MEANCENTRE;
    elseif strcmp(colopt, 'median')
        centre = MEDIANCENTRE;        
    else
        error('Invalid colour centre computation option');
    end
    
    [rows, cols, chan] = size(im);
    if chan ~= 3
        error('Image must be colour');
    end
%%   
    % Convert image to L*a*b* colourspace.  This gives us a colourspace that is
    % nominally perceptually uniform. This allows us to use the euclidean
    % distance between colour coordinates to measure differences between
    % colours.  Note the image becomes double after conversion.  We may want to
    % go to signed shorts to save memory.
%#####################################################################
    im = rgb2lab(im); %执行完这句，图像变成双精度的了

    % Apply median filtering to colour components if mw has been supplied
    % and/or non-zero
    if mw
        if length(mw) == 1
            mw(2) = mw(1);  % Use same filtering for L and chrominance 对L和色度使用相同的滤波
        end
        for n = 1:3
            im(:,:,n) = medfilt2(im(:,:,n), [mw(1) mw(1)]);%medfilt2消除噪声, 中值滤波器
        end
    end
    
    % Nominal spacing between grid elements assuming hexagonal grid
    % 六角网格间的标称间距，说白了就是两个六边形中心间距
    S = sqrt(rows*cols / (k * sqrt(3)/2));
    
    % Get nodes per row allowing a half column margin at one end that alternates
    % from row to row  获取每列的节点, 允许在一端从行到行交替的半列边距
    nodeCols = round(cols/S - 0.5);
    % Given an integer number of nodes per row recompute S 给定整数的每列节点数重新计算 S
    S = cols/(nodeCols + 0.5);

    % Get number of rows of nodes allowing 0.5 row margin top and bottom  获取允许0.5 行边距顶部和底端的节点行数
%     nodeRows = round(rows/(sqrt(3)/2*S));%行节点数
    nodeRows = round(rows/S-0.5);%行节点数_______不明白
%     vSpacing = rows/nodeRows;
    vSpacing = rows/(nodeRows+0.5);

    % Recompute k  重新计算k
    k = nodeRows * nodeCols;
    
    % Allocate memory and initialise clusters, labels and distances.分配内存和初始化中心、标签和距离。
    C = zeros(6,k);          % Cluster centre data  1:3 is mean Lab value,初始化聚类中心
                             % 4:5 is row, col of centre, 6 is No of pixels
    l = -ones(rows, cols);   % Pixel labels.
    d = inf(rows, cols);     % Pixel distances from cluster centres每个像素到聚类中心的距离.inf 表示无穷大
    
    % Initialise clusters on a hexagonal grid
    kk = 1;
    r = vSpacing/2;
 %% 没看   
    for ri = 1:nodeRows
        % Following code alternates the starting column for each row of grid
        % points to obtain a hexagonal pattern. Note S and vSpacing are kept
        % as doubles to prevent errors accumulating across the grid.
        % 下面的代码替换每行栅格点的起始列以获得六角模式。
        % 注意 S 和 vSpacing 保持为双精度, 以防止在网格中累积错误。
        if mod(ri,2),c = S/2; else, c = S;  end%是ri除以2的余数
        
        for ci = 1:nodeCols
            cc = round(c); 
            rr = round(r);
            C(1:5, kk) = [squeeze(im(rr,cc,:)); cc; rr];%
            c = c+S;
            kk = kk+1;
        end
        
        r = r+vSpacing;
    end
 %%   
    % Now perform the clustering.现在执行聚类  10 iterations is suggested but I suspect n
    % could be as small as 2 or even 1
    S = round(S);  % We need S to be an integer from now on
    
    for n = 1:nItr
       for kk = 1:k  % for each cluster每一个聚类中心

           % Get subimage around cluster在聚类周围获得子图像
           rmin = max(C(5,kk)-S, 1);   
           rmax = min(C(5,kk)+S, rows);%C的第5行是聚类中心的行
           cmin = max(C(4,kk)-S, 1);   
           cmax = min(C(4,kk)+S, cols); 
           subim = im(rmin:rmax, cmin:cmax, :);  
           assert(numel(subim) > 0)%判断程序能否/是否需要继续执行
           
           % Compute distances D between C(:,kk) and subimage
           if USEDIST
               D = dist(C(:, kk), subim, rmin, cmin, S, m);%???怎么算的距离
           else
               D = dist2(C(:, kk), subim, rmin, cmin, S, m, eim, We);
           end

           % If any pixel distance from the cluster centre is less than its
           % previous value update its distance and label
           subd =  d(rmin:rmax, cmin:cmax);
           subl =  l(rmin:rmax, cmin:cmax);
           updateMask = D < subd;%？
           subd(updateMask) = D(updateMask);%？
           subl(updateMask) = kk;
           
           d(rmin:rmax, cmin:cmax) = subd;
           l(rmin:rmax, cmin:cmax) = subl;           
       end
       
       % Update cluster centres with mean values
       C(:) = 0;
       for r = 1:rows
           for c = 1:cols
              tmp = [im(r,c,1); im(r,c,2); im(r,c,3); c; r; 1];
              C(:, l(r,c)) = C(:, l(r,c)) + tmp;
           end
       end
       
       % Divide by number of pixels in each superpixel to get mean values
       for kk = 1:k 
           C(1:5,kk) = round(C(1:5,kk)/C(6,kk)); 
       end
       
       % Note the residual error, E, is not calculated because we are using a
       % fixed number of iterations 
    end
    
    % Cleanup small orphaned regions and 'spurs' on each region using
    % morphological opening on each labeled region.  The cleaned up regions are
    % assigned to the nearest cluster. The regions are renumbered and the
    % adjacency matrix regenerated.  This is needed because the cleanup is
    % likely to change the number of labeled regions.
    if seRadius
%#####################################################################
        [l, Am] = mcleanupregions(l, seRadius);
    else
%#####################################################################
        l = makeregionsdistinct(l);
%#####################################################################
        [l, minLabel, maxLabel] = renumberregions(l);
%#####################################################################
        Am = regionadjacency(l);    
    end

    % Recompute the final superpixel attributes and write information into
    % the Sp struct array.
    N = length(Am);
    Sp = struct('L', cell(1,N), 'a', cell(1,N), 'b', cell(1,N), ...
                'stdL', cell(1,N), 'stda', cell(1,N), 'stdb', cell(1,N), ...
                'r', cell(1,N), 'c', cell(1,N), 'N', cell(1,N));
    [X,Y] = meshgrid(1:cols, 1:rows);
    L = im(:,:,1);    
    A = im(:,:,2);    
    B = im(:,:,3);    
    for n = 1:N
        mask = l==n;
        nm = sum(mask(:));
        if centre == MEANCENTRE     
            Sp(n).L = sum(L(mask))/nm;
            Sp(n).a = sum(A(mask))/nm;
            Sp(n).b = sum(B(mask))/nm;
            
        elseif centre == MEDIANCENTRE
            Sp(n).L = median(L(mask));
            Sp(n).a = median(A(mask));
            Sp(n).b = median(B(mask));
        end
        
        Sp(n).r = sum(Y(mask))/nm;
        Sp(n).c = sum(X(mask))/nm;
        
        % Compute standard deviations of the colour components of each super
        % pixel. This can be used by code seeking to merge superpixels into
        % image segments.  Note these are calculated relative to the mean colour
        % component irrespective of the centre being calculated from the mean or
        % median colour component values.
        Sp(n).stdL = std(L(mask));
        Sp(n).stda = std(A(mask));
        Sp(n).stdb = std(B(mask));

        Sp(n).N = nm;  % Record number of pixels in superpixel too.
    end
    
%-- dist -------------------------------------------
%
% Usage:  D = dist(C, im, r1, c1, S, m)
% 
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.  try m in the range [1-40] for L*a*b* space
%
% ?? Might be worth trying the Geometric Mean instead ??
%  Distance = sqrt(dc * ds)
% but having a factor 'm' to play with is probably handy

% This code could be more efficient

function D = dist(C, im, r1, c1, S, m)

    % Squared spatial distance
    %    ds is a fixed 'image' we should be able to exploit this
    %    and use a fixed meshgrid for much of the time somehow...
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);  % x and y dist from cluster centre
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour difference
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    D = sqrt(dc2 + ds2/S^2*m^2);
    
    
    
%--- dist2 ------------------------------------------
%
% Usage:  D = dist2(C, im, r1, c1, S, m, eim)
% 
% Arguments:   C - Cluster being considered
%             im - sub-image surrounding cluster centre
%         r1, c1 - row and column of top left corner of sub image within the
%                  overall image.
%              S - grid spacing
%              m - weighting factor between colour and spatial differences.
%            eim - Edge strength sub-image corresponding to im
%
% Returns:     D - Distance image giving distance of every pixel in the
%                  subimage from the cluster centre
%
% Distance = sqrt( dc^2 + (ds/S)^2*m^2 )
% where:
% dc = sqrt(dl^2 + da^2 + db^2)  % Colour distance
% ds = sqrt(dx^2 + dy^2)         % Spatial distance
%
% m is a weighting factor representing the nominal maximum colour distance
% expected so that one can rank colour similarity relative to distance
% similarity.  try m in the range [1-40] for L*a*b* space
%

function D = dist2(C, im, r1, c1, S, m, eim, We)

    % Squared spatial distance
    %    ds is a fixed 'image' we should be able to exploit this
    %    and use a fixed meshgrid for much of the time somehow...
    [rows, cols, chan] = size(im);
    [x,y] = meshgrid(c1:(c1+cols-1), r1:(r1+rows-1));
    x = x-C(4);
    y = y-C(5);
    ds2 = x.^2 + y.^2;
    
    % Squared colour difference
    for n = 1:3
        im(:,:,n) = (im(:,:,n)-C(n)).^2;
    end
    dc2 = sum(im,3);
    
    % Combine colour and spatial distance measure
    D = sqrt(dc2 + ds2/S^2*m^2);
    
    % for every pixel in the subimage call improfile to the cluster centre
    % and use the largest value as the 'edge distance'
    rCentre = C(5)-r1;   % Cluster centre coords relative to this sub-image
    cCentre = C(4)-c1;
    de = zeros(rows,cols);
    for r = 1:rows
        for c = 1:cols
            v = improfile(eim,[c cCentre], [r rCentre]);
            de(r,c) = max(v);
        end
    end

    % Combine edge distance with weight, We with total Distance.
    D = D + We * de;
    