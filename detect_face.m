function []=detect_face(I)

close all;  %close all window if their open (this is just so I don't go crazy having to do it manually)
Faces = [];
numFaceFound = 0;

I = double(I);

resultI = I;

H = size(I,1); %height
W = size(I,2); %width
R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);

YCbCr = rgb2ycbcr(I);               %WHAT DOES THIS DO?
Y = YCbCr(:,:,1);

minY = min(min(Y));                 %SOUNDS LIKE MAGIC?
maxY = max(max(Y));

Y = 255.0*(Y-minY)./(maxY-minY);
YEye = Y;
Yavg = sum(sum(Y))/(W*H);

T=1; %default T value is 1

if (Yavg<64) %reassign T value
    T=1.4;
elseif (Yavg>192)
    T=0.6;
end

if (T~=1)   %check if T value should change
    RI=R.^T;
    GI=G.^T;
else
    RI=R;
    GI=G;
end


%This part if for Lighing compensation - I DON'T KNOW WHAT THAT IS FOR
C = zeros(H, W, 3);  %new array of the same size as the input image (zeroed)
C(:,:,1) = RI;  %assigning the appropriate T-adjusted colors
C(:,:,2) = GI;
C(:,:,3) = B;   %B stays the same                           %WHY?

YCbCr=rgb2ycbcr(C);
Cr=YCbCr(:,:,3);

%This part is for finding the skin
S = zeros(H,W);
[SkinIndexRow,SkinIndexCol] =find(10<Cr & Cr<45);

for i=1:length(SkinIndexRow)
    S(SkinIndexRow(i),SkinIndexCol(i))=1;
end


%Skin with nose removal
SN=zeros(H,W);
for i=1:H-5
    for j=1:W-5
        localSum=sum(sum(S(i:i+4, j:j+4)));
        SN(i:i+5, j:j+5)=(localSum>12);
    end
end

% figure,imshow(SN);
% title('skin with noise removal');
L = bwlabel(SN,8);
BB  = regionprops(L, 'BoundingBox');
bboxes= cat(1, BB.BoundingBox);
widths=bboxes(:,3);
heights=bboxes(:,4);
hByW=heights./widths;

lenRegions=size(bboxes,1);
foundFaces=zeros(1,lenRegions);

rgb=label2rgb(L);
% figure,imshow(rgb);
% title('face candidates');

for i = 1:lenRegions

    if (hByW(i)>1.75 || hByW(i)<0.75)
        continue;
    end

    if (heights(i)<20 && widths(i)<20)
        continue;
    end

    CurBB=bboxes(i,:);
    XStart=CurBB(1);
    YStart=CurBB(2);
    WCur=CurBB(3);
    HCur=CurBB(4);
    rangeY=int32(YStart):int32(YStart+HCur-1);
    rangeX= int32(XStart):int32(XStart+WCur-1);


    RIC=RI(rangeY, rangeX); %
    GIC=GI(rangeY, rangeX); %creating a new image of the cropped face from the original image.
    BC=B(rangeY, rangeX);   %

    %drawing rectangles on a dup of the original
    figure, imshow(RI/255);
    title('yo-face');

    %# draw a rectangle
    rectangle('Position',[int32(XStart) int32(YStart) int32(WCur-1) int32(HCur-1)], 'LineWidth',2, 'EdgeColor','b');

    M=zeros(HCur, WCur);
    theta=acos( 0.5.*(2.*RIC-GIC-BC) ./ sqrt( (RIC-GIC).*(RIC-GIC) + (RIC-BC).*(GIC-BC) ) );
    theta(isnan(theta))=0;
    thetaMean=mean2(theta);
    [MouthIndexRow,MouthIndexCol] = find(theta < thetaMean / 4);

    for j=1:length(MouthIndexRow)
        M(MouthIndexRow(j),MouthIndexCol(j))=1;
    end

    Hist=zeros(1, HCur);

    for j=1:HCur
        Hist(j)=length(find(M(j,:)==1));
    end

    wMax=find(Hist==max(Hist));
    wMax=wMax(1);

    if (wMax < WCur/6)
        continue;
    end

    % figure, imshow(M);
    % title('Mouth map');
    eyeH = HCur-wMax;
    eyeW = WCur;

    YC = YEye(YStart:YStart+eyeH-1, XStart:XStart+eyeW-1);

    E = zeros(eyeH, eyeW);
    [EyeIndexRow,EyeIndexCol] = find(65<YC & YC<80);
    for j=1:length(EyeIndexRow)
        E(EyeIndexRow(j),EyeIndexCol(j))=1;
    end

    EyeExist=find(Hist > 0.3*wMax);

    if (~(length(EyeExist)>0))
        continue;
    end

    foundFaces(i)=1;
    numFaceFound=numFaceFound+1;

end

    disp('Number of faces found');
    numFaceFound;

    if (numFaceFound > 0)
        disp('Indices of faces found: ');
        ind = find(foundFaces == 1);
        CurBB=bboxes(ind,:);
        CurBB
    else
    close all;
end

end