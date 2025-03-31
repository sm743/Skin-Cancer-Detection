function svmst = train

% lda = waitbar(0,'Db Loading....');
for di=1:1:144
    
    fname = strcat(int2str(di),'.bmp');
    cd database   
       tinp = imread(fname);
    cd ..
image=imresize(tinp,[240 360]);
se = strel('disk',8);
hairs = imbothat(image,se);
%stats
%matrice = logical(rgb2gray(hairs));
replacedImage = roifill(rgb2gray(image), rgb2gray(hairs));
 window=3;
dispmin=0;
dispmax=69;
%scale_factor=5;
scale_factor=floor(255/dispmax); 
span=(window-1)/2;
left=replacedImage;
mriVolumeTranslated = imtranslate(left,[40,30],'OutputView','full');
sizeIn=size(left);
sliceIndex = round(sizeIn(2)/1);
axialSliceOriginal   = left(:,sliceIndex);
axialSliceTranslated = mriVolumeTranslated(:,sliceIndex);



[m,n]=size(left);


leftImage=zeros((m+(2*span)),(n+dispmax+(2*span)));
for i=span+1:1:m+span
    for j=span+1:1:n+span
    leftImage(i,j)=left(i-span,j-span);
    end
end
% right=imread([pathname2,filename2]);
right=double(axialSliceTranslated);
 right=imresize(right,[240 360]);


[k,l]=size(right);
rightImage=zeros((k+(2*span)),(l+dispmax+(2*span)));
for i=span+1:1:k+span
    for j=span+1:1:l+span
    rightImage(i,j)=right(i-span,j-span);
    end
end
[row col]=size(leftImage);
map=zeros(row,col);
finaldisp=99999;

for i=1+span:1:row-span
    for j=1+span:1:col-span-dispmax
               prevbest=99999;
            for disprange=dispmin:1:dispmax
               sad=0;
                for winrow=-span:1:span
                    for wincol=-span:1:span
                        temp = rightImage(i+winrow,j+wincol)-leftImage(i+winrow,j+wincol+disprange);
                        sad=sad+abs(temp);
                      
                    end
                    
                end
                if(sad<prevbest)  
                    finaldisp=disprange;  
                    prevbest=sad;
                end
           
            end
        map(i,j)=finaldisp;
    end
     
end



finalmap=uint8(map);
finalmap=medfilt2(finalmap,[9 9]);

% % figure(1),imshow(mat2gray(finalmap));
finalmap(:,:)=finalmap(:,:)*scale_factor;
finalmap=finalmap(span+1:m+span,span+1:n+span);   
    
     
Dinp = uint8(finalmap);
Dmax = max(Dinp(:));
Dmin = min(Dinp(:));
NumL = Dmax - Dmin;
Gm = graycomatrix(Dinp,'Graylimit',[Dmin Dmax],'Numlevels',NumL);
GProps = graycoprops(Gm);

Fa1 = GProps.Energy;
Fa2 = GProps.Contrast;
Fa3 = GProps.Correlation;
Fa4 = GProps.Homogeneity;
Fa5 = entropy(Gm);

% save Feat2 F2

Feat = [Fa1 ];

 dfeatures(:,di) = [Feat]';
 
 end
% close(lda);

%%%%%Neural network creation and training 

%%%%Assigning target to each class features
Nc = 72; T=1;
save dfeatures dfeatures;
for dfi=1:size(dfeatures,2)
   
    if Nc<1
      T = T+1;
      Nc =71;
      acti(:,dfi) = T; 
    else
      acti(:,dfi) = T;  
      Nc = Nc-1;  
    end
end
       
% actv = ind2vec(acti);   %%%%%Indices to vector creation
%  disp(acti);
svmst =fitcsvm(dfeatures,acti);   %%%%network training

save svmst svmst;

helpdlg('SVM training completed');

return;  

