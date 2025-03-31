function varargout = main(varargin)
% MAIN MATLAB code for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 14-Oct-2019 13:25:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_OpeningFcn, ...
                   'gui_OutputFcn',  @main_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=ones(256,256);
axes(handles.axes1);imshow(a);
axes(handles.axes2);imshow(a);
axes(handles.axes3);imshow(a);
axes(handles.axes4);imshow(a);
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in select.
function select_Callback(hObject, eventdata, handles)
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%load image
cd Inputs
[file path] = uigetfile('*.jpg;*.png;*.bmp','Pick an image file');
image = im2double(imread(file));
image=imresize(image,[240 360]);
cd ..
handles.image = image;
axes(handles.axes1)
imshow(image);
%update handles structure
guidata(hObject, handles);


% --- Executes on button press in seg.
function seg_Callback(hObject, eventdata, handles)
% hObject    handle to seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get hairs using bottomhat filter
image=handles.image;
se = strel('disk',8);
hairs = imbothat(image,se);
%stats
%matrice = logical(rgb2gray(hairs));
replacedImage = roifill(rgb2gray(image), rgb2gray(hairs));
axes(handles.axes2)
imshow(replacedImage);
handles.replacedImage = replacedImage;

%update handles structure
guidata(hObject, handles);

% --- Executes on button press in depth.
function depth_Callback(hObject, eventdata, handles)
% hObject    handle to depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
left=handles.replacedImage

window=3;
dispmin=0;
dispmax=69;
%scale_factor=5;
scale_factor=floor(255/dispmax); 
span=(window-1)/2;
% [filename1,pathname1]= uigetfile('*.png;*.jpg;*.bmp','Select the left image');
% [filename2,pathname2]= uigetfile('*.png;*.jpg;*.bmp','Select the right image');
% left=imread([pathname1,filename1]);
% leftclr=imread([pathname1,filename1]);
%  left=double(left);


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
tic;
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

toc;

finalmap=uint8(map);
finalmap=medfilt2(finalmap,[9 9]);

% % figure(1),imshow(mat2gray(finalmap));
finalmap(:,:)=finalmap(:,:)*scale_factor;
finalmap=finalmap(span+1:m+span,span+1:n+span);
%dlmwrite('C:\phd\pgms\matlabresults\SAD\resized\window13\sad50.txt',finalmap);
% % imwrite(finalmap,'sad50.png');     
axes(handles.axes3)
imshow(finalmap)

handles.finalmap = finalmap;

%update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
finalmap=handles.finalmap;

axes(handles.axes4)

mesh(double(finalmap)); 
% % colormap


% --- Executes on button press in feat.
function feat_Callback(hObject, eventdata, handles)
% hObject    handle to feat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
finalmap=handles.finalmap;
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


Feat = [Fa1];
Feat = Feat';
QFeat = Feat;
save QFeat QFeat

msgbox('Feature Extraction completed')

% --- Executes on button press in train.
function train_Callback(hObject, eventdata, handles)
% hObject    handle to train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
train;

% --- Executes on button press in classify.
function classify_Callback(hObject, eventdata, handles)
% hObject    handle to classify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load QFeat;
load dfeatures;
load svmst;

% qfeatx=max(max(qfeat))
%%%%%%classification
cout = svmclassify(svmst,QFeat);
% cout = vec2ind(cout);

if isequal(cout,1)
msgbox('Malignant')
else isequal(cout,2)
msgbox('Benign');
end