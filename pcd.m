function varargout = pcd(varargin)
% PCD MATLAB code for pcd.fig
%      PCD, by itself, creates a new PCD or raises the existing
%      singleton*.
%
%      H = PCD returns the handle to a new PCD or the handle to
%      the existing singleton*.
%
%      PCD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCD.M with the given input arguments.
%
%      PCD('Property','Value',...) creates a new PCD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pcd_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pcd_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pcd

% Last Modified by GUIDE v2.5 13-May-2015 14:07:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pcd_OpeningFcn, ...
                   'gui_OutputFcn',  @pcd_OutputFcn, ...
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


% --- Executes just before pcd is made visible.
function pcd_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pcd (see VARARGIN)

% Choose default command line output for pcd
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pcd wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pcd_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton8.

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Open_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[namafile, formatfile] = uigetfile({'*.png'}, 'membuka gambar'); %memilih gambar
image = imread([formatfile, namafile]); %membaca gambar
image=uint8(image);
handles.filename = namafile;
guidata(hObject, handles);
axes(handles.axes1); %memilih axes1 sebagai letak gambar yang dimunculkan
imshow(image); %memunculkan gambar


% --------------------------------------------------------------------
function flip_Callback(hObject, eventdata, handles)
% hObject    handle to flip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function horizontal_Callback(hObject, eventdata, handles)
% hObject    handle to horizontal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
s=size(image);
result=zeros(s(2),s(1),s(3)); % 
for i = 1:s(1)
    j=1:s(2);
        k=s(2)-j+1;
        result(i,k,:)=image(i,j,:);
end
image=uint8(result);
imshow(image); 

% --------------------------------------------------------------------
function vertical_Callback(hObject, eventdata, handles)
% hObject    handle to vertical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
s=size(image);
result=zeros(s(2),s(1),s(3)); % 
for j = 1:s(2)
    i=1:s(1);
        k=s(1)-i+1;
        result(k,j,:)=image(i,j,:);
end
image=uint8(result);
imshow(image); 

% --------------------------------------------------------------------
function r90_Callback(hObject, eventdata, handles)
% hObject    handle to r90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
s=size(image);
result=zeros(s(2),s(1),s(3)); % 
for i = 1:s(1)
    j=1:s(2);
        k=s(2)-j+1;
        result(i,k,:)=image(j,i,:);
end
image=uint8(result);
imshow(image); 

% --------------------------------------------------------------------
function r90ccw_Callback(hObject, eventdata, handles)
% hObject    handle to r90ccw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
s=size(image);
result=zeros(s(2),s(1),s(3)); % 
for i = 1:s(1)
    j=1:s(2);
        k=s(2)-j+1;
        result(k,i,:)=image(i,j,:);
end
image=uint8(result);
imshow(image);

% --------------------------------------------------------------------
function R180_Callback(hObject, eventdata, handles)
% hObject    handle to R180 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[h w d] = size(image);
image2 = ones(size(image));
image2 = imrotate(image,180);
image = image2;
imshow(image);



% --- Executes on button press in crop.
function crop_Callback(hObject, eventdata, handles)
% hObject    handle to crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
image = imcrop(image);
imshow(image);


% --- Executes on button press in swirl.
function swirl_Callback(hObject, eventdata, handles)
% hObject    handle to swirl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
s=size(image);
K=30; 
x2=zeros([size(image,1) size(image,2)]);
y2=zeros([size(image,1) size(image,2)]); 
midx=ceil((size(image,1)+1)/2); 
midy=ceil((size(image,2)+1)/2);
for i=1:size(image,1)
    x=i-midx-K; 
    for j=1:size(image,2) %Cartesian to Polar co-ordinates 
        [theta1,rho1]=cart2pol(x,j-midy+K);
        phi=theta1+(rho1/K);
%Polar to Cartesian co-ordinates 
[l,m]=pol2cart(phi,rho1); x2(i,j)=ceil(l)+midx; y2(i,j)=ceil(m)+midy;
    end
end
%The result may produce value lesser than 1 or greater than the image size.
x2=max(x2,1); x2=min(x2,size(image,1));
y2=max(y2,1); y2=min(y2,size(image,2));
        for i=1:size(image,1)
            for j=1:size(image,2)
                result(i,j,:)=image(x2(i,j),y2(i,j),:);
            end
        end
image=uint8(result);
imshow(image);

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
s=size(image);
result=zeros(s(2),s(1),s(3)); % 
for i = 1:s(1)
    j=1:s(2);

        k=s(2)-j+1;
        result(k,i,:)=image(i,j,:);
end
image=uint8(result);
imshow(image);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function zoomout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in cutpaste.
function cutpaste_Callback(hObject, eventdata, handles)
% hObject    handle to cutpaste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[h w d]=size(image);
img = image;
for i=10:60
    for j=10:60
        img(i,j,:)=255;
        img(h-70+i,w-70+j,:)=image(i,j,:); 
    end
end
img=uint8(img);
figure,imshow(img);


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
teks = get(handles.edit1,'String');
x=str2num(teks);
image=imrotate(image,x);
imshow(image);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in histeq.
function histeq_Callback(hObject, eventdata, handles)
% hObject    handle to histeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
image(:,:,1)=histeq(image(:,:,1));
image(:,:,2)=histeq(image(:,:,2));
image(:,:,3)=histeq(image(:,:,3));
imshow(image);

% --- Executes on button press in histo.
function histo_Callback(hObject, eventdata, handles)
% hObject    handle to histo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[h w d]=size(image);
img=zeros(size(image));
for i=1:h
    for j=1:w
        img(i,j,:)=image(min(i+30,h),min(j+20,w),:);
    end
end
img=uint8(img);
figure,imshow(img);


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
% load image


g= rgb2gray(image);
%tampilkan image asal
figure %1
imshow(g)

%lakukan fft & shift hasilnya
fft_g = fft2(g);
fs = fftshift(fft_g);

%tampilkan image hasil fft
figure %2
imshow(log(abs(fft_g)),[])

%kembalikan hasilnya ke image asal
balik = ifft2(fft_g);
figure %3
imshow(uint8(balik))

% -----------------------------------
%coba lihat apa yang terjadi jika hasil fft dimodifikasi

fs1 = fs;
fs1(150:250,1:100) = 0;
%tampilkan image hasil fft
figure %4
imshow(log(abs(fs1)),[])

%kembalikan hasilnya ke image asal
balik = ifft2(ifftshift(fs1));
figure %5
imshow(uint8(balik))

mask = zeros(size(g));
mask(150:250,1:100) = 1;

hasil = mask.*balik;
%tampilkan image hasil fft
figure %6
imshow(log(abs(hasil)),[])

%kembalikan hasilnya ke image asal
balik = ifft2(ifftshift(hasil));
figure %7
imshow(uint8(balik))

% --------------------------------------------------------------------
function histo_rgb_Callback(hObject, eventdata, handles)
% hObject    handle to histo_rgb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
figure, title('Image Histogram');
subplot(2,2,1);
imhist(image(:,:,1));
title('Red');
subplot(2,2,2);
imhist(image(:,:,2));
title('Green');
subplot(2,2,3);
imhist(image(:,:,3));
title('Blue');
hist(reshape(image,[],3),1:max(image(:))); 
colormap([1 0 0; 0 1 0; 0 0 1]);
figure,hist;

% --------------------------------------------------------------------
function histo_gray_Callback(hObject, eventdata, handles)
% hObject    handle to histo_gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
ig = rgb2gray(image);
figure, imshow(ig);
figure,imhist(ig);

% --------------------------------------------------------------------
function zoom_in_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[h w d] = size(image);
image2 = ones([h*2 w*2 d]);
a=1;
for i=1:h 
    b=1;
    for j=1:w
            image2(a,b,:) = image(i,j,:);
            image2(a,b+1,:) = image(i,j,:);
            image2(a+1,b,:) = image(i,j,:);
            image2(a+1,b+1,:) = image(i,j,:);
        b=b+2;
    end
    a=a+2;
end
image2 = uint8(image2);
image3 = ones([h*2 w*2 d]);
for i=1:h
    for j=1:w
        image3(i,j,:) = image(i,j,:);
    end
end
image3 = uint8(image3);
figure, title('Zoom In');
subplot(1,2,1);
h = imshow(image3);
title('Before');
subplot(1,2,2);
h = imshow(image2);
title('After');

% --------------------------------------------------------------------
function zoom_out_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
image2 = imresize(image,0.5);
image3 = zeros(size(image));
sizee = size(image2);
for i=1:sizee(1)
     for j=1:sizee(2)
        image3(i,j,:) = image2(i,j,:);
     end
end
image3 = uint8(image3);
figure, title('Zoom Out');
subplot(1,2,1);
h = imshow(image);
title('Before');
subplot(1,2,2);
h = imshow(image3);
title('After');

% --------------------------------------------------------------------
function mean_filter_Callback(hObject, eventdata, handles)
% hObject    handle to mean_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
img = imresize(image,[512 512]);
img = imnoise(img,'salt & pepper',0.02)
img1 = zeros([518 518 3]);
for i=4:515
    for j=4:515
        img1(i,j,:)=img(i-3,j-3,:);
    end
end
img2 = zeros([512 512 3]);
for m=1:3
    for i=4:515
        for j=4:515
            jlh=0;
            for k=i-1:i+1
                for l=j-1:j+1
                    jlh=jlh+img1(k,l,m);
                end
            end
            img2(i-3,j-3,m)=floor(jlh/9);
        end
    end
end
img2 = uint8(img2);
figure, title('Mean Filtering');
subplot(1,3,1);
h = imshow(img2);
title('mean 3x3');
img2 = zeros([512 512 3]);
for m=1:3
    for i=4:515
        for j=4:515
            jlh=0;
            for k=i-2:i+2
                for l=j-2:j+2
                    jlh=jlh+img1(k,l,m);
                end
            end
            img2(i-3,j-3,m)=floor(jlh/25);
        end
    end
end
img2 = uint8(img2);
subplot(1,3,2);
h = imshow(img2);
title('mean 5x5');
img2 = zeros([512 512 3]);
for m=1:3
    for i=4:515
        for j=4:515
            jlh=0;
            for k=i-3:i+3
                for l=j-3:j+3
                    jlh=jlh+img1(k,l,m);
                end
            end
            img2(i-3,j-3,m)=floor(jlh/49);
        end
    end
end
img2 = uint8(img2);
subplot(1,3,3);
h = imshow(img2);
title('mean 7x7');

% --------------------------------------------------------------------
function modus_filter_Callback(hObject, eventdata, handles)
% hObject    handle to modus_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
img = imresize(image,[512 512]);
img = imnoise(img,'salt & pepper',0.02)
img1 = zeros([518 518 3]);
for i=4:515
    for j=4:515
        img1(i,j,:)=img(i-3,j-3,:);
    end
end
img2 = zeros([512 512 3]);
for m=1:3
    for i=4:515
        for j=4:515
            x = zeros([9 1]);
            counter = 0;
            for k=i-1:i+1
                for l=j-1:j+1
                    counter = counter + 1;
                    x(counter,1) = img1(k,l,m);
                end
            end
            img2(i-3,j-3,m)=mode(x);
        end
    end
end
img2 = uint8(img2);
figure, title('Modus Filtering');
subplot(1,3,1);
h = imshow(img2);
title('modus 3x3');
img2 = zeros([512 512 3]);
for m=1:3
    for i=4:515
        for j=4:515
            x = zeros([25 1]);
            counter = 0;
            for k=i-2:i+2
                for l=j-2:j+2
                    counter = counter + 1;
                    x(counter,1) = img1(k,l,m);
                end
            end
            img2(i-3,j-3,m)=mode(x);
        end
    end
end
img2 = uint8(img2);
subplot(1,3,2);
h = imshow(img2);
title('modus 5x5');
img2 = zeros([512 512 3]);
for m=1:3
    for i=4:514
        for j=4:515
            x = zeros([49 1]);
            counter = 0;
            for k=i-3:i+3
                for l=j-3:j+3
                    counter = counter + 1;
                    x(counter,1) = img1(k,l,m);
                end
            end
            img2(i-3,j-3,m)=mode(x);
        end
    end
end
img2 = uint8(img2);
subplot(1,3,3);
h = imshow(img2);
title('modus 7x7');

% --------------------------------------------------------------------
function median_filter_Callback(hObject, eventdata, handles)
% hObject    handle to median_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
figure, title('Median Filtering');
img = zeros(size(image));
for i=1:3    
    img(:,:,i)=medfilt2(image(:,:,i),[3 3]);
end
img = uint8(img);
subplot(1,3,1);
h = imshow(img);
title('median 3x3');
for i=1:3    
    img(:,:,i)=medfilt2(image(:,:,i),[5 5]);
end
img = uint8(img);
subplot(1,3,2);
h = imshow(img);
title('median 5x5');
for i=1:3    
    img(:,:,i)=medfilt2(image(:,:,i),[7 7]);
end
img = uint8(img);
subplot(1,3,3);
h = imshow(img);
title('median 7x7');

% --------------------------------------------------------------------
function sobel_Callback(hObject, eventdata, handles)
% hObject    handle to sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
img = rgb2gray(image);
img2 = edge(img,'sobel');
figure,imshow(img2);
title('sobel');

% --------------------------------------------------------------------
function prewitt_Callback(hObject, eventdata, handles)
% hObject    handle to prewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
img = rgb2gray(image);
img2 = edge(img,'prewitt');
figure,imshow(img2);
title('prewitt');


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function image_Callback(hObject, eventdata, handles)
% hObject    handle to image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function rotate_Callback(hObject, eventdata, handles)
% hObject    handle to rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function histogram_Callback(hObject, eventdata, handles)
% hObject    handle to histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function highpass_Callback(hObject, eventdata, handles)
% hObject    handle to highpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function lowpass_Callback(hObject, eventdata, handles)
% hObject    handle to lowpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in EROTION.
function EROTION_Callback(hObject, eventdata, handles)
% hObject    handle to EROTION (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
se = strel('square',20);
image = imerode(image,se);
imshow(image);


% --- Executes on button press in DILATION.
function DILATION_Callback(hObject, eventdata, handles)
% hObject    handle to DILATION (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
se = strel('square',20);
image = imdilate(image,se);
imshow(image);


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[h w d]=size(image);
img = zeros(size(image));
thresh = ceil(255/2);
for i=1:h
    for j=1:w
        if((image(i,j,1)>thresh)&&(image(i,j,2)<thresh)&&(image(i,j,3)<thresh))
            img(i,j,:)=image(i,j,:);
        end
    end
end
img = uint8(img);
img2 = im2bw(img,0.01);
a = regionprops(img2,'Image','BoundingBox');
for i=1:size(a)
    box = a(i).BoundingBox;
    img3 = imcrop(img,box);
    figure,imshow(img3);
    s = strcat('Merah ',int2str(i));
    title(s);
end
img = zeros(size(image));
thresh = ceil(255/2);
for i=1:h
    for j=1:w
        if((image(i,j,2)>thresh)&&(image(i,j,1)<thresh)&&(image(i,j,1)<thresh))
            img(i,j,:)=image(i,j,:);
        end
    end
end
img = uint8(img);
img2 = im2bw(img,0.01);
a = regionprops(img2,'Image','BoundingBox');
for i=1:size(a)
    box = a(i).BoundingBox;
    img3 = imcrop(img,box);
    figure,imshow(img3);
    s = strcat('Hijau ',int2str(i));
    title(s);
end
img = zeros(size(image));
thresh = ceil(255/2);
for i=1:h
    for j=1:w
        if((image(i,j,3)>thresh)&&(image(i,j,1)<thresh)&&(image(i,j,2)<thresh))
            img(i,j,:)=image(i,j,:);
        end
    end
end
img = uint8(img);
img2 = im2bw(img,0.01);
a = regionprops(img2,'Image','BoundingBox');
for i=1:size(a)
    box = a(i).BoundingBox;
    img3 = imcrop(img,box);
    figure,imshow(img3);
    s = strcat('Biru ',int2str(i));
    title(s);
end
img = zeros(size(image));
thresh = ceil(255/2);
for i=1:h
    for j=1:w
        if((image(i,j,1)>thresh)&&(image(i,j,3)>thresh)&&(image(i,j,2)<thresh))
            img(i,j,:)=image(i,j,:);
        end
    end
end
img = uint8(img);
img2 = im2bw(img,0.01);
a = regionprops(img2,'Image','BoundingBox');
for i=1:size(a)
    box = a(i).BoundingBox;
    img3 = imcrop(img,box);
    figure,imshow(img3);
    s = strcat('Ungu ',int2str(i));
    title(s);
end



% --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function thresholding_Callback(hObject, eventdata, handles)
% hObject    handle to thresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[h w d]=size(image);
img = zeros(size(image));
thresh = ceil(255/2);
for i=1:h
    for j=1:w
        if((image(i,j,1)>thresh)&&(image(i,j,2)<thresh)&&(image(i,j,3)<thresh))
            img(i,j,:)=image(i,j,:);
        end
    end
end
img = uint8(img);
img2 = im2bw(img,0.01);
a = regionprops(img2,'Image','BoundingBox');
for i=1:size(a)
    box = a(i).BoundingBox;
    img3 = imcrop(img,box);
    figure,imshow(img3);
    s = strcat('Merah ',int2str(i));
    title(s);
end
img = zeros(size(image));
thresh = ceil(255/2);
for i=1:h
    for j=1:w
        if((image(i,j,2)>thresh)&&(image(i,j,1)<thresh)&&(image(i,j,1)<thresh))
            img(i,j,:)=image(i,j,:);
        end
    end
end
img = uint8(img);
img2 = im2bw(img,0.01);
a = regionprops(img2,'Image','BoundingBox');
for i=1:size(a)
    box = a(i).BoundingBox;
    img3 = imcrop(img,box);
    figure,imshow(img3);
    s = strcat('Hijau ',int2str(i));
    title(s);
end
img = zeros(size(image));
thresh = ceil(255/2);
for i=1:h
    for j=1:w
        if((image(i,j,3)>thresh)&&(image(i,j,1)<thresh)&&(image(i,j,2)<thresh))
            img(i,j,:)=image(i,j,:);
        end
    end
end
img = uint8(img);
img2 = im2bw(img,0.01);
a = regionprops(img2,'Image','BoundingBox');
for i=1:size(a)
    box = a(i).BoundingBox;
    img3 = imcrop(img,box);
    figure,imshow(img3);
    s = strcat('Biru ',int2str(i));
    title(s);
end
img = zeros(size(image));
thresh = ceil(255/2);
for i=1:h
    for j=1:w
        if((image(i,j,1)>thresh)&&(image(i,j,3)>thresh)&&(image(i,j,2)<thresh))
            img(i,j,:)=image(i,j,:);
        end
    end
end
img = uint8(img);
img2 = im2bw(img,0.01);
a = regionprops(img2,'Image','BoundingBox');
for i=1:size(a)
    box = a(i).BoundingBox;
    img3 = imcrop(img,box);
    figure,imshow(img3);
    s = strcat('Ungu ',int2str(i));
    title(s);
end

% --------------------------------------------------------------------
function regiongrowing_Callback(hObject, eventdata, handles)
% hObject    handle to regiongrowing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
img1 = zeros(size(image));
img2 = rgb2gray(image);
img3 = img2;
img = zeros(size(image));
[h w d]= size(image);
a = randi([1 h-1],1,3);
b = randi([1 w-1],1,3);
for k=1:3
img1 = zeros(size(image));
x = a(k);
y = b(k);
for i=1:h
    for j=1:w
        if((img3(i,j)>img2(x,y)-10)&&(img3(i,j)<img2(x,y,:)+10))
            img1(i,j,:)=image(x,y,:);
        end
    end
end
img1 = uint8(img1);
figure,imshow(img1);
end


% --------------------------------------------------------------------
function smpling_Callback(hObject, eventdata, handles)
% hObject    handle to smpling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
imgs = cell(4,1);
f = inline('max(x(:))');
for i=1:4
    imgs{i} = imresize(imresize(image,1-(i/5),'box'),i,'box');
end

figure, title('Sampling Image');
for i=1:4
     subplot(2,2,i);
     h = imshow(imgs{i});
     title(num2str(i)); 
end

% --------------------------------------------------------------------
function quantization_Callback(hObject, eventdata, handles)
% hObject    handle to quantization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
figure, title('Sampling Image');
subplot(2,2,1);
h = imshow(image);
title(num2str(1));
subplot(2,2,2);
h = imshow(im2uint16(image));
title('16-bit');
subplot(2,2,3);
X4bit=uint8(4*round(double(image)/4));
h = imshow(X4bit);
title('8-bit');
subplot(2,2,4);
h = imshow(im2bw(image));
title('1-bit');

% --------------------------------------------------------------------
function wrapping_Callback(hObject, eventdata, handles)
% hObject    handle to wrapping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cutnpaste_Callback(hObject, eventdata, handles)
% hObject    handle to cutnpaste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[h w d]=size(image);
img = image;
for i=10:60
    for j=10:60
        img(i,j,:)=255;
        img(h-70+i,w-70+j,:)=image(i,j,:); 
    end
end
img=uint8(img);
figure,imshow(img);

% --------------------------------------------------------------------
function translation_Callback(hObject, eventdata, handles)
% hObject    handle to translation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
[h w d]=size(image);
img=zeros(size(image));
for i=1:h
    for j=1:w
        img(i,j,:)=image(min(i+30,h),min(j+20,w),:);
    end
end
img=uint8(img);
figure,imshow(img);

% --------------------------------------------------------------------
function fft_Callback(hObject, eventdata, handles)
% hObject    handle to fft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
% load image
g= rgb2gray(image);
%tampilkan image asal
figure %1
imshow(g)

%lakukan fft & shift hasilnya
fft_g = fft2(g);
fs = fftshift(fft_g);

%tampilkan image hasil fft
figure %2
imshow(log(abs(fft_g)),[])

%kembalikan hasilnya ke image asal
balik = ifft2(fft_g);
figure %3
imshow(uint8(balik))

% -----------------------------------
%coba lihat apa yang terjadi jika hasil fft dimodifikasi

fs1 = fs;
fs1(150:250,1:100) = 0;
%tampilkan image hasil fft
figure %4
imshow(log(abs(fs1)),[])

%kembalikan hasilnya ke image asal
balik = ifft2(ifftshift(fs1));
figure %5
imshow(uint8(balik))

mask = zeros(size(g));
mask(150:250,1:100) = 1;

hasil = mask.*balik;
%tampilkan image hasil fft
figure %6
imshow(log(abs(hasil)),[])

%kembalikan hasilnya ke image asal
balik = ifft2(ifftshift(hasil));
figure %7
imshow(uint8(balik))

% --------------------------------------------------------------------
function Untitled_34_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function erosion_Callback(hObject, eventdata, handles)
% hObject    handle to erosion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
se = strel('square',20);
image = imerode(image,se);
imshow(image);

% --------------------------------------------------------------------
function dilate_Callback(hObject, eventdata, handles)
% hObject    handle to dilate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
se = strel('square',20);
image = imdilate(image,se);
imshow(image);

% --------------------------------------------------------------------
function thinning_Callback(hObject, eventdata, handles)
% hObject    handle to thinning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
img = rgb2gray(image);
thresh = graythresh(img);
img = im2bw(img,thresh);
figure, title('Thinning');
subplot(1,2,1);
imshow(img);
title('Before');
img = bwmorph(img,'thin');
subplot(1,2,2);
imshow(img);
title('After');


% --------------------------------------------------------------------
function opening_Callback(hObject, eventdata, handles)
% hObject    handle to opening (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
se = strel('disk',50);
image = imopen(image,se);
imshow(image);

% --------------------------------------------------------------------
function closing_Callback(hObject, eventdata, handles)
% hObject    handle to closing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
se = strel('disk',30);
image = imclose(image,se);
imshow(image);


% --------------------------------------------------------------------
function roberts_Callback(hObject, eventdata, handles)
% hObject    handle to roberts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
img = rgb2gray(image);
img2 = edge(img,'roberts');
figure,imshow(img2);
title('Roberts');

% --------------------------------------------------------------------
function laplacian_Callback(hObject, eventdata, handles)
% hObject    handle to laplacian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
img = rgb2gray(image);
img2 = edge(img,'log');
figure,imshow(img2);
title('Laplacian');



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image;
n = str2num(get(handles.edit2,'String'));
image = imrotate(image,n);
imshow(image);
