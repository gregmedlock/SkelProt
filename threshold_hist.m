clear all
close all
clc

img = imread('aTm_S283A_0003.tif');
figure(1)
imhist(img,255);
axis([0 255 0 25E4]);

img_filtered = wiener2(img, [50,50]);
figure(2)
imhist(img_filtered,255);
axis([100 200 0 20E4]);

threshold = graythresh(img);
binary_image = im2bw(img, threshold);
figure(3)
imshow(binary_image)