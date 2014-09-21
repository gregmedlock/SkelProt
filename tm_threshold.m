function newImg = tm_threshold( originalImg )
% Applies Otsu's Method for thresholding to the
% original image after denoising and returns 
% the modified binary image. Input must be a 
% grayscale image data matrix.

% 12x12 median filter chosen, optimizes data loss vs. denoising
denoisedImg = wiener2(originalImg, [9 9]);
threshold = graythresh(denoisedImg);
newImg = im2bw(denoisedImg, threshold);


end

