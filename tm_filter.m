function [ new_img ] = tm_filter(original_img)
% Applies a Weiner filter to the original image and
% returns the filtered image. Input must be a 
% grayscale image data matrix.

[M,N] = size(original_img);
y = zeros(size(original_img));

r = 50;     % Adjust for desired window size

for n = 1+r:N-r
    for m = 1+r:M-r
        % Extract a window of size (2r+1)x(2r+1) around (m,n)
        w = original_img(m+(-r:r),n+(-r:r));

        if (mean(mean(w))) < 50
            y(m,n) = original_img(m,n);
        else
            y(m,n) = 256;
        end
    end
end

new_img = y;
%new_img = weiner2(original_img, [5, 5])

end

