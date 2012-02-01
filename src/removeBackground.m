function Image = removeBackground(Image, diskSize)
    if nargin == 1, diskSize = 8; end
    Image = imtophat(Image, strel('disk', diskSize));
end