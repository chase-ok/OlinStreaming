function Image = removePixelNoise(Image)
    Image = bwareaopen(Image, 4);
end