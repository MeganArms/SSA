function db = convert2double(img)

immax = double(intmax(class(img)));
img = double(img);
db = img./immax;

end