function aperture=circular_aperture(npix, diameter, cent_x, cent_y)
  %{
  Returns a 2D aperture of the desired diameter pixels, centered on (cent_x,cent_y) and on support npix X npix
  %}
outside=0;
centered=true;
if centered == true
  cent_x = (npix+1)/2;
  cent_y = (npix+1)/2;
end

x = ((linspace(1,npix,npix)) - cent_x) / (diameter / 2.);
y = ((linspace(1,npix,npix)) - cent_y) / (diameter / 2.);
xx = repmat(x,npix, 1);
yy = repmat(y,npix, 1)';
rho = (xx.^2 + yy.^2).^(0.5);
aperture = ones(size(rho));
aperture(find(rho>1)) = outside; % this is the aperture mask
end