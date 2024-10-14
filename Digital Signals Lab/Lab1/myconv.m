function myconv(x, nx, y, ny)
   z = conv(x, y); % convolution of x and y
   nz = nx(1)+ny(1):nx(end)+ny(end); % time index for convolution
   stem(nz, z); % display convolution as a function of time index n
   xlabel('Time index n');
   ylabel('Amplitude');
   title('Convolution of two discrete time sequences');
end
