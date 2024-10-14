function myxcorr(x, nx, y, ny)
   [rxy, k] = xcorr(x, y);
   k = k + (min(nx) + min(ny));
   stem(k, rxy);
   xlabel('Time Shift index k');
   ylabel('Amplitude');
   title('Cross-Correlation of two discrete time sequences');
end
