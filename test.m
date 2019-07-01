load count.dat
c = smooth(count(:));
C1 = reshape(c,24,3);


subplot(3,1,1)
plot(count,'.');
hold on
plot(C1,'-');
title('Smooth C1 (All Data)')