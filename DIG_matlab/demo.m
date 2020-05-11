%% Histograms video

for i=1:860

pause(0.1)
drawnow
bar(info.z_mean(i,:)) 

end
 
%% Covariances video

for i=1:860

pause(0.1)
drawnow
imagesc(info.invcov(:,:,i)) 

end 