function ptime(time)
h=waitbar(0, 'SPEX Gratings are Moving');
for i=1:10
    waitbar(i/10,h)
    pause(time/10)
    %java.lang.Thread.sleep(12*1000);
end

delete(h)
