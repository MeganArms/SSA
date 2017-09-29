fit(t,coords(1,:),'poly1')
resids = coords(1,:)-f(t)
err_x = std(resids-mean(resids))
% Instead of fitting ALL the trajectories, subtract the initial position
% from all future positions, average that distance and substract that from
% all displacements from center. Test this.