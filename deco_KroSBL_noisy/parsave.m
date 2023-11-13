function parsave(avg,x,y)
fname = ['./results/noiseless_compare_', num2str(avg),'.mat'];
save(fname,'x','y')
end