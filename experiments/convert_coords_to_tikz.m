function convert_coords_to_tikz(x,y);
%takes two vectors x,y as in "plot" and convert their entries
%to a format that can be pasted in tikzplots

n=length(x);
for i=1:n
    disp(sprintf('(%24.16e,%24.16e)',x(i),y(i)));
end