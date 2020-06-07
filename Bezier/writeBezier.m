function [] = writeBezier(data, filename)
% Write Bezier curve in Chrono::Vehicle format
%
% Radu Serban, 9/28/2107

f = fopen(filename, 'w');
fprintf(f, '%d 9\n', data.n);
for i = 1:data.n
   fprintf(f, '%12.4f %12.4f %12.4f ', data.p(i,1), data.p(i,2), data.p(i,3)); 
   fprintf(f, '%12.4f %12.4f %12.4f ', data.in(i,1), data.in(i,2), data.in(i,3)); 
   fprintf(f, '%12.4f %12.4f %12.4f\n', data.out(i,1), data.out(i,2), data.out(i,3)); 
end
fclose(f);
end