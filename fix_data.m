%%
% Inputs:  real square matrix data with each entry on a single line
% Outputs: real square matrix data that is MATLAB readable
%%
function fix_data(fname,N)
	tempa=dlmread(fname);
	tempa=tempa(:,2);
	temp=zeros(N,N);
	for i=1:N
		temp(i,:)=tempa(1+(i-1)*N:N*i);
	end
	dlmwrite(fname,temp);
endfunction
