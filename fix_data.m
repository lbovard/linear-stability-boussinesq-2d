N=16;
tempa=dlmread("u_0.dat");
tempa=tempa(:,2);
temp=zeros(N,N);
for i=1:N
	temp(i,:)=tempa(1+(i-1)*N:N*i);
end
dlmwrite("u_0.dat",temp);
tempa=dlmread("v_0.dat");
tempa=tempa(:,2);
temp=zeros(N,N);
for i=1:N
	temp(i,:)=tempa(1+(i-1)*N:N*i);
end
dlmwrite("v_0.dat",temp);
tempa=dlmread("uu.dat");
tempa=tempa(:,2);
temp=zeros(N,N);
for i=1:N
	temp(i,:)=tempa(1+(i-1)*N:N*i);
end
dlmwrite("uu.dat",temp);
tempa=dlmread("vv.dat");
tempa=tempa(:,2);
temp=zeros(N,N);
for i=1:N
	temp(i,:)=tempa(1+(i-1)*N:N*i);
end
dlmwrite("vv.dat",temp);
tempa=dlmread("ww.dat");
tempa=tempa(:,2);
temp=zeros(N,N);
for i=1:N
	temp(i,:)=tempa(1+(i-1)*N:N*i);
end
dlmwrite("ww.dat",temp);
