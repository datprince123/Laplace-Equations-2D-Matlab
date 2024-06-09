% Dimensions of the simulation grid in x (xdim) and y (ydim) directions
xdim=100;
ydim=100;
% convergence tolerance
tolerance = 1e-5;
% create vector x, y
x=1:1:xdim+1;
y=1:1:ydim+1;

%--------DIRICHLET BOUNDARY CONDITION-----------------%
% x coordinates for boundary
i=1:1:xdim+1;
% Values of omega tested to find the optimal one
omega=0.5:0.1:1.9;
% Initializing initial temperature matrix used for all values of omega
t_init=zeros(xdim+1,ydim+1);
% A temperature of 100°C is applied on one boundary, the remaining boundaries are going to remain at zero
t_init(i,ydim+1)=100; % right boundary
t_init(1,i)=0; % top boundary
t_init(i,1)=0; % left boundary
t_init(xdim+1,i)=0; % bottom boundary
% Matrix of iterations for various values of omega
iter=zeros(1,length(omega));
% Running for loop to compute number of iterations for different omega
for range=1:1:length(omega)
   % Initializing previous (t_prev) and present (t_now) iterations' temperature matrix
   t_now=t_init;
   t_prev=t_init;
  
   % Giving initial difference between t_now and t_prev to start the iterations
   t_prev(2, ydim)=1;
  
   % Iteration loop
   while(max(max(abs(t_now-t_prev)))>tolerance) % Run this until convergence
      
       % Updating previous iteration matrix as the present iteration matrix to continue iterations
       t_prev=t_now;
      
       % Iteration counter increment for each omega
       iter(range)=iter(range)+1;
      
       % Updating
       for i=2:1:xdim
           for j=2:1:ydim
               t_now(i,j)=t_now(i,j)+omega(range)*(((t_now(i-1,j)+t_now(i+1,j)+t_now(i,j-1)+t_now(i,j+1))/4)-t_now(i,j));
           end
       end
      
       % Movie type colour scaled image plot to see how solution progresses for omega=1.9
       if range==length(omega)
           figure(1);
           imagesc(t_now);
           colorbar;
           colormap(jet);
           title(['\fontsize{7}Temperature distribution at iteration no ',int2str(iter(range)),' for omega=1.9 (Dirichlet BC)'],'Color','k');
           xlabel('x-axis','FontSize',10);
           ylabel('y-axis','FontSize',10);
           set(gca,'FontSize',10);
           getframe;
       end
   end
end
% Plot image
figure(2);
surf(x,y,t_now);
%Plotting alpha v/s iterations
figure(3);
plot(omega,iter);
title(['\fontsize{7}Plot of omega v/s no. of iterations on a grid of ',int2str(xdim),' x ',int2str(ydim), ' Dirichlet BC'],'color','k');
xlabel('omega','FontSize',10);
ylabel('No of iterations','FontSize',10);

%--------MIXED BOUNDARY CONDITION (DIRICHLET AND NEUMANN)-----------------%
% x coordinates for boundary
i=1:1:xdim+1;
% Values of omega tested to find the optimal one
omega=0.5:0.1:1.9;
% Initializing initial temperature matrix used for all values of omega
t_init=zeros(xdim+1,ydim+1);
% Boundary conditions
t_init(i,ydim+1)=50; % right boundary (Dirichlet BC)
t_init(1,i)=100; % top boundary (Dirichlet BC)
t_init(i,1)=75; % left boundary (Dirichlet BC)
% bottom boundary is Neumann BC
% Matrix of iterations for various values of omega
iter=zeros(1,length(omega));
% Running for loop to compute number of iterations for different omega
for range=1:1:length(omega)
   % Initializing previous (t_prev) and present (t_now) iterations' temperature matrix
   t_now=t_init;
   t_prev=t_init;
  
   % Giving initial difference between t_now and t_prev to start the iterations
   t_prev(2, ydim)=1;
  
   % Iteration loop
   while(max(max(abs(t_now-t_prev)))>tolerance) % Run this until convergence
      
       % Updating previous iteration matrix as the present iteration matrix to continue iterations
       t_prev=t_now;
      
       % Iteration counter increment for each omega
       iter(range)=iter(range)+1;
      
       % Updating
       for i=2:1:xdim
           for j=2:1:ydim
               t_now(i,j)=t_now(i,j)+omega(range)*(((t_now(i-1,j)+t_now(i+1,j)+t_now(i,j-1)+t_now(i,j+1))/4)-t_now(i,j));
           end
       end
      
       % Update Neumann BC
       t_now(xdim + 1,:) = (4.*t_now(xdim,:)-t_now(xdim-1,:))/3;
       % Movie type colour scaled image plot to see how solution progresses for omega=1.9
       if range==length(omega)
           figure(4);
           imagesc(t_now);
           colorbar;
           colormap(jet);
           title(['\fontsize{7}Temperature distribution on a at iteration no ',int2str(iter(range)),' for an alpha=1.9 (Mixed BC)'],'Color','k');
           xlabel('x-axis','FontSize',10);
           ylabel('y-axis','FontSize',10);
           set(gca,'FontSize',10);
           getframe;
       end
   end
end
% Plot image
figure(5);
surf(x,y,t_now);
%Plotting alpha v/s iterations
figure(6);
plot(omega,iter);
title(['\fontsize{7}Plot of omega v/s no. of iterations on a grid of ',int2str(xdim),' x ',int2str(ydim), ' Mixed BC'],'color','k');
xlabel('omega','FontSize',10);
ylabel('No of iterations','FontSize',10);