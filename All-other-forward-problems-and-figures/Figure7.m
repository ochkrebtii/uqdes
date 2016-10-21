% MOVIE!!!
% 2D Navier-Stokes pseudo-spectral solver on the torus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navier-Stokes equations in vorticity/stream function formulation on the torus %                                                  %
%                                                                               %
%  Dw/Dt = nu.Laplacian(w)                                                      %
%  Laplacian(psi) = -w                                                          %
%  u = psi_y                                                                    %
%  v =-psi_x                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; 

% set random seed (or comment out to generate one automatically)
RandStream.setGlobalStream(RandStream('mt19937ar','seed',100))


nu                      = 1e-4;                 % low viscosity = high Reynolds number
NX                      = 128;                  % resolution in x
NY                      = 128;                  % resolution in y
dt                      = 1e-3;                 % time step
dx                      = 2*pi/NX;
dy                      = 2*pi/NY;
I                       = sqrt(-1);
TF                      = 800;                 % final time
N                       = floor(TF/dt);         % number of solver steps
B                       = 2;  % draw 6 model realizations of the estimated unknown solution

sr                      = 1;  % subplot rows
sc                      = 2;  %subplot columns

a                       = 0; 
b                       = TF;
t                       = linspace(a,b,N);

% sreen update interval time (NOTE: plotting is usually slow)
TSCREEN                 = 2000;

% define uqdes auxiliary parameters
alpha                   = 1;         % prior precision directly proportional to discretization grid size
lambda                  = 0.5*dt;        % length-scale is proportional to step length

% specify the initial condition type
initial_condition = 'vortices';   % 'vortices' or 'uniform' or 'gaussian'

% Define initial vorticity distribution
switch lower(initial_condition)
    case 'vortices'
        [i,j]       = meshgrid(1:NX,1:NY);
        w           = exp(-((i*dx-pi).^2+(j*dy-pi+pi/4).^2)/(0.2))+exp(-((i*dx-pi).^2+(j*dy-pi-pi/4).^2)/(0.2))-0.5*exp(-((i*dx-pi-pi/4).^2+(j*dy-pi-pi/4).^2)/(0.4));
    case 'uniform'
        w           = random('unif',-1,1,NX,NY);
    case 'gaussian'
        beta        = 5;                        % section 5 NSE paper
        alpha       = 2.2;                      % section 5 NSE paper
        A           = 1;                        % fix this
        w           = 1*randn(NX,NY);           % *beta^2*A^(-alpha);
    otherwise
        disp('Unknown initial conditions !!!');
        return
end


w_hat                   = fft2(w);
w_hat                   = repmat(w_hat,[1,1,B]);

u                       = zeros(NX,NY,B);
v                       = zeros(NX,NY,B);

kx                      = I*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction
ky                      = I*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction

dealias                 = kx<2/3*NX & ky<2/3*NY; % Cutting of frequencies using the 2/3 rule

ksquare_viscous         = kx.^2+ky.^2;        % Laplacian in Fourier space
ksquare_poisson         = ksquare_viscous;
ksquare_poisson(1,1)    = 1;             % fixed Laplacian in Fourier space for Poisson's equation

x = linspace(0,2*pi,NX);
y = linspace(0,2*pi,NY);

kernel     = 'uniform';

switch lower(kernel)
    case 'uniform'
        Ntrim = 1; 
        kern = 'un';
    case 'sqexp'
        Ntrim = N;
        kern = 'se';
    otherwise
        disp('Unknown kernel !!!');
        return
end


% define appropriate kernel convolutions
QQ = str2func(strcat('QQ1d_',kern));
RR = str2func(strcat('RR1d_',kern));
QR = str2func(strcat('QR1d_',kern));
RQ = str2func(strcat('RQ1d_',kern));


Binv  = alpha/RR(t(1),t(1),lambda,a,b);
sapm  = w_hat;

% this is where you would define the forcing function, if any:
%f_hat  = 0;
 f_hat =  fft2(0.005*cos(0.5*repmat(x,NY,1) + 0.5*repmat(y',1,NX)));


%%%%%%%%%%%%%%%%% VIDEO SETUP %%%%%%%%%%%%%%%%%

writerObj = VideoWriter('NavStokes_fast_forced.avi');
open(writerObj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = figure;
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

randnNums       = randn(NX,NY,B,TSCREEN); % Generate random numbers outside of loop

% draw first picture

% warn that the plotting may take some time
disp('Wait while drawing first frame...')

for jj = 1:B
    w = real(ifft2(w_hat(:,:,jj))); % realization
    s = subplot(sr,sc,jj);
    contourf(x,y,w,50); colorbar; shading flat;colormap('jet');
    axis([0,2*pi,0,2*pi])
    axis square
    xlabel('\theta')
    ylabel('\rho')
end

drawnow

for jj = 1:B
    w = real(ifft2(w_hat(:,:,jj))); % realization
    s = subplot(sr,sc,jj);
    contourf(x,y,w,50); colorbar; shading flat;colormap('jet');
    axis([0,2*pi,0,2*pi])
    axis square
    xlabel('\theta')
    ylabel('\rho')
end

drawnow

tic

for k = 1:N-1
    
    tn                          = k*dt;
    qrb                         = QR(t(k+1),t(k),lambda,a,b)*Binv/alpha;
    ss_kp1                      = RR(t(k+1),t(k+1),lambda,a,b)/alpha - RR(t(k+1),t(k),lambda,a,b)*Binv*RR(t(k),t(k+1),lambda,a,b)/(alpha^2);
    %ss_kp1                      = 0;
    sapv                        = QQ(t(k+1),t(k+1),lambda,a,b)/alpha - qrb*RQ(t(k),t(k+1),lambda,a,b)/alpha;
    
    for bb = 1:B
        
        % Compute the stream function and get the velocity and gradient of vorticity
        psi_hat         = -w_hat(:,:,bb)./ksquare_poisson;                  % Solve Poisson's Equation
        u(:,:,bb)       = real(ifft2( ky.*psi_hat));                        % Compute  y derivative of stream function ==> u
        v(:,:,bb)       = real(ifft2(-kx.*psi_hat));                        % Compute -x derivative of stream function ==> v
        w_x             = real(ifft2( kx.*w_hat(:,:,bb) ));                 % Compute  x derivative of vorticity
        w_y             = real(ifft2( ky.*w_hat(:,:,bb) ));                 % Compute  y derivative of vorticity
        
        conv     = u(:,:,bb).*w_x + v(:,:,bb).*w_y;         % evaluate the convective derivative (u,v).grad(w)
        conv_hat = fft2(conv);              % go back to Fourier space
        
        conv_hat = dealias.*conv_hat;   % Spherical dealiasing using the 2/3 rule
        
        sapm(:,:,bb)           = sapm(:,:,bb) + qrb*( nu*ksquare_viscous.*w_hat(:,:,bb) + f_hat - conv_hat );
        
    end
    
    w_hat = sapm + randnNums(:,:,:,mod(k,TSCREEN)+1)*sqrt(sapv);
    
    
    if isnan(sum(sum(real(ifft2(w_hat(:,:,1))))))
        display(['error at iteration ',num2str(k)])
        return
    end
    
    
    % Plotting the vorticity field
    if mod(k,TSCREEN) == 0
        
        % draw some more random numbers
        randnNums       = randn(NX,NY,B,TSCREEN); % Generate random numbers outside of loop
        
        % warn that the plotting may take some time
        disp('Wait while drawing a new frame...')
        
        for jj = 1:B
            % Go back in real space omega in real space for plotting
            w = real(ifft2(w_hat(:,:,jj))); % realization
            s = subplot(sr,sc,jj);
            contourf(x,y,w,50); colorbar; shading flat;colormap('jet');
            axis([0,2*pi,0,2*pi])
            axis square
            xlabel('\theta')
            ylabel('\rho')
        end
        
        drawnow
        
        
        %%%%%%%%%%%%%%%%% VIDEO SETUP %%%%%%%%%%%%%%%%%
        
        frame = getframe(f);
        writeVideo(writerObj,frame);
        writerVideo.FrameRate = 50; %Change the desired frame rate here.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp(['iterations complete: ',num2str(100*k/N),'% itertime: ',num2str(toc/60),' minutes'])
        tic
    end
    
     
end

%%%%%%%%%%%%%%%%% VIDEO SETUP %%%%%%%%%%%%%%%%%

close(writerObj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
