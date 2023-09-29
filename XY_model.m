%% Vortex Dynamics

clc ; clear ; close all

Nx = 100 ;  
Ny = 100 ;
N = Nx * Ny ;
xrange = 0.2*Nx + 1 : 0.8*Nx ; 
yrange = 0.2*Ny + 1 : 0.8*Ny ;
sc = 1 ;                        % grid scale
imsc = 1;                       % image scale
T = 100000 ;
dt = 1 ;
t_plot = 1000;
init_rand = 180 ;
K = + 1 ;
B = 0 ;
eta = 1; 
M = zeros(T , 1);

vorName = sprintf('Vorsim_Nx%03d_K%02d_eta%02d_B%02d' , Nx, 100*K , 100*eta, 100*B);
mkdir(vorName);
cd(vorName);

[X , Y] = meshgrid(1:Nx , 1:Ny);
[Xq , Yq] = meshgrid(1:sc:Nx * imsc , 1:sc:Ny * imsc);

Theta_init = -180 + 2 * init_rand * rand(N,1) ; 
Theta = Theta_init;
Theta_lattice_init = reshape(Theta_init , [Ny , Nx]) ;

M(1) = mean(Theta) ;

theta_time_init = sprintf('thetaBW_field_00001.tif') ;
imshow(rescale(imresize( Theta_lattice_init , 10 ))) ;
imwrite( rescale(imresize( Theta_lattice_init , 10 ) ) , theta_time_init ) ;

Theta_inter = griddata(X , Y , Theta_lattice_init , Xq , Yq);
A = quiver(Xq , Yq , cosd(Theta_inter) , sind(Theta_inter));
flowname = sprintf('Vector_field_00001.tif');
saveas(A , flowname)

s0 = zeros(N , 1);
st = zeros(N , 1);
sb = zeros(N , 1);
sr = zeros(N , 1);
sl = zeros(N , 1);

for i = 2:N-1
    s0(i) = i; 
    st(i) = Ny * floor((i+1)/Ny) + mod(i+1 , Ny) - Ny * ( floor(i / Ny) - floor((i-1)/Ny) );
    sb(i) = Ny * floor((i-1)/Ny) + mod(i-1 , Ny) + Ny * ( floor((i-1)/Ny) - floor((i-2)/Ny)) ;
    if i ~= N - Ny
        sr(i) = mod(i + Ny , N);
    end
    if i ~= Ny
        sl(i) = mod(i - Ny , N);
    end
end

s0(1) = 1 ; s0(N) = N ;
sr(1) = Ny+1 ;  sl(1) = N-Ny+1 ; st(1) = 2 ; sb(1) = Ny ; 
sr(N) = Ny ; sl(N) = N-Ny ; st(N) = N-Ny+1 ; sb(N) = N-1 ;
sl(Ny) = N ; sr(N-Ny) = N ;

HS0 = sparse(s0 , s0 , -1);
HST = sparse(s0 , st , +1) + HS0;
HSB = sparse(s0 , sb , +1) + HS0;
HSR = sparse(s0 , sr , +1) + HS0;
HSL = sparse(s0 , sl , +1) + HS0;

for tt = 2:dt:T
    Theta = Theta + dt * ( K * ( sind(HST * Theta) + sind(HSB * Theta) + ...
        sind(HSR * Theta) + sind(HSL * Theta) ) - B * sind(Theta) + ...
        eta * (rand(N , 1) - rand(N , 1)) );
    Theta = mod(Theta , 360) ;
    Theta_lattice = reshape(Theta , [Ny , Nx]) ;

    if mod(tt , t_plot) == 0
        theta_time = sprintf('thetaBW_field_%05d.tif' , tt) ;
        imwrite( rescale(imresize( Theta_lattice , 10 ) ) , theta_time ) ;

        Theta_inter = griddata(X , Y , Theta_lattice , Xq , Yq);
        A = quiver(Xq , Yq , cosd(Theta_inter) , sind(Theta_inter));
        flowname = sprintf('Vector_field_%05d.tif' , tt);
        saveas(A , flowname)
    end
    
    M(tt) = mean(cosd(Theta));
    
end

Mt = plot(1:T,M,'.');
saveas(Mt , 'MvsT.tif')

