clear all
clc

%%  Complex root finder
% this is a complex root finder based on the Delves & Lynes method


% Author: Rens van Leijden, TU Delft



%% input

f = @(z) z.^4 - 16;
% f1 = @(z) 4*z.^3;
% kk = 2-0.5i;
% kk2 = 0.1-0.2i;
% f = @(z) sqrt(kk2^2-z.^2).*sqrt(z.^2-kk^2).*(z.^15-16)%.*(z.^4-16);
% f = @(z) sqrt(kk2^2-z.^2).*sqrt(z.^2-kk^2)
% f = @(z) sqrt(kk2^2-z.^2).*(z.^4-4);
% f = @(z) z.^15-16
%% input for determinant of elastic layer configuration

c_p = real(1.719084289205639e+03 + 3.322560356224935e+01i);
c_s = real(3.704985415873700e+02 + 1.875674989161015e+01i);
lambda = real(5.114035573342751e+09 + 1.649233839281022e+08i);
mu = real(2.612383107677712e+08 + 2.651871379876653e+07i);

% c_p = 297;
% c_s = 121;
% rho = 1700;
% mu = c_s^2*rho;
% lambda = c_p^2*rho-2*mu;
E = 7e7;
nu = 0.4;
rho = 1700
% lambda = E*nu/((1+nu)*(1-2*nu));
% mu = E/(2*(1+nu));
% c_p = sqrt((lambda+2*mu)/rho)
% c_s = sqrt(mu/rho)
omega = 10*2*pi;
H = 10;
k_p = omega/c_p;
k_s = omega/c_s;


% as matrix from which we take the determinant
% with ux = + psi_y
f = @(k_x) det([-i .* k_x -i .* exp(-i .* sqrt(omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) .* H) .* k_x -i .* sqrt(omega .^ 2 ./ c_s .^ 2 - k_x .^ 2) i .* exp(-i .* sqrt(omega .^ 2 ./ c_s .^ 2 - k_x .^ 2) .* H) .* sqrt(omega .^ 2 ./ c_s .^ 2 - k_x .^ 2); -i .* sqrt(omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) i .* sqrt(omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) .* exp(-i .* sqrt(omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) .* H) i .* k_x i .* exp(-i .* sqrt(omega .^ 2 ./ c_s .^ 2 - k_x .^ 2) .* H) .* k_x; (-omega .^ 2 ./ c_p .^ 2 .* lambda - 2 .* (omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) .* mu) .* exp(-i .* sqrt(omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) .* H) (-2 .* mu - lambda) .* (omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) - k_x .^ 2 .* lambda 2 .* k_x .* exp(-i .* sqrt(omega .^ 2 ./ c_s .^ 2 - k_x .^ 2) .* H) .* sqrt(omega .^ 2 ./ c_s .^ 2 - k_x .^ 2) .* mu -2 .* k_x .* sqrt(omega .^ 2 ./ c_s .^ 2 - k_x .^ 2) .* mu; -2 .* mu .* k_x .* sqrt(omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) .* exp(-i .* sqrt(omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) .* H) 2 .* mu .* k_x .* sqrt(omega .^ 2 ./ c_p .^ 2 - k_x .^ 2) mu .* exp(-i .* sqrt(omega .^ 2 ./ c_s .^ 2 - k_x .^ 2) .* H) .* (2 .* k_x .^ 2 - omega .^ 2 ./ c_s .^ 2) k_x .^ 2 .* mu - (omega .^ 2 ./ c_s .^ 2 - k_x .^ 2) .* mu;]);

% with ux = -psi_y
% f =@(k_x) det([-i * k_x -i * exp(-i * sqrt(omega ^ 2 / c_p ^ 2 - k_x ^ 2) * H) * k_x i * sqrt(omega ^ 2 / c_s ^ 2 - k_x ^ 2) -i * exp(-i * sqrt(omega ^ 2 / c_s ^ 2 - k_x ^ 2) * H) * sqrt(omega ^ 2 / c_s ^ 2 - k_x ^ 2); -i * sqrt(omega ^ 2 / c_p ^ 2 - k_x ^ 2) i * sqrt(omega ^ 2 / c_p ^ 2 - k_x ^ 2) * exp(-i * sqrt(omega ^ 2 / c_p ^ 2 - k_x ^ 2) * H) -i * k_x -i * exp(-i * sqrt(omega ^ 2 / c_s ^ 2 - k_x ^ 2) * H) * k_x; (-omega ^ 2 / c_p ^ 2 * lambda - 2 * (omega ^ 2 / c_p ^ 2 - k_x ^ 2) * mu) * exp(-i * sqrt(omega ^ 2 / c_p ^ 2 - k_x ^ 2) * H) (-2 * mu - lambda) * (omega ^ 2 / c_p ^ 2 - k_x ^ 2) - lambda * k_x ^ 2 -2 * mu * k_x * exp(-i * sqrt(omega ^ 2 / c_s ^ 2 - k_x ^ 2) * H) * sqrt(omega ^ 2 / c_s ^ 2 - k_x ^ 2) 2 * k_x * sqrt(omega ^ 2 / c_s ^ 2 - k_x ^ 2) * mu; -2 * mu * k_x * sqrt(omega ^ 2 / c_p ^ 2 - k_x ^ 2) * exp(-i * sqrt(omega ^ 2 / c_p ^ 2 - k_x ^ 2) * H) 2 * mu * k_x * sqrt(omega ^ 2 / c_p ^ 2 - k_x ^ 2) -mu * exp(-i * sqrt(omega ^ 2 / c_s ^ 2 - k_x ^ 2) * H) * (2 * k_x ^ 2 - omega ^ 2 / c_s ^ 2) -k_x ^ 2 * mu + (omega ^ 2 / c_s ^ 2 - k_x ^ 2) * mu;]);
% 
% % as function directly
% % with ux = -psi_y
f = @(k_x) 32 .* ((((-(0.1e1 ./ 0.16e2) - exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) + sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2))) ./ 16 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) ./ 16 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) ./ 16) .* mu + (-(0.1e1 ./ 0.32e2) - exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) + sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2))) ./ 32 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) ./ 32 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) ./ 32) .* lambda) .* (omega .^ 4) - (-(0.1e1 ./ 0.4e1) - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) ./ 4 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) ./ 4 - exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) + sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2))) ./ 4 + exp(-i .* H .* (sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) + sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2)))) .* (mu .* (c_p .^ 2) + 2 .* (c_s .^ 2) .* (mu + lambda ./ 2)) .* (k_x .^ 2) .* (omega .^ 2) ./ 4 + (-(0.1e1 ./ 0.4e1) - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) ./ 4 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) ./ 4 - exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) + sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2))) ./ 4 + exp(-i .* H .* (sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) + sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2)))) .* (c_p .^ 2) .* (c_s .^ 2) .* mu .* (k_x .^ 4)) .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) + (exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) + exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) - exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) + sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2))) - 1) .* ((omega .^ 4) .* ((0.3e1 ./ 0.4e1) .* mu + lambda ./ 8) - (0.3e1 ./ 0.4e1) .* (mu .* (c_p .^ 2) + (0.4e1 ./ 0.3e1) .* (c_s .^ 2) .* (mu + lambda ./ 4)) .* (k_x .^ 2) .* (omega .^ 2) + (c_p .^ 2) .* (c_s .^ 2) .* (k_x .^ 4) .* mu) .* (k_x .^ 2) ./ 4) .* mu ./ (c_p .^ 2) ./ (c_s .^ 2);
% with ux = + psi_y
% f = @(k_x) 32 .* (sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* (((-exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) + sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2))) ./ 16 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) ./ 16 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) ./ 16 - (0.1e1 ./ 0.16e2)) .* mu + (-exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) + sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2))) ./ 32 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) ./ 32 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) ./ 32 - (0.1e1 ./ 0.32e2)) .* lambda) .* (omega .^ 4) - (-(0.1e1 ./ 0.4e1) - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) ./ 4 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) ./ 4 - exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) + sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2))) ./ 4 + exp(-i .* H .* (sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) + sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2)))) .* (k_x .^ 2) .* (mu .* (c_p .^ 2) + 2 .* (mu + lambda ./ 2) .* (c_s .^ 2)) .* (omega .^ 2) ./ 4 + (-(0.1e1 ./ 0.4e1) - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) ./ 4 - exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) ./ 4 - exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) + sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2))) ./ 4 + exp(-i .* H .* (sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) + sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2)))) .* (k_x .^ 4) .* (c_p .^ 2) .* mu .* (c_s .^ 2)) .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) + (k_x .^ 2) .* (exp(-2.*i .* sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2) .* H) + exp(-2.*i .* sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) .* H) - exp(-2.*i .* H .* (sqrt((-k_x .^ 2 .* c_s .^ 2 + omega .^ 2) ./ c_s .^ 2) + sqrt((-k_x .^ 2 .* c_p .^ 2 + omega .^ 2) ./ c_p .^ 2))) - 1) .* ((omega .^ 4) .* ((0.3e1 ./ 0.4e1) .* mu + lambda ./ 8) - (0.3e1 ./ 0.4e1) .* (k_x .^ 2) .* (mu .* (c_p .^ 2) + (0.4e1 ./ 0.3e1) .* (mu + lambda ./ 4) .* (c_s .^ 2)) .* (omega .^ 2) + (c_p .^ 2) .* (c_s .^ 2) .* (k_x .^ 4) .* mu) ./ 4) .* mu ./ (c_p .^ 2) ./ (c_s .^ 2);
%% plotting for check
dx = 5e-3;
re_min = -max(real(k_s))/0.2;
re_max = 2*max(real(k_s))/0.2;
im_min = -3;
im_max = 50d-6;
% im_max = 0.1;

% meshgrid and final complex value in matrix Z
[X,Y] = meshgrid(re_min:dx:re_max,im_min:dx:im_max);
Z = X + 1i.*Y;

fz = zeros(size(Z));
for ii = 1:size(Z,1)
    for ij = 1:size(Z,2)
        fz(ii,ij) = f(Z(ii,ij));
    end    
end


del = 0;
a   = 1e-6:1e-5:real(k_p);
b   = real(k_p)*imag(k_p)./a;
b0  = (real(k_p)*imag(k_p))*(1+1.e-10)./a;
% 
c   = 1e-6:1e-5:real(k_s);
d   = real(k_s)*imag(k_s)./c;
d0  = (-del+real(k_s)*imag(k_s))*(1-1.e-10)./c;

figure;
contour(real(Z),imag(Z),real(fz),'b--')
hold on 
contour(real(Z),imag(Z),imag(fz),'r-.')
plot(a,b,'--k', 'LineWidth',0.5)
plot(c,d,'--k', 'LineWidth',0.5)
scatter(real(omega/c_p),imag(omega/c_p),'magenta')
scatter(real(omega/c_s),imag(omega/c_s),'cyan')
% surf(X,Y,angle(fz),'EdgeColor','none')
% surf(X,Y,abs(fz),'EdgeColor','none')
ylim([im_min,im_max])
%% plotting for check?
dx = 50e-3;
re_min = -2;
re_max = 2;
im_min = -2;
im_max = 2;

% meshgrid and final complex value in matrix Z
[X,Y] = meshgrid(re_min:dx:re_max,im_min:dx:im_max);
Z = X + 1i.*Y;
% clearvars X Y % delete variables X Y for saving memory

figure
contour(real(Z),imag(Z),real(f(Z)),'b--')
hold on
contour(real(Z),imag(Z),imag(f(Z)),'r-.')
% surf(X,Y,angle(f(Z)),'EdgeColor','none')


%% testing methods

% Z = Vertices_Rectangular(re_min,re_max,im_min,im_max,'Npoints',1e6);
% 
% re_length = re_max-re_min;
% im_length = im_max-im_min;
% % initial circle
% z0 = mean(Z); % middle of the rectangle as z0
% r0 = 0.5*sqrt(re_length^2+im_length^2); % radius of circle equal to diagonal of rectangular / square search domain
% 
% Npoints = length(Z);
% 
% % Args = zeros(Npoints,1);
% % for ii = 1:Npoints-1
% %     Args(ii) = angle(f(Z(ii+1)/f(Z(ii))));
% % end
% % Args(Npoints) = angle(f(Z(1)/f(Z(Npoints))));
% % 
% % Nroots = 1/2/pi*sum(Args)

w_step = min(0.01 * min(re_max - re_min, im_max - im_min), 0.0001);

Nroots = argument_principle(f,re_min, re_max, im_min, im_max, w_step)
Nroots2 = findnumberofroots(f,Z,'option_der','numerical',z0,r0,Npoints);

%% other paper 

% % circle origin and radius
% Npoints = 1e6;
% z0 = 2-2i;
% % z0 = 0-5i;
% % z0 = 5;
% z0 = 0;
% r0 = 1.5;
% % integration parameter
% k = 0:1:Npoints-1;
% 
% Z = z0 + r0*exp(2*pi*1i.*k/Npoints);
% 
% % plot for check
% % figure;
% hold on
% scatter(real(Z),imag(Z))
% hold off
% % title('Search domain circle')
% 
% 
% 
% 
% % n = Npoints;
% % fk = f(Z);
% % c = fft(fk)/n;
% % cp = (1:n-1).*c(2:end);
% % ppzk = n*ifft([cp 0])/r0;
% K = findnumberofroots(f,Z,'option_der','rootsofunity',z0,r0,Npoints)


% s = ifft(ppzk./fk);
% H = hankel(s(2:K+1), s(K+1:2*K));
% H2 = hankel(s(3:K+2), s(K+2:2*K+1));
% w = r0*eig(H2,H)+z0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Circle subdomain - finding and polishing roots at the end

% % search domain
% re_min = -5;
% re_max = 5;
% im_min = -5;
% im_max = 5;
% 
% % plot search domain as rectangular domain inscribing a circle
% Z = Vertices_Rectangular(re_min,re_max,im_min,im_max,1000);
% % figure;
% scatter(real(Z),imag(Z));
% % lengths of domain
% re_length = re_max-re_min;
% im_length = im_max-im_min;
% % initial circle
% z0 = mean(Z); % middle of the rectangle as z0
% r0 = 0.5*sqrt(re_length^2+im_length^2); % radius of circle equal to diagonal of rectangular / square search domain
% k = 0:1:1000-1;
% circ1 = z0 + r0 * exp(2*pi*1i.*k/1000);
% hold on
% h2 = scatter(real(circ1),imag(circ1));
% 
% 
% % amount of zeroes
% Npoints = 1e6;
% Nroots = findnumberofroots(f,Z,'option_der','rootsofunity',z0,r0,Npoints);
% 
% % check if max # of roots are in the initial domain, if so calculate the
% % roots
% Nmax_domain = 2;
% % if Nroots <= Nmax_domain
% %     Npoints = 1e6;
% %     roots_total = findrootsindomain3(f,z0,r0,Npoints,Nroots);
% % % if not then start the while loop
% % else    
%     % intialise the first domain with mid-point, radii
%     Domains = [z0, r0];
%     
%     % domains which will be evaluated further for the roots [r0,z0,Nroots of the domains]
%     SearchDomains = [];    
%     
%     plot = 1;
% 
%     % start the while loop
%     Nfound = 0;
%     while Nfound == 0
%         
%         Domains_new = [];
%         % subdivide the domain(s)
%         for iD = 1:size(Domains,1)
%             new_circles = subdivide_8circles(Domains(iD,1),Domains(iD,2));
%             Domains_new = [Domains_new; ...
%                             new_circles];
% %             % plotting
% %             if plot == 1
% %                 for iCircle = 1:size(new_circles,1)
% %                     plot_circle = new_circles(iCircle,1) + new_circles(iCircle,2)*exp(2*pi*1i.*k/1000);
% %                     scatter(real(plot_circle),imag(plot_circle),'red')
% %                 end
% %             end
%         end
%         
%         Domains = [];
%         % for every domain check the # of roots
%         for iD = 1:size(Domains_new,1)
%             % plotting
%             if plot == 1
%                 delete(h2)
%                 plot_circle = Domains_new(iD,1) + Domains_new(iD,2)*exp(2*pi*1i.*k/1000);
%                 h2 = scatter(real(plot_circle),imag(plot_circle),'red');
%             end
%             % calculate # of roots in the subdomain
%             Nroots_subdomain = findnumberofroots(f,Z,'option_der','rootsofunity',Domains_new(iD,1),Domains_new(iD,2),Npoints);
%             
%             % if the imaginary part is not equal to zero we have crossed a
%             % branch cut and we must subdivide, but if and only if there is
%             % a potential root (thus if real part is 1, 1.5, 2, 2.5, etc.
%             if imag(Nroots_subdomain) ~= 0 && mod(real(Nroots_subdomain),0) - mod(real(Nroots_subdomain),1) ~= 0
%                 Domains = [Domains; ...
%                             Domains_new(iD,:)];                
%             
%             % if we have crossed a branch cut but there are no potential
%             % roots then we don't subdivide (check if imaginary part is not
%             % equal to zero suffices)
%             elseif imag(Nroots_subdomain) ~= 0 
%                 continue
%             
%             % if # of roots in subdomain is smaller or equal to max # of
%             % roots then save domain for root finding for later
%             elseif Nroots_subdomain <= Nmax_domain && Nroots_subdomain > 0 && mod(Nroots_subdomain,1) ~= 0.5
%                 SearchDomains = [SearchDomains; ...
%                                  Domains_new(iD,:), Nroots_subdomain];
%                 if plot == 1
%                     plot_circle = Domains_new(iD,1) + Domains_new(iD,2)*exp(2*pi*1i.*k/1000);
%                     scatter(real(plot_circle),imag(plot_circle),'green');
%                 end
%                 
%             % if not, then check if there are more than Nmax_domain roots
%             % in the domain, then we must subdivide. If not, e.g. equal to
%             % zero, we exclude that subdomain by doing.. nothing (:
%             elseif Nroots_subdomain > Nmax_domain 
%                 Domains = [Domains; ...
%                             Domains_new(iD,:)];
%             
%             % if it gives either NaN or Inf values we know the contour is
%             % over a zero thus we must divide that domain further
%             elseif isnan(Nroots_subdomain) || isinf(Nroots_subdomain)
%                 Domains = [Domains; ...
%                             Domains_new(iD,:)];
%             
%             % if we encounter a branch point we get halve values, when the
%             % Nroots is 0.5 we have no other root, when it is e.g. 1.5
%             % there may be a root so subdivide
%             elseif Nroots_subdomain - mod(Nroots_subdomain,1) ~= 0
%                 Domains = [Domains; ...
%                             Domains_new(iD,:)];
%             end
%         end
%         
%         % end while loop if Domains remains zero in size - SHOULD OPTIMISE
%         % THIS
%         if size(Domains,1) == 0
%             Nfound = 1;
%         end
% 
%     end
% % end    
% 
% %  now use the searchdomains to look for the roots and optimise them
% roots = [];
% for iD = 1:size(SearchDomains,1)
%     roots_found = findrootsindomain3(f,SearchDomains(iD,1),SearchDomains(iD,2),Npoints,SearchDomains(iD,3));
%     roots = [roots; ...
%                     roots_found];
% end
% 
% % remove duplicates
% tol = 1e-6;
% roots = uniquetol([real(roots), imag(roots)],tol, 'ByRows',true);
% 
% % check if the final ones obey the search domain and delete the ones that
% % don't
% roots(find(real(roots) < re_min | real(roots) > re_max | imag(roots) < im_min | imag(roots) > im_max)) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Circle subdomain - finding and polishing roots during the while loop

% search domain
re_min = 0.55;
re_max = 0.565;
im_min = -0.01;
im_max = 0.01;

% plot search domain as rectangular domain inscribing a circle
Z = Vertices_Rectangular(re_min,re_max,im_min,im_max,1000);
% figure;
scatter(real(Z),imag(Z));
% lengths of domain
re_length = re_max-re_min;
im_length = im_max-im_min;
% initial circle
z0 = mean(Z); % middle of the rectangle as z0
r0 = 0.5*sqrt(re_length^2+im_length^2); % radius of circle equal to diagonal of rectangular / square search domain
k = 0:1:1000-1;
circ1 = z0 + r0 * exp(2*pi*1i.*k/1000);
hold on
h2 = scatter(real(circ1),imag(circ1));
h3 = scatter(real(circ1),imag(circ1));

% amount of zeroes
Npoints = 1e6;
Nroots_domain = findnumberofroots(f,Z,'option_der','rootsofunity',z0,r0,Npoints);

% for now lets just say that every time we cross a branch cut we include a
% branch point and thus the total amount of roots must be adjusted
Nroots_domain = real(Nroots_domain) - imag(Nroots_domain)*0.5;

% check if max # of roots are in the initial domain, if so calculate the
% roots
Nmax_domain = 3;
% if Nroots <= Nmax_domain
%     Npoints = 1e6;
%     roots_total = findrootsindomain3(f,z0,r0,Npoints,Nroots);
% % if not then start the while loop
% else    
    % intialise the first domain with mid-point, radii
    Domains = [z0, r0];
    
    % domains which will be evaluated further for the roots [r0,z0,Nroots of the domains]
    SearchDomains = [];    
    
    % option for plotting or not (debugging)
    plot = 1;
    
    % empty array to store the roots (MIGHT THINK ABOUT PREALLOCATING A
    % LARGE  VARIABLE FOR SPEED (NOT SURE IF NECESSARY ATM)
    roots = [];
    % start the while loop
    Nfound = 0;
    while Nfound ~= Nroots_domain
        
        % set / re-set the array for storing the new subdivided search
        % domains
        Domains_new = [];
        % subdivide the domain(s)
        for iD = 1:size(Domains,1)
            new_circles = subdivide_8circles(Domains(iD,1),Domains(iD,2));
            Domains_new = [Domains_new; ...
                            new_circles];
%             % plotting
%             if plot == 1
%                 for iCircle = 1:size(new_circles,1)
%                     plot_circle = new_circles(iCircle,1) + new_circles(iCircle,2)*exp(2*pi*1i.*k/1000);
%                     scatter(real(plot_circle),imag(plot_circle),'red')
%                 end
%             end
        end
        
        % do the same for the old one
        Domains = [];

        % for every domain check the # of roots
        for iD = 1:size(Domains_new,1)
            % plotting
            if plot == 1
                delete(h2)
                plot_circle = Domains_new(iD,1) + Domains_new(iD,2)*exp(2*pi*1i.*k/1000);
                h2 = scatter(real(plot_circle),imag(plot_circle),'red');
            end
            % calculate # of roots in the subdomain
            Nroots_subdomain = findnumberofroots(f,Z,'option_der','rootsofunity',Domains_new(iD,1),Domains_new(iD,2),Npoints);
            
            % if the imaginary part is not equal to zero we have crossed a
            % branch cut and we must subdivide, but if and only if there is
            % a potential root (thus if real part is 1, 1.5, 2, 2.5, etc.
            if imag(Nroots_subdomain) ~= 0 && mod(real(Nroots_subdomain),0) - mod(real(Nroots_subdomain),1) ~= 0
                Domains = [Domains; ...
                            Domains_new(iD,:)];                
            
            % if we have crossed a branch cut but there are no potential
            % roots then we don't subdivide (check if imaginary part is not
            % equal to zero suffices)
            elseif imag(Nroots_subdomain) ~= 0 
                continue
            
            % if # of roots in subdomain is smaller or equal to max # of
            % roots then start the root-finding
            elseif Nroots_subdomain <= Nmax_domain && Nroots_subdomain > 0 && mod(Nroots_subdomain,1) ~= 0.5
                % plot the search domain (for debugging)
                if plot == 1
                    plot_circle = Domains_new(iD,1) + Domains_new(iD,2)*exp(2*pi*1i.*k/1000);
                    scatter(real(plot_circle),imag(plot_circle),'green');
                end
                
                % check whether any previous found roots are within the
                % same domain
                circle_vertices = Domains_new(iD,1) + Domains_new(iD,2)*exp(2*pi*1i.*k/1000);
                
                % set roots_prevfound to an empty array and check for any
                % previous roots
                roots_prevfound = [];
                if any(inpolygon(real(roots),imag(roots),real(circle_vertices),imag(circle_vertices)))
                    % amount of previously founded roots within this domain
                    Nprev_found = length(roots(inpolygon(real(roots),imag(roots),real(circle_vertices),imag(circle_vertices))));
                    % if there are no roots left to be found within the
                    % domain continue with the next iteration
                    if Nroots_subdomain - Nprev_found == 0
                        continue
                    else
                        % otherwise assign the previous found roots to a
                        % vector to use for local deflation
                        roots_prevfound = roots(inpolygon(real(roots),imag(roots),real(circle_vertices),imag(circle_vertices)));
                    end
                end
                
                % find the roots in the sub-domain
                roots_subdomain = findrootsindomain3(f,Domains_new(iD,1),Domains_new(iD,2),Npoints,Nroots_subdomain,roots_prevfound);
                
                % check whether any of the found roots are outside the
                % largest circle (which is larger than the real search
                % domain) and delete those, also make Nfound smaller based
                % on this
                statement = inpolygon(real(roots_subdomain),imag(roots_subdomain),real(circ1),imag(circ1));
                if any(statement) == 0
                    roots_subdomain(statement==0) = [];
                end
                    
                % once finished add the total number of uniquely found
                % roots to Nfound and the total found roots list
                Nfound = Nfound + length(roots_subdomain);
                roots = [roots; ...
                         roots_subdomain];
                
                % plot if wanted as well
                if plot == 1
                    scatter(real(roots_subdomain),imag(roots_subdomain))
                end
                
                % stop the whole for loop of the domain if all roots have
                % been found
                if Nfound == Nroots_domain
                    break
                end

            % if not, then check if there are more than Nmax_domain roots
            % in the domain, then we must subdivide. If not, e.g. equal to
            % zero, we exclude that subdomain by doing.. nothing (:
            elseif Nroots_subdomain > Nmax_domain 
                Domains = [Domains; ...
                            Domains_new(iD,:)];
            
            % if it gives either NaN or Inf values we know the contour is
            % over a zero thus we must divide that domain further
            elseif isnan(Nroots_subdomain) || isinf(Nroots_subdomain)
                Domains = [Domains; ...
                            Domains_new(iD,:)];
            
            % if we encounter a branch point we get halve values, when the
            % Nroots is 0.5 we have no other root, when it is e.g. 1.5
            % there may be a root so subdivide
            elseif Nroots_subdomain - mod(Nroots_subdomain,1) ~= 0
                Domains = [Domains; ...
                            Domains_new(iD,:)];
            end
        end
        
%         % end while loop if Domains remains zero in size - SHOULD OPTIMISE
%         % THIS
%         if size(Domains,1) == 0
%             Nfound = 1;
%         end

    end
% end    


% check if the final ones obey the search domain and delete the ones that
% don't
roots(find(real(roots) < re_min | real(roots) > re_max | imag(roots) < im_min | imag(roots) > im_max)) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Rectangle

% % set the boundaries of the domain and discretise the 4 vertices
% % z = re + i*im
% re_min = -5;
% re_max = 5;
% im_min = 0;
% im_max = 5;
% 
% 
% Z = Vertices_Rectangular(re_min,re_max,im_min,im_max,Npoints);
% 
% % figure;
% % plot(Z,abs(f(Z)))
% 
% % % plot for check
% % figure;
% % scatter(real(Z),imag(Z))
% % title('Search domain rectangle')
% 
% % % construct circle that is around the rectangle/square
% % z0 = mean(Z);
% % x = max(real(Z))-min(real(Z));
% % y = max(imag(Z))-min(imag(Z));
% % r0 = max(x,y);
% % t = linspace(0,1,Npoints);
% % t = t(1:end-1);
% % Z = z0 + r0/2*exp(2*pi*1i.*t);
% % 
% % hold on 
% % scatter(real(Z),imag(Z))
% 
% % amount of zeroes
% N = findnumberofroots(f,Z,'df',f1)
% 
% % finding the roots
% roots_m = findrootsindomain2(f,Z,N,'df',f1)
% 
% 
% 
% 
% 
% % Z = linspace(-5,5,1e5);
% % a = f(Z);
% % b = a.*exp(-floor(log(max(abs(a)))));

%% rectangle subdomain stuff

% re_min = -1;
% re_max = 1;
% im_min = -1;
% im_max = -0.1;


Npoints = 1e5;

% generate first large rectangle
Z = Vertices_Rectangular(re_min,re_max,im_min,im_max);
% plot for check
figure;
scatter(real(Z),imag(Z))
title('Search domain rectangle')


% Nroots = findnumberofroots(f,Z,'df',f1);
w_step = min(0.01 * min(re_max - re_min, im_max - im_min), 0.0001);
Nroots = argument_principle(f,re_min, re_max, im_min, im_max, w_step)

% if original domain only contains Nmax_domain roots then just search for
% them immediately..
Nmax_domain = 1;
if Nroots <= Nmax_domain
%     Roots_total = findrootsindomain2(f,Z,Nroots,'df',d1);
    Roots_total = findrootsindomain2(f,Z,Nroots);
    
else
    % if not then start the while loop
    
    % start with the original large domain [re_min, re_max, im_min, im_max, # total roots in total domain]
    Domains = [re_min,re_max,im_min,im_max, Nroots];

    % value by which to divide the domains
    n_divide = 3;
    
    % domains which will be evaluated further for the roots [# of boxes,    4x boundaries]
    SearchDomains = [];
    
    Nfound = 0;
    while Nfound ~= Nroots

        % loop over each domain, starting of course with the original
        % domain. First determine amount of domains
        nB = size(Domains,1);
        % initialise variable for holding the new Domains already (thus
        % resetting)
        Domains_new = [];
        % now loop over all the current domains
        for inB = 1:nB

            % get boundaries of that domain and calculate the lengths
            re_min = Domains(inB,1);
            re_max = Domains(inB,2);
            im_min = Domains(inB,3);
            im_max = Domains(inB,4);
            length_re = re_max - re_min;
            length_im = im_max - im_min;
            
            % Divide them by arbitrary number n_divide based on the
            % largest side
%             LargestLength = max(length_re,length_im);
%             ind_LargestLength = find([length_re,length_im]==LargestLength);
%             if ind_LargestLength == 1
%                 length_re = length_re/n_divide;
%             else
%                 length_im = length_im/n_divide;
%             end
            % divide all, thus we get only new squares and no rectangles
            length_re = length_re/n_divide;
            length_im = length_im/n_divide;
            
            % determine how many new rectangles there are
            nB_re = (re_max-re_min)/length_re;
            nB_im = (im_max-im_min)/length_im;
            nB_tot = nB_re*nB_im;  
            % cell to not fuck up sizes later on (lazy code rewriting tbh)
            Domains_divided = cell(nB_re,nB_im);
            
            % now create vertices [# real boxes, # imag boxes]
            vertices = cell(nB_re,nB_im); %zeros(nB_re,nB_im,4);
            for inB_re = 1:nB_re
                for inB_im = 1:nB_im
                    % real boundaries
                    Bre_min = re_min + (inB_re-1) * length_re;
                    Bre_max = re_min + inB_re * length_re;
                    % imaginary boundaries
                    Bim_min = im_min + (inB_im-1) * length_im;
                    Bim_max = im_min + inB_im * length_im;
    %                 % store them
                    Domains_divided{inB_re,inB_im} = [Bre_min, Bre_max, Bim_min, Bim_max];
                    % calculate vertices directly
                    vertices{inB_re,inB_im} = Vertices_Rectangular(Bre_min,Bre_max,Bim_min,Bim_max);
                    plot = 1;
                    if plot == 1
                        hold on
                        scatter(real(vertices{inB_re,inB_im}),imag(vertices{inB_re,inB_im}))
                    end    
                end
            end

            % loop over the new boxes and find the number of roots inside each
            % one of them - IMPLEMENT LATER TO SKIP THE BOXES THAT ALREADY HAVE
            % Nmax_domain INSIDE OF THEM
            
            for inB_re = 1:nB_re
                for inB_im = 1:nB_im
                    % original
%                     Nroots_subdomains = findnumberofroots(f,vertices{inB_re,inB_im},'df',f1);
                    % argument principle
                    Bre_min = Domains_divided{inB_re,inB_im}(1);
                    Bre_max = Domains_divided{inB_re,inB_im}(2);
                    Bim_min = Domains_divided{inB_re,inB_im}(3);
                    Bim_max = Domains_divided{inB_re,inB_im}(4);
                    w_step = min(0.01 * min(Bre_max - Bre_min, Bim_max - Bim_min), 0.0001);
                    Nroots_subdomains = argument_principle(f,Bre_min, Bre_max, Bim_min, Bim_max, w_step);
                    % statement to save the domain that satisfies the condition for search later
                    if Nroots_subdomains <= Nmax_domain && Nroots_subdomains ~= 0
                        % store the domains boundaries for saving memory
                        SearchDomains = [SearchDomains; ...
                                        Domains_divided{inB_re,inB_im} Nroots_subdomains];

                        % add to Nfound the amount of roots in this domain
                        Nfound = Nfound + Nroots_subdomains;

                    % statement to check whether we need to keep subdivide current domain 
                    elseif Nroots_subdomains > Nmax_domain 
                        Domains_new = [Domains_new; ...
                                        Domains_divided{inB_re,inB_im} Nroots_subdomains];
                    
                    % statement to check whether we don't have to search further in the current domain    
                    elseif Nroots_subdomains == 0   
                           % do nothing?
                    end
                end
            end
        end
        Domains = Domains_new;
    end
end    

% now loop over Searchdomains to find the roots
Roots_total = [];
for iNsearchdomains = 1:size(SearchDomains,1)
    % amount of roots in specific searchdomain
    Nroots_domain = SearchDomains(iNsearchdomains,end);
    % calculate vertices
    re_min = SearchDomains(iNsearchdomains,1);
    re_max = SearchDomains(iNsearchdomains,2);
    im_min = SearchDomains(iNsearchdomains,3);
    im_max = SearchDomains(iNsearchdomains,4);
    Z = Vertices_Rectangular(re_min,re_max,im_min,im_max);
    % calculate roots
%     Roots = findrootsindomain2(f,Z,Nroots_domain,'df',f1);
    % old method with newton polynomial
%     Roots = findrootsindomain2(f,Z,Nroots_domain);
    % Peter & Kravanja eigenvalue problem
    z0 = mean(Z);
    x = max(real(Z))-min(real(Z));
    y = max(imag(Z))-min(imag(Z));
    r0 = max(x,y);
    Roots = findrootsindomain3(f,z0,r0,Npoints,Nroots_domain)
    % assign roots
    Roots_total = [Roots_total; ...
                    Roots];
end
Roots_total
hold on
scatter(real(Roots_total),imag(Roots_total))
%% old rectangle stuff


% re_min = -7;
% re_max = 5;
% im_min = -7;
% im_max = 5;
% 
% Npoints = 1e5;
% 
% % generate first large rectangle
% Z = Vertices_Rectangular(re_min,re_max,im_min,im_max);
% % plot for check
% figure;
% scatter(real(Z),imag(Z))
% title('Search domain rectangle')
% 
% Nmax_domain = 1; % maximum # of roots in domain before starting to search for the roots themselves
% Nroots = findnumberofroots(f,Z,'df',f1);
% Nfound = 0;
% 
% % if original domain only contains Nmax_domain roots then just search for
% % them immediately..
% % if Nroots == Nmax_domain
% %     roots_m = findrootsindomain2(f,Z,Nroots,'df',d1);
% % else
%     % if not then start the while loop
%     
%     % start with determining the length of the large domain
% 
%     
%     Domains = {[re_min,re_max,im_min,im_max]};
% 
%     % value by which to divide
%     n_divide = 2;
%     
%     % domains which will be evaluated further for the roots [# real boxes, # imag boxes]
%     SearchDomains = [];
% 
% %     while Nfound ~= Nroots
% 
%         % loop over each domain, starting of course with the original
%         % domain. First determine amount of domains
%         [nB_im,nB_re] = size(Domains);
%         
% 
% 
%         % now loop over all the current domains
%         for inB_re = 1:nB_re
%             for inB_im = 1:nB_im
%                 
%                 % get boundaries of that domain and calculate the lengths
%                 re_min = Domains{inB_re,inB_im}(1);
%                 re_max = Domains{inB_re,inB_im}(2);
%                 im_min = Domains{inB_re,inB_im}(3);
%                 im_max = Domains{inB_re,inB_im}(4);
%                 length_re = re_max - re_min;
%                 length_im = im_max - im_min;
%                 
%                 % Divide them by arbitrary number n_divide based on the
%                 % largest side
%                 LargestLength = max(length_re,length_im);
%                 ind_LargestLength = find([length_re,length_im]==LargestLength);
%                 if ind_LargestLength == 1
%                     length_re = length_re/n_divide;
%                 else
%                     length_im = length_im/n_divide;
%                 end
%                 
%                 % determine how many new rectangles there are
%                 nB_re_new = (re_max-re_min)/length_re;
%                 nB_im_new = (im_max-im_min)/length_im;
%                 nB_tot_new = nB_re_new*nB_im_new;  
%                 Domains_divided = cell(nB_re_new,nB_im_new);
%                 
%                 % now create vertices [# real boxes, # imag boxes]
%                 vertices = cell(nB_re_new,nB_im_new); %zeros(nB_re,nB_im,4);
%                 for inB_re_new = 1:nB_re_new
%                     for inB_im_new = 1:nB_im_new
%                         % real boundaries
%                         Bre_min = re_min + (inB_re_new-1) * length_re;
%                         Bre_max = re_min + inB_re_new * length_re;
%                         % imaginary boundaries
%                         Bim_min = im_min + (inB_im_new-1) * length_im;
%                         Bim_max = im_min + inB_im_new * length_im;
%         %                 % store them
%                         Domains_divided{inB_re_new,inB_im_new} = [Bre_min, Bre_max, Bim_min, Bim_max];
%                         % calculate vertices directly
%                         vertices{inB_re_new,inB_im_new} = Vertices_Rectangular(Bre_min,Bre_max,Bim_min,Bim_max);
%                         plot = 1;
%                         if plot == 1
%                             hold on
%                             scatter(real(vertices{inB_re}),imag(vertices{inB_im}))
%                         end    
%                     end
%                 end
% 
%                 % loop over the new boxes and find the number of roots inside each
%                 % one of them - IMPLEMENT LATER TO SKIP THE BOXES THAT ALREADY HAVE
%                 % Nmax_domain INSIDE OF THEM
%                 Nroots_subdomains = zeros(nB_re_new,nB_im_new);
%                 
%                 for inB_re_new = 1:nB_re_new
%                     for inB_im_new = 1:nB_im_new
%                         
%                         Nroots_subdomains(inB_re_new,inB_im_new) = findnumberofroots(f,vertices{inB_re_new,inB_im_new},'df',f1);
%                         
%                         % statement to save the domain that satisfies the condition for search later on and skip over that specific box later on
%                         if Nroots_subdomains(inB_re_new,inB_im_new) <= Nmax_domain && Nroots_subdomains(inB_re_new,inB_im_new) ~= 0
%                             % store the domains boundaries for saving
%                             % memory
%                             SearchDomains = [SearchDomains; ...
%                                             Domains_divided{inB_re_new,inB_im_new}];
% %                             SearchDomains = {SearchDomains; ...
% %                                              {vertices{inB_re_new,inB_im_new}}};
%                         
%                         % statement to check whether we need to keep subdivide current domain 
%                         elseif Nroots_subdomains(inB_re_new,inB_im_new) > Nmax_domain 
%                             Domains_new{inB_re_new,inB_im_new} = Domains_divided{inB_re_new,inB_im_new};
%                         
%                         % statement to check whether we don't have to search further in the current domain    
%                         elseif Nroots_subdomains(inB_re_new,inB_im_new) == 0   
%                                % do nothing?
%                         end
%                     end
%                 end
% 
% 
% 
%             end
%         end
% 
% %  end
        


%% functions to be used

function Rvert = Vertices_Rectangular(min_re,max_re,min_im,max_im,varargin)
    
    % check whether Npoints is given or not
    if any(strcmp(varargin,'Npoints'))
        Npoints = varargin{find(strcmp(varargin,'Npoints'))+1};    
    else
        Npoints = 1e5;
    end

    % define direction, standard is anti-clockwise
    if any(strcmp(varargin,'direction'))
        direction = varargin{find(strcmp(varargin,'direction'))+1};    
    else
        direction = 'anti-clockwise';
    end
    
    
    % vertex 1: top right to top left
    vert1 = linspace(max_re+1i*max_im,min_re+1i*max_im,Npoints);
    % vertex 2: top left to bottom left
    vert2 = linspace(min_re+1i*max_im,min_re+1i*min_im,Npoints);
    % vertex 3: bottom left to bottom right
    vert3 = linspace(min_re+1i*min_im,max_re+1i*min_im,Npoints);
    % vertex 4: bottom right to top right
    vert4 = linspace(max_re+1i*min_im,max_re+1i*max_im,Npoints);
    % everything together
    Rvert = [vert1 vert2 vert3 vert4];

    % if clockwise then flip everything
    if strcmp(direction,'clockwise')
        Rvert = flip(Rvert);
    end

end

%%
% function R = findrootsindomain(f,Z,N,f1,Npoints,option_der,type_optimise)
%     arguments
%         f
%         Z (1,:) % path along to integrate (can be anything) -> not true should change this
%         N (1,1) % number of roots in the subdomain
%         f1 = 0
%         Npoints (1,1) double = 1000 
%         option_der (1,:) char {mustBeMember(option_der,{'symbolic','special1','special2','numerical'})} = 'symbolic'
%         type_optimise (1,:) char {mustBeMember(type_optimise,{'Matlab','Alternate_Newton_Rhapson'})} = 'Matlab'
%     end
% 
%     % if derivative is not given use second option for derivative
%     % determination
%     if isa(f1,'function_handle')
%         % do nothing
%     else
%         option_der = 'numerical';
%     end
%     
%     % calculate N values of s_N
%     s_N = zeros(1,N);  
%     if strcmp(option_der,'symbolic') 
%         
%         for iN = 1:N
%             s_N(iN) = trapz(Z,(Z.^iN) .* f1(Z)./f(Z)) / (2*pi*1i);
%         end
%     
%     elseif strcmp(option_der,'special1')    
%         % or calculate by other way, not using the derivative       # NOT
%         % ROBUST FOR NOW
%         for iN = 1:N
%             a = 0;
%             s_N(iN) = -iN/(2*pi*1i) * trapz(Z,Z.^(iN-1).*log10( ((Z-a).^(-N)) .* f(Z) )) + N*a^iN;
%         end
%     
%     elseif strcmp(option_der,'numerical')
%         % or numerically calculating the derivative
%         dZ = diff(Z);
%         df = diff(f(Z));
%         dfdz = df./dZ;
%         for iN = 1:N
%             s_N(iN) = trapz(Z(1:end-1),(Z(1:end-1).^iN) .* dfdz./f(Z(1:end-1))) / (2*pi*1i);
%         end
%     
%     elseif strcmp(option_der,'special2') % based on subs of w = f(z)    $ NOT WORKING FOR NOW
%         % new integration variable 
%         Z2 = f(Z);
%         for iN = 1:N
%             s_N(iN) = trapz(Z2,(Z.^iN) .* 1./f(Z2)) / (2*pi*1i);
%         end
%     end
%     
%     % calculate the polynomial coefficients (Newtons Identities)
%     e_N = ones(1,N+1);      % thus contains e0, e1, .. eN where eN = e_N(N+1)
%     for iN = 1:N
%         e = 0;
%         for ii = 1:iN
%             e = e + (-(-1)^ii)*e_N(iN+1-ii)*s_N(ii); 
%         end
%         e_N(iN+1) = (1/iN)*e;
%     end
%     
%     % create polynomial
%     p = ones(1,N+1);
%     for iN = 1:N
%         p(iN+1) = (-1)^iN * e_N(iN+1);
%     end
%     
%     
%     % calculate roots of the polynomial
%     roots_p = roots(p)
%     
%     % calculate final roots with polynomial roots as first guesses:
%     R = zeros(1,N);
%     
%     
%     if strcmp(type_optimise,'Matlab')
%         R = fsolve(f,roots_p);
%     elseif strcmp(type_optimise,'Alternate_Newton_Rhapson')
%         for iN = 1:N
%             R(iN) = root_finder_func(f,[min_re,max_re],[min_im,max_im],roots_p(iN));
%         end
%     end
% 
% 
% end

%%
function R = findrootsindomain2(f,Z,N,varargin)
    % define dependent on varargin the input parameters
    
    % check for symbolic derivative of the function
    if any(strcmp(varargin,'df'))
        f1 = varargin{find(strcmp(varargin,'df'))+1};    
        option_der = 'symbolic';     % symbolic derivative has been given
    else
        option_der = 'numerical';    % second choice for derivative of function
    end
    
    % check for which way to calculate derivative
    if any(strcmp(varargin,'option_der'))
        option_der = varargin{find(strcmp(varargin,'option_der'))+1};    % if option is given then use that
    elseif strcmp(option_der,'symbolic')
        % do nothing as symbolic has been chosen either by supplying df or
        % choosing it above
    else    
        option_der = 'numerical';    % second choice for derivative of function
    end
    
    % check for Npoints
    if any(strcmp(varargin,'Npoints'))
        Npoints = varargin{find(strcmp(varargin,'Npoints'))+1};
    else
        Npoints = 1e4; % standard a 10000 points.. must be enough I suppose
    end
    
    % check type to optimise the roots
    if any(strcmp(varargin,'type_optimise'))
       type_optimise = varargin{find(strcmp(varargin,'type_optimise'))+1};
    else
        type_optimise = 'Alternate_Newton_Rhapson'; 
    end
    % check for range_Re and range_Im
    % or maybe not, for now let it depend on Z
    min_re = min(real(Z));
    max_re = max(real(Z));
    min_im = min(imag(Z));
    max_im = max(imag(Z));

    
    
    % calculate N values of s_N
    s_N = zeros(1,N);  
    
    if strcmp(option_der,'symbolic') 
        
        for iN = 1:N
            s_N(iN) = trapz(Z,(Z.^iN) .* f1(Z)./f(Z)) / (2*pi*1i);
        end
    
    elseif strcmp(option_der,'special1')    
        % or calculate by other way, not using the derivative       # NOT
        % ROBUST FOR NOW

        % just take the midpoint, which is the mean of the whole set of Z
        a = mean(Z);
        for iN = 1:N
            s_N(iN) = -iN/(2*pi*1i) * trapz(Z,Z.^(iN-1).*log10( ((Z-a).^(-N)) .* f(Z) )) + N*a^iN;
        end
    
    elseif strcmp(option_der,'numerical')
        % or numerically calculating the derivative
        dfdz = diff(f(Z))./diff(Z);
        for iN = 1:N
            s_N(iN) = trapz(Z(1:end-1),(Z(1:end-1).^iN) .* dfdz./f(Z(1:end-1))) / (2*pi*1i);    % and taking 1 point less because of the function diff 
        end
    
    elseif strcmp(option_der,'special2') % based on subs of w = f(z)    $ NOT WORKING FOR NOW
        % new integration variable 
        Z2 = f(Z);
        for iN = 1:N
            s_N(iN) = trapz(Z2,(Z.^iN) .* 1./f(Z2)) / (2*pi*1i);
        end
    end
    
    
    % calculate the polynomial coefficients (Newtons Identities)
    e_N = ones(1,N+1);      % thus contains e0, e1, .. eN where eN = e_N(N+1)
    for iN = 1:N
        e = 0;
        for ii = 1:iN
            e = e + (-(-1)^ii)*e_N(iN+1-ii)*s_N(ii); 
        end
        e_N(iN+1) = (1/iN)*e;
    end
    
    % create polynomial final coefficients
    p = ones(1,N+1);
    for iN = 1:N
        p(iN+1) = (-1)^iN * e_N(iN+1);
    end
    
    
    % calculate roots of the polynomial
    roots_p = roots(p);
    
    % if N is not equal to 1 then optimise because we have the polynomial
    if N ~= 1
        % calculate final roots with polynomial roots as first guesses:
        R = zeros(N,1);
            
        if strcmp(type_optimise,'Matlab')
            R = fsolve(f,roots_p);
    
        elseif strcmp(type_optimise,'Alternate_Newton_Rhapson')
            for iN = 1:N
                R(iN) = root_finder_func(f,[min_re,max_re],[min_im,max_im],roots_p(iN));
            end
        end
    
    % if it is equal to one the roots is just the root of that polynomial
    % which apparently exact (:
    else
        R = roots_p;
    end
end
%%
function R = findrootsindomain3(f,z0,r0,Npoints,Nroots,roots_prevfound)
    % if no roots_prevfound is given set it to none (to circumvent creating
    % another function)
    if nargin < 6
        roots_prevfound = [];
    end

    % determine Npoints based on a step size, aka circumference divided
    % by the step size. however keep it to a maximum value 
    step_size = 1e-6;
    Npoints = round(2*pi*r0/step_size);
    if Npoints > 1e4
        Npoints = 1e4;
    end    

    % set some values
    n = Npoints;
    K = Nroots;
    
    k = [0:1:Npoints-1];
    Z = z0 + r0*exp(2*pi*1i.*k/Npoints);
    
    % calculate values f_k at points of Z
    fk = f(Z);
    % calculate Taylorcoefficients
    c = fft(fk)/n;
    % Taylor coefficients of fk'
    cp = (1:n-1).*c(2:end);
    % calculate fk'
    ppzk = n*ifft([cp 0])/r0;       % IT SEEMS THAT IF WE ARE CROSSING A BRANCH CUT PPZK GIVES VERY INACCURATE RESULTS
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if length(roots_prevfound) == 0
    % method of Austin & Kravanja:

    % root finding as eigenvalue problem
    % the moments
    s = ifft(ppzk./fk);

    % local deflation if there are any previously found roots - DOESNT WORK
    % FOR this specific method
%     if length(roots_prevfound) ~= 0
%         % change K to the new amount of roots (thus less)
%         K = K - length(roots_prevfound);
%         % loop over # of necessary s values
%         for iK = 2:(2*K+1)
%             s(iK) = s(iK) - sum(((roots_prevfound-z0)).^(iK-1));
% %             s(iK) = s(iK) - sum((roots_prevfound.^(iK-1)-z0)./r0);
%         end
%     end
    
    % Hankel matrices
    H = hankel(s(2:K+1), s(K+1:2*K));
    H2 = hankel(s(3:K+2), s(K+2:2*K+1));
    % scale with r0 and translate by z0     
    w = r0*eig(H2,H)+z0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
    % method of Delves & Lynes
    
%     K = K + length(roots_prevfound);

    % calculate N values of s_N
    s_N = zeros(1,K);  

    % numerically calculating the derivative
    dfdz = diff(f(Z))./diff(Z);
    for iN = 1:K
        s_N(iN) = trapz(Z,(Z.^iN) .* ppzk./fk) / (2*pi*1i);    % and taking 1 point less because of the function diff 
    end

    if length(roots_prevfound) ~= 0
        % change K to the new amount of roots (thus less)
        K = K - length(roots_prevfound);
        % loop over # of necessary s values
        for iK = 1:K
            s_N(iK) = s_N(iK) - sum(roots_prevfound.^(iK));
        end
    end

    % calculate the polynomial coefficients (Newtons Identities)
    e_N = ones(1,K+1);      % thus contains e0, e1, .. eN where eN = e_N(N+1)
    for iN = 1:K
        e = 0;
        for ii = 1:iN
            e = e + (-(-1)^ii)*e_N(iN+1-ii)*s_N(ii); 
        end
        e_N(iN+1) = (1/iN)*e;
    end
    
    % create polynomial final coefficients
    p = ones(1,K+1);
    for iN = 1:K
        p(iN+1) = (-1)^iN * e_N(iN+1);
    end
    
    
    % calculate roots of the polynomial
    w = roots(p);
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % refining the find roots

    % set type_optimise to 'Matlab' for now
    type_optimise = 'Matlab';
%     type_optimise = 'Alternate_Newton_Rhapson';
    % if N is not equal to 1 then optimise because we have the polynomial
%     if Nroots ~= 1
        % calculate final roots with polynomial roots as first guesses:
        R = zeros(Nroots,1);
            
        if strcmp(type_optimise,'Matlab')
            % old, not reliable
%             R = fsolve(f,w);
            % Based on Nelder Mead algorithm
            R = fminsearch(f,w);
    
        elseif strcmp(type_optimise,'Alternate_Newton_Rhapson')
            for iN = 1:Nroots
                R(iN) = root_finder_func(f,[min_re,max_re],[min_im,max_im],w(iN));
            end
        end
    
    % if it is equal to one the roots is just the root of that polynomial
    % which apparently exact (:
%     else
%         R = w;
%     end
% R = w;
end

%%
function Nroots = findnumberofroots(f,Z,varargin)
    
    % check for symbolic derivative of the function
    if any(strcmp(varargin,'df'))
        f1 = varargin{find(strcmp(varargin,'df'))+1};    
        option_der = 'symbolic';     % symbolic derivative has been given
    else
        option_der = 'numerical';    % second choice for derivative of function
    end
    
    % check for which way to calculate derivative
    if any(strcmp(varargin,'option_der'))
        option_der = varargin{find(strcmp(varargin,'option_der'))+1};    % if option is given then use that
    elseif strcmp(option_der,'symbolic')
        % do nothing as symbolic has been chosen either by supplying df or
        % choosing it above
    else    
        option_der = 'numerical';    % second choice for derivative of function
    end
    
    % calculate number of roots
    if strcmp(option_der,'symbolic') 
        Nroots = trapz(Z,f1(Z)./f(Z)) / (2*pi*1i);

    elseif strcmp(option_der,'numerical')
        % or numerically calculating the derivative
        dfdz = diff(f(Z))./diff(Z);

        Nroots = trapz(Z(1:end-1),dfdz./f(Z(1:end-1))) / (2*pi*1i);    % and taking 1 point less because of the function diff 
    elseif strcmp(option_der,'argp')
        % Follow derivation Peter & Kravanja for derivativeless solution
        Npoints = varargin{find(strcmp(varargin,'argp'))+3};
        Args = zeros(Npoints,1)
        for ii = 1:Npoints-1
            Args(ii) = angle(f(Z(ii+1)/f(Z(ii))));
        end
        Args(Npoints) = angle(f(Z(1)/f(Z(Npoints))));
        
        Nroots = 1/2/pi*sum(Args)

    elseif strcmp(option_der,'rootsofunity')
        z0 = varargin{find(strcmp(varargin,'rootsofunity'))+1};
        r0 = varargin{find(strcmp(varargin,'rootsofunity'))+2};
        Npoints = varargin{find(strcmp(varargin,'rootsofunity'))+3};
        
        % determine Npoints based on a step size, aka circumference divided
        % by the step size. however keep it to a maximum value 
        step_size = 1e-6;
        Npoints = round(2*pi*r0/step_size);
        if Npoints > 1e4
            Npoints = 1e4;
        end

        k = 0:1:Npoints-1;
        Z = z0 + r0*exp(2*pi*1i.*k/Npoints);
        n = Npoints;
        fk = f(Z);
        c = fft(fk)/n;
        cp = (1:n-1).*c(2:end);
        ppzk = n*ifft([cp 0])/r0;
        
        % by Austin / Kravanja
        Nroots = (real(mean(Z.*ppzk./fk)));
        % Delves and lynes but new paper approach of derivative
%         Nroots = trapz(Z,ppzk./fk) / (2*pi*1i);
        % fully numerical..
%         Nroots = trapz(Z(1:end-1),(diff(fk)./diff(Z))./fk(1:end-1)) / (2*pi*1i);
    end
    

%     % round to nearest integer
%     Nroots = round(Nroots);
    
    % check for NaN or Inf values
    if isnan(Nroots) || isinf(Nroots)
        return
    end
    
    % IF WE CROSS THE BRANCH CUT ONCE WE GET AN IMAGINARY VALUE OF 1/PI
    % IF WE CROSS THE BRANCH CUT TWICE WE GET AN IMAGINARY VALUE OF 2/PI
    % ETC...

    % check for imaginary value because that means we crossed a
    % branch cut
    if round(imag(Nroots)/(1/pi)) ~= 0
        % check how often 1/pi fits in the imaginary part and round it to
        % get an integer which is the amount of times we cross a branch cut
        N_BC = round(imag(Nroots)/(1/pi));
        
        % if the rounded real part is also equal to zero then we know there
        % aren't any roots (or branch points) so we just stop the search
        if round(real(Nroots)) == 0
            Nroots = 0;
            return
        end
        

%         if mod(real(Nroots),0)-mod(real(Nroots),1) ~= 0 && round(mod(real(Nroots),1)) == 0
%             % when we have whole integer values return also N_BC so we know
%             % we must evade the branch cut (or use it later for something
%             % else)
%             Nroots = round(real(Nroots),2,'significant') + 1i*N_BC;
        
        % check if we have integer values or not and then round to the correct
        % values (i.e. either 0.5, 1.5, 2.5, etc.) and return that value 
        if mod(real(Nroots),0)-mod(real(Nroots),1) ~= 0
            % when we have values of 1.5, 2.5, etc
            Nroots = round(real(Nroots),2,'significant') + 1i*N_BC ;
        else
            % when we have a values of 0.5 (thus one branch point)
            Nroots = round(real(Nroots),1,'significant') + 1i*N_BC;
        end

        return
    end

    % check whether the rounded value goes to zero, if so, make it zero and
    % return that value
    if round(Nroots) == 0
        Nroots = 0;
        return
    end

    % if the above is not true, round to nearest integer
    Nroots = round(real(Nroots));
end
%%

function new_circles = subdivide_8circles(z0,r0)

    % function to divide initial circle into smaller circles
    con_circle = [z0, r0/2];
    r0_8 = 5*r0/12;
    k_8 = (0:1:8-1).';
    z0_8 = z0 + (3/4)*r0*exp(2*pi*1i.*k_8/8);
    new_circles = [con_circle; ...
                    z0_8, ones(8,1).*r0_8];
end

%%

function n_sing = argument_principle(f,wr_min, wr_max, wi_min, wi_max, w_step)
    % Initialize constants
    branch_limit = 1e-10;
    step_limit = 1e-25;
    corner_limit = 1e-16;

    % Initialize variables
    n_sing_sucks = true;
    step = w_step;
    n_sing = 0;
    n_sing_exact = 0;
    n_sing_previous = 0;
    calc_counter = 0;
    wr_all = [wr_min,wr_max,wi_min,wi_max];
    while n_sing_sucks
        calc_counter = calc_counter + 1;
        n_count = 0;
        step_orig = step;
        denom_cmplxA = zeros(4,1);
        % Corner points
        % 1st corner
        wcA = wr_min + 1i * wi_max;
        denom_cmplxA(1) = f(wcA); 
        denom_realC1 = real(denom_cmplxA(1));
        denom_imagC1 = imag(denom_cmplxA(1));

        % 2nd corner
        wcA = wr_max + 1i * wi_max;
        denom_cmplxA(2) = f(wcA);
        denom_realC2 = real(denom_cmplxA(2));
        denom_imagC2 = imag(denom_cmplxA(2));

        % 3rd corner
        wcA = wr_max + 1i * wi_min;
        denom_cmplxA(3) = f(wcA);
        denom_realC3 = real(denom_cmplxA(3));
        denom_imagC3 = imag(denom_cmplxA(3));

        % 4th corner
        wcA = wr_min + 1i * wi_min;
        denom_cmplxA(4) = f(wcA);
        denom_realC4 = real(denom_cmplxA(4));
        denom_imagC4 = imag(denom_cmplxA(4));

        % 1ST PART: THE HORIZONTAL PATH FROM THE 1ST TO THE 2ND CORNER
        denom_realA = real(denom_cmplxA(1));
        denom_imagA = imag(denom_cmplxA(1));
        
        wr = wr_min;
        while wr < wr_max
            if wr + step < wr_max
                wrB = wr + step;
                wcB = wrB + 1i * wi_max; % first point
                denom_cmplxB = f(wcB);
                denom_realB = real(denom_cmplxB);
                denom_imagB = imag(denom_cmplxB);
            else
                wrB = wr_max;
                step = wr_max - wr;
                denom_realB = denom_realC2;
                denom_imagB = denom_imagC2;
            end
        
            if denom_realA * denom_realB < 0 && denom_imagA * denom_imagB < 0
                step = 0.3 * step;
                if step < step_limit
                    disp(['kr stuck at ', num2str(wcB)]);
                    disp(['det ', num2str(denom_cmplxB)]);
                    disp(['real(kr)*imag(kr) ', num2str(real(wcB) * imag(wcB))]);
                    
%                     n_sing = 2001;
%                     return;
                end
            else
                if denom_realA > 0 && denom_imagA < 0 && denom_realB > 0 && denom_imagB > 0
                    n_count = n_count + 1;
                end
                if denom_realA > 0 && denom_imagA > 0 && denom_realB < 0 && denom_imagB > 0
                    n_count = n_count + 1;
                end
                if denom_realA < 0 && denom_imagA > 0 && denom_realB < 0 && denom_imagB < 0
                    n_count = n_count + 1;
                end
                if denom_realA < 0 && denom_imagA < 0 && denom_realB > 0 && denom_imagB < 0
                    n_count = n_count + 1;
                end
        
                if denom_realA > 0 && denom_imagA > 0 && denom_realB > 0 && denom_imagB < 0
                    n_count = n_count - 1;
                end
                if denom_realA < 0 && denom_imagA > 0 && denom_realB > 0 && denom_imagB > 0
                    n_count = n_count - 1;
                end
                if denom_realA < 0 && denom_imagA < 0 && denom_realB < 0 && denom_imagB > 0
                    n_count = n_count - 1;
                end
                if denom_realA > 0 && denom_imagA < 0 && denom_realB < 0 && denom_imagB < 0
                    n_count = n_count - 1;
                end
        
                denom_realA = denom_realB;
                denom_imagA = denom_imagB;
                step = step_orig;
                wr = wrB;
            end
        end
        % 2nd part: The vertical path from the 2nd to the 3rd corner
        denom_realA = real(denom_cmplxA(2));
        denom_imagA = imag(denom_cmplxA(2));
        wi = wi_max;
        
        while wi > wi_min
            if wi - step > wi_min
                wiB = wi - step;
                wcB = wr_max + 1i * wiB;
                denom_cmplxB = f(wcB);
                denom_realB = real(denom_cmplxB);
                denom_imagB = imag(denom_cmplxB);
            else
                wiB = wi_min;
                step = wi - wi_min;
                denom_realB = denom_realC3;
                denom_imagB = denom_imagC3;
            end
        
            if denom_realA * denom_realB < 0 && denom_imagA * denom_imagB < 0
                step = 0.3 * step;
                if step < step_limit
                    disp(['kr stuck at ', num2str(wcB)]);
                    disp(['det ', num2str(denom_cmplxB)]);
                    disp(['real(kr)*imag(kr) ', num2str(real(wcB) * imag(wcB))]);
        
%                     n_sing = 2002;
%                     return;
                end
            else
                if denom_realA > 0 && denom_imagA < 0 && denom_realB > 0 && denom_imagB > 0
                    n_count = n_count + 1;
                end
                if denom_realA > 0 && denom_imagA > 0 && denom_realB < 0 && denom_imagB > 0
                    n_count = n_count + 1;
                end
                if denom_realA < 0 && denom_imagA > 0 && denom_realB < 0 && denom_imagB < 0
                    n_count = n_count + 1;
                end
                if denom_realA < 0 && denom_imagA < 0 && denom_realB > 0 && denom_imagB < 0
                    n_count = n_count + 1;
                end
        
                if denom_realA > 0 && denom_imagA > 0 && denom_realB > 0 && denom_imagB < 0
                    n_count = n_count - 1;
                end
                if denom_realA < 0 && denom_imagA > 0 && denom_realB > 0 && denom_imagB > 0
                    n_count = n_count - 1;
                end
                if denom_realA < 0 && denom_imagA < 0 && denom_realB < 0 && denom_imagB > 0
                    n_count = n_count - 1;
                end
                if denom_realA > 0 && denom_imagA < 0 && denom_realB < 0 && denom_imagB < 0
                    n_count = n_count - 1;
                end
                % Other conditions for counting
        
                denom_realA = denom_realB;
                denom_imagA = denom_imagB;
                step = step_orig;
                wi = wiB;
            end
        end
        % 3rd part: The horizontal path from the 3rd to the 4th corner
        denom_realA = real(denom_cmplxA(3));
        denom_imagA = imag(denom_cmplxA(3));
        wr = wr_max;
        
        while wr > wr_min
            if wr - step > wr_min
                wrB = wr - step;
                wcB = wrB + 1i * wi_min;
        
                denom_cmplxB = f(wcB);
                denom_realB = real(denom_cmplxB);
                denom_imagB = imag(denom_cmplxB);
            else
                wrB = wr_min;
                step = wr - wr_min;
                denom_realB = denom_realC4;
                denom_imagB = denom_imagC4;
            end
        
            if denom_realA * denom_realB < 0 && denom_imagA * denom_imagB < 0
                step = 0.3 * step;
               if step < step_limit
                    disp(['kr stuck at ', num2str(wcB)]);
                    disp(['det ', num2str(denom_cmplxB)]);
                    disp(['real(kr)*imag(kr) ', num2str(real(wcB) * imag(wcB))]);
        
%                     n_sing = 2003;
%                     return;
                end
            else
                if denom_realA > 0 && denom_imagA < 0 && denom_realB > 0 && denom_imagB > 0
                    n_count = n_count + 1;
                end
                if denom_realA > 0 && denom_imagA > 0 && denom_realB < 0 && denom_imagB > 0
                    n_count = n_count + 1;
                end
                if denom_realA < 0 && denom_imagA > 0 && denom_realB < 0 && denom_imagB < 0
                    n_count = n_count + 1;
                end
                if denom_realA < 0 && denom_imagA < 0 && denom_realB > 0 && denom_imagB < 0
                    n_count = n_count + 1;
                end
        
                if denom_realA > 0 && denom_imagA > 0 && denom_realB > 0 && denom_imagB < 0
                    n_count = n_count - 1;
                end
                if denom_realA < 0 && denom_imagA > 0 && denom_realB > 0 && denom_imagB > 0
                    n_count = n_count - 1;
                end
                if denom_realA < 0 && denom_imagA < 0 && denom_realB < 0 && denom_imagB > 0
                    n_count = n_count - 1;
                end
                if denom_realA > 0 && denom_imagA < 0 && denom_realB < 0 && denom_imagB < 0
                    n_count = n_count - 1;
                end
        
                denom_realA = denom_realB;
                denom_imagA = denom_imagB;
                step = step_orig;
                wr = wrB;
            end
        end
        % 4th part: The vertical path from the 4th to the 1st corner
        denom_realA = real(denom_cmplxA(4));
        denom_imagA = imag(denom_cmplxA(4));
        wi = wi_min;
        
        while wi < wi_max
            if wi + step < wi_max
                wiB = wi + step;
                wcB = wr_min + 1i * wiB;
                denom_cmplxB = f(wcB);
        
                denom_realB = real(denom_cmplxB);
                denom_imagB = imag(denom_cmplxB);
            else
                wiB = wi_max;
                step = wi_max - wi;
                denom_realB = denom_realC1;
                denom_imagB = denom_imagC1;
            end
        
            if denom_realA * denom_realB < 0 && denom_imagA * denom_imagB < 0
                step = 0.3 * step;
                if step < step_limit
                    disp(['kr stuck at ', num2str(wcB)]);
                    disp(['det ', num2str(denom_cmplxB)]);
                    disp(['real(kr)*imag(kr) ', num2str(real(wcB) * imag(wcB))]);
        
%                     n_sing = 2004;
%                     return;
                end
            else
                if denom_realA > 0 && denom_imagA < 0 && denom_realB > 0 && denom_imagB > 0
                    n_count = n_count + 1;
                end
                if denom_realA > 0 && denom_imagA > 0 && denom_realB < 0 && denom_imagB > 0
                    n_count = n_count + 1;
                end
                if denom_realA < 0 && denom_imagA > 0 && denom_realB < 0 && denom_imagB < 0
                    n_count = n_count + 1;
                end
                if denom_realA < 0 && denom_imagA < 0 && denom_realB > 0 && denom_imagB < 0
                    n_count = n_count + 1;
                end
        
                if denom_realA > 0 && denom_imagA > 0 && denom_realB > 0 && denom_imagB < 0
                    n_count = n_count - 1;
                end
                if denom_realA < 0 && denom_imagA > 0 && denom_realB > 0 && denom_imagB > 0
                    n_count = n_count - 1;
                end
                if denom_realA < 0 && denom_imagA < 0 && denom_realB < 0 && denom_imagB > 0
                    n_count = n_count - 1;
                end
                if denom_realA > 0 && denom_imagA < 0 && denom_realB < 0 && denom_imagB < 0
                    n_count = n_count - 1;
                end
        
                denom_realA = denom_realB;
                denom_imagA = denom_imagB;
                step = step_orig;
                wi = wiB;
            end
        end

        % Calculate n_sing and n_sing_exact
        n_sing = abs(n_count / 4);
        n_sing_exact = abs(n_count / 4);

        % Check convergence criteria
        if calc_counter == 1
            step = step * 0.5;
            n_sing_previous = n_sing;
        elseif n_sing ~= n_sing_previous
            step = step * 0.5;
            n_sing_previous = n_sing;
        else
            n_sing_sucks = false;
        end
    end
end
