function Roots = ComplexRootFinder_circular(f,SearchDomain)
% Root finder that finds complex roots with circular subdomains

%% input
% f --> anonymous function with one variable
% SearchDomain = [re_min, re_max, im_min, im_max]


%% main part

%% Circle subdomain - finding and polishing roots during the while loop


% plot search domain as rectangular domain inscribing a circle
Z = Vertices_Rectangular(re_min,re_max,im_min,im_max,1000);
% figure;
% scatter(real(Z),imag(Z));
% lengths of domain
re_length = re_max-re_min;
im_length = im_max-im_min;
% initial circle
z0 = mean(Z); % middle of the rectangle as z0
r0 = 0.5*sqrt(re_length^2+im_length^2); % radius of circle equal to diagonal of rectangular / square search domain
k = 0:1:1000-1;
circ1 = z0 + r0 * exp(2*pi*1i.*k/1000);

plot = 1
if plot == 1
    figure
    hold on
    h2 = scatter(real(circ1),imag(circ1));
    h3 = scatter(real(circ1),imag(circ1));
end

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
    plot = 0;
    
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






%% extra functions needed

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
            
%         if strcmp(type_optimise,'Matlab')
            R = fsolve(f,w);
    
%         elseif strcmp(type_optimise,'Alternate_Newton_Rhapson')
%             for iN = 1:Nroots
%                 R(iN) = root_finder_func(f,[min_re,max_re],[min_im,max_im],w(iN));
%             end
%         end
    
    % if it is equal to one the roots is just the root of that polynomial
    % which apparently exact (:
%     else
%         R = w;
%     end

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
%         n = Npoints;
        fk = f(Z);
%         c = fft(fk)/n;
%         cp = (1:n-1).*c(2:end);
%         ppzk = n*ifft([cp 0])/r0;
        
        % by Austin / Kravanja
%         Nroots = (real(mean(Z.*ppzk./fk)));
        % Delves and lynes but new paper approach of derivative
%         Nroots = trapz(Z,ppzk./fk) / (2*pi*1i);
        % fully numerical..
        Nroots = trapz(Z(1:end-1),(diff(fk)./diff(Z))./fk(1:end-1)) / (2*pi*1i);
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





end