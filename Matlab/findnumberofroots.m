function Nroots = findnumberofroots(f,Z,varargin)
    
%% check for symbolic derivative of the function
    if any(strcmp(varargin,'df'))
        f1 = varargin{find(strcmp(varargin,'df'))+1};    
        option_der = 'symbolic';     % symbolic derivative has been given
    else
        option_der = 'numerical';    % second choice for derivative of function
    end
    
%% check for which way to calculate derivative
    if any(strcmp(varargin,'option_der'))
        option_der = varargin{find(strcmp(varargin,'option_der'))+1};    % if option is given then use that
    elseif strcmp(option_der,'symbolic')
        % do nothing as symbolic has been chosen either by supplying df or
        % choosing it above
    else    
        option_der = 'numerical';    % second choice for derivative of function
    end
    
%% calculate number of roots
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
        
        Nroots = 1/2/pi*sum(Args);

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
%%

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