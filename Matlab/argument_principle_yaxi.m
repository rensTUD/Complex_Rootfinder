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