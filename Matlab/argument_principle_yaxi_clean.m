function n_sing = Argument_Principle(f, wr_min, wr_max, wi_min, wi_max, w_step)
% Argument_Principle: Find the number of roots of a complex analytic function within a rectangle
% in the complex plane using the Argument Principle method.

% Inputs:
% - f(z): determinant depending on a complex variable 'z'
% - wr_min, wr_max: Real part boundaries
% - wi_min, wi_max: Imaginary part boundaries
% - w_step: Initial step size

% Output:
% - n_sing: Number of singularities (zeros) within the rectangle

% Note: The function f(z) must be defined externally.

% Constants
branch_limit = 1e-10; % Not used in this code but included for completeness
step_limit = 1e-12;   % Minimum allowable step size
corner_limit = 1e-16; % Tolerance for proximity to singularities

% Initialize variables
calc_counter = 0;
n_sing_sucks = true;
step = w_step;

% Initialize variables for the corner points
% 1st corner
wcA = wr_min + 1i * wi_max;
denom_cmplxA = f(wcA);
denom_realC1 = real(denom_cmplxA);
denom_imagC1 = imag(denom_cmplxA);

% 2nd corner
wcA = wr_max + 1i * wi_max;
denom_cmplxA = f(wcA);
denom_realC2 = real(denom_cmplxA);
denom_imagC2 = imag(denom_cmplxA);

% 3rd corner
wcA = wr_max + 1i * wi_min;
denom_cmplxA = f(wcA);
denom_realC3 = real(denom_cmplxA);
denom_imagC3 = imag(denom_cmplxA);

% 4th corner
wcA = wr_min + 1i * wi_min;
denom_cmplxA = f(wcA);
denom_realC4 = real(denom_cmplxA);
denom_imagC4 = imag(denom_cmplxA);

% Start Argument Principle calculation
while n_sing_sucks
    calc_counter = calc_counter + 1;
    n_count = 0;
    step_orig = step;
    
    % 1ST PART: HORIZONTAL PATH FROM 1ST TO 2ND CORNER
    denom_realA = denom_realC1;
    denom_imagA = denom_imagC1;
    wr = wr_min;
    while wr < wr_max
        if wr + step < wr_max
            wrB = wr + step;
            wcB = wrB + 1i * wi_max;
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
%                 n_sing = 2001;
                return
            end
        else
            n_count = update_n_count(n_count, denom_realA, denom_imagA, denom_realB, denom_imagB);
            denom_realA = denom_realB;
            denom_imagA = denom_imagB;
            step = step_orig;
            wr = wrB;
        end
    end

    % 2ND PART: VERTICAL PATH FROM 2ND TO 3RD CORNER
    denom_realA = denom_realC2;
    denom_imagA = denom_imagC2;
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
%                 n_sing = 2002;
                return
            end
        else
            n_count = update_n_count(n_count, denom_realA, denom_imagA, denom_realB, denom_imagB);
            denom_realA = denom_realB;
            denom_imagA = denom_imagB;
            step = step_orig;
            wi = wiB;
        end
    end

    % 3RD PART: HORIZONTAL PATH FROM 3RD TO 4TH CORNER
    denom_realA = denom_realC3;
    denom_imagA = denom_imagC3;
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
%                 n_sing = 2003;
                return
            end
        else
            n_count = update_n_count(n_count, denom_realA, denom_imagA, denom_realB, denom_imagB);
            denom_realA = denom_realB;
            denom_imagA = denom_imagB;
            step = step_orig;
            wr = wrB;
        end
    end

    % 4TH PART: VERTICAL PATH FROM 4TH TO 1ST CORNER
    denom_realA = denom_realC4;
    denom_imagA = denom_imagC4;
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
%                 n_sing = 2004;
                return
            end
        else
            n_count = update_n_count(n_count, denom_realA, denom_imagA, denom_realB, denom_imagB);
            denom_realA = denom_realB;
            denom_imagA = denom_imagB;
            step = step_orig;
            wi = wiB;
        end
    end

    % Final calculation
    n_sing = abs(n_count / 4);
    n_sing_exact = abs(n_count / 4); % Unused variable but kept for consistency

    if calc_counter == 1
        step = step * 0.5;
        n_sing_previous = n_sing;
    elseif n_sing ~= n_sing_previous
        step = step * 0.5;
        n_sing_previous = n_sing;
    else
        n_sing_sucks = false;
    end
end % End of while loop

end % End of function

% Helper function to update n_count based on conditions
function n_count = update_n_count(n_count, denom_realA, denom_imagA, denom_realB, denom_imagB)
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
end
