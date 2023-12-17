%Hassan Jaber - 330214 - Project number 5 - Group B
function I = Gauss_Legendre_Composite(f, a, b, c, d, nx, ny)
    % nx: number of subintervals in x
    % ny: number of subintervals in y

    % Gauss-Legendre quadrature points and weights
    p = [-0.5773502691896257, 0.5773502691896257];
    w = [1.0000000000000000, 1.0000000000000000];

    % Initialize result
    I = 0;

    % Step sizes for x and y
    hx = (b - a) / nx;
    hy = (d - c) / ny;

    % Loop over subintervals in x
    for i = 1:nx
        x_start = a + (i - 1) * hx;
        x_end = a + i * hx;

        % Loop over subintervals in y for each x interval
        for j = 1:ny
            y_start = c + (j - 1) * hy;
            y_end = c + j * hy;

            % Coordinate linear transformation for x
            xp = 0.5 * ((x_end - x_start) * p + (x_end + x_start));

            % Coordinate linear transformation for y
            yp = 0.5 * ((y_end - y_start) * p + (y_end + y_start));

            % Initialize result for current subinterval
            I1 = 0;

            % Gauss-Legendre quadrature for x variable
            for k = 1:length(w)
                % Substitution of x with xp(k)
                ex = f(xp(k), yp);

                % Calculating the integral for x
                I1 = I1 + ex * w(k) * (x_end - x_start) / 2;
            end

            % Gauss-Legendre quadrature for y variable
            for l = 1:length(w)
                % Substitution of y with yp(l)
                ey = I1 * w(l) * (y_end - y_start) / 2;

                % Summing along the x-axis
                I = I + sum(ey);
            end
        end
    end

    % Dividing the final answer by two, In the double integration process, 
    % I have nested loops over subintervals in both the x and y directions. 
    % For each subinterval in x, I iterate over subintervals in y, and 
    % within each y subinterval, I perform Gauss-Legendre quadrature over points in x.
    % To avoid double-counting the contributions from each subinterval, the final result is divided by 2.
    I = I /2;

end
