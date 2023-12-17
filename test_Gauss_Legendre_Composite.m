
function test_Gauss_Legendre_Composite(f, a, b, c, d, nx, ny)
    
    % Expected result (you should replace this with the actual expected result)
    expected_result = integral2(f, a, b, c, d);

    % Actual result from Gauss_Legendre_Composite
    approximated_result = Gauss_Legendre_Composite(f, a, b, c, d, nx, ny);

    % Calculate error
    error = abs(approximated_result - expected_result);

    % Display the results
    disp(['Expected Result: ', num2str(expected_result)]);
    disp(['Actual Result: ', num2str(approximated_result)]);

    % Display errors
    disp(['Error: ', num2str(error)]);
    

    
end
