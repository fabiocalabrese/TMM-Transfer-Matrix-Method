% Dati originali (lambda in µm, n)
lambda_data = [0.2952, 0.3024, 0.3100, 0.3179, 0.3263, 0.3351, 0.3444, 0.3542, ...
               0.3647, 0.3757, 0.3875, 0.3999, 0.4133, 0.4275, 0.4428, 0.4592, ...
               0.4769, 0.4959, 0.5166, 0.5391, 0.5636, 0.5904, 0.6199, 0.6525, ...
               0.6888, 0.7293, 0.7749, 0.8266];  % µm

n_data = [3.810, 3.692, 3.601, 3.538, 3.501, 3.485, 3.495, 3.531, ...
          3.596, 3.709, 3.938, 4.373, 4.509, 5.052, 4.959, 4.694, ...
          4.492, 4.333, 4.205, 4.100, 4.013, 3.940, 3.878, 3.826, ...
          3.785, 3.742, 3.700, 3.666];

% Convertiamo lambda da µm a nm
lambda_data_nm = lambda_data * 1000;

% Nuovo vettore di interpolazione
lambda_interp = 300:5:800;

% Interpolazione (pchip = spline a pendenza continua)
n_interp = interp1(lambda_data_nm, n_data, lambda_interp, 'pchip', 'extrap');

% Plot
plot(lambda_interp, n_interp, 'b', 'LineWidth', 2);
xlabel('\lambda (nm)');
ylabel('n');
title('Interpolazione dell''indice di rifrazione di GaAs');
grid on;
