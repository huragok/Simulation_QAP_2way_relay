function h = plot_mapping(X, map, figureName, flagBinary)

Q = length(map); % Size of constellation
Nbps = round(log2(Q));

x_max = max(real(X));
x_min = min(real(X));
y_max = max(imag(X));
y_min = min(imag(X));

d_x = abs(real(X(:)) * ones(1, Q) - ones(Q, 1) * real(X(:))');
delta_x_min = min(d_x(d_x > 1e-5));
d_y = abs(imag(X(:)) * ones(1, Q) - ones(Q, 1) * imag(X(:))');
delta_y_min = min(d_y(d_y > 1e-5));

h = figure('name', figureName);
for q = 1 : Q
    plot(real(X(map(q))), imag(X(map(q))), 'bo', 'LineWidth', 2), hold on;
    if (flagBinary)
        text(real(X(map(q))) + 0.2 * delta_x_min, imag(X(map(q))) + 0.2 * delta_y_min, dec2bin(q - 1, Nbps), 'fontsize', 12);
    else
        text(real(X(map(q))) + 0.2 * delta_x_min, imag(X(map(q))) + 0.2 * delta_y_min, num2str(q - 1), 'fontsize', 12);
    end
end

axis equal;
xlim([x_min - 0.4 * delta_x_min, x_max + 0.4 * delta_x_min]), ylim([y_min - 0.4 * delta_y_min, y_max + 0.4 * delta_y_min]);
grid on
set(gca, 'fontsize', 16);
xlabel('real'), ylabel('image')