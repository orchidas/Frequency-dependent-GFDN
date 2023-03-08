%% test real-time filter design
% PRODUCES FIGURE 2 IN PAPER
close all; clear all; clc;

phi =  (0:0.25:1) * (pi/4);
N = 2;
order = 1; % Sebastian: order 1 seems enough
filter_coeff = zeros(N,N,order+1);
nbins = 1024;
w = linspace(0,pi,nbins);
fs = 44100;
alpha = ones(N, N, order);
beta =(0:1:1) * pi/2;

phi_selection = [2 5];

%% get the colors right
color_gradient = zeros(length(beta), 3, length(phi));

orange = [255,153,51]/255;
yellow = [255, 204, 153]/255;
color_gradient(:,:,1) = create_color_gradient(orange, yellow, length(beta));

red = [204 0 0]/255;
pink = [255, 192, 203]/255;
color_gradient(:,:,2) = create_color_gradient(red, pink, length(beta));

green = [51, 102, 0]/255;
light_green = [153, 255, 153]/255;
color_gradient(:,:,3) = create_color_gradient(green, light_green, length(beta));

purple = [102, 0, 204]/255;
light_purple = [204,153,255]/255;
color_gradient(:,:,4) = create_color_gradient(purple, light_purple, length(beta));

blue = [0 0 1];
light_blue = [91, 207, 244] / 255;
color_gradient(:,:,5) = create_color_gradient(blue, light_blue, length(beta));

lgdstr = {};
l = {};

%% compute filter coefficients
zeroes = zeros(N,N,order,length(phi),length(beta));
for i = phi_selection %1:length(phi)
    for m = 1:length(beta)

        % basic 2D fundamental orthonormal matrix
        rotation_2D = [cos(phi(i)),sin(phi(i)); -sin(phi(i)), cos(phi(i))];
        %build complete set of orthonormal idempotents
        P = build_complete_set_orthonormal_idempotents(rotation_2D);

        %sign matrix, alpha, whose modulus is 1
        alpha(:,:,order+1) = exp(1j*beta(m));

        % P_1  + alpha P_2 z - i have checked, this does the right thing
        filter_coeff(:,:,:,i,m) = P.* alpha;

        % always true
        % isParaunitary(filter_coeff(:,:,:,i,m))

    end
end

%% setup figures
fig = figure('Units','inches', 'Position',[0 0 3.25 3.3],'PaperPositionMode','auto');
for j = 1:N
    for k = 1:N
        b = sub2ind([N,N],j,k);
        ax(j,k) = subplot(N,N,b); hold on;
        [hz(j,k), hp(j,k), ht(j,k)] = zplane([], []);
        title('')
    end
end

for i = phi_selection %1:length(phi)
    lgdstr{i} = sprintf('$\\phi = %2.2f$', round(phi(i)/(pi/4),3));
end

%% plot filter coefficients
for j = 1:N
    for k = 1:N
        for i = 1:length(phi)
            for m = 1:length(beta)

                b = squeeze(filter_coeff(j,k,:,i,m));

                plot_handle(j,k,i,m) = plot(ax(j,k),real(b),imag(b),'.',...
                    'Color', color_gradient(m,:,i),...
                    'MarkerSize', 5);
            end
        end
    end
end
%% legend
Lgnd = legend(squeeze(plot_handle(1,1,phi_selection,1)),lgdstr{phi_selection},'Interpreter','latex');
Lgnd.NumColumns = length(phi);
Lgnd.Position(1) = 0.4;
Lgnd.Position(2) = 0.95;

saveas(gcf,'../figures/filter_coefficients_2x2.png')


%% setup figures
fig = figure('Units','inches', 'Position',[0 0 3.25 3.3]*1.1,'PaperPositionMode','auto');

%% plot magnitude response
f = w / pi * fs/2;

for j = 1:N
    for k = 1:N
        for i = phi_selection %1:length(phi)
            for m = 1:length(beta)
                b = squeeze(filter_coeff(j,k,:,i,m));
                H(j,k,:,i,m) = freqz(b, 1, w);
            end
        end
    end
end

for i = phi_selection %1:length(phi)
    for m = 1:length(beta)
        [plotAxes,plot_handle(:,:,i,m)] = plotImpulseResponseMatrix(f,mag2db(abs(H(:,:,:,i,m))),...
            'Color', color_gradient(m,:,i),...
            'MarkerSize', 5, ...
            'xlabel','Frequency (Hz)','ylabel','Magnitude (dB)',...
            'xlim',[200, fs/2],'ylim',[-60, 10]);
    end
end
set(plotAxes,'XScale','log');
set(plotAxes,'xtick',[100 1000 10000]);

%% legend
Lgnd = legend(squeeze(plot_handle(1,1,phi_selection,1)),lgdstr{phi_selection},'Interpreter','latex');
Lgnd.NumColumns = length(phi);
Lgnd.Position(1) = 0.25;
Lgnd.Position(2) = 0.95;

exportgraphics(gcf,'../figures/filter_magnitude_respose_2x2.pdf','BackgroundColor','none','ContentType','vector')

%% returns closed form magnitude response for particular values of phi and beta
function [H_exp] = closed_form_magnitude_response(phi, beta, w, i, j)

if (i == j)
    H_exp = sqrt(cos(phi)^4 + sin(phi)^4 + (2*(cos(phi)^2)*(sin(phi)^2))...
        .*cos(w - beta));
else
    H_exp = sqrt(sin(2*phi).*(1-cos(w-beta)));
end

end

function [P] = build_complete_set_orthonormal_idempotents(M)

% P - NxNxN set of idempotent matrices
% M - NxN unitary matrix (square) from whose elements we will form the idempotents

N = size(M,1);
P = zeros(N,N,N);

for  k = 1:N
    col_M = M(:,k);
    P(:,:,k) = col_M*col_M';
end

if ~check_orthonormal_idempotents(P)
    error('Set of matrices do not constitute a complete set of orthonormal idempotents');
end
end


function [flag] = check_orthonormal_idempotents(P)

% given a set of matrices P, check if they form a complete set of symmetric
% orthonormal idempotents
N = size(P,1);
tol = 1e-10;
flag = true;

for i = 1:N
    cur_mat = P(:,:,i);
    for j = 1:N
        other_mat = P(:,:,j);
        if i == j
            expected = cur_mat;
        else
            expected = zeros(N,N);
        end

        if any(abs(sum(cur_mat*other_mat - expected)) > tol)
            flag = false;
            return;
        end
    end
end

end


function [color_gradient] = create_color_gradient(col1, col2, nterms)
color_gradient = [linspace(col1(1),col2(1),nterms)',...
    linspace(col1(2),col2(2),nterms)',linspace(col1(3),col2(3),nterms)'];
end

