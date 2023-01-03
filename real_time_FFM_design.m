%% test real-time filter design
% PRODUCES FIGURE 2 IN PAPER

close all; clear all;
phi =  (0:0.25:1) * (pi/4);
N = 2;
order = 2;
filter_coeff = zeros(N,N,order+1);
nbins = 1024;
w = linspace(0,pi,nbins);
fs = 44100;
alpha = ones(N, N, order);
beta =(0:0.1:1) * pi/2;


%% get the colors right

color_gradient = zeros(length(beta), 3, length(phi));
blue = [0 0 1];
light_blue = [91, 207, 244] / 255;
color_gradient(:,:,5) = create_color_gradient(blue, light_blue, length(beta));

red = [204 0 0]/255;
pink = [255, 192, 203]/255;
color_gradient(:,:,2) = create_color_gradient(red, pink, length(beta));

green = [51, 102, 0]/255;
light_green = [153, 255, 153]/255; 
color_gradient(:,:,3) = create_color_gradient(green, light_green, length(beta));

purple = [102, 0, 204]/255;
light_purple = [204,153,255]/255;
color_gradient(:,:,4) = create_color_gradient(purple, light_purple, length(beta));

orange = [255,153,51]/255;
yellow = [255, 204, 153]/255;
color_gradient(:,:,1) = create_color_gradient(orange, yellow, length(beta));

lgdstr = {};
l = {};

fig = figure('Units','inches', 'Position',[0 0 6.58 4.2],'PaperPositionMode','auto');


%% run loop and plot locus

for i = 1:length(phi)
    for m = 1:length(beta)
    
        % basic 2D fundamental orthonormal matrix
        rotation_2D = [cos(phi(i)),sin(phi(i)); -sin(phi(i)), cos(phi(i))];
        %build complete set of orthonormal idempotents
        P = build_complete_set_orthonormal_idempotents(rotation_2D);


        %sign matrix, alpha, whose modulus is 1
        for j = 1:N
             for k = 1:N
                 alpha(j,k,order) = exp(1j*beta(m));
            end
        end

        % P_1  + alpha P_2 z - i hav checked, this does the right thing
        filter_coeff(:,:,(order+1)-N+1:end) = P.*alpha;


        %plot zeros of the top diagonal

        count = 1;
        for j = 1:N
           for k = 1:N
                b = reshape(filter_coeff(j,k,:), [order+1,1]);
                fig; subplot(N,N,count); 
                [hz1, hp1, ht1] = zplane(b, zeros(order+1));hold on;
                if (m ==1)
                    l{i} = hz1;
                end
                set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
                set(findobj(hz1, 'Type', 'line'), 'Color', color_gradient(m,:,i));
                set(findobj(hz1, 'Type', 'line'), 'MarkerFaceColor', color_gradient(m,:,i));
                set(findobj(hz1, 'Type', 'line'), 'MarkerSize', 3);
                xlim([-1.1, 1.1]); ylim([-1.1,1.1]);
                axis square; 
                drawnow;

%               % uncomment to plot magnitude response instead               
%                 if (i == length(phi))
%                     
%                     [H] = freqz(b, 1, w);
%                     H_exp = closed_form_magnitude_response(phi(i), beta(m), w, j, k);
%                     
%                     figure(1); subplot(N,N,count);
%                     semilogx(w * fs/2, mag2db(abs(H)), 'Color', color_gradient(m,:,i)); grid on; hold on;          
%                     xlim([200, fs/2]); ylim([-60, 10]);
%                     xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
%                     axis square; 
%                  
%                     figure(2); subplot(N,N,count);
%                     semilogx(w * fs/2, mag2db(abs(H_exp)), 'Color', color_gradient(m,:,i)); grid on; hold on;          
%                     xlim([200, fs/2]); ylim([-60, 10]);
%                     xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
%                     axis square;
%                    
%                     
% %                     set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
%                     drawnow;
%                 end

                count = count + 1;

           end
           lgdstr{i} = strcat('$\bar{\phi} = $', num2str(round(phi(i)/(pi/4),3)));
        end
    end
    
  
end

%% plot legend and save

% figure(1); hold off; sgtitle('Calculated magnitude response from zeros');
% figure(2); hold off; sgtitle('Closed form magnitude response');

fig; hold off; 
Lgnd = legend([l{1}, l{2}, l{3}, l{4}, l{5}],lgdstr);
Lgnd.Position(1) = 0.43;
Lgnd.Position(2) = 0.45;
Lgnd.Interpreter = 'latex';

print('figures/real_time_paraunitary_zeros.eps', '-depsc');


%% 

%% returns closed form magnitude response for particular values of phi and beta

function [H_exp] = closed_form_magnitude_response(phi, beta, w, i, j)
    
    if (i == j)
        H_exp = sqrt(cos(phi)^4 + sin(phi)^4 + (2*(cos(phi)^2)*(sin(phi)^2)).*cos(w - beta));
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
    color_gradient = [linspace(col1(1),col2(1),nterms)', linspace(col1(2),col2(2),nterms)'...
        , linspace(col1(3),col2(3),nterms)'];
end

