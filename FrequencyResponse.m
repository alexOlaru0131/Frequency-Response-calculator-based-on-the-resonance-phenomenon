%%
%

clc
clear
close all

data = readmatrix('scope_71.csv','NumHeaderLines',2);

t = data(:,1);
u = data(:,2);
y1 = data(:,3);
y2 = data(:,4);

[max_y, i_h_y] = findpeaks(y1);
[min_y, i_l_y] = findpeaks(-y1);
[max_u, i_h_u] = findpeaks(u);
[min_u, i_l_u] = findpeaks(-u);

[max_y_total, y_max] = max(max_y);
[max_u_total, u_max] = max(max_u);
[min_u_total, u_min] = min(min_u);

eMPN_v = 1;
tol_p = 5;
tol_lr = 30;

for i = -tol_p:1:tol_p
    for j = -tol_p:1:tol_p
        for k = -tol_lr:1:tol_lr
            eMPN = 1;
            Mr = (abs(y1(i_h_y(y_max+i)+k)) - abs(y1(i_l_y(y_max+j+1)+k)))/(abs(max_u(u_max)) - abs(min_u(u_min)));
            if Mr > 1
                tita = roots([-4*Mr^2 0 4*Mr^2 0 -1]);
                if(length(tita) ~= 0)
                    tita = tita(4);
                else 
                    tita = 1;
                end
                
                    for l = -tol_lr:1:tol_lr
                        Tr = t(i_h_y(y_max+1)+i+l) - t(i_h_y(y_max)+j+k);
                        wr = (2*pi)/Tr;
                        wn = wr / sqrt(1 - 2*tita^2);
            
                        if(wr ~= Inf)
                            K = mean(y1)/mean(u);
                
                            A = [ 0 1 ; -wn^2 -2*tita*wn];
                            B = [ 0 ; K * wn^2];
                            C = [1 0];
                            D = 0;
                
                            y_m = lsim(A,B,C,D , u, t, [y1(1) (y1(1)-y1(2))/(t(1)-t(2))]);
                
                            J = norm(y1 - y_m) / sqrt(length(y1));
                            eMPN = norm(y1 - y_m) / norm(y1 - mean(y1));
                            if(eMPN < eMPN_v) 
                                eMPN_v = eMPN;
                                K_v = K;
                                wn_v = wn;
                                tita_v = tita;
                                y_m_v = y_m;
                
                                i1 = i_h_y(y_max+i+j)+k;
                                i2 = i_l_y(y_max+i+j)+k;
                            end
                        end
                    end
            end
        end
    end
end

H_y1 = tf(K_v*wn_v^2, [1 2*tita_v*wn_v wn_v^2]);

logspaceCounter = 0;
min_i = min([length(i_h_y), length(i_l_y), length(i_h_u), length(i_l_u)]);
while true
    if wn >= 10^logspaceCounter && wn <= 10^(logspaceCounter+1)
        break
    else
        logspaceCounter = logspaceCounter + 1;
    end
end

w = logspace(logspaceCounter, logspaceCounter+1, min_i);
[M, Ph] = bode(K_v*wn_v^2, [1 2*tita_v*wn_v wn_v^2], w);

for i = 2:1:min_i-1
    Mc(i) = 20*log10((max_y(i) + min_y(i))/(max_u(i) + min_u(i)));
end

for i = 1:1:min_i-1
        dt = t(i_l_y(i))-t(i_h_y(i+1));
        t01 = t(i_l_y(i))...
             + (mean(y1)-y1(i_l_y(i)+1))/(y1(i_h_y(i)-1)-y1(i_l_y(i)+1))*dt;
        t02 = t(i_h_y(i+1))...
             + (mean(y1)-y1(i_l_y(i+1)+1))/(y1(i_h_y(i+1)-1)-y1(i_l_y(i+1)+1))*dt;
        wc(i) = pi/(t02-t01);
end

for i = 1:1:min_i-length(i_h_u)+length(i_h_y)
    if length(i_h_y) < length(i_h_u)
        phic(i) = -rad2deg(wc(i)*(t(i_h_y(i))-t(i_h_u(i+length(i_h_u)-length(i_h_y)))));
    else
        phic(i) = -rad2deg(wc(i)*(t(i_h_y(i))-t(i_h_u(i+length(i_h_y)-length(i_h_u)))));
    end
end

n_afis = 0;
subplot(211)
semilogx(w, 20*log10(M)), hold on
for i=1:1:min_i-1
    for k=1:1:min_i-1
                if (sqrt((wc(i) - w(k))^2 + (abs(Mc(i) - 20*log10(M(k))))^2) < 10e+2 && abs(Mc(i) - 20*log10(M(k))) < 1)
                    semilogx(wc(i), Mc(i),'x')
                    n_afis = n_afis + 1;
                    Mc_afis(n_afis) = Mc(i);
                    wc_afis(n_afis) = wc(i);
                end
    end
end

grid, shg, xlim([10^logspaceCounter 10^(logspaceCounter+1)]), ylim(20*log10([min(M) max(M)]))

nr_afis = 1;
wc_f_afis = [];
phic_afis = [];
subplot(212), semilogx(w, Ph), hold on
for i = 1:1:min_i-1
    for j = 1:1:min_i-1
        if (sqrt((wc(i) - w(j))^2 + (phic(i) - Ph(j))^2) < 10e+3 && abs(phic(i) - Ph(j)) < 2)
            wc_f_afis(nr_afis) = wc(i);
            phic_afis(nr_afis) = phic(i);
            semilogx(wc_f_afis(nr_afis),phic_afis(nr_afis),'x')
            nr_afis = nr_afis + 1;
        end
    end
end

grid, shg, xlim([10^logspaceCounter 10^(logspaceCounter+1)]), ylim([min(Ph) max(Ph)])

%