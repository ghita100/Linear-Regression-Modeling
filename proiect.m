clear all
close all
clc
load('proj_fit_03.mat')
x1 = id.X{1};
x2 = id.X{2};
yid = id.Y;
yvalidare = val.Y;
xvalidare1 = val.X{1};
xvalidare2 = val.X{2};
dim1 = id.dims(1);
dim2 = val.dims(1); 
m=30;
mse_val =zeros(1,m);
mse_identif =zeros(1,m);
for g=1:m
    cont = 1;
    nrtermeni =(g+2)*(g+1)/2;
    fi = zeros(length(yid), nrtermeni); 
for i = 1:length(x1)
    for j = 1:length(x2)
        a = 0;b = 0;
        for k = 1:nrtermeni
            fi(cont, k) = (x1(i) .^ (a - b)) .* (x2(j) .^ b);
            if (a == b && a < g)
                a = b + 1;
                b = 0;
            else
                b = b + 1;
            end
        end
        cont = cont + 1;
    end
end
k = 1;
%convertim matricea de indentificare intr un vector
y_idconvertitdinmatrice = zeros(dim1^2, 1);
for i = 1:length(x1)
    for j = 1:length(x2)
        y_idconvertitdinmatrice(k, 1) = yid(i, j);
        k = k + 1;
    end
end
teta = fi \ y_idconvertitdinmatrice;
Yid_aprox= fi*teta;
k=1;
%rezultatul modelului pentru datele de identificare
y_idaproxmat =zeros(dim1,dim1);
for i = 1:length(x1)
    for j = 1:length(x2)
        y_idaproxmat(i,j) = Yid_aprox(k);
        k = k + 1;
    end
end
eid = zeros(length(Yid_aprox), 1);
for k = 1:length(Yid_aprox)
    eid(k) = (Yid_aprox(k) - y_idconvertitdinmatrice(k))^2;
end
s = sum(eid);
MSE_id = s / length(y_idconvertitdinmatrice);
    mse_identif(g) =MSE_id;
fival = zeros(length(yvalidare), nrtermeni);
cont = 1;
for i = 1:length(xvalidare1)
    for j = 1:length(xvalidare2)
        a = 0;b = 0;
        for k = 1:nrtermeni
            fival(cont, k) = (xvalidare1(i) .^ (a - b)) .* (xvalidare2(j) .^ b);
            if (a == b && a < g)
                a = a + 1;
                b = 0;
            else
                b = b + 1;
            end
        end
        cont = cont + 1;
    end
end
Yarpox = fival * teta;
k = 1; 
y_aproxmat = zeros(dim2, dim2); 
for i = 1:length(xvalidare1)
    for j = 1:length(xvalidare2)
        y_aproxmat(i, j) = Yarpox(k);
        k = k + 1;
    end
end
k = 1; 
y_valconvertitdinmatrice = zeros(dim2^2, 1);
for i = 1:length(xvalidare1)
    for j = 1:length(xvalidare2)
        y_valconvertitdinmatrice(k, 1) = yvalidare(i, j);
        k = k + 1;
    end
end
e = zeros(length(y_valconvertitdinmatrice), 1);
for k = 1:length(y_valconvertitdinmatrice)
    e(k) = (Yarpox(k) - y_valconvertitdinmatrice(k))^2;
end
s = sum(e);
MSE_val = s / length(y_valconvertitdinmatrice);
    mse_val(g) =MSE_val;
if (g == 11)
figure
subplot(211)
mesh(xvalidare1, xvalidare2, y_aproxmat)
title('Aproximare date de validare la MSE cel mai mic')
xlabel('X1')
ylabel('X2')
zlabel('Y')
subplot(212)
mesh(xvalidare1, xvalidare2, yvalidare)
title('Date de validare')
end
%am plotat m cel mai optim
end
[min_mseid, min_gradid] = min(mse_identif);
[min_mse, min_grad] = min(mse_val);
figure;
subplot(211)
plot(1:m, mse_val, '-o')
xlabel('Gradul polinomului -> g')
ylabel('MSE_validare ')
title('MSE_validare vs Gradul Polinomului')
grid on
subplot(212)
plot(1:m, mse_identif, '-o')
xlabel('Gradul polinomului -> g')
ylabel('MSE_identificare ')
title('MSE_identificare vs Gradul Polinomului')
grid on
fprintf('MSE_val minim: %.2f, obținut la gradul: %d\n', min_mse, min_grad);
fprintf('MSE_id minim: %.2f, obținut la gradul: %d\n', min_mseid, min_gradid);