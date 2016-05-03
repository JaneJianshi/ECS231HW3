%------ bcsstk15 ------%
bcsstk=spconvert(load('bcsstk15.mtx'));
bc=bcsstk+bcsstk';

x_true = rand(size(bc, 2), 1); 
b = bc * x_true; 
x_0 = zeros(size(bc,2), 1); 

[bcx_msd, bc_itersd, bc_msd] = mysd(bc, b, x_0, 0.000001, 100000);
[bcx_mcg, bc_itercg, bc_mcg] = mycg(bc, b, x_0, 0.000001, 100000);
 
%matlab function

[bcx,flag,tempres,bc_iter, bc_res] = cgs(bc,b,0.000001, 100000);
bc_res = bc_res/norm(b); 



plot(0:bc_itersd, bc_msd, 'ro-');
hold on ;
plot(0:bc_itercg, bc_mcg, 'cx-');
hold on;
plot(0:bc_iter, bc_res, 'bs-');
hold off;

lgd = legend ('mysd', 'mycg', 'Matlab');
set(lgd,'FontSize', 16);  
ylabel('Relative Residual') ;
xlabel('Iteration');
title('bcsstk15 matrix','FontSize', 20);
xlim([0 100]);
ylim([0 0.5]);



%------ nos3 ------%
nos=spconvert(load('nos3.mtx'));
nos = nos + nos';

x_true = rand(size(nos, 2), 1); 
b = nos * x_true; 
x_0 = zeros(size(nos,2), 1); 

[nosx_msd, nos_itersd, nos_msd] = mysd(nos, b, x_0, 0.000001, 100000);
[nosx_mcg, nos_itercg, nos_mcg] = mycg(nos, b, x_0, 0.000001, 100000);
 
%matlab function

[nosx,flag,tempres,nos_iter, nos_res] = cgs(nos,b, 0.000001, 100000);
nos_res = nos_res/norm(b); 



plot(0:nos_itersd, nos_msd, 'ro-');
hold on ;
plot(0:nos_itercg, nos_mcg, 'cx-');
hold on;
plot(0:nos_iter, nos_res, 'bs-');
hold off;

lgd = legend ('mysd', 'mycg', 'Matlab');
set(lgd,'FontSize', 16);  
ylabel('Relative Residual') ;
xlabel('Iteration');
title('NOS3 matrix','FontSize', 20);
%xlim([0 30]);
ylim([0 0.5]);


%------ WEST0479 ------%
west=spconvert(load('west0479.mtx'));


x_true = rand(size(west, 2), 1); 
b = west * x_true; 
x_0 = zeros(size(west,2), 1);

[westx_mmr, west_itermr, west_mmr] = mymr(west, b, x_0, 0.001, 5000);
[westx_mgm20, west_itergm20, west_mgm20] = mygmres(west, b, x_0, 0.001, 5000, 20);
[westx_mgm50, west_itergm50, west_mgm50] = mygmres(west, b, x_0, 0.001, 5000, 50);
[westx_mgm300, west_itergm300, west_mgm300] = mygmres(west, b, x_0, 0.001, 5000, 300);
%matlab function

[westx,flag,tempres,west_iter20, west_res20] = gmres(west,b, 20, 0.001, 5000);
west_res20 = west_res20/norm(b); 
index20 = 1:20:(west_iter20(1) * west_iter20(2)+1);
[westx,flag,tempres,west_iter50, west_res50] = gmres(west,b, 50, 0.001, 5000);
west_res50 = west_res50/norm(b);
index50 = 1:50:(west_iter50(1) * west_iter50(2)+1);
[westx,flag,tempres,west_iter300, west_res300] = gmres(west,b, 300, 0.001, 5000);
west_res300 = west_res300/norm(b);
index300 = [1 west_iter300(2)+1] ;




figure;
subplot(1,2,1);
plot(0:west_itermr, west_mmr, 'co-');
hold on ;
plot(0:west_itergm20, west_mgm20, 'r+-');
hold on;
plot(0:west_iter20(1), west_res20(index20), 'bs-');
hold off;
lgd = legend ('mymr', 'mygmres(m = 20)', 'Matlab(m = 20)'); 
  
ylabel('Relative Residual') ;
xlabel('Iteration');
title('west0479 matrix');
xlim([0 50]);
ylim([0 1]);

subplot(1,2,2);
plot(0:west_itermr, west_mmr, 'co-');
hold on ;
plot(0:west_itergm50, west_mgm50, 'r+-');
hold on;
plot(0:west_iter50(1), west_res50(index50), 'bs-');
hold off;
lgd = legend ('mymr', 'mygmres(m = 50)', 'Matlab(m = 50)'); 
set(lgd,'FontSize', 3);  
ylabel('Relative Residual') ;
xlabel('Iteration');
title('west0479 matrix');
xlim([0 50]);
ylim([0 1]);

figure;
subplot(1,2,1);
plot(0:west_itermr, west_mmr, 'co-');
hold on ;
plot(0:west_itergm300, west_mgm300, 'r+-');
hold on;
plot(0:1, west_res300(index300), 'bs-');
hold off;
lgd = legend ('mymr', 'mygmres(m = 300)', 'Matlab(m = 300)'); 
set(lgd,'FontSize', 3);  
ylabel('Relative Residual') ;
xlabel('Iteration');
title('west0479 matrix');
xlim([0 50]);
ylim([0 1]);

subplot(1, 2, 2);
plot(0:west_iter300(2), west_res300, 'bs-');


  ylabel('Relative Residual') ;
  xlabel('Inner Iteration');
  
  title('west0479 matrix with gmres(restart = 300)');
  xlim([0 80]);
  ylim([0 0.7]);




%------ mahindas ------%
mah=spconvert(load('mahindas.mtx'));

x_true = rand(size(mah, 2), 1); 
b = mah * x_true; 
x_0 = zeros(size(mah,2), 1); 
restrt = 4; 
[mahx_mmr, mah_itermr, mah_mmr] = mymr(mah, b, x_0, 0.001, 5000);
[mahx_mgm, mah_itergm, mah_mgm] = mygmres(mah, b, x_0, 0.001, 5000, restrt);
 
%matlab function

[mahx,flag,tempres,mah_iter, mah_res] = gmres(mah,b, restrt, 0.001, 5000);
mah_res = mah_res/norm(b); 



  index = 1:restrt:(mah_iter(1) * mah_iter(2)+1);




plot(0:mah_itermr, mah_mmr, 'ro-');
hold on ;
plot(0:mah_itergm, mah_mgm, 'cx-');
hold on;
plot(0:mah_iter(1), mah_res(index), 'bs-');
hold off;

lgd = legend ('mymr', 'mygmres(m = 4)', 'Matlab');
set(lgd,'FontSize', 16);
ylabel('Relative Residual') ;
xlabel('Iteration');
title('MAHINDAS matrix','FontSize', 20 );
xlim([0 50]);
ylim([0 1]);