% infA = [98 99; 97 98; 96 97];
% supA=[100 101; 99 100; 98 99];
% infb=[1;2;3];
% supb=[2;3;4];

infA = [2 0; 4 6; 0 2];
supA = [6 6; 8 8; 4 4];
infb = [-1; 3; -2];
supb = [11; 17; 10];

infA2 = [2 4 0; 0 6 2];
supA2 = [6 8 4; 6 8 4];
infb2 = [-7; -2];
supb2 = [35; 22];

infA2_ = [4 0; 6 2];
supA2_ = [8 4; 8 4];
infb2_ = [-7; -2];
supb2_ = [35; 22];

% EqnTol2D(infA, supA, infb, supb)
% infA = [98 99];
% supA=[100 101];
% infb=[1];
% supb=[2];

% % infA = [ 2 1 1; -1 15 3 ];
% % infA2 = [ 2 -1 ;1 15; 1 3 ];
% % supA = [ 4 5 4; 5 19 10];
% % supA2 = [ 4 5; 5 19; 4 10];
% % infb = [-2; 2];
% % supb = [1; 7];
% % 
% % infb2 = [-4.5; -5.5; 0.5];
% % supb2 = [-1.5; -0.5; 4.5];
% % OrientPoints = 1;
% % transparency = 1;
% % 
% EqnTolR3(infA2,supA2,infb2,supb2,0,1);
% scatter3(-0.4, 0.3, 0, 'm*');
figure()
% plot3([-1.0637 0.0897], [-0.2567 -0.2567], [-0.5767 -0.5767], 'g');
% hold on
% plot3([-1.0637 0.0897], [0.8967 0.8967], [-0.5767 -0.5767], 'g');
% plot3([-1.0637 0.0897], [-0.2567 -0.2567], [0.5767 0.5767], 'g');
% plot3([-1.0637 0.0897], [0.8967 0.8967], [0.5767 0.5767], 'g');
% plot3([-1.0637 -1.0637], [-0.2567  0.8967], [-0.5767 -0.5767], 'm');
% plot3([-1.0637 -1.0637], [-0.2567  0.8967], [0.5767 0.5767], 'm');
% plot3([0.0897 0.0897], [-0.2567  0.8967], [-0.5767 -0.5767], 'm');
% plot3([0.0897 0.0897], [-0.2567  0.8967], [0.5767 0.5767], 'm');
% plot3([-1.0637 -1.0637], [-0.2567 -0.2567], [-0.5767 0.5767], 'm');
% plot3([-1.0637 -1.0637], [0.8967 0.8967], [-0.5767 0.5767], 'm');
% plot3([0.0897 0.0897], [-0.2567 -0.2567], [-0.5767 0.5767], 'm');
% plot3([0.0897 0.0897], [0.8967 0.8967], [-0.5767 0.5767], 'm');
% hold off
% EqnTolR2(infA,supA,infb,supb)
% hold on
% plot(0.54285, 0.48571, 'r.')
% text(0.54285, 0.58571, 'argmaxTol')
% hold off
x = [-0.278, -0.278, 3.135, 3.135, -0.278];
y = [-1.708, 1.708, 1.708, -1.708, -1.708];
EqnTolR2(infA2_, supA2_, infb2_, supb2_)
hold on
plot(1.42857034, 1.63777828e-06, 'r.')
plot(x,y, 'k')
text(1.42857034, 0.6, 'argmaxTol')
hold off


