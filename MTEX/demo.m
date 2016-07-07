%% qpyTex MTEX Demonstration File

hexSym = crystalSymmetry('6/mmm', [3.23 3.23 5.15], 'X||a');
dataDir = '/home/chris/Python/qpyTex/dat/';
cd '/home/chris/Python/qpyTex/MTEX'
%%
n = 0;
dn = 1;
indx = [4125:4160];
step = 5;
for ind = indx
    pname = [dataDir sprintf('%05d',ind) '/'];
    if (n == 0)
        measPF = getPFfromFiles(pname, hexSym);
        n = 1;
    else
        nextPF = getPFfromFiles(pname,hexSym);
        measPF = union(measPF,rotate(getPFfromFiles(pname, hexSym),axis2quat(xvector, step * n * dn * degree)));
        n = n + 1;
    end
    fprintf('%2d / %2d\n',n,length(indx));
end
measPF = normalize(measPF);
%%
myPF = loadSample(dataDir, hexSym, 4089, 0);
figure(1), plot(myPF, 'contourf', [0:5])

myPF2 = loadSample(dataDir, hexSym, 4125, 0);
figure(2), plot(myPF, 'contourf', [0:5])

figure(3), plot(myPF + myPF2, 'contourf', [0:5])
%%


rot1 = axis2quat(xvector, (65) * degree);
rot2 = axis2quat(yvector, 90 * degree);
rot3 = axis2quat(xvector, 90 * degree);
figure(1), plot(rotate(measPF, rot2 * rot1),'contourf');
%%
figure(2),plotPDF(rotate(myODF,rot1),{Miller(0,0,2,hexSym), Miller(1,0,0,hexSym)},'contourf',0:5)
%%
figure(3), plot(rotate(measPF, rot1),'contourf',[0 5]);
%%
rotation_180 = axis2quat(xvector, 180 * degree);
figure(1),plot(measPF,'MarkerSize',4,'colorrange',[0 5], 'contourf')
%%
figure(2),plot(rotate(pf2,rotation_180),'MarkerSize',4,'colorrange',[0 5],'contourf')
%%
 odf1 = calcODF(pf1({1});
%%
% figure(3), plotPDF(odf2,pf2.h,'antipodal','contourf',[0:5])
%% Rotate to bring sample into proper orientation
rot1 = axis2quat(xvector, (65) * degree);
rot2 = axis2quat(yvector, 90 * degree);
figure(4), plot(rotate(rotate(pf2,rot1),rot2),'colorrange', [0 5])
%%
rot1 = axis2quat(xvector, (65) * degree);

figure(2), plot(rotate(measPF,rot1),'contourf',0:5)
%% Locate Area of Maximum Intensity and Plot on Figure
[theta,rho] = polar(measPF.allR{1});
PFmax = measPF({1}).max;
PFmaxInd = find(measPF({1}).allI{1}==PFmax);
maxTheta =theta(PFmaxInd);
maxRho = rho(PFmaxInd);
v = vector3d('polar',maxTheta,maxRho);
hold off, plot(measPF({1}),'contour'), hold on
plot(v)

%% Calculate Kearns factors
[f1, f2, f3] = Kearns(measPF({1}));