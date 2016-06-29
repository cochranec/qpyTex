%% qpyTex MTEX Demonstration File

hexSym = crystalSymmetry('6/mmm', [3.23 3.23 5.15], 'X||a');
dataDir = '/home/chris/Python/qpyTex/dat/';
cd '/home/chris/Python/qpyTex/MTEX'
%%
n = 0;
dn = 1;
indx = [4089:4124];
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
myODF = calcODF(normalize(measPF),'resolution', 10 * degree);
%%

rot1 = axis2quat(xvector, (65) * degree);
rot2 = axis2quat(yvector, 90 * degree);
rot3 = axis2quat(xvector, 90 * degree);
figure(1), plot(rotate(measPF, rot2 * rot1),'contourf');
figure(2),plotPDF(rotate(odf2,rot2 * rot1),{Miller(0,0,2,hexSym), Miller(1,0,0,hexSym)})
%%
figure(3), plot(rotate(normalize(measPF), rot1),'colorrange',[0 5]);