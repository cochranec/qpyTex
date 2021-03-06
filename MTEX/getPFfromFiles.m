function myPF = getPFfromFiles(fname, mSym)
    myDir = dir([fname '*']);
    mH = {};
    fnames = {};
    for i=1:numel(myDir)
        if ~strcmp(myDir(i).name(1),'.') && ~strcmp(myDir(i).name(end-3:end),'.png')
            h = str2num(myDir(i).name(end-2));
            k = str2num(myDir(i).name(end-1));
            l = str2num(myDir(i).name(end)); 
            mH{i} = Miller(h,k,l,mSym);
            fnames{i} = [fname myDir(i).name];
        end
    end
    mH = mH(~cellfun('isempty',mH));
    fnames = fnames(~cellfun('isempty',fnames));
    pf = loadPoleFigure(fnames, mH,mSym,'ColumnNames', {'Polar Angle' 'Azimuth Angle' 'Intensity'},'InfoLevel',0);
    myPF = pf;
end