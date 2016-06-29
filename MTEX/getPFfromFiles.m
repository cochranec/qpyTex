function myPF = getPFfromFiles(fname, mSym)
    myDir = dir([fname '*']);
    mH = {};
    fnames = {};
    for i=1:numel(myDir)
        if ~strcmp(myDir(i).name(1),'.') && ~strcmp(myDir(i).name(end-3:end),'.png') && ~strcmp(myDir(i).name(end-3:end),'.002')
            h = str2num(myDir(i).name(end-2));
            k = str2num(myDir(i).name(end-1)); %#ok<*ST2NM>
            l = str2num(myDir(i).name(end)); 
            mH{i} = Miller(h,k,-h-k,l,mSym);
            fnames{i} = [fname myDir(i).name];
        end
    end
    mH = mH(~cellfun('isempty',mH));
    fnames = fnames(~cellfun('isempty',fnames));
    pf = loadPoleFigure(fnames, mH,mSym,'ColumnNames', {'Polar Angle' 'Azimuth Angle' 'Intensity'},'infolvl',0);
%     for i=1:numel(pf)
%         cpf = pf(i);
%         cpf(cpf.isOutlier) = [];
%         pf(i) = cpf;
%     end

    myPF = pf;
end