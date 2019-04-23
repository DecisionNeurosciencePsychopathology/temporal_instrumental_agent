function a=getexploreid(x)
    [~,indx] = size(strsplit(x,filesep));
    indr = strsplit(x,filesep);
    a=char(indr(indx-1));
end