function [mt, r, tcspc, head, mthigh, num ] = SpcReader(name, cnts)

if ~(name(end-2:end)=='spc')
    disp('Not a t3r-file!');
    return
end
fin = fopen(name,'r');
if (fin==-1)
    errordlg('Cannot open specified file. Please try again.');
    return
end
fseek(fin, 0,1);
head.NCounts = ftell(fin)/6;
if nargin>1
    fseek(fin, 0,-1);
    fseek(fin, 6*(cnts(1)), -1);
    [tmpx num] = fread(fin, cnts(2)/6, 'ubit48');
    tmpx=uint64(tmpx);
end
if num>0
    tcspc = 4095 - uint32(bitand(tmpx, hex2dec('fff')));
    invalid = uint32(bitand(bitshift(tmpx, -12), 1));
    mtov = double(bitand(bitshift(tmpx, -13), 1));
    gap = uint32(bitand(bitshift(tmpx, -14), 1));
    mt = uint32(bitand(bitshift(tmpx, -32), hex2dec('ffff'))) + uint32(bitand(tmpx, hex2dec('ff0000')));
    r = uint32(bitand(bitshift(tmpx, -24), hex2dec('ff')));

    mtovcnt = uint64(cumsum(mtov));
    mt = (mt + bitshift(uint32(bitand(mtovcnt,hex2dec('ff'))),24));
    mthigh = double(bitshift(mtovcnt, -8))*2^32;
%     mt = mt(invalid==0);

end
head.NChannels = 4096;
head.GlobClock = 50;
head.Resolution = 0.008;
head.Ident = 'SPC-630    ';
fclose(fin);
