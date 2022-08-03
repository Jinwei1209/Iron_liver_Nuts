function kspace = myScanArchiveRead_T1wT2wQSM(filename, nframe, necho_QSM)

if (nargin<1)
    filename = dir('.\*.h5');
    filename = filename(end).name;
end
archive = GERecon('Archive.Load', filename);

xRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
yRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres - 1;
stop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.stop_rcv;
start = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.start_rcv;
necho = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nechoes;
nChannels = stop - start + 1;

% Keep track of the current pass
pass = 1;
zRes = archive.SlicesPerPass(pass);
idx_all = 0;
T = nframe*(necho_QSM+2);

% Allocate K-space
% kspace = complex(zeros(xRes, yRes, nChannels, zRes, necho));
kspace = complex(zeros(xRes, yRes, nChannels, zRes, necho_QSM+2));
for i = 1:archive.ControlCount
    
    if (mod(i, 100) == 0)
        disp(i);
    end
    
    control = GERecon('Archive.Next', archive);
    
    % Sort only programmable packets in range
    if(control.opcode == 1 && ...
            control.viewNum > 0 && ...
            control.viewNum <= yRes && ...
            control.sliceNum < zRes)
        
        idx_all = idx_all + 1;
        idx = mod(idx_all, T);
        if idx == 0
            idx = T;
        end
        if (idx <= nframe)
            kspace(:,control.viewNum,:,control.sliceNum+1,necho_QSM+1) = control.Data;
        else
            if (idx> nframe && idx <= nframe*(necho_QSM+1))
                kspace(:,control.viewNum,:,control.sliceNum+1,control.echoNum+1) = control.Data;
            else
                kspace(:,control.viewNum,:,control.sliceNum+1,necho_QSM+2) = control.Data;
            end
        end
    end
end
% kspace = ifft(kspace,[],4);
kspace = permute(kspace, [1, 2, 4, 3, 5]);

end