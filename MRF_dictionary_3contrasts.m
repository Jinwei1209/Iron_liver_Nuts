function D = MRF_dictionary_3contrasts(T1, T2, T2s, alpha, IR_eff, R, TR, TEs, nframe)

c2g = @(x) x;
g2c = @(x) x;

%% setup TR and TD/TE
Nrep = 5;   % number of repetations
necho = length(TEs);
TI = 0.0191 - 0.01606 + 0.299;  % + 0.299
ncontrast_T1 = 1;
ncontrast_mGRE = 1;
ncontrast_T2 = 1;
Twait = 1.0;  % 1
T2echo = 0.0825 + 0.00;  % additinal 0.02?
ncontrast = ncontrast_T1 + ncontrast_mGRE + ncontrast_T2;
TD = TR - TEs(end);
% TD_mGRE = TD + 0.00388*8;
TD_mGRE = 0.01898 - TEs(end);

if isempty(T2s)
    T2s = [0, 0];
end

%% Create a Vector with T1 and T2
if isempty(IR_eff) == 0
    T1 = repmat(T1,  [length(T2), 1, length(T2s), length(alpha)]);
    T2 = repmat(T2', [1, size(T1,2), size(T1,3), length(alpha)]);
    T2s = repmat(reshape(T2s, [1 1 length(T2s)]), [size(T2,1), size(T1,2), 1, length(alpha)]);
    alpha  = repmat(reshape(alpha, [1 1 1 length(alpha)]), [size(T2,1), size(T1,2), size(T2s, 3), 1]);
    alpha = alpha(:); T1 = T1(:); T2 = T2(:); T2s = T2s(:);
    
    num_TR_eff = length(IR_eff);
    IR_eff = repmat(IR_eff, [length(alpha), 1]);
    alpha = repmat(alpha, [num_TR_eff, 1]);
    T1 = repmat(T1, [num_TR_eff, 1]);
    T2 = repmat(T2, [num_TR_eff, 1]);
    T2s = repmat(T2s, [num_TR_eff, 1]);
    IR_eff = IR_eff(:);
    
    T12sidx = (T1>=T2) .* (T2>=T2s);
    T1 = T1(T12sidx==1);
    T2 = T2(T12sidx==1);
    T2s  =  T2s(T12sidx==1);
    alpha  =  alpha(T12sidx==1);
    IR_eff  =  IR_eff(T12sidx==1);
    T12_dict = [T1, T2, T2s, alpha, IR_eff];
else
    T1 = repmat(T1, [length(T2), 1]);
    T2 = repmat(T2', [1 size(T1,2)]);
    T1 = T1(:); T2 = T2(:);
    T12_dict = [T1, T2];
end
T1 = c2g(T1.');
T2 = c2g(T2.');
T2s = c2g(T2s.');
alpha = c2g(alpha.');
IR_eff = c2g(IR_eff.');

if T2s(1) == 0
    T2s = T2;
end


%% Bloch Simulation in the complex SO(2) group
Rot = @(alpha) [cos(alpha) -sin(alpha);
    sin(alpha)  cos(alpha)];
Rot_reshape = @ (Rot) permute(cat(3, Rot(:, 1:length(alpha)), Rot(:, length(alpha)+1:end)), [3, 2, 1]);

x = c2g(zeros(ncontrast,length(T1)));
% z = c2g(zeros(ncontrast,length(T1)));

M = complex(zeros(2, length(T1)));
M = c2g(M);
M(2, :) = 1;
for rep = 1:Nrep
    if rep == Nrep
        ip = 1;
    end
    disp(rep);
    %% IR prep
    % 180 degree tip down (considering inversion efficiency)
    M(2, :) = M(2, :) .* IR_eff;
%         Mtmp = Rot(slice_profile(is)*pi) * real(M);
%         M = Mtmp + 1i*imag(M);
    % gradient spoiler
    M(1, :) = 0;
    % IR time relaxation
    M(1,:) = M(1,:) .* exp(-TI./T2);
    M(2,:) = ones(1,size(M,2)) + (M(2,:)-ones(1,size(M,2))) .* exp(-TI./T1);
    
    for t = 1:ncontrast_T1
        for i = 1:nframe
            % Apply RF-pulse (that now acts only on the real part)
            Mtmp = permute(squeeze(sum(Rot_reshape(Rot(alpha)) .* real(M), 1)), [2, 1]);
            M = Mtmp + 1i*imag(M);

            for j = 1:necho
                % Relaxation
                if j == 1
                    delta_te = TEs(j);
                else
                    delta_te = TEs(j)-TEs(j-1);
                end
                M(1,:) = M(1,:) .* exp(-delta_te./T2s);
                M(2,:) = ones(1,size(M,2)) + (M(2,:)-ones(1,size(M,2))) .* exp(-delta_te./T1);

                % store the signal at the echo time
                if rep == Nrep && i == nframe && j == 1
                    x(ip,:) = M(1,:);
    %                 z(ip,:) = M(2,:);
                    ip = ip + 1;
                end
            end

            % gradient spoiler
            M(1, :) = 0;

            % Relaxation
            M(1,:) = M(1,:) .* exp(-TD./T2);
            M(2,:) = ones(1,size(M,2)) + (M(2,:)-ones(1,size(M,2))) .* exp(-TD./T1);
        end
    end
    
    %% mGRE
    for t = 1:ncontrast_mGRE
        for i = 1:nframe
            % Apply RF-pulse (that now acts only on the real part)
            Mtmp = permute(squeeze(sum(Rot_reshape(Rot(alpha*1.875)) .* real(M), 1)), [2, 1]);
            M = Mtmp + 1i*imag(M);

            for j = 1:necho
                % Relaxation
                if j == 1
                    delta_te = TEs(j);
                else
                    delta_te = TEs(j)-TEs(j-1);
                end
                M(1,:) = M(1,:) .* exp(-delta_te./T2s);
                M(2,:) = ones(1,size(M,2)) + (M(2,:)-ones(1,size(M,2))) .* exp(-delta_te./T1);

                % store the signal at the echo time
                if rep == Nrep && i == nframe && j == 1
                    x(ip,:) = M(1,:);
    %                 z(ip,:) = M(2,:);
                    ip = ip + 1;
                end
            end

            % gradient spoiler
            M(1, :) = 0;

            % Relaxation
            M(1,:) = M(1,:) .* exp(-TD_mGRE./T2);
            M(2,:) = ones(1,size(M,2)) + (M(2,:)-ones(1,size(M,2))) .* exp(-TD_mGRE./T1);
        end
    end

    % wait time Relaxation
    M(1,:) = M(1,:) .* exp(-Twait./T2);
    M(2,:) = ones(1,size(M,2)) + (M(2,:)-ones(1,size(M,2))) .* exp(-Twait./T1);
    
    %% T2prep
    % 90 degree flip
    Mtmp = Rot(pi/2) * real(M);
    M = Mtmp + 1i*imag(M);
    % T2echo relaxation
    M(1,:) = M(1,:) .* exp(-T2echo./T2);
    % 90 degree flip back
    Mtmp = Rot(-pi/2) * real(M);
    M = Mtmp + 1i*imag(M);
    
    for t = 1:ncontrast_T2
        for i = 1:nframe   %nframe or nframe-1? something wrong in the sequence
            % Apply RF-pulse (that now acts only on the real part)
            Mtmp = permute(squeeze(sum(Rot_reshape(Rot(alpha*0.75)) .* real(M), 1)), [2, 1]);
            M = Mtmp + 1i*imag(M);

            for j = 1:necho
                % Relaxation
                if j == 1
                    delta_te = TEs(j);
                else
                    delta_te = TEs(j)-TEs(j-1);
                end
                M(1,:) = M(1,:) .* exp(-delta_te./T2s);
                M(2,:) = ones(1,size(M,2)) + (M(2,:)-ones(1,size(M,2))) .* exp(-delta_te./T1);

                % store the signal at the echo time
                if rep == Nrep && i == 1 && j == 1  % i == 1 or i == 2? something wrong in the sequence
                    x(ip,:) = M(1,:);
    %                 z(ip,:) = M(2,:);
                    ip = ip + 1;
                end
            end

            % gradient spoiler
            M(1, :) = 0;

            % Relaxation
            M(1,:) = M(1,:) .* exp(-TD./T2);
            M(2,:) = ones(1,size(M,2)) + (M(2,:)-ones(1,size(M,2))) .* exp(-TD./T1);

        end
    end

%     % wait time Relaxation
%     M(1,:) = M(1,:) .* exp(-Twait./T2);
%     M(2,:) = ones(1,size(M,2)) + (M(2,:)-ones(1,size(M,2))) .* exp(-Twait./T1);
end

if ~isempty(R)
    x_ = x;
    % SVD of the dictionary
    [u_,s_,~]=svd(g2c(x_), 'econ');
    u_ = u_(:,1:R);
    x_ = u_'*x_;
    D.u_    = u_;
    D.s_    = diag(s_);
    
    sos_dict_ = l2_norm(x_, 1);
    x_ = (x_./repmat(sos_dict_, [size(x_,1) 1]));
    D.normalization_ = g2c(sos_dict_);
    D.magnetization_ = g2c(x_);
    
else
    D.u_ = [];
end

sos_dict = l2_norm(x, 1);
x = (x./repmat(sos_dict, [size(x,1) 1]));

D.magnetization = x;
% D.z             = g2c(z);
D.normalization = sos_dict;
D.lookup_table  = T12_dict;
D.TR            = TR;
D.TEs            = TEs;
D.parameter{1}  = 'T1';
D.parameter{2}  = 'T2';
if length(alpha) > 1
    D.parameter{3}  = 'alpha';
end

% reset(gpuDevice(2));


end