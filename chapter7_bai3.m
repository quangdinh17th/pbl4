clear;
Nf = 10;
L = 2; 
max_delay = 1.25e-5; 
decay_base = 1;
Nc = 16; 
T = 8 * max_delay;
t_step = (T / Nc) / 16;
f_delta = 1 / T;
t_vector = 0:t_step:(T - t_step);
Ns = length(t_vector);
GI = 1/4;
Ns_in_GI = ceil(Ns * GI);
Ns_total = Ns + Ns_in_GI;
Eb = 1;
EbN0dBvector = 0:3:18;
subcarrier_matrix = zeros(Nc, Ns);
F = zeros(1, Nc);
datasymbols_in_OFDMframe = zeros(Nf, Nc);
estimated_data_symbols_in_OFDMframe = zeros(Nf, Nc);
for k = 0:(Nc - 1)
    subcarrier = 1/sqrt(T) * exp(1i * 2 * pi * k * f_delta * t_vector);
    subcarrier_matrix(k + 1, :) = subcarrier;
end
pilot_sk = ones(1, Nc);
xt = zeros(1, Ns);
for k = 0:(Nc - 1)
    s_k = pilot_sk(k + 1);
    xt = xt + s_k * subcarrier_matrix(k + 1, :);
end
xt_tail = xt((Ns - Ns_in_GI + 1):Ns);
pilot_OFDM_symbol = [xt_tail xt];
for snr_i = 1:length(EbN0dBvector)
    EbN0dB = EbN0dBvector(snr_i);
    EbN0 = 10^(EbN0dB / 10);
    N0 = Eb / EbN0;
    vn = N0 / (2 * t_step);
    bitcnt = 0;
    errcnt = 0;
    while errcnt < 1000
        OFDM_frame = [];
        for m = 1:Nf
            datasymbols_in_OFDM_symbol = sign(rand(1, Nc) - 0.5) + 1i * sign(rand(1, Nc) - 0.5);
            datasymbols_in_OFDMframe(m, :) = datasymbols_in_OFDM_symbol;
            xt = zeros(1, Ns);
            for k = 0:(Nc - 1)
                s_k = datasymbols_in_OFDM_symbol(k + 1);
                xt = xt + s_k * subcarrier_matrix(k + 1, :);
            end
            xt_tail = xt((Ns - Ns_in_GI + 1):Ns);
            xt = [xt_tail xt];
            OFDM_frame = [OFDM_frame xt];
        end
        % Rayleigh fading channel impulse response
        ht = sqrt(1/2) * (randn(1, length(OFDM_frame)) + 1i * randn(1, length(OFDM_frame)));
        OFDM_frame_after_ht = conv(OFDM_frame, ht);
        frame_sample_length = length(OFDM_frame_after_ht);
        noise = sqrt(vn) * (randn(1, frame_sample_length) + 1i * randn(1, frame_sample_length));
        rt_frame = OFDM_frame_after_ht + noise;
        rt_pilot = conv(pilot_OFDM_symbol, ht);
        for k = 0:(Nc - 1)
            D = t_step * sum(rt_pilot(Ns_in_GI + (1:Ns)) .* conj(subcarrier_matrix(k + 1, :)) / sqrt(T));
            F(k + 1) = D / pilot_sk(k + 1);
        end
        for m = 1:Nf
            first_index_of_mth_OFDMsymbol = (m - 1) * Ns_total + Ns_in_GI + 1;
            mth_OFDM_symbol_in_rt = rt_frame(first_index_of_mth_OFDMsymbol + (0:Ns - 1));
            for k = 0:(Nc - 1)
                D = t_step * sum(mth_OFDM_symbol_in_rt .* conj(subcarrier_matrix(k + 1, :)) / sqrt(T));
                Dc = D / F(k + 1); 
                estimated_data_symbols_in_OFDMframe(m, k + 1) = sign(real(Dc)) + 1i * sign(imag(Dc));
            end
        end
        Ierrs = sum(sum(real(datasymbols_in_OFDMframe) ~= real(estimated_data_symbols_in_OFDMframe)));
        Qerrs = sum(sum(imag(datasymbols_in_OFDMframe) ~= imag(estimated_data_symbols_in_OFDMframe)));
        errcnt = errcnt + (Ierrs + Qerrs);
        bitcnt = bitcnt + Nc * Nf * 2;
    end
    BER(snr_i) = errcnt / bitcnt;
    BERtheory(snr_i) = 1/2 - EbN0^(1/2) / (2 * (EbN0 + 1)^(1/2));
end

figure
semilogy(EbN0dBvector, BER, 'b')
hold on
semilogy(EbN0dBvector, BERtheory, 'r')
grid
legend('BER simulation', 'BER theory (Rayleigh fading)')
