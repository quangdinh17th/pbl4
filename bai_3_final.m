clear; clc;
Nf = 10;
L = 3; 
max_delay = 1.25e-5; 
decay_base = 1;
Nc = 64; 
T = 8*max_delay;
t_step = (T/Nc)/16;
f_delta = 1/T;
t_vector = 0:t_step:(T-t_step); %=t_step*(0:Nc-1)
Ns = length(t_vector); %Number of samples per OFDM symbol(T seconds) before inserting Cyclic Preﬁx
GI = 1/4;
Ns_in_GI = ceil(Ns*GI);
Ns_total = Ns+Ns_in_GI;
Eb = 1;
EbN0dBvector = 0:3:18;
for k = 0:(Nc-1)
    subcarrier = 1/sqrt(T)*exp(j*2*pi*k*f_delta*t_vector);
    subcarrier_matrix(k+1,:) = subcarrier;
end
pilot_sk = ones(1,Nc);
xt = zeros(1,Ns);
for k = 0:(Nc-1)
    s_k = pilot_sk(k+1);
    xt = xt+s_k*subcarrier_matrix(k+1,:);
end
xt_tail = xt((Ns-Ns_in_GI+1):Ns);
pilot_OFDM_symbol = [xt_tail xt];
for snr_i = 1:length(EbN0dBvector)
    EbN0dB = EbN0dBvector(snr_i);
    EbN0 = 10^(EbN0dB/10);
    N0 = Eb/EbN0;
    vn = N0/(2*t_step); % Refer to Section 1.C and 1.D of Chapter 21.
    bitcnt = 0; errcnt=0;
        while errcnt < 1000 %Reduce 1000 to a smaller value if the simulation takes long time. The simulation error will increase instead.
            OFDM_frame = [];
        for m=1:Nf
            datasymbols_in_OFDM_symbol = sign(rand(1,Nc)-0.5)+j*sign(rand(1,Nc)-0.5);
            datasymbols_in_OFDMframe(m,:) = datasymbols_in_OFDM_symbol;
            xt = zeros(1,Ns);
            for k = 0:(Nc-1)
                s_k = datasymbols_in_OFDM_symbol(k+1);
                xt = xt+s_k*subcarrier_matrix(k+1,:);
            end
            xt_tail = xt((Ns-Ns_in_GI+1):Ns);
            xt = [xt_tail xt];
            OFDM_frame = [OFDM_frame xt];
        end
ht = ht_mp_ch(max_delay,L,decay_base,t_step);
OFDM_frame_after_ht = conv(OFDM_frame,ht);
frame_sample_length = length(OFDM_frame_after_ht);
noise = sqrt(vn)*(randn(1,frame_sample_length)+j*randn(1,frame_sample_length));
rt_frame = OFDM_frame_after_ht+noise;
rt_pilot = conv(pilot_OFDM_symbol,ht);
for k = 0:(Nc-1)
    D = t_step*sum(rt_pilot(Ns_in_GI+(1:Ns)).*conj(subcarrier_matrix(k+1,:))/sqrt(T));
    F(k+1) = D/pilot_sk(k+1);
end
for m = 1:Nf
    first_index_of_mth_OFDMsymbol=(m-1)*Ns_total+ Ns_in_GI + 1;
    mth_OFDM_symbol_in_rt = rt_frame(first_index_of_mth_OFDMsymbol+(0:Ns-1));
for k = 0:(Nc-1)
D = t_step*sum(mth_OFDM_symbol_in_rt.*conj(subcarrier_matrix(k+1,:))/sqrt(T));
Dc = D/F(k+1); 
estimated_data_symbols_in_OFDMframe(m,k+1) = sign(real(Dc))+j*sign(imag(Dc));
end
end
Ierrs = sum(sum(real(datasymbols_in_OFDMframe) ~= real(estimated_data_symbols_in_OFDMframe)));
Qerrs = sum(sum(imag(datasymbols_in_OFDMframe) ~= imag(estimated_data_symbols_in_OFDMframe)));
errcnt = errcnt+(Ierrs+Qerrs);
bitcnt = bitcnt+Nc*Nf*2;
end
BER(snr_i) = errcnt/bitcnt
%BERtheory(snr_i)=qfunc(sqrt(2*EbN0))
BERtheory(snr_i) = 1/2 - EbN0^(1/2)/(2*(EbN0 + 1)^(1/2))
end
figure

semilogy(EbN0dBvector, BER,'b')
hold on
semilogy(EbN0dBvector, BERtheory,'r')
grid
xlabel('Eb/N0 (dB)');
ylabel('BER')
legend('BER simulation','BER theory (Rayleigh fading)')