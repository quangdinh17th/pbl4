function impulse_response=ht_mp_ch(max_delay,L,decay_base,t_step)
t_vector=0:t_step:max_delay;
mp_tmp=0*(t_vector);
path_delays=[0 sort(rand(1,L-1)*max_delay)];
impulse_positions=floor(path_delays/t_step);
mp_tmp(impulse_positions+1)=exp(j*2*pi*rand(1,L));
mp_tmp=mp_tmp.*(decay_base.^(t_vector/max_delay));
impulse_response=mp_tmp/sqrt(sum(abs(mp_tmp).^2));