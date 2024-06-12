%xóa không gian làm việc và cửa sổ lệnh, sau đó đọc tệp âm thanh 't.wav'
clear; clc;
% Đọc file âm thanh
[audio, fs] = audioread('t.wav');

% Tạo bộ lọc thông thấp thiết kế bộ lọc thông thấp Butterworth bậc 6 với tần số cắt (fc) là 1000 Hz
fc = 1000; % Tần số cắt
[b, a] = butter(6, fc/(fs/2), 'low');

% Áp dụng bộ lọc thông thấp cho tín hiệu âm thanh 
filtered_audio = filter(b, a, audio);

% Tạo FFT cho tín hiệu gốc và tín hiệu lọc. Biến đổi Fourier nhanh (FFT) được tính toán cho cả tín hiệu âm thanh gốc và tín hiệu âm thanh được lọc. Hàm abs tính toán độ lớn của kết quả FFT, đưa ra biểu diễn miền tần số của tín hiệu.
fft_original = abs(fft(audio));
fft_filtered = abs(fft(filtered_audio));

% Tạo vector thời gian để vẽ đồ thị cho tín hiệu âm thanh
t = (0:length(audio)-1)/fs;

% Vẽ đồ thị
subplot(2,2,1);
plot(t, audio);
title('Tín hiệu âm thanh gốc (Analog)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,2,2);
plot(t, filtered_audio);
title('Tín hiệu sau khi lọc (Analog)');
xlabel('Time (s)');
ylabel('Amplitude');
% Tạo vectơ tần số để vẽ đồ thị
frequencies = linspace(0, fs, length(audio));
subplot(2,2,3);
plot(frequencies(1:length(audio)/2), fft_original(1:length(audio)/2));
title('Biểu đồ phổ của tín hiệu gốc (Digital)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(2,2,4);
plot(frequencies(1:length(audio)/2), fft_filtered(1:length(audio)/2));
title('Biểu đồ phổ của tín hiệu sau khi lọc (Digital)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

