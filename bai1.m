% Dòng này đặt hệ số lấy mẫu quá mức thành 4, nghĩa là mỗi bit trong tín hiệu gốc sẽ được mở rộng thành 4 mẫu trong tín hiệu lấy mẫu quá mức.
overSampling_Factor=4;
% Xác định chuổi bit đầu vào là [1 1 1 0]
input_bit = [1 1 1 0];
% Lấy mẫu lại chuỗi bit đầu vào. Chức năng lấy mẫu tăng tốc độ lấy mẫu theo hệ số lấy mẫu quá mức. Cụ thể, nó chèn các số 0 (overSampling_Factor - 1) giữa mỗi bit của chuỗi đầu vào. Chuỗi kết quả input_bit_os sẽ là: [1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0]
input_bit_os = upsample(input_bit,overSampling_Factor); 
% xác định đáp ứng xung của bộ lọc hình chữ nhật. Bộ lọc có độ dài là 10, với 4 giá trị đầu tiên là 1/sqrt(overSampling_Factor) và 6 giá trị còn lại là số 0. pt= [0.5 0.5 0.5 0.5 0 0 0 0 0 0]
pt = [ones(1,overSampling_Factor) 0 0 0 0 0 0]/sqrt(overSampling_Factor);
% Hàm tích chập kết hợp tín hiệu đầu vào được lấy mẫu quá mức với bộ lọc hình chữ nhật. Tích chập là một phép toán được sử dụng để áp dụng bộ lọc cho tín hiệu đầu vào, giúp làm mịn hoặc lấy trung bình tín hiệu một cách hiệu quả trên số lượng mẫu được chỉ định.
output_of_rect_filter = conv(input_bit_os, pt);
% tín hiệu đầu ra của bộ lọc hình chữ nhật sử dụng biểu đồ gốc, hiển thị các điểm dữ liệu rời rạc. Biểu đồ bao gồm tiêu đề và nhãn cho trục x (mẫu) và trục y (biên độ).
stem(output_of_rect_filter);
title('Output of Rectangular Filter at Tx side')
xlabel('Samples')
ylabel('Amplitude')
