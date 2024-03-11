建议使用python2.7运行，若使用python3.5，可能会有一些int和float转换的问题。
在main.py底下调用main函数，输入为文件路径。
需要在parameter.py修改当前配置，如采样点数、chirp数、cfar阈值等。
当前处理流程为，对输入数据（range fft）做速度维fft，随后进行cfar检测，对cfar检测的所有点做dbf，并保存dbf图。对所有点做完dbf后，会画rdm + 0速cfar门限图 + 前
27个点（如果有那么多点的话，若超出27个点无法显示）的点云数据，包含距离、速度索引值、水平角、俯仰角。另外还会保留一个存有所有点点云信息的txt。
