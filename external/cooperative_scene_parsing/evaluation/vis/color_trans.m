function color = color_trans(ori_color)


color = ori_color;

color = [hex2dec(color(2:3)) hex2dec(color(4:5)) hex2dec(color(6:7))]/255;

color = color * 0.6;

end