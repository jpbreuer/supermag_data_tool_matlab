function jd2000 = jd2000_new(iy, im, id, ut)
% Julian day 2000
% Input
%  year  = calendar year  [yyyy]
%  month = calendar month [1 - 12]
%  day   = calendar day   [1 - 31]
%  ut    = Universal Time [0 - 24]hours
% Output
%  jd2000 = Julian day 2000

y = iy;
m = im;

if (m <= 2)
   y = y - 1;
   m = m + 12;
end

a = fix(y / 100);
b = 2 - a + floor(a / 4);

jd = fix(365.25 * y) + fix(30.6001 * (m + 1));
% jd2000 = jd + id + b + 1720994.5 + ut/24 - 2451544.5;
jd2000 = jd + id + b - 730550 + ut/24;