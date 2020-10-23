gcc fast_tri.c -o fast_tri -O3
set /A range = 2
set /A iter = 10000
fast_tri dummy %range% %iter%
python unit-test-print.py %range% %iter%
