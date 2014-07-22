#
#
#

CC = clang

OPENCV = -I/usr/include/opencv -lopencv_highgui -lopencv_core -lopencv_imgproc

#
#
#

output:
	$(CC) rnt.c tme.c -lm -lrt $(OPENCV) -DUSE_RGB -O3 -o rnt
	$(CC) lrn.c tme.c -lm -lrt $(OPENCV) -O3 -o lrn
	$(CC) wrp.c -lm -lrt $(OPENCV) -O3 -o wrp