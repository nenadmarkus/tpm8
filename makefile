#
#
#

OPENCV = -I/usr/include/opencv -lopencv_highgui -lopencv_core -lopencv_imgproc

#
#
#

output:
	$(CC) rnt.c -lm -lrt $(OPENCV) -O3 -o rnt
	$(CC) wrp.c -lm -lrt $(OPENCV) -O3 -o wrp
	$(CC) lrn.c -lm -lrt $(OPENCV) -O3 -o lrn