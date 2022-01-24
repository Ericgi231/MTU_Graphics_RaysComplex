#standard build
ray3: main.c stb_image_write.h
	gcc -g -Wall -Wextra -std=gnu99 -o ray3 main.c -lm

#build fancy version with anti aliasing and multiple image output
fancy: mainFancy.c stb_image_write.h
	gcc -g -Wall -Wextra -std=gnu99 -o ray3Fancy mainFancy.c -lm

#build multiple images into gif
gif:
	convert -delay 0 -loop 0 *.png zoom.gif

#clean files
.PHONY: clean
clean:
	rm -f ray3 ray3Fancy *.png zoom.gif ray3.zip

#zip for submit
submit:
	zip ray3.zip main.c mainFancy.c makefile stb_image_write.h