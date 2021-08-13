all:
	g++ -O3 -std=c++14 -Wall mandelbrot.cpp -I/home/hyun/iRRAM/iRRAMx/include -L/home/hyun/iRRAM/iRRAMx/lib -liRRAMx -liRRAM -lmpfr -lgmp -o mandelbrot

test:
	g++ -std=c++14 -O3 -Wall mandelbrot.cpp -liRRAM -lmpfr -lgmp -o mandelbrot

clean:
	rm -rf mandelbrot
