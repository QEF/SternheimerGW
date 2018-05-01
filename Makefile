all:
	make -C util
	make -C data
	make -C algo
	make -C phys
	make -C user 

test:
	make -C util test
	make -C data test
	make -C algo test
	make -C phys test
	make -C user test 

clean:
	make -C util clean
	make -C data clean
	make -C algo clean
	make -C phys clean
	make -C user clean
	
.PHONY: all test clean
