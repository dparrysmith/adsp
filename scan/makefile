scp_new: scp_new.o nrr_new.o
	cc -o scp_new nrr_new.o scp_new.o -pthread -lm

nrr_new.o: nrr_new.c
	cc -c nrr_new.c -lm -o nrr_new.o 

scp_new.o: scp_new.c
	cc -c scp_new.c -pthread -o scp_new.o 

