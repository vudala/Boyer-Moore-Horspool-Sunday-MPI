# Eduardo Vudala Senoski GRR20195689

FLAGS=-g -O3

LIBS=-lm

CC=mpicc

RM=rm -f

EXEC=dna

all: $(EXEC)

$(EXEC):
	$(CC) $(FLAGS) dna.c utils.c -o $(EXEC) $(LIBS)

run: all
	mpirun --hostfile hosts_file.txt -np 8 ./$(EXEC)

clean:
	@$(RM) dna.out

purge: clean
	@$(RM) $(EXEC)