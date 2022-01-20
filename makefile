CC = g++

all: compile run_test

compile: main.cpp exffunction.cpp stdafx.h
	${CC} stopwatch.cpp exffunction.cpp -o ExpForce

run_test: ExpForce
	./ExpForce fb_full.txt 1
#	./ExpForce test_edgelist.txt 1
	
