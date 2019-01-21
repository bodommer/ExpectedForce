CC = g++

all: compile run_test

compile: main.cpp exffunction.cpp stdafx.h
	${CC} main.cpp exffunction.cpp -o ExpForce

run_test: ExpForce
	./ExpForce fb_full
	
