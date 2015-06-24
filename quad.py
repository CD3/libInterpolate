#!/usr/bin/env python2
# encoding: utf-8

xmin = -10.
xmax = 10.
N = 200
dx = (xmax - xmin)/(N-1)

f = open("test.txt", "w")

for i in range(N):
    x = xmin+i*dx
    y = x*x

    f.write(str(x)+"   "+str(y)+"\n")
