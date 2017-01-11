#! /usr/bin/python

import dpath
import subprocess, re, tempfile

import os, shutil, random

import numpy


try:
  os.mkdir('scratch')
except:
  pass


def run( cmd, workdir = './' ):
  output = ""
  with tempfile.TemporaryFile() as f:
    status = subprocess.call( "cd %s; %s" % (workdir,cmd), stdout=f, stderr=f, shell=True )
    f.seek(0)
    output += f.read()

  return status,output

def tail(text, n=10):
  return '\n'.join( text.split('\n')[-n:-1] )

def close(a,b,tol=0.001):
  return abs(a-b) < (a+b)*tol


def test_linear_interp():

  with open("tests.linear_interp.input-data.txt","w") as f:
    f.write("0.0 10.0\n")
    f.write("1.0 30.0\n")
    f.write("2.0 50.0\n")
    f.write("3.0 70.0\n")
    f.write("4.0 90.0\n")

  with open("tests.linear_interp.x-values.txt","w") as f:
    f.write("0.5 0.0\n")
    f.write("1.5 1.0\n")
    f.write("2.5 2.0\n")
    f.write("3.5 3.0\n")

  status,output = run( "./interp-cli --method=linear tests.linear_interp.input-data.txt tests.linear_interp.x-values.txt tests.linear_interp.output.txt" )

  if status == 0:
    print tail(output)
  else:
    print output
  
  assert status == 0
  assert re.search("ERROR", output) is None


  data = numpy.loadtxt("tests.linear_interp.output.txt")

  assert close( data[0][0] ,0.5  )
  assert close( data[0][1] ,20.0 )
  assert close( data[1][0] ,1.5  )
  assert close( data[1][1] ,40.0 )
  assert close( data[2][0] ,2.5  )
  assert close( data[2][1] ,60.0 )
  assert close( data[3][0] ,3.5  )
  assert close( data[3][1] ,80.0 )




