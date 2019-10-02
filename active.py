#!/usr/bin/python
# -*- coding: utf-8 -*-

# coding=utf-8
# run active task for /pa/ and /pe/ stimuli continuum
# Adrielle C. Santana  2019

import sys, os
reload(sys)
sys.setdefaultencoding('utf8')

import easygui
import re
import subprocess
import numpy
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import random
from random import shuffle, randint
import wave
import contextlib
import pyxid
import time
from pygame import mixer
import termios, fcntl
import select
import pickle

### Enter file name
file_name=easygui.enterbox('Enter dat file name for the subject: ')

### Read file and identificate middle continuum stimuli
ro.r('load(file=file.path("Identification/data/'+file_name+'"))')
pa_amb=ro.r("pa.amb")[0]
amb=ro.r("amb")[0]
pe_amb=ro.r("pe.amb")[0]

#Amount of repetitions of each one of the 5 stimuli selected      
rep=200

#Total number of blocks
nblocks=5

#Index of the stimuli
idx=[1,pa_amb,amb,pe_amb,200]*rep

#Total number of trials
ntrials=len(idx)

for pos in range(len(idx)):
    if ((pos+1)%2==0):
        idx[pos]=-idx[pos]

shuffle(idx)

#Duration of the wav file. Use any wav file once their duration are the same

fname = "/home/adrielle/Desktop/Piloto/Formantes/stim-1.wav"
with contextlib.closing(wave.open(fname,'r')) as f:
    frames = f.getnframes()
    rate = f.getframerate()
    duration = (frames / float(rate))*1000

    
mixer.pre_init (rate, -16, 2, 2048)
mixer.init ()

###########################################################################
# Configuration to detect button pressed
fd = sys.stdin.fileno()
newattr = termios.tcgetattr(fd)
newattr[3] = newattr[3] & ~termios.ICANON
newattr[3] = newattr[3] & ~termios.ECHO
termios.tcsetattr(fd, termios.TCSANOW, newattr)

oldterm = termios.tcgetattr(fd)
oldflags = fcntl.fcntl(fd, fcntl.F_GETFL)
fcntl.fcntl(fd, fcntl.F_SETFL, oldflags | os.O_NONBLOCK)
##########################################################################

count=0
rt=[]
btn=[]
stim=[]
    
for b in range(nblocks):
    
    for t in range(rep):
        if (int(idx[count])<0):
            sound_file = "/home/adrielle/Desktop/Piloto/Formantes/stim%g.wav" % int(idx[count])
        else:
            sound_file = "/home/adrielle/Desktop/Piloto/Formantes/stim+%g.wav" % int(idx[count])

        sound = mixer.Sound (sound_file)
        sound.play ()

        t0 = time.time()


        while True:
              inp, outp, err = select.select([sys.stdin], [], [])
              k = sys.stdin.read()
            
              tm=time.time () - t0
              
              if (tm>2.5):
                  button = ""
                  break
              elif k == 'a':
                  button = "pa"
                  break
              elif k == 'd':
                  button = "pe"
                  break
              else:
                  button = ""
                  break

        rt=rt+[tm]
        btn=btn+[button]
        stim=stim+[int(idx[count])]

        count=count+1
        
        #delay
        if (tm>1):
            time.sleep((900+(600*random.random()))/1000)
        else:
            time.sleep((1-tm)+((500+600*random.random())/1000))

    if (b<nblocks):
        raw_input("Pressione ENTER para continuar...\n")
        time.sleep((900+(400*random.random()))/1000)

all={'idx': idx, 'rt': rt, 'button': btn, 'stimuli': stim}

name="active_"+(file_name[0:(len(file_name)-3)])+"pickle"
myfile = open(name, 'wb')
pickle.dump(all, myfile)
myfile.close()

# Reset the terminal:
termios.tcsetattr(fd, termios.TCSAFLUSH, oldterm)
fcntl.fcntl(fd, fcntl.F_SETFL, oldflags)

#########################################################
# To load data, use the following commands

# myfile = open(name, 'rb')
# all = pickle.load(myfile)  # variables come out in the order you put them in
# myfile.close()
