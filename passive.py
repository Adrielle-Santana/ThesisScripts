#!/usr/bin/python
# -*- coding: utf-8 -*-

# coding=utf-8
# run identification task for /pa/ and /pe/ stimuli continuum

import easygui
import sys, os
reload(sys)
sys.setdefaultencoding('utf8')
import numpy
import scipy as sp
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import random
from random import shuffle, randint
import wave
import contextlib
import time
from pygame import mixer
import pickle

### Enter file name
file_name=easygui.enterbox('Enter dat file name for the subject: ')

### Read file and identificate middle continuum stimuli
ro.r('load(file=file.path("Identification/data/'+file_name+'"))')
pa_amb=ro.r("pa.amb")[0]
amb=ro.r("amb")[0]
pe_amb=ro.r("pe.amb")[0]

#Amount of repetitions of each one of the stimuli selected      
rep=350

#Index of the stimuli
idx=[1,pa_amb,amb,pe_amb,200]*rep

#Total number of trials
ntrials=len(idx)

for pos in range(len(idx)):
    if ((pos+1)%2==0):
        idx[pos]=-idx[pos]

shuffle(idx)

###########################################
name="passive_"+(file_name[0:(len(file_name)-3)])+"pickle"
myfile = open(name, 'wb')
pickle.dump(idx, myfile)
myfile.close()
###########################################

#Duration of the wav file. Use any wav file once their duration are the same

fname = "/home/adrielle/Desktop/Piloto/Formantes/stim-1.wav"
with contextlib.closing(wave.open(fname,'r')) as f:
    frames = f.getnframes()
    rate = f.getframerate()
    duration = (frames / float(rate))*1000

    
mixer.pre_init (rate, -16, 2, 4096)
mixer.init ()

count=0

for t in range(ntrials):
    if (int(idx[count])<0):
       sound_file = "/home/adrielle/Desktop/Piloto/Formantes/stim%g.wav" % int(idx[count])
    else:
       sound_file = "/home/adrielle/Desktop/Piloto/Formantes/stim+%g.wav" % int(idx[count])

    sound = mixer.Sound (sound_file)
    sound.play ()

    count=count+1

    time.sleep(1.0+(300*random.random())/1000)
