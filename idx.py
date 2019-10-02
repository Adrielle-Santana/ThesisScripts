#!/usr/bin/python
# -*- coding: utf-8 -*-

# coding=utf-8
# converts python pickle file from the  identification task for /pa/ and /pe/ stimuli continuum in a MATLAB mat file

import pickle
import scipy.io as sio


file=open("passive_felipe.pickle",'rb')
idx=pickle.load(file)
file.close()

adict={}
adict['idx']=idx
sio.savemat('idx_passive.mat', adict)

#######################################################################

file=open("active_felipe.pickle",'rb')
idx=pickle.load(file)
file.close()

adict={}
adict['idx']=idx['idx']
sio.savemat('idx_active.mat', adict)

adict={}
adict['rt']=idx['rt']
sio.savemat('rt.mat', adict)
