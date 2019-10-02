#!/usr/bin/python
# -*- coding: utf-8 -*-

# coding=utf-8

import sys
reload(sys)
sys.setdefaultencoding('utf8')

### Identification of ba-pa continuum
import expyriment as exp
import numpy
import scipy as sp
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from random import shuffle, randint
import wave
import contextlib
import os
  
from expyriment import control, stimuli, io, design, misc

### Directories
prog_dir = os.path.dirname (os.path.realpath (__file__))
sound_dir = os.path.join ('/home','adrielle','Desktop','Piloto','Formantes')
exp.io.defaults.datafile_directory = os.path.join (prog_dir, 'data')
exp.io.defaults.eventfile_directory = os.path.join (prog_dir, 'events')

# settings
design.defaults.experiment_background_colour = misc.constants.C_BLACK

#Total number of trials
ntrials=200  

#Number of trials per Block
TpB=20

#Sequencies for jitter
seq1=numpy.linspace(0,500,10)
seq2=numpy.linspace(0,300,10)

#Number of blocks - choose ntrials and TpB so that nblocks can be round
nblocks=int(ntrials/TpB)

#Index of the stimuli
idx=numpy.linspace(1,200,ntrials)

for pos in range(len(idx)):
    if ((pos+1)%2==0):
        idx[pos]=-idx[pos]

shuffle(idx)
count=0

#Duration of the wav file

fname = sound_dir+"/stim-1.wav"
with contextlib.closing(wave.open(fname,'r')) as f:
    frames = f.getnframes()
    rate = f.getframerate()
    duration = (frames / float(rate))*1000

curve = exp.design.Experiment(name="Experimento")
exp.control.initialize(curve)

def bloco(nblc, sound):

   #Definição dos blocos
   bloco = exp.design.Block(name="Bloco %02d" % nblc)
   trial = exp.design.Trial()
   stim = exp.stimuli.TextLine(text="pa                pe", text_size=90, text_colour=(255,255,255), background_colour=(0,0,0))
   stim.preload()
   trial.add_stimulus(stim)

   stim = exp.stimuli.TextLine(text="pe                pa", text_size=90, text_colour=(255,255,255), background_colour=(0,0,0))
   stim.preload()
   trial.add_stimulus(stim)

   if (sound<0):
       sound= sound*(-1)
       stim = exp.stimuli.Audio(sound_dir+"/stim-%g.wav" % sound)
       stim.preload()
       trial.add_stimulus(stim)
   else:
       stim = exp.stimuli.Audio(sound_dir+"/stim+%g.wav" % sound)
       stim.preload()
       trial.add_stimulus(stim)

   stim = exp.stimuli.TextLine(text="pressione barra de espaço para continuar", text_size=40, text_colour=(255,255,255), background_colour=(0,0,0))
   stim.preload()
   trial.add_stimulus(stim)

   bloco.add_trial(trial)
   curve.add_block(bloco)

try:
        cruz=exp.stimuli.FixCross(size=(100, 100), line_width=3, colour=(255,255,255))
        cruz.preload()

        exp.control.start(skip_ready_screen=True)

        curve.data_variable_names = ['block.name','stim', 'tecla', 'rt','st']
        for blc in range(nblocks):
           
            for it in range(TpB):
               
                #Jitter
                jitter1=shuffle(seq1)
                jitter1=round(seq1[0])
                jitter2=shuffle(seq2)
                jitter2=round(seq2[0])


                bloco(blc, int(round(idx[count])))
                cruz.present()
                curve.clock.wait(1000+jitter1)

                #Toca o áudio
                curve.blocks[blc].trials[0].stimuli[2].present()

                #Aguarda um tempo antes da tela de resposta
                curve.clock.wait(duration+jitter2)

                #Randomize the query window
                st=randint(0, 1)
                curve.blocks[blc].trials[0].stimuli[st].present()
                
                #Tempo para a apresentação do áudio
                #Pessoa tem 2.5s para pressionar a tecla (esse tempo inclui o tempo do áudio aqui) #K_RCTRL
                tecla, rt = curve.keyboard.wait(keys=[exp.misc.constants.K_LCTRL, exp.misc.constants.K_LALT, exp.misc.constants.K_q], duration=2500)
                #Desmarque o comando a seguir se desejar que o programa espere até o fim da duração estabelecida em duration para seguir
                #curve.clock.reset_stopwatch()
                #curve.clock.wait(3000 - teste.clock.stopwatch_time)

                #Guarda os dados num arquivo
                curve.data.add([curve.blocks[blc].name, int(round(idx[count])), tecla, rt, st])
                count=count+1

                #key to exit at any moment
                if (tecla==113):
                  exp.control.end(goodbye_text="Saindo do experimento!", goodbye_delay=1000)
                  
                #Espera pressionar barra de espaço para novo bloco
                if (it == TpB-1) and (blc != (nblocks-1)):
                        curve.blocks[blc].trials[0].stimuli[3].present()
                        curve.keyboard.wait(exp.misc.constants.K_SPACE)
                
                if (it != TpB-1):        
                        curve.remove_block(blc)
                        
        exp.control.end(goodbye_text="Fim do experimento. Obrigada!", goodbye_delay=1000)

except Exception as e:
        print (str (e))
        exp.control.end(goodbye_text="Deu Ruim! :-(")

### Plot the psicometric curve

ro.r ('source ("plot-curve.r")')
ro.r ('plot.curve ("data/' + curve.data.filename +'")')
ro.r('dev.off()')
os.system ('evince data/' + curve.data.filename + ".pdf")
