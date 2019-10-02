### Sampling frequency of EEG (in Hz)
fs.eeg <- 5000

### Indicates for which subject you want to compute regression
subjects <- c (1,3,4,5,6,7,9,10,11)

nb.subjects <- length(subjects)

### Number of points for the DWT
nb.points <- 2048

### Number of different stims (considering pos and neg)
nstim <- 10

### Amount of groups by stimulus (normal or inverted)
groups <- 7 ###sqrt(175) quant. de trials por estímulo

