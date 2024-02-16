# https://github.com/tanzirhasan/ECG_Simulation

# main distribution functions -----

# ECG (Electrocardiogram) measurements are fundamental in diagnosing and monitoring heart health, providing insights into the heart's electrical activity through various waveforms. Each component of the ECG waveform represents a different phase of the cardiac cycle. Here's a brief description of the values used to define ECG measurements:
#   
#   - **P Wave**: Represents atrial depolarization, which is the electrical activity that leads to the contraction of the atria. It's the first small upward deflection on the ECG.
# 
# - **Q Wave**: The first negative (downward) deflection after the P wave, but preceding the R wave. It represents the initial negative deflection produced by ventricular depolarization.
# 
# - **R Wave**: The first upward deflection following the P wave, representing the rapid depolarization of the right and left ventricles. The R wave is usually the most prominent part of the ECG trace.
# 
# - **S Wave**: The first downward deflection immediately following the R wave, representing the final depolarization of the ventricles, specifically at the base of the heart.
# 
# - **T Wave**: Represents ventricular repolarization, which is when the ventricles recover from the state of depolarization. The T wave is usually a modest upward waveform following the S wave.
# 
# - **U Wave**: A small wave that may follow the T wave but is not always present. Its origin is not completely understood, but it may represent the repolarization of the Purkinje fibers.
# 
# - **RR Interval**: The time between two successive R waves. It's a direct measure of the heart rate. The RR interval is used to calculate the heart rate by determining the time between beats. It is especially useful in detecting heart disorders for several reasons:
#   - **Arrhythmia Detection**: Variability in the RR interval can indicate irregular heartbeats or arrhythmias.
# - **Heart Rate Variability (HRV)**: Analysis of the variations in the RR intervals over time can provide insights into autonomic nervous system function and is linked to various cardiac and non-cardiac diseases.
# - **Cardiac Health**: Changes in the RR interval length or pattern can indicate stress on the heart or potential cardiac disorders, making it a valuable tool for monitoring overall heart health and function.
# 
# By analyzing these components, clinicians can assess the heart's electrical activity, diagnose various heart conditions such as arrhythmias, blockages, and other disorders, and monitor the effectiveness of treatments.

library(ggplot2)


ecg_wave <- function(x, c){
  
  p_wave <- function(x,c) {
    l=c# half heartbeat cycle.eg.s/heartbeat=1 cycle
    a=0.25#amplitude of p wave
    x=x+(l/1.8)# x is the wave starting point shifted.
    b=3# 2l/b is duration of p wave.
    n=100 # fourier series levels, the bigher the more accurate
    p1=0 # baseline
    p2=0 # p wave
    # fourier series to creat p wave
    for (i in 1:n){
      harm1<-(((sin((pi/(2*b))*(b-(2*i))))/(b- (2*i))+(sin((pi/(2*b))*(b+(2*i))))
               /(b+(2*i)))*(2/pi))*cos((i*pi*x)/l)
      p2<-p2+harm1
    }
    pwav1=p1+p2
    pwav=a*pwav1
    return(pwav)
  }
  
  ####### q wave
  q_wave<-function (x,c){
    l=c
    x=x+l/6
    a=0.03# amplitude 
    b=16#duration
    n=100
    q1=0 #(a/(2*b))*(2-b)
    q2=0
    for (i in 1:n){
      harm5=(((2*b*a)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l);
      q2=q2+harm5}
    qwav=-1*(q1+q2)
    return (qwav)
  }
  
  
  #### qrs wave
  qrs_wave<-function(x,c){
    l=c
    a=1
    b=5
    n=100
    qrs1=0 #(a/(2*b))*(2-b)
    qrs2=0
    for (i in 1:n){
      harm=(((2*b*a)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l);
      qrs2=qrs2+harm}
    qrswav=qrs1+qrs2
    return(qrswav)
  }
  
  ##### s wave
  s_wave<-function(x,c){
    l=c
    x=x-l/6
    a=0.25
    b=15
    n=100;
    s1=0 # (a/(2*b))*(2-b);
    s2=0;
    for (i in 1:n){
      harm3=(((2*b*a)/(i*i*pi*pi))*(1-cos((i*pi)/b)))*cos((i*pi*x)/l);
      s2=s2+harm3}
    swav=-1*(s1+s2)
    return(swav)
  }
  
  
  ##### t wave
  t_wave<-function(x, c){
    l=c
    a=0.35
    x=x-l/1.8
    b=7
    n=100
    t1=0 # 1/l
    t2=0
    for (i in 1:n){
      harm2 = (((sin((pi/(2*b))*(b-(2*i))))/(b-(2*i))+(sin((pi/(2*b))*(b+(2*i))))/(b+(2*i)))*(2/pi))*cos((i*pi*x)/l);             
      t2 = t2+harm2
    }
    twav1 = t1+t2
    twav = a*twav1
    return (twav)
  }
  
  p_wave(x, c) + qrs_wave(x, c) + t_wave(x, c) + s_wave(x, c) + q_wave(x, c)
}


fit_ecg_waves_to_beats <- function(Q, sampling_rate=125){
  
  PP <- diff(c(0, Q))  # re-create P as "PP": might just want to use P
  
  x <- seq(0, sum(PP), by=1/sampling_rate)
  x_start = 0
  
  ecg <- numeric(length(x))
  
  for (i in 1:length(PP)) {
    x_end = x_start + PP[i]
    if (i == length(PP)) {
      x_used = subset(x, x>=x_start)
    } else {
      x_used = subset(x, x>=x_start & x<x_end)
    }
    
    next_wave <- ecg_wave(x_used - x_start, PP[i]/2)
    
    # ecg <- c(ecg, next_wave)
    ecg[x %in% x_used] <- next_wave
    
    x_start <- x_end
  }
  
  return(ecg)
}



# variability functions ----
# Function to simulate RR intervals with variability
simulate_rr_intervals <- function(T, baseline_rr, variability) {
  # T is the total number of heartbeats
  # baseline_rr is the baseline RR interval (in seconds)
  # variability is the maximum deviation from the baseline_rr
  rr_intervals <- baseline_rr + runif(T, -variability, variability)
  return(rr_intervals)
}

# Function to convert RR intervals to cumsum of beats for ECG simulation
rr_intervals_to_beats <- function(rr_intervals) {
  # Convert RR intervals from seconds to number of beats
  beats <- cumsum(rr_intervals * SAMPLING_RATE) # Assuming SAMPLING_RATE is in Hz
  return(beats)
}

# simulate groups ----

# Function to simulate RR intervals for a given group
simulate_rr_intervals_group <- function(N, baseline_rr, variability, group_label) {
  rr_intervals <- baseline_rr + runif(N, -variability, variability)
  group <- rep(group_label, N)
  data.frame(RR_Interval = rr_intervals, Group = group)
}
