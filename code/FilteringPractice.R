freqOfInterest = 133 # Hz

per1 = c(rep(0, 400 -2), 1, -1) 
y = rep(per1, freqOfInterest)
timesteps = seq(0, 1, length.out = length(y))
sampFreq = length(y) / max(timesteps)
sampFreq
dev.off()

plot(y, x = timesteps, type = 'l')


# filtSize should be wider than # points between peaks of interest
# but not wider than about 4-5 X the period 
# in this case, the period is 500 pts
# so filter should be between 500 and 2000
# I'll shoot for 1000 pts


periodForFOI = sampFreq / freqOfInterest

filtSize = 3 * periodForFOI


filttt =  dnorm(x = seq(-4, 4, length.out =  filtSize))

plot(y, x = timesteps, type = 'l')
lines(y = stats::filter(x=y, filt = filttt, circular = TRUE), col = 'red', 
      x = timesteps)



spectro(y, f = sampFreq, osc = TRUE)
dev.off()
seewave::spec(y, f = 1000)

filteredData = stats::filter(x=y, filt = filttt, circular = TRUE)
filteredData = filteredData - mean(filteredData)
plot(y = filteredData, col = 'red', 
      x = timesteps, type = 'l')


dev.off()
spp = seewave::spec(filteredData, f = sampFreq, flim = c(0, 1), wl = 2**12)

# returns frequency info
spp[,1][which.max(spp[,2])] * 1000

spectro(filteredData, f = sampFreq, osc = TRUE, scale = FALSE)
            


spikySignalSmoother = function(spikySignal = y, sampFreq = sampFreq, freqOfInterest = 3){
  # Makes a smooth waveform with the same frequency as
  # the spiky signal
  #
  # Args:
  #   spikySignal: Spiky signal to be smoothed
  #   sampFreq: The sampling frequency of the spiky signal
  #   freqOfInterest: A rough guess at the minimum frequency you're interested in resolving.
  #     If you guess too below the frequency of interest, then you won't get a smooth waveform. 
  #     If you guess too high, then you also won't get a smooth waveform.
  #
  # Returns:
  #   A smooth waveform with the same length as the original spiky signal
  
  periodForFOI = sampFreq / freqOfInterest
  filtSize = 4 * periodForFOI
  
  # defineFilter
  gaussianFilter =  dnorm(x = seq(-4, 4, length.out =  filtSize))
  filteredData = stats::filter(x=(spikySignal), filt = gaussianFilter, circular = TRUE)
  
  # rescale
  scaledFilteredData <- (filteredData-min(filteredData))/(max(filteredData)-min(filteredData)) * 2 - 1
  
  # warning if returned waveform is not smooth
  deriv <- diff(scaledFilteredData)
  if(max(abs(deriv)) > 0.5){
    warning("Returned wave may not be smooth. Try increasing freqOfInterest to be higher than signal frequency")
  }
  else if(sum(deriv == 0) / length(deriv) > 0.1 ){
    warning("Returned wave may not be smooth. Try decreasing freqOfInterest")
  }
  
  return(scaledFilteredData)
}

dev.off()
plot(y, x= timesteps, type = "l", xlim = c(0, 0.1))

filtY = spikySignalSmoother(y, sampFreq, 200)

lines(x = timesteps, y = filtY, col = 'red')

plot(x = timesteps, y = filtY, col = 'red', type = 'l', xlim = c(0, 0.1))

spectro(filtY, f = sampFreq, osc = TRUE, scale = FALSE)
dev.off()
spp = seewave::spec(filtY, f = sampFreq, flim = c(0, 1), wl = 2**12)

# returns frequency info
spp[,1][which.max(spp[,2])] * 1000




spikySignalSmoother_rollingMean = function(spikySignal = y, sampFreq = sampFreq, freqOfInterest = 3){
  # Makes a smooth waveform with the same frequency as
  # the spiky signal
  #
  # Args:
  #   spikySignal: Spiky signal to be smoothed
  #   sampFreq: The sampling frequency of the spiky signal
  #   freqOfInterest: A rough guess at the minimum frequency you're interested in resolving.
  #     If you guess too below the frequency of interest, then you won't get a smooth waveform. 
  #     If you guess too high, then you also won't get a smooth waveform.
  #
  # Returns:
  #   A smooth waveform with the same length as the original spiky signal
  
  periodForFOI = sampFreq / freqOfInterest
  filtSize = 4 * periodForFOI
  
  # defineFilter
  b2 <- zoo::rollmean(abs(spikySignal), k = periodForFOI )
  filteredData = stats::filter(x=(spikySignal), filt = gaussianFilter, circular = TRUE)
  
  # rescale
  scaledFilteredData <- (filteredData-min(filteredData))/(max(filteredData)-min(filteredData)) * 2 - 1
  
  return(scaledFilteredData)
}

