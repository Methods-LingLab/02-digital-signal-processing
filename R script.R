library(magrittr)

# download function
dl_from_dropbox <- function(x, key) {
  bin <- RCurl::getBinaryURL(paste0("https://dl.dropboxusercontent.com/s/", key, "/", x),
                             ssl.verifypeer = FALSE)
  con <- file(x, open = "wb")
  writeBin(bin, con)
  close(con)
}

# set audio player
sound::setWavPlayer('/bin/aplay')


# load first data set
dl_from_dropbox('EMA-EMG_data.Rda', 'gvpe7yr7rugx3cv')
load('EMA-EMG_data.Rda')

# play audio
audio <- sound::as.Sample(EMA.EMG$audio, EMA.EMG$SR$audio)
sound::play(audio)

# get sampling rates
aud.sr <- EMA.EMG$SR$audio
ema.sr <- EMA.EMG$SR$EMA

# get time points for word
t1 <- round(EMA.EMG$segments$start[16] * aud.sr)
t2 <- round(EMA.EMG$segments$end[20] * aud.sr)

t1; t2


# plot word signal
word <- EMA.EMG$segments$word[16]

plot(EMA.EMG$audio[t1:t2], type='l', 
     main=paste("Audio signal for word:",word))


# get time points for all segments of word
t1 <- EMA.EMG$segments$start[16]
t2 <- EMA.EMG$segments$start[17]
t3 <- EMA.EMG$segments$start[18]
t4 <- EMA.EMG$segments$start[19]
t5 <- EMA.EMG$segments$start[20]
t6 <- EMA.EMG$segments$end[20]


# get samples 
aud.seg <- round(t1 * aud.sr):round(t6 * aud.sr)
ema.seg <- round(t1 * ema.sr):round(t6 * ema.sr)

length(aud.seg); length(ema.seg)


# example of numeric sequence
seq(t1, t6, length.out=15)


# get data and time series
audio       <- data.frame(data=numeric(length(aud.seg)), time=numeric(length(aud.seg)))
audio$data  <- EMA.EMG$audio[aud.seg]
audio$time  <- seq(t1, t6, length.out=length(audio$data))

jawrot      <- data.frame(data=numeric(length(ema.seg)), time=numeric(length(ema.seg)))
jawrot$data <- EMA.EMG$EMA$jawrot[ema.seg]
jawrot$time <- seq(t1, t6, length.out=length(jawrot$data))


# function for plotting vertical line
vline <- function(x = 0, color = 'black') {
  list(
    type = 'line', y0 = 0, y1 = 1, yref = 'paper', x0 = x, x1 = x, 
    line = list(color = color, dash = 'dash')
  )
}

# create audio and jaw rotation sub-plots
p1 <- plotly::plot_ly(audio, x = ~time, y = ~data, type='contour') %>% 
  plotly::add_lines(name = "audio") %>% 
  plotly::layout(shapes = list(vline(t1),vline(t2),vline(t3),vline(t4),vline(t5),vline(t6)),
                 yaxis = list(zeroline=F, title=NA), 
                 xaxis = list(title = 'Time (s)'), 
                 title = word)

p2 <- plotly::plot_ly(jawrot,x=~time,y=~data,type='contour') %>% 
  plotly::add_lines(name = "rotation") %>% 
  plotly::layout(shapes = list(vline(t1),vline(t2),vline(t3),vline(t4),vline(t5),vline(t6)),
                 yaxis = list(zeroline=F, title=NA), 
                 xaxis = list(title = 'Time (s)'))

# create plot
plotly::subplot(p1, p2, nrows=2, shareX=T, shareY=T)



# load second data set
dl_from_dropbox('MRI_data.Rda', 'yf27qm0x5okp1sy')
load('MRI_data.Rda')


# look at spectrogram of noisy audio
sr <- MRI.data$SR$audio

phonTools::spectrogram(MRI.data$data3$audio, fs=sr)

# adjust spectrogram
phonTools::spectrogram(MRI.data$data3$audio, fs=sr, colors=F, dynamicrange=65, maxfreq=8000)

# listen to noisy audio
audio <- sound::as.Sample(MRI.data$data3$audio, sr)
sound::play(audio)


# FFT of noise
fftlen <- 1024

noise <- MRI.data$data3$audio[1:fftlen]

plot(noise, type='l')

fft0 <- fft(noise)
mag0 <- Mod(fft0)


# plot spectrum of noise
plot(mag0[1:(fftlen/2)], type='l')

hertz <- seq(1, sr/2, length.out=fftlen/2)

plot(hertz, mag0[1:(fftlen/2)], type='l', xlab="Frequency (Hz)", ylab="Amplitude")


# grab a section of audio to demonstrate spectral subtraction
chunk <- MRI.data$data3$audio[7580:8603]

phonTools::spectrogram(chunk, fs=sr, color=F, dynamicrange=65, maxfreq=sr/2)

fft1 <- fft(chunk)
mag1 <- Mod(fft1)

plot(hertz, mag1[1:(fftlen/2)], type='l')


# compare noisy audio and baseline noise
plot(hertz, mag1[1:(fftlen/2)], type='l', col='red')
lines(hertz, mag0[1:(fftlen/2)], col='blue')


# spectral subtraction
mag2 <- mag1 - mag0

plot(hertz, mag2[1:(fftlen/2)], type='l')
abline(a=0, b=0, col='red')

mag2[mag2<0] <- 0

plot(hertz, mag2[1:(fftlen/2)], type='l')


# inverse FFT
phs1  <- Arg(fft1)
j     <- sqrt(as.complex(-1))

dn        <- mag2 * exp(j*phs1)
denoised  <- Re(fft(dn,inverse=T)/fftlen)

# look at cleaned IFFT waveform
phonTools::spectrogram(denoised, fs=sr, color=F, dynamicrange=65, maxfreq=sr/2)


# let's do it for the whole audio signal!
N <- length(MRI.data$data3$audio)

cleaned <- numeric(length=N)

# overlapping moving window
overlap <- fftlen/2

# Hann window/kernel
hannwin <- ( 0.5 - (0.5 * cos(2*pi*(0:(fftlen-1))/(fftlen-1))) )
hannwin <- phonTools::windowfunc(fftlen, type="hann")

plot(hannwin, type='l')

# half-Hann kernels
hann1 <- hannwin[(1 + fftlen/2):fftlen]
hann2 <- hannwin[1:(fftlen/2)]

# additive property of two half-Hann kernels
plot(hann1, type='l', col='red')
lines(hann2, col='blue')
lines(hann1 + hann2, col='purple')

chunk <- MRI.data$data3$audio[1:512]

plot(hann1*chunk, type='l', col='red')
lines(hann2*chunk, col='blue')


# how many overlapping windows will fit?
bins <- floor(N/(fftlen-overlap))


# moving window over whole audio signal
for (bin in 1:bins) {
  chunk <- MRI.data$data3$audio[(1+(bin-1)*(fftlen-overlap)):(fftlen+(bin-1)*(fftlen-overlap))]
  
  fft1 <- fft(chunk)
  mag1 <- Mod(fft1)
  phs1 <- Arg(fft1)
  
  mag2 <- mag1 - mag0
  
  #mag2[mag2<0] <- 0
  
  dn        <- mag2 * exp(j*phs1)
  denoised  <- Re(fft(dn,inverse=T)/fftlen)
  
  overlap1  <- cleaned[(1+(bin-1)*(fftlen-overlap)):(overlap+(bin-1)*(fftlen-overlap))]
  
  overlap2  <- denoised[1:overlap]
  
  cleaned[(1+(bin-1)*(fftlen-overlap)):(overlap+(bin-1)*(fftlen-overlap))] <- overlap1*hann1 + overlap2*hann2
  
  cleaned[(1+overlap+(bin-1)*(fftlen-overlap)):(fftlen+(bin-1)*(fftlen-overlap))] <- denoised[(1+overlap):fftlen]
}

cleaned[is.na(cleaned)] <- 0


# view spectrogram and listen to cleaned audio
phonTools::spectrogram(cleaned, fs=sr, color=F, dynamicrange=65, maxfreq=sr/2)

audio <- sound::as.Sample(cleaned, sr)
sound::play(audio)





## DAY 2
dl_from_dropbox('MRI_data.Rda', 'yf27qm0x5okp1sy')
load('MRI_data.Rda')


# inspect MRI and audio data for the first word
MRI.data$data1$word

dim(MRI.data$data1$images)

sr    <- MRI.data$SR$audio
audio <- sound::as.Sample(MRI.data$data1$audio, sr)
sound::play(audio)


# inspect MRI and audio data for the second word
MRI.data$data2$word

dim(MRI.data$data2$images)

audio <- sound::as.Sample(MRI.data$data2$audio, sr)
sound::play(audio)


# view MR image
grayscale <- gray(seq(0, 1, length.out=256))

image(MRI.data$data1$images[,,1], col=grayscale)


# images of velum up and velum down
image1 <- MRI.data$data1$images[,,1]
image2 <- MRI.data$data1$images[,,19]

image(image1, col=grayscale)
image(image2, col=grayscale)

# difference image
image(image1-image2, col=grayscale)

# make an ROI around VP port
coords <- locator()

# create polygon object from ROI coordinates
polygon(as.vector(coords$x), as.vector(coords$y), border='red', lwd=2)

poly <- sp::Polygon(coords)

# create mask
poly.mask <- sp::SpatialPolygons(list(sp::Polygons(list(poly),"mask")))


# view (incorrect) masked image
image(
  raster::as.matrix(
    raster::mask(
      raster::raster(MRI.data$data1$images[,,1]), poly.mask))
  , col=grayscale)

# flip the mask around
masked <- 
  raster::as.matrix(
    raster::mask(
      raster::flip(
        raster::raster(
          t(MRI.data$data1$images[,,1])), 2), poly.mask))

masked <- MRI.data$data1$images[,,1] %>% 
  t() %>%
  raster::raster() %>% 
  raster::flip(2) %>% 
  raster::mask(poly.mask) %>%
  raster::as.matrix()


# view (correct) masked image
image(
  raster::as.matrix(
    raster::flip(
      raster::raster(
        t(masked)), 1))
  , col=grayscale)



# get number of images in first word
imgnum <- dim(MRI.data$data1$images)[3]

# pre-allocate vector for ROI intensity measure
mr.dat <- vector(length=imgnum)

# get average intensity in ROI for each image
for (img in 1:imgnum) {
  masked <- raster::as.matrix(
    raster::mask(
      raster::flip(
        raster::raster(
          t(MRI.data$data1$images[,,img])), 2), poly.mask))
  
  mr.dat[img] <- mean(masked, na.rm=T)
}


# plot the intensity signal
plot(mr.dat, type='l')


# add time points of the start and end of target vowel
abline(v=round(MRI.data$data1$vstart*MRI.data$SR$MRI), col="blue")
abline(v=round(MRI.data$data1$vend*MRI.data$SR$MRI), col="red")



# similar thing, but now with PCA!
image1 <- MRI.data$data1$images[,,1]
image2 <- MRI.data$data1$images[,,19]

image(image1-image2, col=grayscale)


# other options to look at total image variance...
avgimg <- apply(MRI.data$data1$images, 1:2, mean)
image(avgimg, col=grayscale)

varimg <- apply(MRI.data$data1$images, 1:2, var)
image(varimg, col=grayscale)

# ...or variance over time
dists <- apply(MRI.data$data1$images, 3, function(x) sqrt(sum((x - avgimg)^2)))
plot(dists, type='l')


# larger ROI over whole region of velum movement
image(image1-image2, col=grayscale)

coords <- locator()

polygon(as.vector(coords$x), as.vector(coords$y), border='red', lwd=2)


# combine images for one PCA model
MRI.data$combined <- c()

MRI.data$combined$images <- array(c(MRI.data$data1$images, MRI.data$data2$images), 
                                  dim=c(dim(MRI.data$data1$images)[1:2], 
                                        sum(c(dim(MRI.data$data1$images)[3],
                                              dim(MRI.data$data2$images)[3]))
                                  ))


imgnum <- dim(MRI.data$combined$images)[3]

# check that the number of images makes sense
imgnum

dim(MRI.data$data1$images)[3]; dim(MRI.data$data2$images)[3]



# create polygon mask rom ROI coordinates
poly <- sp::Polygon(coords)
poly.mask <- sp::SpatialPolygons(list(sp::Polygons(list(poly), "mask")))


# how many voxels are inside the ROI?
test <- raster::as.matrix(
  raster::mask(
    raster::flip(
      raster::raster(
        t(MRI.data$combined$images[,,1])), 2), poly.mask))

voxels <- sum(!is.na(test))

voxels


# preallocate matrix for PCA
mr.dat <- matrix(data=NA,
                 nrow=imgnum,
                 ncol=voxels)

# create a vector of the masked pixels for each image
for (img in 1:imgnum) {
  masked <- raster::as.matrix(
    raster::mask(
      raster::flip(
        raster::raster(
          t(MRI.data$combined$images[,,img])), 2), poly.mask))
  
  mr.dat[img,] <- as.vector(masked[!is.na(masked)])
}

# too many dimensions!
dim(mr.dat)


# we need to down-sample the image resolution...
imfact  <- 10

pca.dat <- matrix(nrow=imgnum, 
                  ncol=round(voxels/imfact))

dim(pca.dat)


# down-sample the images, and vectorize the masked data
for (frame in 1:imgnum){
  this.frame  <- mr.dat[frame,]
  new.frame   <- signal::resample(this.frame, ncol(pca.dat)/voxels)
  pca.dat[frame, ] <- as.vector(new.frame)
}


# PCA model
pca.mod <- prcomp(pca.dat)

summary(pca.mod)


# plot PC1 scores
plot(pca.mod$x[,1], type='l')


# make velum movement signal from PC1 scores, scale
MRI.data$combined$velum <- (pca.mod$x[,1]-min(pca.mod$x[,1])) / (max(pca.mod$x[,1])-min(pca.mod$x[,1]))

plot(MRI.data$combined$velum, type='l')


# create frame offset for the two sets of word images
offset <- (dim(MRI.data$data1$images)[3]+1) / MRI.data$SR$MRI

# create time points including offsetted times for second word
MRI.data$combined$vstarts <- c(MRI.data$data1$vstart, MRI.data$data2$vstart + offset)
MRI.data$combined$vends   <- c(MRI.data$data1$vend, MRI.data$data2$vend + offset)

# plot velum movement signal
plot(MRI.data$combined$velum, type='l')

# add target vowel time points
abline(v=round(offset*MRI.data$SR$MRI), lty=2)
abline(v=round(MRI.data$combined$vstarts * MRI.data$SR$MRI), col='blue')
abline(v=round(MRI.data$combined$vends * MRI.data$SR$MRI), col='red')