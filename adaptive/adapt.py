#------------------------------------------------------------------------------
# IIR filter test program.
#
# Harri Salokorpi, 2004
#------------------------------------------------------------------------------

import math
import wave
import array

#------------------------------------------------------------------------------
# Class definitions
#------------------------------------------------------------------------------


# Test frequency class
class Freq:
    # Constructor
    def __init__(self, freq_):
        self.freq = freq_
        self.maxval = 0
        return

    # This function seeks maximal value
    def collect(self, sample):
        self.maxval = max(abs(sample), self.maxval)
        return


# IIR SecondOrderSection class
class SecondOrderSection:
    # Constructor
    def __init__(self, a1_, a2_, b0_, b1_, b2_):
        self.a1, self.a2 = a1_, a2_
        self.b0, self.b1, self.b2 = b0_, b1_, b2_
        self.buf = [0, 0, 0]
        return

    # Filter sample (givin bits-parameter enables the integer mode)
    def filter(self, sample, bits=None):
        # Feedback
        if bits == None:
            # Floatpoint mode
            sample = sample + self.a1 * self.buf[1]
            sample = sample + self.a2 * self.buf[2]
        else:
            # Integer mode
            sample = sample + (self.a1 * self.buf[1]) / (2**(bits - 1))
            sample = sample + (self.a2 * self.buf[2]) / (2**(bits - 1))

        self.buf[0] = sample

        # Forward
        if bits == None:
            # Floatpoint mode
            sample = self.b0 * self.buf[0]
            sample = sample + self.b1 * self.buf[1]
            sample = sample + self.b2 * self.buf[2]
        else:
            # Integer mode
            sample = (self.b0 * self.buf[0]) / (2**(bits - 1))
            sample = sample + (self.b1 * self.buf[1]) / (2**(bits - 1))
            sample = sample + (self.b2 * self.buf[2]) / (2**(bits - 1))

        # Shift data forward
        self.buf[1], self.buf[2] = self.buf[0], self.buf[1]

        return sample


# IIR sine generator
class IIRSineGenerator:
    # Constructor
    def __init__(self, samplerate_, freq_=1000, gain_=1.2):
        self.samplerate = samplerate_

        # Create IIR filter that generates sine
        w = 2 * math.pi * freq_ / samplerate_
        a1 = 2 * math.cos(w)
        a2 = -1
        b0 = 0
        b1 = gain_ * math.sin(w)
        b2 = 0
        self.iir = SecondOrderSection(a1, a2, b0, b1, b2)
        return

    # This function generates sine sample
    def generate(self):
        return self.iir.filter(0.125) - 0.38


# Sine generator
class SineGenerator:
    # Constructor
    def __init__(self, samplerate_):
        self.samplerate = samplerate_
        self.angle = 0.0
        pass

    # This function generates sine sample
    def generate(self, freq):
        # Generate new angle
        anglechange = 2.0 * math.pi * \
         ( float( freq.freq ) / float( self.samplerate ) )
        self.angle = self.angle + anglechange

        # Generate sine
        sample = math.sin(self.angle)
        return sample


# Adaptive filter class
class AdaptiveFilter:
    # Constructor
    def __init__(self, blocksize, roc):
        self.blocksize = blocksize
        self.roc = roc
        self.coeffs = [0 for x in range(blocksize)]
        self.nbuf = [0 for x in range(2 * blocksize)]
        return

    # Filter noise from signal with noise (swn)
    # Note: inputs are lists with equal amount of items
    def filter(self, swn, noise):
        outbuf = []
        # Check sequence sizes
        if len(swn) != len(noise) and len(swn) != self.blocksize:
            raise Exception, "Invalid incoming data blocksize"

        # Move noise data forward (our buffer is 2x blocksize)
        for i in range(self.blocksize):
            self.nbuf[i] = self.nbuf[i + self.blocksize]
            self.nbuf[i + self.blocksize] = noise[i]

        # Process block
        for i in range(self.blocksize):
            n = 0
            # Compute adaptive filter output
            for j in range(self.blocksize):
                n += self.coeffs[j] * self.nbuf[self.blocksize + i - j]

            # Compute the error
            e = swn[i] - n

            outbuf.append(e)

            # Update coeffs
            for j in range(self.blocksize):
                self.coeffs[j] = self.coeffs[j] + \
                 2 * self.roc * e * self.nbuf[self.blocksize + i - j ]

        return outbuf


#------------------------------------------------------------------------------
# Data tables and general settings
#------------------------------------------------------------------------------

# Sample rate
samplerate = 16000.0

# How many samples are run through the filter / test frequency band
scount = 5000

# roc (rate of convergence)
roc = 0.01

# Blocksize
blocksize = 32

# Bit count for integer filtering
bits = 16

# Generate test frequency table (from 100Hz to 3600 Hz in 100 Hz steps)
freqtable = map(Freq, range(100, 8000, 100))

# Create sine generator
sinegen = SineGenerator(samplerate)
iirsinegen = IIRSineGenerator(samplerate)

# Adaptive filter
af = AdaptiveFilter(blocksize, roc)

#------------------------------------------------------------------------------
# Main functions
#------------------------------------------------------------------------------


# Main function
def main():

    # ---------- Floating point -------------

    waveout = wave.open("float.wav", "w")
    waveout.setnchannels(1)
    waveout.setsampwidth(2)
    waveout.setframerate(samplerate)
    waveout.setcomptype("NONE", "Uncompressed")
    sbuf = array.array('h')

    # Initialize noise and signal+noise buffers
    nbuf = [0 for x in range(blocksize)]
    nsbuf = [0 for x in range(blocksize)]

    # Test all frequency bands (float)
    for freq in freqtable:
        # Run scount samples with same frequence through filter
        for i in range(int(scount / blocksize)):

            # Generate signal + noise
            for j in range(blocksize):
                signal = iirsinegen.generate()
                noise = sinegen.generate(freq)
                nbuf[j] = noise / 2.0
                nsbuf[j] = (noise + signal) / 2.0

            # Run block through adaptive filter
            out = af.filter(nsbuf, nbuf)

            # Append data output buffer
            for j in range(blocksize):
                sbuf.append(int(out[j] * (2**(bits - 1))))

    waveout.writeframes(sbuf.tostring())
    waveout.close()
    print("Wrote file float.wav")
    return


main()
