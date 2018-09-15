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
        anglechange = 2.0 * math.pi * (
            float(freq.freq) / float(self.samplerate))
        self.angle = self.angle + anglechange

        # Generate sine
        sample = math.sin(self.angle)
        return sample


#------------------------------------------------------------------------------
# Data tables and general settings
#------------------------------------------------------------------------------

# Sample rate
samplerate = 16000.0

# How many samples are run through the filter / test frequency band
scount = 3000

# Bit count for integer filtering
bits = 16

# Incoming sample scaling
prescale = 7.7336

# Generate test frequency table (from 100Hz to 3600 Hz in 100 Hz steps)
freqtable = map(Freq, range(100, 5000, 100))

# Second order sections with coefficients (floating point)
sections = [
    SecondOrderSection(1.7598, -0.8026, 0.1239, -0.1107, 0.1239),
    SecondOrderSection(1.8028, -0.9445, 0.0146, -0.0247, 0.0146)
]

# Second order sections with coefficients (fixed point)
intsections = [
    SecondOrderSection(57664, -26300, 4058, -3626, 4058),
    SecondOrderSection(59074, -30951, 480, -811, 480)
]

# Create sine generator
sinegen = SineGenerator(samplerate)

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

    # Test all frequency bands (float)
    for freq in freqtable:

        # Run scount samples with same frequence through filter
        for i in range(scount):
            sample = prescale * sinegen.generate(freq)

            # Run sample through filter sections
            for sos in sections:
                sample = sos.filter(sample)

            freq.collect(sample)
            sbuf.append(int(sample * (2**(bits - 1))))

    waveout.writeframes(sbuf.tostring())
    waveout.close()

    # Print float results
    print "-----------------------------------------------------------------"
    print "                 Floating point results                          "
    print "-----------------------------------------------------------------"
    for freq in freqtable:
        print 'Freq: %4d, Maxvalue: %3f, Atten: %-3.1f dB' % \
         (freq.freq, freq.maxval, 20 * math.log10( freq.maxval ) )

    # ---------- Fixed point -------------
    waveout = wave.open("integer.wav", "w")
    waveout.setnchannels(1)
    waveout.setsampwidth(2)
    waveout.setframerate(samplerate)
    waveout.setcomptype("NONE", "Uncompressed")
    sbuf = array.array('h')

    # Test all frequency bands (integer)
    for freq in freqtable:

        # Run scount samples with same frequence through filter
        for i in range(scount):
            sample = prescale * sinegen.generate(freq)
            sample = int(sample * (2**(bits - 1)))

            # Run sample through filter sections
            for sos in intsections:
                sample = sos.filter(sample, bits)

            freq.collect(sample)
            sbuf.append(sample)

    waveout.writeframes(sbuf.tostring())
    waveout.close()

    # Print int results
    print "-----------------------------------------------------------------"
    print "                   Fixed point results                           "
    print "-----------------------------------------------------------------"
    for freq in freqtable:
        print 'Freq: %4d, Maxvalue: %8d, Atten: %2.1f dB' % \
         (freq.freq, freq.maxval, 20 * math.log10( float( freq.maxval ) / 32768.0 ) )

    return


if __name__ == "__main__":
    main()
