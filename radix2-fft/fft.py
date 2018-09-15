#!/usr/bin/python

import sys
import math

#-------------------------------------------------------------------------------
# This file contains Decimation-In-Frequency implementation of Radix-2 FFT
# 19.9.2004 Harri Salokorpi <harri.salokorpi@iki.fi>
#
# List comprehension explained:
# 	http://docs.python.org/tut/node7.html#SECTION007140000000000000000
#-------------------------------------------------------------------------------

# Default FFT length
default_fft_len = 8


#-------------------------------------------------------------------------------
# Progress printing
#-------------------------------------------------------------------------------
class progress:
    """Progress listener class"""

    # FFT start format string
    fft_start_format = ["Length: %(length)i => Stages: %(stage)i"]

    # Stage start format string
    stage_start_format = [
        "'-Groups at stage %(stage)i: %(groups)i (Groupsize=%(groupsize)s)"
    ]

    # Group start format string
    group_start_format = [
        "| '-Bflies at stage %(stage)i", "|   '-Group      : %(group)i",
        "|   '-Butterflies: %(bflycount)i"
    ]

    # Butterfly start format
    bfly_start_format = [
        "|     '- Initial values:", "|     |  '- a  = %(a)s",
        "|     |  '- b  = %(b)s", "|     |  '- k  = %(k)s",
        "|     |  '- i1 = %(i1)s", "|     |  '- i1 = %(i2)s"
    ]

    # Butterfly info format
    bfly_info_format = [
        "|     '- Results:", "|     |  '- a: %(orig_a)s => %(new_a)s",
        "|     |  '- b: %(orig_b)s => %(new_b)s"
    ]

    # Printing helper
    def p(self, format, values):
        for string in format:
            print string % values
        return

    # FFT started
    def fft_start(self, length, stages):
        values = {"length": length, "stage": stages}
        self.p(self.fft_start_format, values)
        return

    # FFT ended
    def fft_end(self, length, stages):
        print "'- FFT END"
        return

    # Stage started
    def stage_start(self, stage, groupcount, groupsize):
        values = {"stage": stage, "groups": groupcount, "groupsize": groupsize}
        self.p(self.stage_start_format, values)
        return

    # Stage ended
    def stage_end(self, stage, groupcount):
        return

    # Group started
    def group_start(self, stage, group, bflycount):
        values = {"stage": stage, "group": group, "bflycount": bflycount}
        self.p(self.group_start_format, values)
        return

    # Group ended
    def group_end(self, stage, group, butterflycount):
        return

    # Butterfly started
    def butterfly_start(self, stage, group, butterfly, a, b, k, i1, i2):
        values = {"a": a, "b": b, "k": k, "i1": i1, "i2": i2}
        self.p(self.bfly_start_format, values)
        return

    # Butterfly info
    def butterfly_info(self, stage, group, butterfly, a, b, x, y):
        values = {
            "orig_a": complex2string(a),
            "orig_b": complex2string(b),
            "new_a": complex2string(x),
            "new_b": complex2string(y)
        }
        self.p(self.bfly_info_format, values)
        return

    # Butterfly info
    def butterfly_end(self, stage, group, butterfly, a, b, index_a, index_b):
        return


def complex2string(c):
    """Converts complex number to nicely formatted string"""
    return "(%(real).2f, j%(complex).2f)" % {"real": c.real, "complex": c.imag}


#-------------------------------------------------------------------------------
# FFT functions
#-------------------------------------------------------------------------------


#
# Initialize FFT (calculate butterflies)
#
def fft_init(length, sign=1):
    """
		Calculate butterfly coefficients
		@param length:		FFT lenght
		return:				butterfly array
	"""
    # calculate butterfly coefficients with list comprehension
    return [
        complex(
            math.cos(-2.0 * sign * math.pi * float(x) / float(length)),
            -math.sin(2.0 * sign * math.pi * float(x) / float(length)))
        for x in range(0, length / 2)
    ]


#
# Reorder FFT sequence
# Borrowed from: http://www.python.org/topics/scicomp/recipes_in_python.html
#
def fft_reorder(x):
    """Reorder FFT elements"""
    N, x = len(x), x[:]

    for i in range(N):
        k, b, a = 0, N >> 1, 1
        while b >= a:
            if b & i: k = k | a
            if a & i: k = k | b
            b, a = b >> 1, a << 1
        if i < k:  # important not to swap back
            x[i], x[k] = x[k], x[i]

    return x


#
# Process butterfly
#
def fft_butterfly(a, b, k, butterflies):
    """
		Process single butterfly
		@param a:			Complex number
		@param b: 			Complex number
		@param k: 			butterfly coeff index
		@param butterflies: butterfly coeff array
		@return x,y:		a+b and (a-b)*W (where W is butterfly coeff)
	"""
    x = a + b
    y = (a - b) * butterflies[k]
    return x, y


#
# Calculate FFT transform
#
def fft_transform(progress, sequence, butterflies):
    # Check sequence length
    x = 0
    while 2**x < len(sequence):
        x += 1
    if 2**x != len(sequence):
        raise "Sequence length is not power of two!"

    # Sequence length
    length = len(sequence)

    # calculate stage count
    stages = int(math.log(length, 2))

    progress.fft_start(length, stages)

    # initialize buffer to zero with list comprehension
    fft_buffer = [complex(sequence[x]) for x in range(0, length)]

    for stage in range(0, stages):
        groups = 1 << stage
        groupsize = length >> (stage + 1)

        progress.stage_start(stage, groups, groupsize)

        for group in range(0, groups):

            # Butterfly count. Actual formula is FFT_LEN / 2**(stage+1)
            bflies = length >> stage + 1

            progress.group_start(stage, group, bflies)

            for bfly in range(0, bflies):
                k = bfly << stage
                i1 = group * groupsize * 2 + bfly
                i2 = group * groupsize * 2 + bfly + groupsize
                a = complex(fft_buffer[i1])
                b = complex(fft_buffer[i2])

                progress.butterfly_start(stage, group, bflies, a, b, k, i1, i2)

                # Process butterfly
                x, y = fft_butterfly(a, b, k, butterflies)

                progress.butterfly_info(stage, group, bflies, a, b, x, y)

                fft_buffer[i1] = x
                fft_buffer[i2] = y

    progress.fft_end(length, stages)

    # Unscramble the sequence
    fft_buffer = fft_reorder(fft_buffer)

    return fft_buffer


#
# Calculate FFT for sequence
#
def fft(progress, sequence):
    """
		Calculate FFT.
		@param sequence:	Sample sequence to transform
		@param butterflies:	butterfly coefficients
		@return:			FFT result
	"""
    # Calculate butterflies
    butterflies = fft_init(len(sequence))
    fft = fft_transform(progress, sequence, butterflies)
    return fft


#
# Calculate IFFT for sequence
#
def ifft(progress, sequence):
    butterflies = fft_init(len(sequence), -1)
    fft = fft_transform(progress, sequence, butterflies)
    return [fft[x] / len(fft) for x in range(0, len(fft) - 1)]


#
# Main function
#
def main():
    """
		Main method
	"""

    if len(sys.argv) > 1:
        fftsize = int(sys.argv[1])
    else:
        fftsize = default_fft_len

    print "Decimation-In-Frequency Radix-2 FFT"
    print "Harri Salokorpi <harri.salokorpi@iki.fi"
    print ""
    print "Usage: python fft-python.py <fft size>"
    print "  fft size defaults to %(size)s" % {"size": default_fft_len}
    print ""
    print "Used FFT length:", fftsize

    # Initialize test sequence
    sequence = [float(x) for x in range(1, fftsize + 1)]

    # Initialize progress monitor
    p = progress()

    # Calculate fft
    fft_result = fft(p, sequence)

    print ""
    print "FFT results:"
    for index, value in enumerate(fft_result):
        print index, ":", complex2string(value)

    print "Inverse FFT:"
    print ""

    ifft_result = ifft(p, fft_result)

    print ""
    print "FFT results:"
    for index, value in enumerate(ifft_result):
        print index, ":", value.real

    return


if __name__ == "__main__":
    main()
