Project to change the suppression matrix of a calibration systems for XBPM readings.

XBPMs provide photocurrents readings to calculate X-ray positions. Nevertheless, due to its own characteristics and biases in mechanical and electronic structures, the calculations of beam position is not accurate. This project is an attempt to calculate the correction matrix by an annealing method. The matrix acts upon the values read by the XBPM blades and corrects them to reduce the biasing in the linear region, namely, the central area between the four blades. After correction, a simple linear scaling to the real system dimension is applied. The final corrected values might provide reasonably accurate values for the X-ray position.

The program may start from a matrix with suppressions equal to 1 or from a previously calculated matrix (by methods as linear analysis of central the region).

This program is under GPL.
