Project to change the suppression matrix of a calibration systems for XBPM readings.

XBPMs provide photocurrents readings to calculate X-ray positions. Nevertheless, due to its own characteristics and biases in mechanical and electronic structures, the calculations of beam position is not accurate. This project is an attempt to calculate the correction matrix by a Monte Carlo method. The matrix acts upon the values read by the XBPM blades and corrects them to reduce the biasing of positions in the linear region, namely, the central area between the four blades.

The random walk procedures change the element of the matrix by small steps, trying to reduce the differences between calculated positions and their nominal values. At each step, a linear scaling by least square method is necessary. The final corrected values might provide reasonably accurate values for the X-ray position.

The program may start from a matrix with suppressions equal to 1 or from a previously calculated matrix (by methods as linear analysis of central the region).

This program is under GPL.
