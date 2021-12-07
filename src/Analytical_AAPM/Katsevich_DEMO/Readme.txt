To do a simulation, you only need to execute the function FunctionKat.

FunctionKat has three parameters:
(1) COEF: this is the filtering coeficent Tau, it can be specified between 0 and 1.
(2) FilteringMode: Specify the filtering mode.
                   0: Single filtering operator paramized by Tau
                   1: Dual filtering operators parametized by Tau and 1-Tau
                   2: Dual symmetric operators parametized by Tau at the top and 1-Tau at the bottom
(3) Phantom: 0: Disk Phantom
             1: Head Phantom

Note: For standard Katsevich's algorithm, you can set COEF=0.5 and FilteringMode=0.

Othes simulation parameters are set between line 19 and 36 in "FunctionKat". The default
parameters will take about 1 hour to do one simulation on my computer.