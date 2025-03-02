# Frequency-Response-calculator-based-on-the-resonance-phenomenon

This projects shows how you can calculate the `frequency response (Bode plots)` based on the output and input signals. The system has `2 poles and no zeros` and has the following transfer function:

![alt text](transferFunction.png)

The input and output signals are the following:

![alt text](Signals.png)

The `input signal (blue)` is a `sine wave` with variable frequency that has a `continuous component` and an `alternative component`, and also the noise is present.
The `output signal (red)` is also a `sine wave` with same characteristics as the input signal, but after it's start it shows the `resonance phenomenon` that modifies the modulus up to a maximum that we will call `resonance modulus ($M_r$)`.

First, we need to calculate $\omega_n$, $zeta$, and $K$. I calculated them based on the following formulas:

![alt text](Screenshot_1.png) , which means that `K is the ratio between the mean value of y and the mean value of u`

![alt text](zeta.png)

![alt text](naturalOscillations.png) , which means that we can calculate the `value of natural oscillations` $\omega_n$ based on the `value of the resonance oscillations $\omega_r$ and damping factor $\zeta$` 
