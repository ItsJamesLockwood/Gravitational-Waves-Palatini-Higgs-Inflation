## 21/11

### Declaring new variables: 	
  when using real(dl)  to initialize a new variable, this must apparently be done before any variables are assigned values.

### HLATTICE begins then stalls on startup:
  for the tanh^4 potential, the executable ran the first step (i.e. initialising time) and didn’t proceed further. After troubleshooting, it turns out a parameter was a few orders of magnitude too high, thus the COSH function returned values beyond its allowed limit and the code stalled. This was resolved by correctly defining the parameter (in this case: xi – minimal coupling).

### parameter definition: 
  parameters cannot be defined using simple variables (e.g. the real(dl) lambda was multiplying a parameter in the definition of a second parameter, thus raising an error)
