# GPU ensemble comparisons

Use identical equations, parameter arrays, output requirements, and explicit
precision across libraries. Compare controllers at measured achieved error;
equal tolerances do not establish equal accuracy. Keep the original PI-controlled
`GPUTsit5` alongside opt-in controller alternatives.

State whether timings include transfers, setup, and result materialization.
Use documented public solver APIs and fail if a requested GPU engine falls back
to another implementation. Do not label host-to-host measurements kernel timings.
