# fBlockKit report

## Summary

Subject: work\g09_h2co_ts.out
Counts: Warning 1

## Findings

### 1. [Warning] A transition-state optimisation without a frequency verification of the final geometry

- Suggested action: The optimised structure has not been verified as a saddle point. Run a final frequency analysis on the optimised geometry (expect exactly one imaginary mode, at the same level of theory as the optimisation), or add Freq to the OptTS job so that the frequency step follows the optimisation. The tutorial's recommended procedure computes an exact Hessian first and then runs OptTS reading it (Calc_Hess true in the %geom block, or inhess read with a Hessian file from a previous calculation). A converged OptTS run on its own is not evidence: the optimisation converges to a stationary point without saying of which order it is.
- Rule: D4-TS-NO-FREQUENCY-VERIFICATION
- Evidence:
  - [Manual] "By adding the !Freq simple input keyword aside !OptTS, a frequency calculation after the TS optimization is requested. This is useful to verify that the optimization yielded a transition state that is characterized by presence of exactly one imaginary frequency." -- "Run a direct !OptTS calculation requesting a prior exact Hessian calculation (recommended)."
    - Source: ORCA 6.1 tutorial, Reaction path / TS optimization (https://www.faccts.de/docs/orca/6.1/tutorials/react/tsopt.html)
  - [Measured] Two fixtures started from one OptTS guess: the job without Freq ("! r2SCAN-3c OptTS") converges (HURRAY, OPTIMIZATION RUN DONE) and contains no VIBRATIONAL FREQUENCIES block at all, while the same guess with Freq ("! r2SCAN-3c OptTS Freq") converges and the frequency block after the optimisation reports one imaginary mode at -90.48 cm^-1.
    - Source: Fixtures fixtures/orca/fhh_optts_nofreq.out and fixtures/orca/fhh_optts_freq.out (ORCA 6.1.1)

## Provenance

- [Manual] ORCA 6.1 tutorial, Reaction path / TS optimization
- [Measured] Fixtures fixtures/orca/fhh_optts_nofreq.out and fixtures/orca/fhh_optts_freq.out (ORCA 6.1.1)
