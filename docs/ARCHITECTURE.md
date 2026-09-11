# Architecture

The package is organized by physics responsibility rather than notebook execution order.

Data flow:

halo -> interaction/Wilson coefficients -> nuclear response -> recoil rate -> detector response -> backgrounds/statistics -> sensitivity

Design rules:

- Each layer should depend only on information required for its own physics task.
- Wilson coefficients are the native interaction-level representation; standard SI O1 scattering is a special case.
- Halo models expose velocity distributions and velocity moments rather than only the mean inverse speed.
- Targets are isotope-aware and nuclear responses are kept separate from particle-response coefficients.
- Detector response is separate from recoil physics.
- Background models are separate from signal models and should expose provenance for empirical inputs.
- Statistical methods consume expected signal/background spectra or counts, not WIMP-specific functions.
- Validation tests are added alongside each implemented module.
- The original notebook remains a historical/reference calculation; production physics should live under src/.
