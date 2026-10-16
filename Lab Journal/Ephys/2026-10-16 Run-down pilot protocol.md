# Run-down / stability pilot (aCSF only, no blockers)

Aim: see how long a cell stays usable at one fixed holding potential, and whether Rs, Rm, Cm and sEPSC/sIPSC frequency and amplitude drift over ~25 min. Sets the Rs-drift exclusion rule and tells me whether the randomised holding order is needed.

## Setup
- Bath: aCSF only (no DNQX, AP5, gabazine, TTX)
- Cells: 3-4, ideally from at least 2 slices. Tick cells off below.
- Holding: -50 (same list as the holding pilot). Check this is the value I type in Clampex, not the true value after LJP (typed = true + 15).
- Recording: 2 min gap-free each, 10 recordings, about 25 min per cell
- Manual seal test between every recording. Write down Rs, Rm, Cm each time.
- Start criteria: Rs stable, Rm good, cell looks healthy. Otherwise discard before starting.

## Cell ___  (slice ___, date ___, DIV ___, sex ___)

Manual check 1   Rs ___  Rm ___  Cm ___

1. -50

Manual check 2   Rs ___  Rm ___  Cm ___

2. -50

Manual check 3   Rs ___  Rm ___  Cm ___

3. -50

Manual check 4   Rs ___  Rm ___  Cm ___

4. -50

Manual check 5   Rs ___  Rm ___  Cm ___

5. -50

Manual check 6   Rs ___  Rm ___  Cm ___

6. -50

Manual check 7   Rs ___  Rm ___  Cm ___

7. -50

Manual check 8   Rs ___  Rm ___  Cm ___

8. -50

Manual check 9   Rs ___  Rm ___  Cm ___

9. -50

Manual check 10   Rs ___  Rm ___  Cm ___

10. -50

Manual check 11   Rs ___  Rm ___  Cm ___

## Notes during recording
- Stop and note the time if Rs changes by more than about 20% from Manual check 1, or if the holding current drifts a lot
- Note any solution or flow change, bubbles, or noise bursts

## Analysis
- Plot Rs, Rm, Cm against time per cell
- Per 2 min recording: sEPSC and sIPSC frequency and amplitude
- Compare recordings 1-2 with 9-10. Run-down if frequency or amplitude drops more than about 20%
- Write down the time point where the cell stops being usable (Rs drift over 20% or events clearly declining). Use that as the cutoff for future experiments.
- If run-down is small over 25 min, the randomised holding order is less critical. If it is large, keep the balanced randomisation.
