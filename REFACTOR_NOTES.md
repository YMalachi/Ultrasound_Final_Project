# Refactoring Notes

## Objective

Refactor the original ultrasound course project into a clean,
documented portfolio repository without changing its scientific results.

## Baseline Environment

- MATLAB version: 24.2.0.2740171 (R2024b) Update 1
- Operating system: Windows 11
- Date tested: 12/06/2026

## Script Status

| Script | Runs successfully? | Output produced | Errors or warnings |
|---|---:|---|---|
| Question_1.m | Yes | 6 figures | No |
| Question_2.m | Not tested | | |
| Question_3.m | Not tested | | |
| test.m | Not tested | | |
| test2.m | Not tested | | |
| testtest.m | Not tested | | |
| three_b_backup.m | Not tested | | |
| ultrasoud_projQ3.m | Not tested | | |

### Question 1 baseline

- Script completed: Yes
- Number of figures: 6
- Concave transducer geometry displayed: Yes
- Time-domain signals displayed: Yes
- Frequency-domain signals displayed: Yes
- 2D pressure field displayed: Yes
- Lateral cuts displayed: Yes
- Errors: 0
- Warnings: 0

## Refactoring Rules

1. Do not change working numerical behavior unintentionally.
2. Make one logical change per commit.
3. Run the relevant script after every meaningful change.
4. Do not delete experimental files until the authoritative version is identified.
5. Use descriptive English names for scripts, functions, and variables.