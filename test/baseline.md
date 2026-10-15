# Baseline verification record

Date: 2026-08-05

Command intended by the rebuild plan:

```text
Rscript test/run_tests.R
```

Result: **setup blocker**. Neither `Rscript` nor `rig` was available on `PATH`
in the Windows project shell. No analytical code or tests were changed to
bypass the missing runtime.

The repository was clean at baseline on commit `4ec97ac`.

## Toolchain follow-up

The blocker was resolved on the same date by installing `rig` 0.8.1 and R
4.6.1. The unchanged baseline suite then passed:

```text
distance-matrix: .............
gam-imputation: ...........
load-and-prepdata: ........
DONE
```
