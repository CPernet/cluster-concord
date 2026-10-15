# Windows R toolchain

The validated development toolchain is:

- `rig` 0.8.1
- R 4.6.1 (2026-06-24)
- `renv` 1.2.4
- sequential execution (`workers = 1`) by default

`renv.lock` records the runtime, test, and report dependencies. The project
library under `renv/library` is intentionally ignored and must not be committed.

## Restore

From PowerShell in the repository root:

```powershell
& 'C:\Program Files\rig\rig.exe' list
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' -e "renv::restore(prompt = FALSE)"
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' -e "renv::status()"
```

The `.Rprofile` activates `renv` automatically. A version mismatch reported by
`renv` should be resolved by selecting/installing R 4.6.1 with `rig`, not by
silently updating the lockfile.

## Verify

Run the fast deterministic suite:

```powershell
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' test\run_tests.R
```

Run the opt-in synthetic full-pipeline validation:

```powershell
$env:RUN_SLOW_TESTS = 'true'
& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe' test\run_tests.R
```

Production analysis retains `resamples = 1000`. Synthetic known-truth validation
sets `resamples = 0`; permutation and bootstrap behavior is checked separately
by the deterministic resampling tests.

## Parallel resampling

`workers > 1` starts local PSOCK worker processes. Worker initialization copies
the main process's `.libPaths()`, sources all canonical `src/R` modules using
absolute paths, and verifies required functions and package namespaces before
executing an iteration. Consequently, the repository checkout and restored
`renv` library must remain readable at the same paths for the duration of the
run. On a batch system, request the desired CPUs from the scheduler first and
set `workers` no higher than that allocation.

The two-worker regression cases require localhost socket connections. If a
container or security sandbox blocks worker sockets, run the suite in the normal
host or batch-job environment rather than treating the socket denial as an
analytical failure.
