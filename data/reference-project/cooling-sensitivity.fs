; Run from the repository root:
; frankensim --json run data/reference-project/cooling-sensitivity.fs sensitivity.db
;
; Reuse the existing native pick-freeze experiment: 16 rows, 64 physical solves.
; The program's seed 7 is the physical seed. The referenced study explicitly
; declares sampling seed 29. Neither seed, the solver nor the laws are changed.
; Add :budget N to cooling.study to limit this invocation's evaluations; zero
; retains an empty resumable study. A stop returns its ordinary native run ID:
; frankensim --json study --resume <study-run-id> sensitivity.db
; Whole-program resume and aggregate program wall metering are not implied.
(study "cooling-global-sensitivity"
  (seed 0x7)
  (versions (constellation :lock "2026-07"))
  (budget (wall 120s) (mem 64MiB))
  (capability :cores 1 :mem 64MiB :wall 120s :ops (cooling.*))
  (let project (cooling.project "cooling-reference.fsim"))
  (cooling.study project :source "cooling-sensitivity.fsim"))
