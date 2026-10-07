(study "marquee-bracket-2d-journey-b"
  (seed 0x539)
  (versions (constellation :lock "2026-07"))
  (budget (wall 300s) (mem 512MiB))
  (capability :cores 1 :mem 512MiB :wall 300s :ops (study.*))
  (study.run "bracket-2d.fsim" :hash "5cdb3f5e8cedba44f4cedb80c28ff36dcb1e92098fb2c5f6ea5359090252907d"))
