(study "heatsink-fan-journey-a"
  (seed 0x7)
  (versions (constellation :lock "2026-07"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1 :mem 64MiB :wall 60s :ops (cooling.*))
  (let project (cooling.project "heatsink-fan.fsim"
    :hash "eb3c83dfbcc91cb7a2def04f89206b56802ff1770c6ac4405ffdccd275b01c07"))
  (cooling.import project :sources ("heatsink.stl") :unit "m" :max-hole-edges 0)
  (cooling.run project :materials ("../../data/reference-project/aa6061.fsmcdpk")))
