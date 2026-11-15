(study "heatsink-fan-journey-a"
  (seed 0x7)
  (versions (constellation :lock "2026-07"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1 :mem 64MiB :wall 60s :ops (cooling.*))
  (let project (cooling.project "heatsink-fan.fsim"
    :hash "063b2b2c44b48974f594ad33b898427137a48b5b722708067aef8398e9291848"))
  (cooling.import project :sources ("heatsink.stl") :unit "m" :max-hole-edges 0)
  (cooling.run project :materials ("../../data/reference-project/aa6061.fsmcdpk")))
