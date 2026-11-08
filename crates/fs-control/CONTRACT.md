# CONTRACT: fs-control

## Fixed-period synchronous current regulation

`fs-control` is an L3 sampled policy runtime. Its first production owner,
`synchronous::SynchronousCurrentController`, executes authored d/q PI gains and
emits actual averaged three-phase high-leg duty fractions. Its only runtime
dependency is deterministic `fs-math`; `fs-phs` is a test-only plant dependency.

Every alpha/beta and d/q current, flux and voltage is power invariant. Positive
phase order is a,b,c; alpha is `sqrt(2/3)*(a-(b+c)/2)`, beta is `(b-c)/sqrt(2)`.
Positive q is counterclockwise from d. The sampled angle is relative MECHANICAL
rotor-to-stator angle, multiplied by the authored positive pole-pair count.
Electrical speed is `p*(omega_rotor-omega_stator)` using both absolute shaft
speeds. There is no hidden 3/2 factor or fixed-stator assumption. The authored
current bound is a power-invariant d/q vector magnitude, not a phase RMS rating.

The requested current vector is radially limited to that bound. Torque requests
must first pass through the machine's original constitutive torque/current
allocator, which owns current/torque authority and any saturation curves. The
controller does not substitute a second nominal torque or magnetic law. Actual
measured flux supplies speed voltage `[-omega_e*lambda_q,+omega_e*lambda_d]`;
`R_winding*i_reference` supplies resistive feedforward. Inverter channel drops
remain plant-owned and can be rejected by the PI feedback through actual current.

For each declared sample period Ts, the exact digital rule is:

- `e = i_reference-i_measured`;
- `z_predict = z + Ts*Ki*e`;
- `u_raw = R*i_reference + speed_voltage + Kp*e + z_predict`;
- physical modulation maps `u_raw` to `u_command`;
- `z_next = z_predict + gamma*(u_command-u_raw)`, where
  `gamma = Ts/(Taw+Ts)` per axis, and `Taw=0` means immediate tracking.

With error, feedforward and saturated commanded voltage held fixed, the tracking
pole is `1-gamma` in `[0,1)`, even when Ts is larger than Taw. This is a scalar
anti-windup recurrence property, not a certificate for a changing hexagon or
coupled closed-loop plant. Unsaturated operation retains the ordinary full
integral increment. Integrator coordinates remain d/q as the rotor moves.

Physical modulation uses the complete two-level bridge hexagon: balanced phase
voltages are feasible exactly when their maximum-minus-minimum span does not
exceed sampled DC voltage. A single radial voltage scale and centered common-mode
offset produce all three duties in `[0,1]`; individual phases are not clipped.
The returned ideal bridge voltage is reconstructed from those actual duties and
the measured bus. Zero bus emits `[0.5,0.5,0.5]` and zero differential voltage
with active tracking; negative bus refuses. Future bus evolution and device
voltage drops can change the actual held-interval terminal voltage. The caller
holds the returned duties for the full Ts, including any plant substeps.

`CurrentControllerState` is caller owned. A sample returns one
`CurrentControlProposal`, containing candidate state, physical duties, requested
and admitted currents, actual transformed measurements and voltage-limit reports.
The caller adopts `next_state` only after the corresponding plant interval is
accepted. Refusal, cancellation or a rejected plant interval therefore leaves
accepted controller state intact; retry from identical inputs is bit-identical
on the same ISA. Cancellation is polled before arithmetic and before publication.
The fixed-size, fixed-work runtime allocates no heap memory and contains no unsafe
code. Numerical overflows, invalid domains and unrepresentable sampled gains
return typed errors; there is no hidden repair or partial output publication.

Focused tests independently reconstruct phase/bus power and the physical voltage
hexagon, including headroom beyond the inscribed circle; verify actual-flux
feedforward and both-shaft speed invariance; check admitted current and visible
measured overcurrent; solve the saturated tracking recurrence analytically; and
exercise zero-bus, cancellation, refusal and deterministic retry. Native finite-
bus electromechanical integration supplies the closed-loop product fixture.

`tests/synchronous_drive.rs` closes the real controller around the original
bounded saturated d/q sources, finite 24 V / 0.3 F bus, asymmetric channel losses
and both finite shaft inertias. Authored 100 microsecond samples command +1.2,
-1.2 and zero N m for 50 ms each through the machine-owned torque allocator.
The fixture checks settled torque, exercised voltage saturation, net regeneration
with rising bus voltage, original-source magnetic/kinetic/capacitor energy,
positive physical heat and equal/opposite shaft impulses. A physical negative-bus
candidate refusal and an after-stamp cancellation leave accepted PI/plant state
and duties unchanged; the successful retry matches a fresh native solve exactly.

This slice does not implement gain synthesis, automatic MTPA or field weakening,
observers, multirate clocks, delay/quantization/fault models, PWM
events, switching losses or protection hardware. Current-reference limiting does
not guarantee that a transient physical current remains inside that bound. It
does not claim closed-loop stability, HIL timing, experimental motor validation
or completion of the broader control-runtime/E4 beads. No feature flags apply.


## Relative-shaft speed, power limits and thermal derating

`speed::SpeedController` regulates rotor-minus-stator speed into signed electromagnetic torque. Authored positive/negative torque caps and separate sampled motoring/regeneration power caps use torque times relative speed in all four quadrants. Exactly zero relative speed retains finite starting torque without a denominator floor. Demand predicts the fixed-period PI contribution; finish tracks the ORIGINAL machine allocator's achievable torque using `gamma=Ts/(Taw+Ts)`. Torque and power endpoints are checked independently with 32 machine-epsilon relative slack on each endpoint's own magnitude; zero capabilities remain exact. Both outer and inner PI proposals are accepted only with their successful native plant interval.

`thermal::ThermalDeratingPolicy` applies authored full-output/trip/reset temperatures to 1..=64 original measured components. The minimum linear component factor scales all four directional capabilities without averaging temperatures. Any trip latches zero requested torque until an explicit reset with every channel at or below its reset temperature. The policy has fixed bounded sample work, allocates nothing, preserves the original cap asymmetry and proposes its latch state transactionally. It changes command capability, never physical energy, current, motion or cooling.

Pure tests independently cover signed quadrants/common rotation, analytic PI/tracking, exact disabled capabilities, thermal bands, latch/reset/cancellation and the thermal-to-speed composition. The native battery-synchronous fixture connects the original torque allocator and FOC controller, checking load recovery, regenerative chemical recharge, independently reconstructed energy, both shaft impulses and refused/cancelled interval retry.

These are sampled command limits. They do not guarantee transient current, battery power or between-sample thermal bounds. Zero torque demand is not electrical isolation, commanded braking or hardware protection. Allocator anti-windup does not certify voltage feasibility; gain synthesis, automatic MTPA/field weakening, scheduling and closed-loop stability remain separate capabilities.
