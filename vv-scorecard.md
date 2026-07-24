# FrankenSim public V&V scorecard

schema: 1
corpus_schema: 3
corpus_authority: seeded
corpus_digest: 8c522e2cfbe268f56fcf116e5af7ccf2d712af6b616fe274c80c7c06b7847a5d
adversarial_registry: 275940f3404a8d22984bed11df014472b12b112d9d13cd01aa4892145858539d
datasets: 28
datasets_by_level: A=19 B=5 C=4 D=0 E=0
datasets_by_axis: numerical-verification=19 cross-code-agreement=5 controlled-experimental-validation=4 blind-predictive-validation=0 field-monitoring=0 transferability-across-regimes=0 independent-reproduction=0
ledgered_run_records: 0
executed_adversarial_challenges: 0/8
false_acceptance_total: NO-DATA (0 executed challenges)
interval_coverage: NO-DATA (empirical interval-coverage machinery (e07) is not live; nominal coverage is never extrapolated)
external_axes: cross-code-agreement, controlled-experimental-validation, blind-predictive-validation, field-monitoring

This scorecard is a deterministic projection of the registered validation corpus and the supplied ledgered run results. It grants no authority: a cell with data reports outcome arithmetic; a cell without data reports NO-DATA, never zero. Corpus claim caps remain in force regardless of anything shown here.

## Known gaps

- qoi=center-temperature-rise regime=slab-thickness-m in [0.1, 0.1] m external_datasets=0
- qoi=component-peak-temperature regime=(no dataset registered) external_datasets=0
- qoi=fin-efficiency regime=m-times-l in [1, 1] 1 external_datasets=0
- qoi=normalized-temperature-excess regime=biot-number in [0, 0.1] 1; normalized-time in [1, 1] 1 external_datasets=0
- qoi=nusselt-number regime=reynolds-number in [1, 2300] 1 external_datasets=0
- qoi=observed-l2-order regime=element-degree in [1, 1] 1; mesh-size-m in [0.0625, 0.25] m external_datasets=0
- qoi=observed-l2-order regime=element-degree in [2, 2] 1; mesh-size-m in [0.0625, 0.25] m external_datasets=0
- qoi=outward-heat-flux regime=biot-number in [1, 1] 1 external_datasets=0
- qoi=outward-heat-flux regime=slab-thickness-m in [0.2, 0.2] m external_datasets=0
- qoi=probe-temperature regime=rectangle-aspect-ratio in [2, 2] 1 external_datasets=0
- qoi=temperature-nonuniformity regime=(no dataset registered) external_datasets=0
- qoi=thermal-conductance regime=radius-ratio in [2, 2] 1 external_datasets=0
- qoi=thermal-interface-resistance regime=(no dataset registered) external_datasets=0
- qoi=thermal-resistance regime=interface-resistance-k-per-w in [0.1, 0.1] kg^-1·m^-2·s^3·K external_datasets=0
- qoi=view-factor-12 regime=gap-to-extent-ratio in [0, 0] 1 external_datasets=0

## Per-QoI/regime cells

Axis order in the `axes` column: numerical-verification / cross-code-agreement / controlled-experimental-validation / blind-predictive-validation / field-monitoring / transferability-across-regimes / independent-reproduction.

| qoi | regime | refs | axes | external | prediction error | envelope | coverage | false acceptance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| average-nusselt-number | mass-flux in [500, 750] kg·m^-2·s^-1 | 1 | 0/0/1/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| average-nusselt-number | reynolds-number in [810, 3800] 1 | 1 | 0/0/1/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| center-temperature-rise | slab-thickness-m in [0.1, 0.1] m | 1 | 1/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| convective-thermal-resistance | reynolds-number in [810, 3800] 1 | 1 | 0/0/1/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| effective-heat-flux | mass-flux in [1000, 1200] kg·m^-2·s^-1; inlet-subcooling in [10, 20] K; wall-superheat in [-15, 6] K | 1 | 0/0/1/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| fin-efficiency | m-times-l in [1, 1] 1 | 1 | 1/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| friction-factor | mass-flux in [500, 750] kg·m^-2·s^-1 | 1 | 0/0/1/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| hotspot_thermal_margin | reference_cost_work_units in [250, 250] 1 | 1 | 0/1/0/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| normalized-temperature-excess | biot-number in [0, 0.1] 1; normalized-time in [1, 1] 1 | 1 | 1/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| nusselt-number | reynolds-number in [1, 2300] 1 | 2 | 2/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| observed-l2-order | element-degree in [1, 1] 1; mesh-size-m in [0.0625, 0.25] m | 5 | 5/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| observed-l2-order | element-degree in [2, 2] 1; mesh-size-m in [0.0625, 0.25] m | 2 | 2/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| outward-heat-flux | biot-number in [1, 1] 1 | 1 | 1/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| outward-heat-flux | slab-thickness-m in [0.2, 0.2] m | 1 | 1/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| pressure-drop | mass-flux in [1000, 1200] kg·m^-2·s^-1; inlet-subcooling in [10, 20] K; wall-superheat in [-15, 6] K | 1 | 0/0/1/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| pressure-drop | reynolds-number in [810, 3800] 1 | 1 | 0/0/1/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| probe-temperature | rectangle-aspect-ratio in [2, 2] 1 | 1 | 1/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| probe-temperature-k | same-discretization-tet-count in [576, 576] 1 | 1 | 0/1/0/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| probe-temperature-k | same-discretization-tet-count in [720, 720] 1 | 1 | 0/1/0/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| probe-temperature-k | same-discretization-tet-count in [768, 768] 1 | 1 | 0/1/0/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| probe-temperature-k | same-discretization-tet-count in [960, 960] 1 | 1 | 0/1/0/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| surge-front-position-z | t_star in [0.5, 2] 1 | 1 | 0/0/1/0/0/0/0 | 1 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| thermal-conductance | radius-ratio in [2, 2] 1 | 2 | 2/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| thermal-resistance | interface-resistance-k-per-w in [0.1, 0.1] kg^-1·m^-2·s^3·K | 1 | 1/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |
| view-factor-12 | gap-to-extent-ratio in [0, 0] 1 | 1 | 1/0/0/0/0/0/0 | 0 | NO-DATA | NO-DATA | NO-DATA | NO-DATA |

## Regime limitations

- `average-nusselt-number` references apply only within mass-flux in [500, 750] kg·m^-2·s^-1; outside this declared context the corpus asserts no claim.
- `average-nusselt-number` references apply only within reynolds-number in [810, 3800] 1; outside this declared context the corpus asserts no claim.
- `center-temperature-rise` references apply only within slab-thickness-m in [0.1, 0.1] m; outside this declared context the corpus asserts no claim.
- `convective-thermal-resistance` references apply only within reynolds-number in [810, 3800] 1; outside this declared context the corpus asserts no claim.
- `effective-heat-flux` references apply only within mass-flux in [1000, 1200] kg·m^-2·s^-1; inlet-subcooling in [10, 20] K; wall-superheat in [-15, 6] K; outside this declared context the corpus asserts no claim.
- `fin-efficiency` references apply only within m-times-l in [1, 1] 1; outside this declared context the corpus asserts no claim.
- `friction-factor` references apply only within mass-flux in [500, 750] kg·m^-2·s^-1; outside this declared context the corpus asserts no claim.
- `hotspot_thermal_margin` references apply only within reference_cost_work_units in [250, 250] 1; outside this declared context the corpus asserts no claim.
- `normalized-temperature-excess` references apply only within biot-number in [0, 0.1] 1; normalized-time in [1, 1] 1; outside this declared context the corpus asserts no claim.
- `nusselt-number` references apply only within reynolds-number in [1, 2300] 1; outside this declared context the corpus asserts no claim.
- `observed-l2-order` references apply only within element-degree in [1, 1] 1; mesh-size-m in [0.0625, 0.25] m; outside this declared context the corpus asserts no claim.
- `observed-l2-order` references apply only within element-degree in [2, 2] 1; mesh-size-m in [0.0625, 0.25] m; outside this declared context the corpus asserts no claim.
- `outward-heat-flux` references apply only within biot-number in [1, 1] 1; outside this declared context the corpus asserts no claim.
- `outward-heat-flux` references apply only within slab-thickness-m in [0.2, 0.2] m; outside this declared context the corpus asserts no claim.
- `pressure-drop` references apply only within mass-flux in [1000, 1200] kg·m^-2·s^-1; inlet-subcooling in [10, 20] K; wall-superheat in [-15, 6] K; outside this declared context the corpus asserts no claim.
- `pressure-drop` references apply only within reynolds-number in [810, 3800] 1; outside this declared context the corpus asserts no claim.
- `probe-temperature` references apply only within rectangle-aspect-ratio in [2, 2] 1; outside this declared context the corpus asserts no claim.
- `probe-temperature-k` references apply only within same-discretization-tet-count in [576, 576] 1; outside this declared context the corpus asserts no claim.
- `probe-temperature-k` references apply only within same-discretization-tet-count in [720, 720] 1; outside this declared context the corpus asserts no claim.
- `probe-temperature-k` references apply only within same-discretization-tet-count in [768, 768] 1; outside this declared context the corpus asserts no claim.
- `probe-temperature-k` references apply only within same-discretization-tet-count in [960, 960] 1; outside this declared context the corpus asserts no claim.
- `surge-front-position-z` references apply only within t_star in [0.5, 2] 1; outside this declared context the corpus asserts no claim.
- `thermal-conductance` references apply only within radius-ratio in [2, 2] 1; outside this declared context the corpus asserts no claim.
- `thermal-resistance` references apply only within interface-resistance-k-per-w in [0.1, 0.1] kg^-1·m^-2·s^3·K; outside this declared context the corpus asserts no claim.
- `view-factor-12` references apply only within gap-to-extent-ratio in [0, 0] 1; outside this declared context the corpus asserts no claim.

## Adversarial regime limitations

schema: 1
registry: 275940f3404a8d22984bed11df014472b12b112d9d13cd01aa4892145858539d
false_acceptance_count: 0

| case | regime | attacked assumption | evidence | assessment | dominant uncertainty | limitation |
| --- | --- | --- | --- | --- | --- | --- |
| biot-extremes-lumped-breakdown | biot-number-validity-boundary | lumped-temperature | retained:thermal-a-lumped-transient | NO-DATA | NO-DATA | Lumped-capacitance predictions are restricted to their declared small-Biot domain; high-Biot cases require spatial resolution or an explicit demotion. |
| contact-dominated-two-layer-stack | series-contact-resistance-dominant | known-contact-resistance | retained:thermal-a-contact-series | NO-DATA | NO-DATA | Perfect-contact or silently defaulted interface laws are outside scope; the interface card and its uncertainty must dominate the reported budget when appropriate. |
| fan-stall-multiple-operating-points | fan-stall-negative-slope | stable-fan-operating-point | NO-DATA:frankensim-extreal-program-f85xj.4.5 | NO-DATA | NO-DATA | A single steady fan-curve intersection is not admissible through stall or hysteresis; unresolved operating-point multiplicity requires refusal or demotion. |
| material-lot-property-variability | multi-lot-as-manufactured-material-state | fixed-material-properties | NO-DATA:frankensim-extreal-program-f85xj.4.5 | NO-DATA | NO-DATA | Nominal handbook values do not establish an as-manufactured lot; unresolved lot variability must remain in the material-property budget. |
| natural-convection-cavity-reversal | buoyancy-dominated-enclosure-cavity | forced-convection-closure | NO-DATA:frankensim-extreal-program-f85xj.4.3 | NO-DATA | NO-DATA | Forced-convection correlations do not cover buoyancy-driven reversal; use an independently retained cavity reference or refuse the regime. |
| radiation-dominated-low-flow-enclosure | low-flow-high-emissivity-enclosure | convection-dominates | retained:thermal-a-parallel-plate-view-factor | NO-DATA | NO-DATA | Convection-only predictions are outside scope when radiative exchange is material; emissivity, view-factor, and nonlinear-radiation discrepancy must remain visible. |
| recirculation-behind-strip-fins | forced-air-strip-fin-re-810-3800 | attached-flow | retained:pires-fonseca-2024-flat-strip-fins | NO-DATA | NO-DATA | Attached-flow correlations are not validated inside separated fin wakes; predict within the retained experimental envelope or demote for flow-topology uncertainty. |
| uncertain-blockable-vent-leakage | sealed-to-leaky-enclosure-transition | known-vent-leakage | NO-DATA:frankensim-extreal-program-f85xj.4.5 | NO-DATA | NO-DATA | Nominal vent geometry cannot stand in for as-built leakage; unknown leakage area must remain a boundary-condition uncertainty. |

identity: bf23ffb556c9e058ae51a587e0783902682070ad8be92fc7c283ed645645f308
