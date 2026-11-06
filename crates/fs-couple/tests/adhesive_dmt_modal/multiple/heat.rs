//! Full native junction work and an independent shared-cell caloric oracle.
//! Constants are analytical fixtures, not measured material properties.
use super::*;
use fs_couple::friction_heat::adhesive::dmt_modal::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use std::cell::Cell;
use std::collections::BTreeMap;

const INDENT: f64 = 0.2 * Z;
const LENGTH: f64 = 1e-6;
const CAP: [f64; 3] = [1e-16, 1.5e-16, 2e-16];
const K: [f64; 3] = [0.08, 0.12, 0.20];
const THETA: [f64; 3] = [6.0, 0.0, -6.0];
const SLIP: [[f64; 2]; 2] = [[1.4 * Z, -0.8 * Z], [-1.1 * Z, 0.4 * Z]];
struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > -300.0 {
            Ok(300.0 + h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        from: f64,
        to: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let t = self.temperature_k(from)?;
        self.temperature_k(to)?;
        Ok(Some(((to - from) / t).ln_1p()))
    }
}
static CALORIC: Caloric = Caloric;

fn thermal(area: f64) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..3 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: CAP[i],
                specific_h_j_kg: THETA[i],
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, 2.0 * LENGTH * i as f64],
                    end_m: [0.0, 0.0, 2.0 * LENGTH * (i + 1) as f64],
                    volume_m3: 2.0 * LENGTH * area,
                    surface_area_m2: 2.0 * area + 4.0 * LENGTH * (PI * area).sqrt(),
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: [20, 10][i],
                a: CellId(i as u64),
                b: CellId(i as u64 + 1),
                area_m2: area,
                conductance_w_k: series_conductance_w_k(area, [LENGTH; 2], [K[i], K[i + 1]])
                    .unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0 / f64::from(RATE),
            energy_tolerance_j: 1e-27,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn routes(t: &DepositionSimulation<'_>) -> [FrictionHeatRoute; 2] {
    [
        FrictionHeatRoute::from_link(t, 20, [K[0], K[1]]).unwrap(),
        FrictionHeatRoute::from_link(t, 10, [K[1], K[2]]).unwrap(),
    ]
}
fn bindings(r: &[FrictionHeatRoute; 2], reverse: bool) -> [DmtModalHeatContact<'_>; 2] {
    (if reverse { [1, 0] } else { [0, 1] }).map(|i| DmtModalHeatContact {
        contact_index: i,
        route: &r[i],
        temperature_guess_k: 300.0,
    })
}
fn heat_controls() -> DmtModalHeatControls {
    DmtModalHeatControls {
        maximum_fractional_area_change: 1e-7,
        maximum_relative_clock_error: 1e-12,
        thermal: DmtModalThermalControls {
            maximum_cells: 3,
            maximum_interfaces: 2,
            maximum_unknowns: 5,
            maximum_trials: 256,
            maximum_newton_iterations: 20,
            difference_step_j: 1e-20,
            temperature_scale_j_per_k: 1e-16,
            heat: FrictionHeatTolerances {
                cell_j: 1e-26,
                thermal_j: 1e-26,
                interface_j: 1e-26,
                total_j: 1e-26,
            },
            combined_energy_tolerance_j: 1e-24,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-27,
                relative_tolerance: 1e-11,
                linear_restart: 5,
                forcing_minimum: 1e-6,
                ..Default::default()
            },
            audit_entropy: true,
        },
    }
}
fn source(
    _: &DepositionSimulation<'_>,
    ids: &[u64],
    temperatures: &[f64],
) -> Result<DmtModalHeatSource<[u64; 2]>, &'static str> {
    assert_eq!(ids, &[10, 20]);
    assert_eq!(temperatures.len(), 2);
    Ok(DmtModalHeatSource {
        interfaces: ids
            .iter()
            .map(|&link_id| DmtModalHeatConductivity {
                link_id,
                conductivity_w_m_k: if link_id == 10 {
                    [K[1], K[2]]
                } else {
                    [K[0], K[1]]
                },
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: [ids[0], ids[1]],
        },
    })
}
fn owner<'a>(
    n: [&'a SourceDmtSphereContact; 2],
    s: [&'a SourceAdhesiveShear; 2],
    d: f64,
) -> MultipleAdhesiveDmtModalSystem<'a> {
    MultipleAdhesiveDmtModalSystem::new(
        rich_network(),
        vec![
            contact(n[0], s[0], [d, 0.0, 0.0], columns()),
            contact(n[1], s[1], [d, 0.0, 0.0], second_columns()),
        ],
        controls(),
        &gate(),
    )
    .unwrap()
}

// Manufacture prescribed rates through the COMPLETE original response. The
// free native increment comes from the real network, retaining both memories.
fn rates(m: &MultipleAdhesiveDmtModalSystem<'_>, endpoint: f64) -> [[f64; 3]; 2] {
    let mut free = rich_network();
    let before = native_states(&free);
    free.step(&[0.0; 6]).unwrap();
    let delta: Vec<_> = native_states(&free)
        .iter()
        .zip(before)
        .map(|(a, b)| a.displacement_m_sqrt_kg - b.displacement_m_sqrt_kg)
        .collect();
    let mut forces = [[0.0; 3]; 2];
    for i in 0..2 {
        let n = m.normal(i).unwrap();
        let start = m.state(i).unwrap().position_m()[0];
        let mean_indent = if endpoint >= 0.0 {
            0.5 * (start + endpoint)
        } else {
            0.5 * start * start / (start - endpoint)
        };
        let capacity = m.shear(i).unwrap().law().kinetic_shear_strength_pa()
            * PI
            * n.effective_radius_m()
            * mean_indent;
        let length = SLIP[i][0].hypot(SLIP[i][1]);
        forces[i] = [
            if endpoint == start {
                n.evaluate(start).unwrap().normal_force_n
            } else {
                n.reversible_work_between(start, endpoint).unwrap() / (endpoint - start)
            },
            capacity * SLIP[i][0] / length,
            capacity * SLIP[i][1] / length,
        ];
    }
    core::array::from_fn(|i| {
        core::array::from_fn(|axis| {
            let response: f64 = (0..6)
                .map(|j| m.compliance_m_per_n()[(3 * i + axis) * 6 + j] * forces[j / 3][j % 3])
                .sum();
            let prediction: f64 = m.work_columns(i).unwrap()[axis]
                .iter()
                .zip(&delta)
                .map(|(b, q)| b * q)
                .sum();
            let change = if axis == 0 {
                endpoint - m.state(i).unwrap().position_m()[0]
            } else {
                SLIP[i][axis - 1]
            };
            (change + response - prediction) / m.sample_period_s()
        })
    })
}

// Eliminate both massless faces, then the two outer cells. This independent
// three-cell backward-Euler Schur solution uses original half-face conductances
// and actual mechanical loss; it never calls a production thermal residual.
fn caloric_reference(area: f64, dt: f64, heat: [f64; 2]) -> ([f64; 3], [f64; 2]) {
    let g = K.map(|k| k * area / LENGTH);
    let w = [
        dt * g[0] * g[1] / (g[0] + g[1]),
        dt * g[1] * g[2] / (g[1] + g[2]),
    ];
    let a = [CAP[0] + w[0], CAP[1] + w[0] + w[1], CAP[2] + w[1]];
    let b = [
        CAP[0] * THETA[0] + g[0] / (g[0] + g[1]) * heat[0],
        CAP[1] * THETA[1] + g[1] / (g[0] + g[1]) * heat[0] + g[1] / (g[1] + g[2]) * heat[1],
        CAP[2] * THETA[2] + g[2] / (g[1] + g[2]) * heat[1],
    ];
    let middle = (b[1] + w[0] * b[0] / a[0] + w[1] * b[2] / a[2])
        / (a[1] - w[0] * w[0] / a[0] - w[1] * w[1] / a[2]);
    let body = [
        (b[0] + w[0] * middle) / a[0],
        middle,
        (b[2] + w[1] * middle) / a[2],
    ];
    let face = core::array::from_fn(|i| {
        300.0 + (g[i] * body[i] + g[i + 1] * body[i + 1] + heat[i] / dt) / (g[i] + g[i + 1])
    });
    (body, face)
}

#[test]
fn native_path_and_shared_cell_heat_match_independent_mechanical_and_caloric_references() {
    let (n0, s0) = sources(0.01, 1e8, 6e7);
    let (n1, s1) = sources(0.014, 8e7, 4e7);
    let make = || owner([&n0, &n1], [&s0, &s1], INDENT);
    let mut reference = make();
    let rates = rates(&reference, INDENT);
    let mechanical = reference.step(&[0.0; 6], &rates).unwrap().clone();
    let area = n0.evaluate(INDENT).unwrap().contact_area_m2;
    let loss = core::array::from_fn(|i| mechanical.contacts[i].friction_dissipation_j);
    let (body, face) = caloric_reference(area, reference.sample_period_s(), loss);
    for i in 0..2 {
        close(
            mechanical.contacts[i].position_after_m[0],
            INDENT,
            2e-21,
            2e-8,
        );
        close(
            loss[i],
            [6e7, 4e7][i] * area * SLIP[i][0].hypot(SLIP[i][1]),
            1e-26,
            2e-8,
        );
        assert!(loss[i] > 0.0);
    }
    for reverse in [false, true] {
        let mut m = make();
        let mut t = thermal(area);
        let r = routes(&t);
        let report = m
            .step_with_friction_heat(
                &mut t,
                &bindings(&r, reverse),
                &[0.0; 6],
                &rates,
                heat_controls(),
                &[0.0; 3],
                source,
                &gate(),
            )
            .unwrap();
        assert_eq!(report.mechanical, mechanical);
        assert_eq!(bits(m.network()), bits(reference.network()));
        assert_eq!(m.network().last_frame(), reference.network().last_frame());
        assert_eq!(report.source, [10, 20]);
        assert_eq!(
            report
                .footprints
                .iter()
                .map(|f| f.contact_index)
                .collect::<Vec<_>>(),
            [1, 0]
        );
        assert_eq!(t.time_s(), m.sample_period_s());
        assert!(m.network().material_memory().iter().any(|x| *x != 0.0));
        assert!(m.network().radiation_state().iter().any(|x| *x != 0.0));
        for i in 0..3 {
            close(
                t.cells()[&CellId(i as u64)].specific_h_j_kg,
                body[i],
                1e-8,
                2e-9,
            );
        }
        for (j, i) in [1, 0].into_iter().enumerate() {
            let h = &report.heat.interfaces[j];
            assert_eq!(h.link_id, [20, 10][i]);
            close(h.interface_temperature_k, face[i], 1e-8, 2e-10);
            assert_eq!(h.dissipated_j, loss[i]);
            for side in 0..2 {
                let g = K[i + side] * area / LENGTH;
                close(h.conductances_w_k[side], g, 1e-24, 2e-13);
                close(
                    h.body_heat_j[side],
                    m.sample_period_s() * g * (face[i] - 300.0 - body[i + side]),
                    1e-25,
                    2e-9,
                );
            }
        }
        assert!(
            report.combined_energy_residual_j.abs()
                <= heat_controls().thermal.combined_energy_tolerance_j
        );
        assert!(
            report.absolute_stage_energy_residual_j
                <= heat_controls().thermal.combined_energy_tolerance_j
        );
    }
}

#[test]
fn frozen_area_and_mapping_guards_precede_heat_and_zero_slip_adds_no_heat() {
    let (n, s) = sources(0.01, 1e8, 6e7);
    let make = |d| owner([&n; 2], [&s; 2], d);
    let area = n.evaluate(INDENT).unwrap().contact_area_m2;
    let calls = Cell::new(0);
    for (initial, area_scale, endpoint, allowance) in [
        (INDENT, 1.0, 1.05 * INDENT, 0.01),
        (INDENT, 1.1, INDENT, 0.1),
        (0.0, 1.0, INDENT, 0.1),
        (-0.4 * Z, 1.0, INDENT, 0.1),
        (INDENT, 1.0, -0.1 * Z, 0.9),
    ] {
        let mut m = make(initial);
        let rates = if initial > 0.0 {
            rates(&m, endpoint)
        } else {
            [[0.0; 3]; 2]
        };
        let mut t = thermal(area * area_scale);
        let r = routes(&t);
        let before = bits(m.network());
        let cells = t.cells().clone();
        let mut c = heat_controls();
        c.maximum_fractional_area_change = allowance;
        assert!(
            m.step_with_friction_heat(
                &mut t,
                &bindings(&r, false),
                &[0.0; 6],
                &rates,
                c,
                &[0.0; 3],
                |t, ids, temps| {
                    calls.set(calls.get() + 1);
                    source(t, ids, temps)
                },
                &gate()
            )
            .is_err()
        );
        assert_eq!(calls.get(), 0);
        assert_eq!(bits(m.network()), before);
        assert_eq!(t.cells(), &cells);
        assert_eq!(t.time_s(), 0.0);
        assert!(m.last_frame().is_none());
    }
    let mut m = make(INDENT);
    let rates = rates(&m, 1.05 * INDENT);
    let mut t = thermal(area);
    let r = routes(&t);
    let mut c = heat_controls();
    c.maximum_fractional_area_change = 0.1;
    let report = m
        .step_with_friction_heat(
            &mut t,
            &bindings(&r, false),
            &[0.0; 6],
            &rates,
            c,
            &[0.0; 3],
            source,
            &gate(),
        )
        .unwrap();
    for f in report.footprints {
        close(f.fractional_area_change, 0.05, 2e-8, 2e-8);
    }
    // Diagonal native response and balanced normal forces retain exact zero
    // tangential motion. Ordinary hot-to-cold transport remains signed and real.
    let mut m = MultipleAdhesiveDmtModalSystem::new(
        network([0.0; 3], false, RATE),
        vec![
            contact(&n, &s, [INDENT, 0.0, 0.0], columns()),
            contact(&n, &s, [INDENT, 0.0, 0.0], columns()),
        ],
        controls(),
        &gate(),
    )
    .unwrap();
    let external: Vec<_> = columns()[0]
        .iter()
        .map(|b| 2.0 * b * n.evaluate(INDENT).unwrap().normal_force_n)
        .collect();
    let mut t = thermal(area);
    let r = routes(&t);
    let report = m
        .step_with_friction_heat(
            &mut t,
            &bindings(&r, false),
            &external,
            &[[0.0; 3]; 2],
            heat_controls(),
            &[0.0; 3],
            source,
            &gate(),
        )
        .unwrap();
    assert_eq!(report.mechanical.friction_dissipation_j, 0.0);
    assert_eq!(report.heat.dissipated_j, 0.0);
    assert!(
        report
            .heat
            .interfaces
            .iter()
            .all(|h| h.body_heat_j[0] < 0.0 && h.body_heat_j[1] > 0.0)
    );
    let (body, _) = caloric_reference(area, m.sample_period_s(), [0.0; 2]);
    for i in 0..3 {
        close(
            t.cells()[&CellId(i as u64)].specific_h_j_kg,
            body[i],
            1e-8,
            2e-9,
        );
    }
}

#[test]
fn thermal_source_receiver_cancel_and_budgets_rollback_to_exact_retry() {
    let (n, s) = sources(0.01, 1e8, 6e7);
    let make = || owner([&n; 2], [&s; 2], INDENT);
    let area = n.evaluate(INDENT).unwrap().contact_area_m2;
    let mut m = make();
    let mut t = thermal(area);
    let r = routes(&t);
    let b = bindings(&r, true);
    let rates = rates(&m, INDENT);
    let native = bits(m.network());
    let states = [m.state(0), m.state(1)];
    let cells = t.cells().clone();
    let links = t.links().clone();
    let unchanged = |m: &MultipleAdhesiveDmtModalSystem<'_>, t: &DepositionSimulation<'_>| {
        assert_eq!(bits(m.network()), native);
        assert_eq!([m.state(0), m.state(1)], states);
        assert!(m.last_frame().is_none());
        assert_eq!(t.time_s(), 0.0);
        assert_eq!(t.cells(), &cells);
        assert_eq!(t.links(), &links);
    };
    assert!(
        m.step_with_friction_heat(
            &mut t,
            &b,
            &[0.0; 6],
            &rates,
            heat_controls(),
            &[0.0; 3],
            |_, _, _| Err::<DmtModalHeatSource<[u64; 2]>, _>("original thermal source refused"),
            &gate()
        )
        .is_err()
    );
    unchanged(&m, &t);
    for trials in [false, true] {
        let mut c = heat_controls();
        if trials {
            c.thermal.maximum_trials = 1;
        } else {
            c.thermal.maximum_cells = 2;
        }
        assert!(
            m.step_with_friction_heat(&mut t, &b, &[0.0; 6], &rates, c, &[0.0; 3], source, &gate())
                .is_err()
        );
        unchanged(&m, &t);
    }
    for cancel in [false, true] {
        let cancelled = gate();
        let refused = m.step_with_friction_heat_admission(
            &mut t,
            &b,
            &[0.0; 6],
            &rates,
            heat_controls(),
            &[0.0; 3],
            source,
            &cancelled,
            |report, candidate| {
                assert_eq!(report.mechanical.sample, 1);
                assert_eq!(candidate.time_s(), 1.0 / f64::from(RATE));
                assert_ne!(candidate.cells(), &cells);
                if cancel {
                    cancelled.request();
                    Ok(())
                } else {
                    Err(DmtModalHeatError::Invalid("receiver refused"))
                }
            },
        );
        assert!(refused.is_err());
        unchanged(&m, &t);
    }
    // Duplicate mechanical coverage and duplicate thermal ownership refuse
    // before resolving temperature-dependent material inputs.
    for bad in [
        [b[0], b[0]],
        [
            DmtModalHeatContact {
                contact_index: 0,
                ..b[0]
            },
            b[1],
        ],
    ] {
        assert!(
            m.step_with_friction_heat(
                &mut t,
                &bad,
                &[0.0; 6],
                &rates,
                heat_controls(),
                &[0.0; 3],
                |_, _, _| -> Result<DmtModalHeatSource<[u64; 2]>, &'static str> {
                    panic!("invalid mapping reached heat")
                },
                &gate()
            )
            .is_err()
        );
        unchanged(&m, &t);
    }
    let mut replay = make();
    let mut replay_thermal = thermal(area);
    let a = m
        .step_with_friction_heat(
            &mut t,
            &b,
            &[0.0; 6],
            &rates,
            heat_controls(),
            &[0.0; 3],
            source,
            &gate(),
        )
        .unwrap();
    let z = replay
        .step_with_friction_heat(
            &mut replay_thermal,
            &bindings(&r, false),
            &[0.0; 6],
            &rates,
            heat_controls(),
            &[0.0; 3],
            source,
            &gate(),
        )
        .unwrap();
    assert_eq!(a.mechanical, z.mechanical);
    assert_eq!(bits(m.network()), bits(replay.network()));
    assert_eq!(t.time_s().to_bits(), replay_thermal.time_s().to_bits());
    assert_eq!(t.cells(), replay_thermal.cells());
    assert_eq!(a.footprints, z.footprints);
    assert_eq!(a.source, z.source);
    let native = bits(m.network());
    let mut wrong_clock = thermal(area);
    assert!(
        m.step_with_friction_heat(
            &mut wrong_clock,
            &b,
            &[0.0; 6],
            &rates,
            heat_controls(),
            &[0.0; 3],
            source,
            &gate()
        )
        .is_err()
    );
    assert_eq!(bits(m.network()), native);
    assert_eq!(wrong_clock.time_s(), 0.0);
}
