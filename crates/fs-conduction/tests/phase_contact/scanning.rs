//! G1/G4: moving deposition and original phase-dependent interface feedback.
//! The two native tetrahedra stay on distinct latent plateaus. The oracle
//! integrates the declared P1 source and solves one independent contact
//! quadratic; it does not call a production transport solver or tangent.

#[path = "scanning/history_controls.rs"]
mod history_controls;

use super::*;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{
    EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun,
};
use fs_conduction::bc::history::march::enthalpy::scanning::{
    ScanTransport, ScanningEnthalpyRun,
};
use fs_conduction::bc::history::march::moving::{MovingHeatBudget, MovingHeatSource};
use fs_qty::QtyAny;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const POWER_W: f64 = 2.0;
const OLD: [f64; 8] = [25.0, 25.0, 25.0, 30.0, 60.0, 60.0, 60.0, 65.0];

fn scan_history<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(
        cx,
        ThermalBoundaryBuilder::new(mesh).adiabatic_remainder(),
        &Environment::earth_lab(),
        Vec::new(),
        0.0,
    )
    .unwrap()
}

fn scan_base<'h, 'm, 'c>(
    cx: &Cx<'_>,
    boundary: &'h ThermalBoundaryHistory<'m>,
    conductivity: &'h ConductivityModel,
    curves: &'c [EquilibriumEnthalpyPhaseCurve; 2],
    max_steps: usize,
) -> EnthalpyHistoryRun<'h, 'm, 'c> {
    EnthalpyHistoryRun::new(
        cx,
        ThermalHistoryProblem {
            boundary,
            material: conductivity,
            element_materials: None,
            interfaces: None,
            source: HeatSourceHistory::Fixed(ScalarField::Uniform(0.0)),
        },
        EnthalpyHistoryMaterials {
            table: &[
                ReferenceEnthalpyMaterial {
                    curve: &curves[0],
                    reference_density_kg_m3: 1.0,
                },
                ReferenceEnthalpyMaterial {
                    curve: &curves[1],
                    reference_density_kg_m3: 2.0,
                },
            ],
            of_element: &[0, 1],
            phase_conductivity: None,
        },
        &OLD,
        EnthalpyHistoryBudget {
            spatial: EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
            max_steps,
        },
    )
    .unwrap()
}

fn heater(z: TimeSignal) -> MovingHeatSource {
    MovingHeatSource::new(
        [
            TimeSignal::Constant(QtyAny::new(0.0, LENGTH)),
            TimeSignal::Constant(QtyAny::new(0.0, LENGTH)),
            z,
        ],
        [0.5, 0.5, 1.5],
        TimeSignal::Constant(QtyAny::new(POWER_W, COMPONENT_POWER_DIMS)),
        0.0,
        MovingHeatBudget {
            max_vertices: 8,
            max_elements: 2,
            min_active_vertices: 4,
            power_relative_tolerance: 1e-12,
        },
    )
    .unwrap()
}

// Only the two coincident origins and the two opposite tips see this footprint.
// Each nodal volume is 1/24 m^3. A P1 source q contributes
// (sum(q) + q_i)/(5*rho) J/(kg*s), including the off-diagonal source load.
// On the plateaus, R=1+0.02*mean(hA)+0.01*mean(hB), Q=50/R W.
// Contact adds 8*dt*Q to A's face nodes and -4*dt*Q to B's face nodes:
// Q*(R_predictor + 0.12*dt*Q)=50. No temperature increment encodes this feedback.
fn endpoint(old: &[f64; 8], dt: f64, center_z: f64) -> ([f64; 8], f64) {
    assert!((-0.25..=0.25).contains(&center_z));
    let origin = (1.0 - (center_z / 1.5).powi(2)).powi(2);
    let lower = (1.0 - ((1.0 + center_z) / 1.5).powi(2)).powi(2);
    let upper = (1.0 - ((1.0 - center_z) / 1.5).powi(2)).powi(2);
    let integral = (2.0 * origin + lower + upper) / 24.0;
    let source = [origin, 0.0, 0.0, lower, origin, 0.0, 0.0, upper]
        .map(|weight| POWER_W * weight / integral);
    let sums = [
        source[..4].iter().sum::<f64>(),
        source[4..].iter().sum::<f64>(),
    ];
    let predictor: [f64; 8] = std::array::from_fn(|i| {
        let side = usize::from(i >= 4);
        let density = if side == 0 { 1.0 } else { 2.0 };
        old[i] + dt * (sums[side] + source[i]) / (5.0 * density)
    });
    let resistance = 1.0
        + 0.02 * predictor[..3].iter().sum::<f64>() / 3.0
        + 0.01 * predictor[4..7].iter().sum::<f64>() / 3.0;
    let heat_w = 100.0 / (resistance + (resistance * resistance + 24.0 * dt).sqrt());
    let h = std::array::from_fn(|i| {
        let contact = if i % 4 == 3 {
            0.0
        } else if i < 4 {
            8.0 * heat_w
        } else {
            -4.0 * heat_w
        };
        predictor[i] + dt * contact
    });
    (h, heat_w)
}

#[test]
fn g1_moving_source_keeps_latent_contact_fraction_feedback_and_exact_p1_energy() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = scan_history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source_card = source(&cards, 1.0, false, false, 1.0);
        let query_point = point();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let contact_laws = contacts(
            cx, &mesh, &initial, &cards, &curves, &source_card, &query_point,
            sides(&mesh, 4, 2, 0.0), 0, false,
        );
        let path = TimeSignal::Ramp {
            t_start: 0.0,
            t_end: 0.1,
            from: QtyAny::new(-0.25, LENGTH),
            to: QtyAny::new(0.25, LENGTH),
        };
        let mut scan = ScanningEnthalpyRun::new(
            cx, scan_base(cx, &boundary, &conductivity, &curves, 4),
            vec![heater(path)], 1, None,
        )
        .unwrap()
        .with_contacts_newton(cx, contact_laws, config())
        .unwrap();

        let mut old = OLD;
        let mut previous = 0.0;
        for (ordinal, end) in [0.025, 0.05, 0.075, 0.1].into_iter().enumerate() {
            let dt = end - previous;
            let center_z = -0.25 + 5.0 * end;
            let (expected, expected_heat) = endpoint(&old, dt, center_z);
            let report = scan.advance_to(cx, end, config().step).unwrap();
            let ScanTransport::Contact(contact) = &report.transport else {
                panic!("phase-dependent contact was omitted from the moving-source solve");
            };
            let solved = report.solution();
            for i in 0..8 {
                assert!((0.0..100.0).contains(&expected[i]));
                near(solved.specific_enthalpy_j_kg[i], expected[i], 2e-8);
                near(solved.liquid_mass_fraction[i], expected[i] / 100.0, 2e-10);
            }
            assert_eq!(solved.temperature, [400.0, 400.0, 400.0, 400.0, 500.0, 500.0, 500.0, 500.0]);
            near(report.deposition[0].center_m[2], center_z, 1e-14);
            assert_eq!(report.deposition[0].active_vertices, 4);
            near(report.deposition[0].delivered_power_w, POWER_W, 1e-12);
            near(solved.source_w, POWER_W, 1e-12);
            near(solved.stored_energy_change_j, dt * POWER_W, 1e-9);
            near(contact.physical.energy_residual_j, 0.0, config().step.energy_tolerance_j);
            assert!(contact.physical.residual_norm_j <= config().residual_tolerance_j);

            let traces = contact.physical_contacts.material_trace_states();
            assert!(!traces.is_empty());
            for trace in traces {
                near(trace.liquid_mass_fractions[0], expected[..3].iter().sum::<f64>() / 300.0, 2e-10);
                near(trace.liquid_mass_fractions[1], expected[4..7].iter().sum::<f64>() / 300.0, 2e-10);
            }
            let flux = contact.physical_contacts.frozen_interfaces().fluxes(&solved.temperature).unwrap();
            near(-flux[0].heat_rate_a_to_b_w, expected_heat, 2e-8);
            let frozen_heat = 50.0 / (1.0
                + 0.02 * old[..3].iter().sum::<f64>() / 3.0
                + 0.01 * old[4..7].iter().sum::<f64>() / 3.0);
            assert!((expected_heat - frozen_heat).abs() > 1e-3);
            assert_eq!(scan.conduction().accepted_steps(), ordinal + 1);
            assert_eq!(scan.conduction().time_s(), end);
            old = expected;
            previous = end;
        }
    });
}

#[test]
fn g4_latent_scan_contact_policy_and_source_domain_refusals_preserve_resume() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = scan_history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [common::bulk("cold"), common::bulk("hot")];
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source_card = source(&cards, 1.0, false, false, 0.35);
        let query_point = point();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let path = TimeSignal::Table {
            times: vec![0.0, 0.1, 0.2],
            values: vec![-0.25, 0.25, 100.0],
            dims: LENGTH,
            interp: Interp::Hold,
        };
        let create = || {
            let contact_laws = contacts(
                cx, &mesh, &initial, &cards, &curves, &source_card, &query_point,
                sides(&mesh, 4, 2, 0.0), 0, false,
            );
            ScanningEnthalpyRun::new(
                cx, scan_base(cx, &boundary, &conductivity, &curves, 2),
                vec![heater(path.clone())], 1, None,
            )
            .unwrap()
            .with_contacts_newton(cx, contact_laws, config())
            .unwrap()
        };
        let mut scan = create();
        let before_t = scan.conduction().temperature().to_vec();
        let before_f = scan.conduction().liquid_mass_fraction().to_vec();
        let mut changed = config().step;
        changed.energy_tolerance_j *= 2.0;
        // At 0.1 s the physical endpoint leaves the original fa <= 0.35 claim.
        // At 0.2 s the nonzero footprint has no mesh support at all.
        assert!(endpoint(&OLD, 0.1, 0.25).0[..3].iter().all(|h| *h > 35.0));
        for (end, policy) in [(0.025, changed), (0.1, config().step), (0.2, config().step)] {
            assert!(scan.advance_to(cx, end, policy).is_err());
            assert_eq!(scan.conduction().specific_enthalpy_j_kg(), &OLD);
            assert_eq!(scan.conduction().temperature(), before_t.as_slice());
            assert_eq!(scan.conduction().liquid_mass_fraction(), before_f.as_slice());
            assert_eq!(scan.conduction().time_s(), 0.0);
            assert_eq!(scan.conduction().accepted_steps(), 0);
        }
        let accepted = scan.advance_to(cx, 0.025, config().step).unwrap();
        let (expected, _) = endpoint(&OLD, 0.025, -0.25);
        for (actual, expected) in accepted.solution().specific_enthalpy_j_kg.iter().zip(expected) {
            near(*actual, expected, 2e-8);
        }
        let replay = create().advance_to(cx, 0.025, config().step).unwrap();
        assert_eq!(accepted.solution().specific_enthalpy_j_kg, replay.solution().specific_enthalpy_j_kg);
        assert_eq!(scan.conduction().accepted_steps(), 1);
        assert_eq!(scan.conduction().time_s(), 0.025);
    });
}
