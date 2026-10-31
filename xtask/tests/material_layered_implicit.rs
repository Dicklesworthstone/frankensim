//! Exact layered-body energy Jacobians share the original implicit heat residual.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedPack,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::thermal::{
    conduction::LayerConductionConfig,
    transport::{
        stack::{
            LayerStackSourceQuery,
            transient::{
                LayeredThermalIds, PreparedLayeredThermalBody, TransientLayerSourceQuery,
                implicit::{
                    LayeredThermalLinearization, LayeredThermalSolveControls,
                    LayeredThermalSolveError, LayeredThermalStep, LayeredThermalWorkspace,
                },
            },
        },
        storage::{HeatStorageAmount, HeatStorageComponentQuery},
    },
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    Boundary, CellId, DepositionBudget, DepositionError, DepositionSimulation, Event, MaterialId,
    Reservoir, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use std::{path::Path, process::Command};

const CONFIG: LayerConductionConfig = LayerConductionConfig {
    max_iterations: 128,
    absolute_heat_tolerance_w: 1e-12,
    relative_heat_tolerance: 1e-13,
};
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const CVOL: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const IDS: LayeredThermalIds = LayeredThermalIds {
    first_cell: 11,
    first_material: 21,
    first_link: 31,
    body_id: 7,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, tolerance {tolerance:.3e}"
    );
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn keys() -> [PropertyKey; 2] {
    [
        PropertyKey::with_quantity(
            "thermal-conductivity",
            quantity(QuantityKind::ThermalConductivity),
        ),
        PropertyKey::new("volumetric-heat-capacity", CVOL),
    ]
}
fn state(name: &str) -> MaterialStateId {
    MaterialStateId {
        chemistry: name.into(),
        phase: "solid".into(),
        process: "original fixed source context".into(),
        revision: 0,
    }
}
fn analytic_card(name: &str, conductivity: [f64; 2], capacity: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, dims, [base, slope]) in [
        (keys()[0].clone(), Dims([1, 1, -3, -1, 0, 0]), conductivity),
        (keys()[1].clone(), CVOL, capacity),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [280.0, 300.0, 340.0]
                        .map(|t| (t, base + slope * (t - 300.0)))
                        .to_vec(),
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    280.0,
                    340.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored affine test law".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(state(name), claims, vec![]).unwrap()
}
fn point(card: &MaterialCard, key: &PropertyKey, t: f64) -> QueryPoint {
    let (_, claim) = card
        .claims()
        .claims_ordered()
        .find(|(_, c)| &c.key == key)
        .unwrap();
    let mut point = QueryPoint::new();
    for (axis, &(lo, hi)) in claim.validity.bounds() {
        let value = if axis == "temperature" {
            t
        } else {
            assert_eq!(lo, hi);
            lo
        };
        point = point
            .with_quantity(axis, claim.validity.axis_quantities()[axis], value)
            .unwrap();
    }
    point
}
fn prepare(
    cards: &[MaterialCard],
    lengths: &[f64],
    counts: &[usize],
    contacts: &[f64],
    area: f64,
    span: [f64; 2],
    datum: f64,
) -> PreparedLayeredThermalBody {
    let keys = keys();
    let lower: Vec<_> = cards
        .iter()
        .map(|card| keys.each_ref().map(|key| point(card, key, span[0])))
        .collect();
    let upper: Vec<_> = cards
        .iter()
        .map(|card| keys.each_ref().map(|key| point(card, key, span[1])))
        .collect();
    let queries: Vec<_> = cards
        .iter()
        .enumerate()
        .map(|(i, card)| TransientLayerSourceQuery {
            conduction: LayerStackSourceQuery {
                id: 8 + i as u64,
                card,
                conductivity: &keys[0],
                thickness_m: lengths[i],
                temperature_axis: "temperature",
                lower: &lower[i][0],
                upper: &upper[i][0],
                selection: SINGLE,
            },
            capacity: HeatStorageComponentQuery {
                id: 8 + i as u64,
                card,
                capacity: &keys[1],
                amount: HeatStorageAmount::ReferenceVolumeM3(area * lengths[i]),
                temperature_axis: "temperature",
                lower: &lower[i][1],
                upper: &upper[i][1],
                selection: SINGLE,
            },
            mass_kg: 0.1,
            cell_count: counts[i],
        })
        .collect();
    PreparedLayeredThermalBody::try_new(&queries, contacts, area, datum).unwrap()
}
fn timeline<'a>(
    events: Vec<TimedEvent>,
    materials: &[(MaterialId, &'a dyn ThermalLaw)],
) -> DepositionSimulation<'a> {
    DepositionSimulation::new(
        events,
        materials,
        DepositionBudget {
            max_cells: 32,
            max_links: 32,
            max_events: 128,
            max_substeps: 1024,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn dense(j: &LayeredThermalLinearization) -> Vec<f64> {
    let n = j.dimension();
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        matrix[i * n + i] = j.diagonal()[i];
        if i + 1 < n {
            matrix[(i + 1) * n + i] = j.lower_diagonal()[i];
            matrix[i * n + i + 1] = j.upper_diagonal()[i];
        }
    }
    matrix
}
fn action(
    j: &LayeredThermalLinearization,
    which: usize,
    x: &[f64],
    out: &mut [f64],
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    match which {
        0 => j.apply_into(x, out, checkpoint),
        1 => j.apply_transpose_into(x, out, checkpoint),
        _ => j.apply_inverse_diagonal_into(x, out, checkpoint),
    }
}

#[test]
fn g1_constant_layer_energy_jacobian_includes_exact_exposed_and_reservoir_loss() {
    let cards = [
        analytic_card("constant a", [0.4, 0.0], [2e6, 0.0]),
        analytic_card("constant b", [0.2, 0.0], [1e6, 0.0]),
        analytic_card("constant c", [0.8, 0.0], [3e6, 0.0]),
    ];
    let body = prepare(
        &cards,
        &[0.02, 0.03, 0.01],
        &[1; 3],
        &[0.03, 0.04],
        0.02,
        [280.0, 340.0],
        300.0,
    );
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, IDS).unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    let initial_t = [310.0, 300.0, 290.0];
    let mut events = bound
        .birth_events(0.0, &initial_t, CONFIG, || Ok(()))
        .unwrap();
    let boundaries = [
        Boundary {
            ambient_k: 295.0,
            htc_w_m2_k: 10.0,
            source_w: 3.0,
            reservoir: Some(Reservoir {
                area_m2: 0.005,
                conductance_w_k: 0.2,
                temperature_k: 315.0,
            }),
        },
        Boundary {
            ambient_k: 280.0,
            htc_w_m2_k: 7.0,
            source_w: -1.5,
            reservoir: None,
        },
        Boundary {
            ambient_k: 305.0,
            htc_w_m2_k: 5.0,
            source_w: 2.5,
            reservoir: None,
        },
    ];
    for (i, boundary) in boundaries.into_iter().enumerate() {
        let Event::Birth(cell) = &mut events[i].event else {
            panic!("birth order")
        };
        cell.boundary = boundary;
    }
    let initial = timeline(events, &materials);
    for (id, area) in bound.cell_ids().into_iter().zip([0.015, 0.0, 0.02]) {
        close(initial.exposed_area_m2(id).unwrap(), area, 1e-17);
    }
    let dt = 50.0;
    let step = LayeredThermalStep::new(&bound, &initial, dt, CONFIG, || Ok(())).unwrap();
    assert_eq!(step.cell_ids(), &bound.cell_ids());
    assert_eq!(step.duration_s(), dt);
    assert_eq!(step.target_s(), dt);
    let increments = [16.0, -12.0, 6.0];
    let trial = step.evaluate(&increments, || Ok(())).unwrap();
    let j = trial.payload();
    let capacities = [800.0, 600.0, 600.0];
    let t: [f64; 3] = std::array::from_fn(|i| initial_t[i] + increments[i] / capacities[i]);
    let g = [
        0.02 / (0.01 / 0.4 + 0.03 + 0.015 / 0.2),
        0.02 / (0.015 / 0.2 + 0.04 + 0.005 / 0.8),
    ];
    let q = [g[0] * (t[0] - t[1]), g[1] * (t[1] - t[2])];
    let external = [
        3.0 + 0.15 * (295.0 - t[0]) + 0.2 * (315.0 - t[0]),
        -1.5,
        2.5 + 0.1 * (305.0 - t[2]),
    ];
    let watts = [
        external[0] - q[0],
        external[1] + q[0] - q[1],
        external[2] + q[1],
    ];
    let diag = [
        1.0 + dt * (g[0] + 0.35) / capacities[0],
        1.0 + dt * (g[0] + g[1]) / capacities[1],
        1.0 + dt * (g[1] + 0.1) / capacities[2],
    ];
    for i in 0..3 {
        close(j.diagonal()[i], diag[i], 1e-13);
        close(j.cell_temperatures_k()[i], t[i], 1e-12);
        close(j.heat_capacities_j_k()[i], capacities[i], 1e-10);
        close(trial.residual_j()[i], increments[i] - dt * watts[i], 2e-9);
    }
    for i in 0..2 {
        close(j.lower_diagonal()[i], -dt * g[i] / capacities[i], 1e-13);
        close(j.upper_diagonal()[i], -dt * g[i] / capacities[i + 1], 1e-13);
        close(j.link_responses()[i].power_a_to_b_w, q[i], 1e-11);
    }
    close(
        trial.audit().external_heat_j,
        dt * external.iter().sum::<f64>(),
        1e-9,
    );
    let matrix = dense(j);
    for col in 0..3 {
        let mut plus = increments;
        let mut minus = increments;
        plus[col] += 0.1;
        minus[col] -= 0.1;
        let a = step.evaluate(&plus, || Ok(())).unwrap();
        let b = step.evaluate(&minus, || Ok(())).unwrap();
        for row in 0..3 {
            close(
                (a.residual_j()[row] - b.residual_j()[row]) / 0.2,
                matrix[row * 3 + col],
                3e-9,
            );
        }
    }
}

#[test]
fn g1_g4_nonlinear_energy_operator_transpose_and_diagonal_actions_are_atomic() {
    let k = [[0.2, 0.002], [0.5, -0.002], [0.12, 0.001]];
    let c = [[1e6, 10000.0], [2e6, -5000.0], [0.8e6, 7000.0]];
    let cards: [MaterialCard; 3] =
        std::array::from_fn(|i| analytic_card(&format!("variable layer {i}"), k[i], c[i]));
    let body = prepare(
        &cards,
        &[0.02, 0.03, 0.015],
        &[2, 1, 2],
        &[0.02, 0.04],
        0.02,
        [280.0, 340.0],
        300.0,
    );
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, IDS).unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    let initial_t = [325.0, 317.0, 307.0, 299.0, 289.0];
    let initial = timeline(
        bound
            .birth_events(0.0, &initial_t, CONFIG, || Ok(()))
            .unwrap(),
        &materials,
    );
    let step = LayeredThermalStep::new(&bound, &initial, 75.0, CONFIG, || Ok(())).unwrap();
    let increments = [12.0, -9.0, 8.0, -4.0, 3.0];
    let trial = step.evaluate(&increments, || Ok(())).unwrap();
    let j = trial.payload();
    let layer = [0, 0, 1, 2, 2];
    let volume = [0.0002, 0.0002, 0.0006, 0.00015, 0.00015];
    for i in 0..5 {
        let [base, slope] = c[layer[i]];
        let x = initial_t[i] - 300.0;
        let e = volume[i] * (base * x + 0.5 * slope * x * x) + increments[i];
        let a = volume[i] * base;
        let b = volume[i] * slope;
        let temperature = 300.0 + 2.0 * e / (a + (a * a + 2.0 * b * e).sqrt());
        close(j.cell_temperatures_k()[i], temperature, 2e-12);
        close(
            j.heat_capacities_j_k()[i],
            volume[i] * (base + slope * (temperature - 300.0)),
            1e-9,
        );
    }
    let matrix = dense(j);
    for col in 0..5 {
        let mut plus = increments;
        let mut minus = increments;
        plus[col] += 0.02;
        minus[col] -= 0.02;
        let a = step.evaluate(&plus, || Ok(())).unwrap();
        let b = step.evaluate(&minus, || Ok(())).unwrap();
        for row in 0..5 {
            close(
                (a.residual_j()[row] - b.residual_j()[row]) / 0.04,
                matrix[row * 5 + col],
                3e-7,
            );
        }
    }
    assert!(
        j.lower_diagonal()
            .iter()
            .zip(j.upper_diagonal())
            .any(|(a, b)| (a - b).abs() > 0.01)
    );
    let x = [1.25, -2.5, 3.75, 0.125, -0.75];
    let y = [-0.5, 4.25, 1.0, -3.5, 2.0];
    let mut jx = [0.0; 5];
    let mut jty = [0.0; 5];
    j.apply_into(&x, &mut jx, || Ok(())).unwrap();
    j.apply_transpose_into(&y, &mut jty, || Ok(())).unwrap();
    close(
        y.iter().zip(jx).map(|(a, b)| a * b).sum(),
        x.iter().zip(jty).map(|(a, b)| a * b).sum(),
        2e-12,
    );
    for i in 0..5 {
        close(
            jx[i],
            (0..5).map(|col| matrix[i * 5 + col] * x[col]).sum(),
            2e-14,
        );
        close(
            jty[i],
            (0..5).map(|row| matrix[row * 5 + i] * y[row]).sum(),
            2e-14,
        );
    }
    let mut inverse = [0.0; 5];
    j.apply_inverse_diagonal_into(&x, &mut inverse, || Ok(()))
        .unwrap();
    for i in 0..5 {
        close(inverse[i], x[i] / j.diagonal()[i], 1e-14);
    }
    for which in 0..3 {
        let sentinel = [91.0, 92.0, 93.0, 94.0, 95.0];
        let mut output = sentinel;
        assert!(action(j, which, &x[..4], &mut output, || Ok(())).is_err());
        assert_eq!(output, sentinel);
        assert!(action(j, which, &x, &mut output[..4], || Ok(())).is_err());
        assert_eq!(output, sentinel);
        for invalid in [f64::NAN, f64::INFINITY] {
            let mut bad = x;
            bad[2] = invalid;
            assert!(action(j, which, &bad, &mut output, || Ok(())).is_err());
            assert_eq!(output, sentinel);
        }
        if which < 2 {
            assert!(action(j, which, &[f64::MAX; 5], &mut output, || Ok(())).is_err());
            assert_eq!(output, sentinel);
        }
        let mut checkpoints = 0;
        action(j, which, &x, &mut [0.0; 5], || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
        for stop in 1..=checkpoints {
            let mut called = 0;
            let result = action(j, which, &x, &mut output, || {
                called += 1;
                if called == stop {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            });
            assert!(matches!(result, Err(DepositionError::Cancelled)));
            assert_eq!(output, sentinel);
        }
    }
}

fn compiled(slug: &str) -> MaterialCard {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "layered-implicit-{}-{nonce}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    MaterialCard::assemble(state(slug), pack.claims().clone(), vec![]).unwrap()
}

#[test]
fn g1_g3_real_fsri_wall_uses_sparse_production_solve_and_original_heat_owner() {
    let cards = [
        compiled("fsri-gypsum-wallboard-unconditioned-2021"),
        compiled("fsri-mineral-wool-unconditioned-2021"),
        compiled("fsri-plywood-unconditioned-2021"),
    ];
    let substituted = [cards[0].clone(), cards[2].clone(), cards[2].clone()];
    let initial_t = [310.0, 308.0, 305.0, 302.0, 297.0, 294.0];
    let mut outcomes = Vec::new();
    for selection in [&cards, &substituted] {
        let body = prepare(
            selection,
            &[0.012, 0.04, 0.018],
            &[2; 3],
            &[0.02, 0.03],
            0.01,
            [288.15, 313.15],
            293.15,
        );
        let laws = body.thermal_laws(true).unwrap();
        let bound = body.bind(&laws, IDS).unwrap();
        let materials: Vec<_> = bound.materials().into_iter().collect();
        let initial = timeline(
            bound
                .birth_events(0.0, &initial_t, CONFIG, || Ok(()))
                .unwrap(),
            &materials,
        );
        let step = LayeredThermalStep::new(&bound, &initial, 3000.0, CONFIG, || Ok(())).unwrap();
        let mut workspace = LayeredThermalWorkspace::new(6).unwrap();
        assert!(
            step.evaluate(&[0.0; 6], || Ok(()))
                .unwrap()
                .payload()
                .diagonal()
                .iter()
                .any(|d| *d > 10.0)
        );
        let controls = LayeredThermalSolveControls {
            maximum_newton_iterations: 12,
            cell_tolerance_j: 1e-7,
            total_tolerance_j: 1e-7,
            ..LayeredThermalSolveControls::default()
        };
        let (accepted, linearization, report) = step
            .solve_with_entropy(controls, &mut workspace, || Ok(()))
            .unwrap();
        assert!(report.newton_iterations > 0 && report.newton_iterations < 12);
        assert!(report.maximum_cell_residual_j < 1e-7);
        let entropy = report.entropy.unwrap();
        assert!(entropy.thermal_storage_change_j_k > 0.0);
        assert!(entropy.conduction_generation_j_k > 0.0);
        assert!(entropy.temporal_defect_j_k > 0.0);
        close(accepted.energy_j(), initial.energy_j(), 1e-8);
        assert_eq!(report.energy.external_heat_j, 0.0);
        assert_eq!(linearization.dimension(), 6);
        assert_eq!(linearization.link_responses().len(), 5);
        let temperatures: Vec<_> = bound
            .cell_ids()
            .into_iter()
            .map(|id| accepted.temperature_k(id).unwrap())
            .collect();
        assert!(temperatures.iter().all(|t| *t > 294.0 && *t < 310.0));
        assert!(temperatures[0] < 309.0 && temperatures[5] > 295.0);
        close(accepted.mass_kg(), initial.mass_kg(), 0.0);
        assert_eq!(initial.time_s(), 0.0);
        outcomes.push(temperatures);
    }
    assert!(
        outcomes[0]
            .iter()
            .zip(&outcomes[1])
            .any(|(a, b)| (a - b).abs() > 0.1)
    );
}

#[test]
fn g1_g3_sparse_layer_solve_matches_heterogeneous_closed_form_and_replays() {
    let cards = [
        analytic_card("high capacity", [0.4, 0.0], [2e6, 0.0]),
        analytic_card("low capacity", [0.2, 0.0], [1e6, 0.0]),
    ];
    let body = prepare(
        &cards,
        &[0.02, 0.03],
        &[1, 1],
        &[0.03],
        0.02,
        [280.0, 340.0],
        300.0,
    );
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, IDS).unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    let initial = timeline(
        bound
            .birth_events(0.0, &[330.0, 290.0], CONFIG, || Ok(()))
            .unwrap(),
        &materials,
    );
    let initial_cells = initial.cells().clone();
    let dt = 100_000.0;
    let conductance = 0.02 / (0.01 / 0.4 + 0.03 + 0.015 / 0.2);
    let capacities = [800.0, 600.0];
    let difference = 40.0 / (1.0 + dt * conductance * (1.0 / capacities[0] + 1.0 / capacities[1]));
    let heat = dt * conductance * difference;
    let expected = [330.0 - heat / capacities[0], 290.0 + heat / capacities[1]];
    let step = LayeredThermalStep::new(&bound, &initial, dt, CONFIG, || Ok(())).unwrap();
    let controls = LayeredThermalSolveControls::default();
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    assert_eq!(workspace.cells(), 2);
    let (next, tangent, report) = step
        .solve_with_entropy(controls, &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(report.newton_iterations, 1);
    assert_eq!(report.residual_evaluations, 3);
    assert_eq!(report.backtracks, 0);
    for (id, t) in bound.cell_ids().into_iter().zip(expected) {
        close(next.temperature_k(id).unwrap(), t, 3e-12);
    }
    close(next.energy_j(), initial.energy_j(), 1e-8);
    assert_eq!(report.energy.external_heat_j, 0.0);
    let entropy = report.entropy.unwrap();
    assert!(entropy.thermal_storage_change_j_k > 0.0);
    assert!(entropy.conduction_generation_j_k > 0.0);
    assert!(entropy.temporal_defect_j_k > 0.0);
    assert_eq!(next.time_s(), dt);
    let (replay, repeated, replay_report) = step
        .solve_with_entropy(controls, &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(replay.cells(), next.cells());
    assert_eq!(repeated, tangent);
    assert_eq!(replay_report, report);
    let (without_entropy, _, ordinary) = step.solve(controls, &mut workspace, || Ok(())).unwrap();
    assert_eq!(without_entropy.cells(), next.cells());
    assert_eq!(ordinary.entropy, None);
    assert_eq!(ordinary.residual_evaluations, 2);
    assert_eq!(initial.cells(), &initial_cells);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g1_sparse_layer_solve_recovers_first_order_diffusion_decay_under_time_refinement() {
    let cards = [analytic_card(
        "diffusive relaxation",
        [0.4, 0.0],
        [2e6, 0.0],
    )];
    let body = prepare(&cards, &[0.04], &[2], &[], 0.02, [280.0, 340.0], 300.0);
    let laws = body.thermal_laws(false).unwrap();
    let bound = body.bind(&laws, IDS).unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    let initial = timeline(
        bound
            .birth_events(0.0, &[330.0, 290.0], CONFIG, || Ok(()))
            .unwrap(),
        &materials,
    );
    let target = 2000.0;
    // Each C = 800 J/K and G = 0.4 W/K, so the exact semidiscrete
    // temperature difference decays at lambda = 2G/C = 0.001 /s.
    let expected_difference = 40.0 * (-0.001_f64 * target).exp();
    let mut errors = Vec::new();
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    for count in [8, 16, 32] {
        let mut state = initial.clone();
        for index in 1..=count {
            let step = LayeredThermalStep::new(
                &bound,
                &state,
                target * index as f64 / count as f64,
                CONFIG,
                || Ok(()),
            )
            .unwrap();
            state = step
                .solve(
                    LayeredThermalSolveControls::default(),
                    &mut workspace,
                    || Ok(()),
                )
                .unwrap()
                .0;
        }
        let difference =
            state.temperature_k(CellId(11)).unwrap() - state.temperature_k(CellId(12)).unwrap();
        errors.push((difference - expected_difference).abs());
        close(state.energy_j(), initial.energy_j(), 1e-8);
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            ratio > 1.8 && ratio < 2.2,
            "first-order error ratio {ratio}: {errors:?}"
        );
    }
}

#[test]
fn g0_g4_sparse_layer_solve_refuses_work_domain_and_cancellation_atomically() {
    let cards = [analytic_card("nonlinear", [0.3, 0.003], [1e6, 15000.0])];
    let body = prepare(&cards, &[0.04], &[2], &[], 0.02, [280.0, 340.0], 300.0);
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, IDS).unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    let initial = timeline(
        bound
            .birth_events(0.0, &[335.0, 285.0], CONFIG, || Ok(()))
            .unwrap(),
        &materials,
    );
    let before = initial.cells().clone();
    let step = LayeredThermalStep::new(&bound, &initial, 3000.0, CONFIG, || Ok(())).unwrap();
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let controls = LayeredThermalSolveControls::default();
    for invalid in [
        LayeredThermalSolveControls {
            maximum_newton_iterations: 0,
            ..controls
        },
        LayeredThermalSolveControls {
            maximum_residual_evaluations: 0,
            ..controls
        },
        LayeredThermalSolveControls {
            cell_tolerance_j: f64::NAN,
            ..controls
        },
        LayeredThermalSolveControls {
            total_tolerance_j: -1.0,
            ..controls
        },
    ] {
        assert!(matches!(
            step.solve(invalid, &mut workspace, || Ok(())),
            Err(LayeredThermalSolveError::Thermal(
                DepositionError::InvalidInput(_)
            ))
        ));
    }
    assert!(matches!(
        step.solve(
            LayeredThermalSolveControls {
                maximum_residual_evaluations: 1,
                ..controls
            },
            &mut workspace,
            || Ok(())
        ),
        Err(LayeredThermalSolveError::Thermal(DepositionError::Budget(
            _
        )))
    ));
    assert!(matches!(
        step.solve(
            LayeredThermalSolveControls {
                maximum_newton_iterations: 1,
                ..controls
            },
            &mut workspace,
            || Ok(())
        ),
        Err(LayeredThermalSolveError::NotConverged { .. })
    ));
    assert!(LayeredThermalWorkspace::new(0).is_err());
    assert!(LayeredThermalWorkspace::new(513).is_err());
    assert!(LayeredThermalWorkspace::new(512).is_ok());
    assert!(
        step.solve(
            controls,
            &mut LayeredThermalWorkspace::new(1).unwrap(),
            || Ok(())
        )
        .is_err()
    );
    let mut calls = 0;
    let (next, _, report) = step
        .solve_with_entropy(controls, &mut workspace, || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    assert!(report.newton_iterations > 1);
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        let result = step.solve_with_entropy(controls, &mut workspace, || {
            count += 1;
            if count == stop {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            result,
            Err(LayeredThermalSolveError::Thermal(
                DepositionError::Cancelled
            ))
        ));
    }
    // Final entropy evaluation consumes the SAME work budget. No accepted
    // thermal state is published if only that last required query cannot run.
    assert!(matches!(
        step.solve_with_entropy(
            LayeredThermalSolveControls {
                maximum_residual_evaluations: report.residual_evaluations - 1,
                ..controls
            },
            &mut workspace,
            || Ok(())
        ),
        Err(LayeredThermalSolveError::Thermal(DepositionError::Budget(
            _
        )))
    ));
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
    assert!(next.temperature_k(CellId(11)).unwrap() < 335.0);
    close(next.energy_j(), initial.energy_j(), 1e-8);

    // No endpoint in the source domain can store this supplied heat. A line
    // search may retry in-domain candidates, but must never extrapolate a
    // capacity/conductivity chart or publish a clipped temperature.
    let mut events = bound
        .birth_events(0.0, &[335.0, 285.0], CONFIG, || Ok(()))
        .unwrap();
    let Event::Birth(cell) = &mut events[0].event else {
        panic!("birth order")
    };
    cell.boundary.source_w = 1e6;
    let heated = timeline(events, &materials);
    let original = heated.cells().clone();
    let impossible = LayeredThermalStep::new(&bound, &heated, 3000.0, CONFIG, || Ok(())).unwrap();
    assert!(
        impossible
            .solve(controls, &mut workspace, || Ok(()))
            .is_err()
    );
    assert_eq!(heated.cells(), &original);
    assert_eq!(heated.time_s(), 0.0);
}

#[test]
fn g0_g4_dedicated_step_refuses_wrong_scope_events_controls_and_stale_bindings() {
    let cards = [analytic_card("dedicated body", [0.3, 0.002], [1e6, 5000.0])];
    let body = prepare(&cards, &[0.04], &[2], &[], 0.02, [280.0, 340.0], 300.0);
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, IDS).unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    let events = bound
        .birth_events(0.0, &[320.0, 290.0], CONFIG, || Ok(()))
        .unwrap();
    let initial = timeline(events.clone(), &materials);
    let snapshot = (
        initial.cells().clone(),
        initial.links().clone(),
        initial.energy_j(),
    );
    for target in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(LayeredThermalStep::new(&bound, &initial, target, CONFIG, || Ok(())).is_err());
    }
    for config in [
        LayerConductionConfig {
            max_iterations: 0,
            ..CONFIG
        },
        LayerConductionConfig {
            absolute_heat_tolerance_w: -1.0,
            ..CONFIG
        },
        LayerConductionConfig {
            relative_heat_tolerance: f64::NAN,
            ..CONFIG
        },
    ] {
        assert!(LayeredThermalStep::new(&bound, &initial, 10.0, config, || Ok(())).is_err());
    }
    assert!(matches!(
        LayeredThermalStep::new(&bound, &initial, 10.0, CONFIG, || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    ));
    let lookalike_laws = body.thermal_laws(true).unwrap();
    let stale = body.bind(&lookalike_laws, IDS).unwrap();
    assert!(LayeredThermalStep::new(&stale, &initial, 10.0, CONFIG, || Ok(())).is_err());
    let mut extra_events = events.clone();
    let Event::Birth(mut extra_cell) = events[0].event.clone() else {
        panic!("birth order")
    };
    extra_cell.id = CellId(99);
    extra_events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Birth(extra_cell),
    });
    let extra = timeline(extra_events, &materials);
    assert!(LayeredThermalStep::new(&bound, &extra, 10.0, CONFIG, || Ok(())).is_err());
    let mut extra_links = events.clone();
    extra_links.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 99,
            a: CellId(11),
            b: CellId(12),
            area_m2: 0.005,
            conductance_w_k: 0.1,
        }),
    });
    let extra = timeline(extra_links, &materials);
    assert!(LayeredThermalStep::new(&bound, &extra, 10.0, CONFIG, || Ok(())).is_err());
    let mut scheduled = events.clone();
    scheduled.push(TimedEvent {
        time_s: 5.0,
        event: Event::Boundary(
            CellId(11),
            Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 1.0,
                reservoir: None,
                source_w: 1.0,
            },
        ),
    });
    let pending = timeline(scheduled, &materials);
    assert!(LayeredThermalStep::new(&bound, &pending, 10.0, CONFIG, || Ok(())).is_err());
    assert!(LayeredThermalStep::new(&bound, &pending, 5.0, CONFIG, || Ok(())).is_ok());
    let step = LayeredThermalStep::new(&bound, &initial, 10.0, CONFIG, || Ok(())).unwrap();
    for bad in [
        vec![],
        vec![0.0],
        vec![0.0; 3],
        vec![f64::NAN, 0.0],
        vec![f64::INFINITY, 0.0],
        vec![1e12, 0.0],
    ] {
        assert!(step.evaluate(&bad, || Ok(())).is_err());
    }
    let mut calls = 0;
    assert!(matches!(
        step.evaluate(&[0.0; 2], || {
            calls += 1;
            if calls == 4 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(calls, 4);
    assert!(step.evaluate(&[0.0; 2], || Ok(())).is_ok());
    // Zero physical heat can still have an unrepresentable derivative. Refuse
    // the overflowing Jacobian instead of publishing a misleading finite row.
    let single_body = prepare(&cards, &[0.04], &[1], &[], 0.02, [280.0, 340.0], 300.0);
    let single_laws = single_body.thermal_laws(false).unwrap();
    let single_bound = single_body.bind(&single_laws, IDS).unwrap();
    let single_materials: Vec<_> = single_bound.materials().into_iter().collect();
    let mut single_events = single_bound
        .birth_events(0.0, &[300.0], CONFIG, || Ok(()))
        .unwrap();
    let Event::Birth(cell) = &mut single_events[0].event else {
        panic!("birth order")
    };
    cell.boundary.htc_w_m2_k = 1e6;
    let single = timeline(single_events, &single_materials);
    let enormous =
        LayeredThermalStep::new(&single_bound, &single, f64::MAX, CONFIG, || Ok(())).unwrap();
    assert!(enormous.evaluate(&[0.0], || Ok(())).is_err());
    assert_eq!(initial.cells(), &snapshot.0);
    assert_eq!(initial.links(), &snapshot.1);
    assert_eq!(initial.energy_j(), snapshot.2);
    assert_eq!(initial.time_s(), 0.0);
}
