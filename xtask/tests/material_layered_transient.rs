//! Spatial source-backed storage uses the original conservative thermal timeline.
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
                LayeredThermalIds, MAX_LAYERED_THERMAL_CELLS, PreparedLayeredThermalBody,
                TransientLayerSourceQuery,
            },
        },
        storage::{HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw},
    },
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    CellId, DepositionBudget, DepositionError, DepositionSimulation, Event, MaterialId, StepConfig,
    ThermalLaw, TimedEvent, backward_euler::BackwardEulerStep,
};
use std::{f64::consts::PI, path::Path, process::Command};

const CONFIG: LayerConductionConfig = LayerConductionConfig {
    max_iterations: 128,
    absolute_heat_tolerance_w: 1e-11,
    relative_heat_tolerance: 1e-12,
};
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const CVOL: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
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
        process: "fixed original source context".into(),
        revision: 0,
    }
}
fn analytic_card(name: &str, conductivity: [f64; 2], capacity: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, dims, [base, slope]) in [
        (keys()[0].clone(), Dims([1, 1, -3, -1, 0, 0]), conductivity),
        (keys()[1].clone(), CVOL, capacity),
        (
            PropertyKey::new("specific-heat-capacity", CP),
            CP,
            [1000.0, 0.0],
        ),
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
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        280.0,
                        340.0,
                    )
                    .with_quantity(
                        "fixture-context",
                        QuantitySpec::dimensional(Dims::NONE),
                        1.0,
                        2.0,
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
    for (axis, &(low, high)) in claim.validity.bounds() {
        let value = if axis == "temperature" {
            t
        } else {
            assert!(low == high || axis == "fixture-context");
            low
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
    masses: &[f64],
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
            mass_kg: masses[i],
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
            max_cells: 512,
            max_links: 512,
            max_events: 1024,
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

#[test]
fn g1_insulated_cosine_mode_separates_spatial_and_backward_euler_refinement() {
    let cards = [analytic_card("homogeneous slab", [0.8, 0.0], [2e6, 0.0])];
    let (length, area, duration, amplitude) = (0.08, 0.05, 400.0, 8.0);
    let continuum_rate = 0.8 / 2e6 * (PI / length).powi(2);
    let mut previous_spatial_error = f64::INFINITY;
    for n in [2, 4, 8, 16] {
        let body = prepare(
            &cards,
            &[length],
            &[n],
            &[1.0],
            &[],
            area,
            [280.0, 340.0],
            300.0,
        );
        assert_eq!(body.cell_count(), n);
        assert_eq!(body.layer_cell_counts(), &[n]);
        assert_eq!(body.storages().len(), 1);
        close(body.total_mass_kg(), 1.0, 0.0);
        close(body.total_volume_m3(), area * length, 0.0);
        let laws = body.thermal_laws(true).unwrap();
        let bound = body.bind(&laws, IDS).unwrap();
        let materials: Vec<_> = bound.materials().into_iter().collect();
        let mode: Vec<_> = (0..n)
            .map(|i| (PI * (i as f64 + 0.5) / n as f64).cos())
            .collect();
        let initial_t: Vec<_> = mode.iter().map(|v| 300.0 + amplitude * v).collect();
        let initial = timeline(
            bound
                .birth_events(0.0, &initial_t, CONFIG, || Ok(()))
                .unwrap(),
            &materials,
        );
        let dx = length / n as f64;
        let cell_capacity = 2e6 * area * dx;
        let rate = 4.0 * (0.8 / 2e6) * (PI / (2.0 * n as f64)).sin().powi(2) / dx.powi(2);
        for cell in initial.cells().values() {
            close(cell.mass_kg, 1.0 / n as f64, 0.0);
            close(cell.geometry.volume_m3, area * dx, 2e-19);
            close(cell.geometry.end_m[0] - cell.geometry.start_m[0], dx, 2e-17);
            close(cell.geometry.surface_area_m2, 2.0 * area, 0.0);
            assert_eq!(cell.boundary.htc_w_m2_k, 0.0);
        }
        for link in initial.links().values() {
            close(link.conductance_w_k, 0.8 * area / dx, 1e-10);
        }
        let mut previous_temporal_error = f64::INFINITY;
        for steps in [1, 2, 4, 8] {
            let dt = duration / steps as f64;
            let mut current = initial.clone();
            let mut old_amplitude = amplitude;
            for step_index in 1..=steps {
                let next_amplitude = old_amplitude / (1.0 + dt * rate);
                let increments: Vec<_> = mode
                    .iter()
                    .map(|v| cell_capacity * (next_amplitude - old_amplitude) * v)
                    .collect();
                let step =
                    BackwardEulerStep::new(&current, step_index as f64 * dt, || Ok(())).unwrap();
                let trial = step
                    .evaluate(
                        &increments,
                        |candidate| bound.step_coefficients(candidate, CONFIG, || Ok(())),
                        || Ok(()),
                    )
                    .unwrap();
                assert!(trial.maximum_cell_residual_j() < 1e-7);
                let (next, _, audit) = trial.accept(1e-7, 1e-7, || Ok(())).unwrap();
                close(audit.external_heat_j, 0.0, 0.0);
                close(next.energy_j(), initial.energy_j(), 1e-7);
                current = next;
                old_amplitude = next_amplitude;
            }
            for (id, shape) in bound.cell_ids().into_iter().zip(&mode) {
                close(
                    current.temperature_k(id).unwrap(),
                    300.0 + old_amplitude * shape,
                    2e-11,
                );
            }
            let observed =
                (current.temperature_k(CellId(IDS.first_cell)).unwrap() - 300.0) / mode[0];
            let temporal_error = (observed - amplitude * (-rate * duration).exp()).abs();
            assert!(temporal_error < previous_temporal_error / 1.7);
            previous_temporal_error = temporal_error;
            if steps == 8 {
                // Same time discretization, exact continuum spatial eigenvalue.
                let reference = amplitude / (1.0 + dt * continuum_rate).powi(steps);
                let spatial_error = (observed - reference).abs();
                assert!(spatial_error < previous_spatial_error / 3.3);
                previous_spatial_error = spatial_error;
            }
        }
        assert_eq!(initial.time_s(), 0.0);
    }
}

#[test]
fn g1_g3_nonlinear_half_cells_contacts_storage_and_entropy_match_manufactured_solution() {
    let k = [[0.2, 0.002], [0.5, -0.003], [0.12, 0.001]];
    let capacity = [[1e6, 10000.0], [2e6, -5000.0], [0.8e6, 7000.0]];
    let cards: [MaterialCard; 3] =
        std::array::from_fn(|i| analytic_card(&format!("affine layer {i}"), k[i], capacity[i]));
    let final_t = [320.0, 305.0, 295.0];
    let contacts = [0.02, 0.04];
    let area = 0.01;
    let primitive =
        |c: [f64; 2], lo: f64, hi: f64| (hi - lo) * (c[0] + c[1] * (0.5 * (lo + hi) - 300.0));
    // The middle half-width is fixed. These prescribed middle-layer face
    // temperatures determine two different fluxes, contact jumps and outer widths.
    let half_middle = 0.01;
    let flux = [
        primitive(k[1], 305.0, 307.0) / half_middle,
        primitive(k[1], 302.0, 305.0) / half_middle,
    ];
    let faces = [
        [307.0 + contacts[0] * flux[0], 307.0],
        [302.0, 302.0 - contacts[1] * flux[1]],
    ];
    let lengths = [
        2.0 * primitive(k[0], faces[0][0], final_t[0]) / flux[0],
        2.0 * half_middle,
        2.0 * primitive(k[2], final_t[2], faces[1][1]) / flux[1],
    ];
    let body = prepare(
        &cards,
        &lengths,
        &[1; 3],
        &[0.5, 0.2, 0.05],
        &contacts,
        area,
        [280.0, 340.0],
        300.0,
    );
    let laws = body.thermal_laws(true).unwrap();
    let bound = body.bind(&laws, IDS).unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    let target = timeline(
        bound
            .birth_events(0.0, &final_t, CONFIG, || Ok(()))
            .unwrap(),
        &materials,
    );
    let coefficients = bound.step_coefficients(&target, CONFIG, || Ok(())).unwrap();
    for i in 0..2 {
        let response = &coefficients.payload[i];
        close(response.power_a_to_b_w, area * flux[i], 1e-10);
        close(response.interface_temperature_k[0], faces[i][0], 2e-9);
        close(response.interface_temperature_k[1], faces[i][1], 2e-9);
        let at = |layer: usize, t: f64| k[layer][0] + k[layer][1] * (t - 300.0);
        let r = 0.5 * lengths[i] / at(i, faces[i][0])
            + contacts[i]
            + 0.5 * lengths[i + 1] / at(i + 1, faces[i][1]);
        close(
            response.endpoint_derivative_w_k[0],
            area / r * at(i, final_t[i]) / at(i, faces[i][0]),
            1e-11,
        );
        close(
            response.endpoint_derivative_w_k[1],
            -area / r * at(i + 1, final_t[i + 1]) / at(i + 1, faces[i][1]),
            1e-11,
        );
        for side in 0..2 {
            let mut values = [0.0; 2];
            for (j, offset) in [-0.001, 0.001].into_iter().enumerate() {
                let mut t = final_t;
                t[i + side] += offset;
                let shifted = timeline(
                    bound.birth_events(0.0, &t, CONFIG, || Ok(())).unwrap(),
                    &materials,
                );
                values[j] = bound
                    .step_coefficients(&shifted, CONFIG, || Ok(()))
                    .unwrap()
                    .payload[i]
                    .power_a_to_b_w;
            }
            close(
                response.endpoint_derivative_w_k[side],
                (values[1] - values[0]) / 0.002,
                2e-8,
            );
        }
    }
    let duration = 10.0;
    let q = flux.map(|j| area * j);
    let increments = [-duration * q[0], duration * (q[0] - q[1]), duration * q[1]];
    let initial_t: [f64; 3] = std::array::from_fn(|i| {
        let v = area * lengths[i];
        let x = final_t[i] - 300.0;
        let energy = v * (capacity[i][0] * x + 0.5 * capacity[i][1] * x * x) - increments[i];
        let c = v * capacity[i][0];
        let slope = v * capacity[i][1];
        300.0 + 2.0 * energy / (c + (c * c + 2.0 * slope * energy).sqrt())
    });
    let initial = timeline(
        bound
            .birth_events(0.0, &initial_t, CONFIG, || Ok(()))
            .unwrap(),
        &materials,
    );
    let step = BackwardEulerStep::new(&initial, duration, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_entropy(
            &increments,
            |candidate| bound.step_coefficients(candidate, CONFIG, || Ok(())),
            || Ok(()),
        )
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-7);
    let storage_entropy: f64 = (0..3)
        .map(|i| {
            area * lengths[i]
                * ((capacity[i][0] - 300.0 * capacity[i][1])
                    * ((final_t[i] - initial_t[i]) / initial_t[i]).ln_1p()
                    + capacity[i][1] * (final_t[i] - initial_t[i]))
        })
        .sum();
    let generated: f64 = (0..2)
        .map(|i| duration * q[i] * (1.0 / final_t[i + 1] - 1.0 / final_t[i]))
        .sum();
    let entropy = trial.entropy_audit().unwrap();
    close(entropy.thermal_storage_change_j_k, storage_entropy, 1e-10);
    close(entropy.conduction_generation_j_k, generated, 1e-10);
    assert!(entropy.thermal_storage_change_j_k > 0.0);
    assert!(entropy.temporal_defect_j_k > 0.0);
    let (accepted, _, audit) = trial.accept(1e-7, 1e-7, || Ok(())).unwrap();
    for (id, expected) in bound.cell_ids().into_iter().zip(final_t) {
        close(accepted.temperature_k(id).unwrap(), expected, 2e-12);
    }
    close(accepted.energy_j(), initial.energy_j(), 1e-8);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(accepted.mass_kg(), initial.mass_kg());
    assert_eq!(initial.time_s(), 0.0);
}

fn compiled(slug: &str) -> MaterialCard {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "layered-transient-{}-{nonce}.fsmatpk",
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
fn g1_g3_real_fsri_cells_store_heat_and_material_substitution_changes_transient() {
    let cards = [
        compiled("fsri-gypsum-wallboard-unconditioned-2021"),
        compiled("fsri-mineral-wool-unconditioned-2021"),
        compiled("fsri-plywood-unconditioned-2021"),
    ];
    let substituted = [cards[0].clone(), cards[2].clone(), cards[2].clone()];
    let initial_t = [310.0, 308.0, 305.0, 302.0, 297.0, 294.0];
    let mut outcomes = Vec::new();
    let mut identities = Vec::new();
    for selection in [&cards, &substituted] {
        // Geometric volume is explicit; the supplied reference masses do not
        // convert FSRI's native volumetric capacity into an inferred density.
        let body = prepare(
            selection,
            &[0.012, 0.04, 0.018],
            &[2; 3],
            &[0.1; 3],
            &[0.02, 0.03],
            0.01,
            [288.15, 313.15],
            293.15,
        );
        identities.push(body.identity());
        assert_eq!(body.temperature_span_k(), [288.15, 313.15]);
        assert_eq!(body.storages().len(), 3);
        for (i, storage) in body.storages().iter().enumerate() {
            assert_eq!(
                storage.components()[0].amount(),
                HeatStorageAmount::ReferenceVolumeM3(0.01 * [0.012, 0.04, 0.018][i] / 2.0)
            );
            for line in [
                body.conduction().layers()[i].source_line(),
                storage.components()[0].source_line(),
            ] {
                for (_, answer) in &line.knots {
                    selection[i]
                        .claims()
                        .verify_receipt(&answer.receipt)
                        .unwrap();
                }
            }
        }
        let laws = body.thermal_laws(true).unwrap();
        let bound = body.bind(&laws, IDS).unwrap();
        let materials: Vec<_> = bound.materials().into_iter().collect();
        let initial = timeline(
            bound
                .birth_events(0.0, &initial_t, CONFIG, || Ok(()))
                .unwrap(),
            &materials,
        );
        let initial_cells = initial.cells().clone();
        let step = BackwardEulerStep::new(&initial, 1.0, || Ok(())).unwrap();
        // Small fixed dt gives a contractive test-local fixed-point iteration.
        // Only the original owner's residual admits the endpoint; this is not
        // a new production solver or a convergence claim for arbitrary dt.
        let mut increments = vec![0.0; 6];
        let mut converged = false;
        for _ in 0..32 {
            let trial = step
                .evaluate(
                    &increments,
                    |candidate| bound.step_coefficients(candidate, CONFIG, || Ok(())),
                    || Ok(()),
                )
                .unwrap();
            if trial.maximum_cell_residual_j() < 1e-9 {
                converged = true;
                break;
            }
            for (increment, residual) in increments.iter_mut().zip(trial.residual_j()) {
                *increment -= residual;
            }
        }
        assert!(converged);
        let trial = step
            .evaluate_with_entropy(
                &increments,
                |candidate| bound.step_coefficients(candidate, CONFIG, || Ok(())),
                || Ok(()),
            )
            .unwrap();
        assert!(trial.entropy_audit().unwrap().thermal_storage_change_j_k > 0.0);
        assert!(trial.entropy_audit().unwrap().conduction_generation_j_k > 0.0);
        let (accepted, response, audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
        assert_eq!(response.len(), 5);
        close(accepted.energy_j(), initial.energy_j(), 1e-8);
        assert_eq!(audit.external_heat_j, 0.0);
        let temperatures: Vec<_> = bound
            .cell_ids()
            .into_iter()
            .map(|id| accepted.temperature_k(id).unwrap())
            .collect();
        assert!(temperatures[0] < initial_t[0]);
        assert!(temperatures[5] > initial_t[5]);
        assert!(temperatures.windows(2).all(|w| w[0] > w[1]));
        assert!(
            temperatures
                .iter()
                .zip(initial_t)
                .all(|(a, b)| (a - b).abs() > 1e-7)
        );
        assert_eq!(initial.cells(), &initial_cells);
        assert_eq!(initial.time_s(), 0.0);
        outcomes.push(temperatures);
    }
    assert_ne!(identities[0], identities[1]);
    assert!(
        outcomes[0]
            .iter()
            .zip(&outcomes[1])
            .any(|(a, b)| (a - b).abs() > 1e-3)
    );
}

#[test]
fn g0_g4_bad_sources_geometry_bindings_and_work_refuse_without_changing_timeline() {
    let card = analytic_card("admitted layer", [0.3, 0.002], [2e6, 0.0]);
    let other_card = analytic_card("different source card", [0.3, 0.002], [2e6, 0.0]);
    let keys = keys();
    let lower = point(&card, &keys[0], 280.0);
    let upper = point(&card, &keys[0], 340.0);
    let q = TransientLayerSourceQuery {
        conduction: LayerStackSourceQuery {
            id: 8,
            card: &card,
            conductivity: &keys[0],
            thickness_m: 0.04,
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selection: SINGLE,
        },
        capacity: HeatStorageComponentQuery {
            id: 8,
            card: &card,
            capacity: &keys[1],
            amount: HeatStorageAmount::ReferenceVolumeM3(0.0004),
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selection: SINGLE,
        },
        mass_kg: 1.0,
        cell_count: 2,
    };
    let prepare_one = |query| PreparedLayeredThermalBody::try_new(&[query], &[], 0.01, 300.0);
    assert!(PreparedLayeredThermalBody::try_new(&[], &[], 0.01, 300.0).is_err());
    for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(PreparedLayeredThermalBody::try_new(&[q], &[], invalid, 300.0).is_err());
        assert!(
            prepare_one(TransientLayerSourceQuery {
                mass_kg: invalid,
                ..q
            })
            .is_err()
        );
    }
    for count in [0, MAX_LAYERED_THERMAL_CELLS + 1, usize::MAX] {
        assert!(
            prepare_one(TransientLayerSourceQuery {
                cell_count: count,
                ..q
            })
            .is_err()
        );
    }
    assert!(
        prepare_one(TransientLayerSourceQuery {
            conduction: LayerStackSourceQuery {
                thickness_m: 0.0,
                ..q.conduction
            },
            ..q
        })
        .is_err()
    );
    assert!(
        prepare_one(TransientLayerSourceQuery {
            capacity: HeatStorageComponentQuery {
                amount: HeatStorageAmount::ReferenceVolumeM3(0.0008),
                ..q.capacity
            },
            ..q
        })
        .is_err()
    );
    assert!(
        prepare_one(TransientLayerSourceQuery {
            capacity: HeatStorageComponentQuery {
                card: &other_card,
                ..q.capacity
            },
            ..q
        })
        .is_err()
    );
    let cp = PropertyKey::new("specific-heat-capacity", CP);
    let mass_based = TransientLayerSourceQuery {
        capacity: HeatStorageComponentQuery {
            capacity: &cp,
            amount: HeatStorageAmount::MassKg(1.0),
            ..q.capacity
        },
        ..q
    };
    assert!(prepare_one(mass_based).is_ok());
    assert!(
        prepare_one(TransientLayerSourceQuery {
            mass_kg: 2.0,
            ..mass_based
        })
        .is_err()
    );
    let changed_lower = lower
        .clone()
        .with_quantity(
            "fixture-context",
            QuantitySpec::dimensional(Dims::NONE),
            2.0,
        )
        .unwrap();
    let changed_upper = upper
        .clone()
        .with_quantity(
            "fixture-context",
            QuantitySpec::dimensional(Dims::NONE),
            2.0,
        )
        .unwrap();
    assert!(
        prepare_one(TransientLayerSourceQuery {
            capacity: HeatStorageComponentQuery {
                lower: &changed_lower,
                upper: &changed_upper,
                ..q.capacity
            },
            ..q
        })
        .is_err()
    );
    assert!(PreparedLayeredThermalBody::try_new(&[q], &[], 0.01, 279.0).is_err());
    let second = TransientLayerSourceQuery {
        conduction: LayerStackSourceQuery {
            id: 9,
            thickness_m: 0.02,
            ..q.conduction
        },
        capacity: HeatStorageComponentQuery {
            id: 9,
            amount: HeatStorageAmount::ReferenceVolumeM3(0.0002),
            ..q.capacity
        },
        mass_kg: 0.5,
        ..q
    };
    let body = PreparedLayeredThermalBody::try_new(&[q, second], &[0.02], 0.01, 300.0).unwrap();
    let laws = body.thermal_laws(true).unwrap();
    assert!(body.bind(&laws[..1], IDS).is_err());
    for ids in [
        LayeredThermalIds {
            first_cell: u64::MAX,
            ..IDS
        },
        LayeredThermalIds {
            first_material: u64::MAX,
            ..IDS
        },
        LayeredThermalIds {
            first_link: u64::MAX,
            ..IDS
        },
    ] {
        assert!(body.bind(&laws, ids).is_err());
    }
    let mut wrong_mass_laws = laws.clone();
    wrong_mass_laws[0] = HeatStorageThermalLaw::try_new(&body.storages()[0], 0.6, 300.0).unwrap();
    assert!(body.bind(&wrong_mass_laws, IDS).is_err());
    let bound = body.bind(&laws, IDS).unwrap();
    let materials: Vec<_> = bound.materials().into_iter().collect();
    let temperatures = [330.0, 320.0, 310.0, 290.0];
    for bad in [vec![300.0], vec![279.0; 4], vec![f64::NAN; 4]] {
        assert!(bound.birth_events(0.0, &bad, CONFIG, || Ok(())).is_err());
    }
    assert!(
        bound
            .birth_events(-1.0, &temperatures, CONFIG, || Ok(()))
            .is_err()
    );
    let mut calls = 0;
    assert!(matches!(
        bound.birth_events(0.0, &temperatures, CONFIG, || {
            calls += 1;
            if calls == 3 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(calls, 3);
    let events = bound
        .birth_events(0.0, &temperatures, CONFIG, || Ok(()))
        .unwrap();
    let initial = timeline(events.clone(), &materials);
    let snapshot = (
        initial.cells().clone(),
        initial.links().clone(),
        initial.energy_j(),
    );
    let lookalike_laws = body.thermal_laws(true).unwrap();
    let wrong_binding = body.bind(&lookalike_laws, IDS).unwrap();
    assert!(
        wrong_binding
            .step_coefficients(&initial, CONFIG, || Ok(()))
            .is_err()
    );
    for fault in 0..7 {
        let mut altered = events.clone();
        if fault < 5 {
            let Event::Birth(cell) = &mut altered[0].event else {
                panic!("birth order")
            };
            match fault {
                0 => cell.mass_kg *= 1.01,
                1 => cell.geometry.volume_m3 *= 1.01,
                2 => cell.geometry.end_m[0] += 0.001,
                3 => cell.material_id = MaterialId(IDS.first_material + 1),
                _ => cell.body_id += 1,
            }
        } else {
            let Event::Contact(link) = &mut altered[4].event else {
                panic!("contact order")
            };
            if fault == 5 {
                std::mem::swap(&mut link.a, &mut link.b);
            } else {
                link.area_m2 *= 0.5;
            }
        }
        let bad = timeline(altered, &materials);
        assert!(
            bound.step_coefficients(&bad, CONFIG, || Ok(())).is_err(),
            "binding fault {fault}"
        );
    }
    assert!(matches!(
        bound.step_coefficients(
            &initial,
            LayerConductionConfig {
                max_iterations: 1,
                ..CONFIG
            },
            || Ok(())
        ),
        Err(DepositionError::Budget(_))
    ));
    assert!(matches!(
        bound.step_coefficients(&initial, CONFIG, || Err(DepositionError::Cancelled)),
        Err(DepositionError::Cancelled)
    ));
    let step = BackwardEulerStep::new(&initial, 1.0, || Ok(())).unwrap();
    assert!(
        step.evaluate(
            &[1e12, 0.0, 0.0, 0.0],
            |candidate| bound.step_coefficients(candidate, CONFIG, || Ok(())),
            || Ok(())
        )
        .is_err()
    );
    assert!(
        step.evaluate(
            &[0.0; 4],
            |candidate| bound
                .step_coefficients(candidate, CONFIG, || Err(DepositionError::Cancelled)),
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(initial.cells(), &snapshot.0);
    assert_eq!(initial.links(), &snapshot.1);
    assert_eq!(initial.energy_j(), snapshot.2);
    assert_eq!(initial.time_s(), 0.0);
    assert!(bound.step_coefficients(&initial, CONFIG, || Ok(())).is_ok());
}
