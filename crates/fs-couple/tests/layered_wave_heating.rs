#![cfg(feature = "layered-wave-heating")]
//! Authored constitutive cards and independent optical/caloric oracles.
//! These tests establish implementation behavior, not material qualification.

#[path = "layered_wave_heating/thermo_optical.rs"]
mod thermo_optical;

#[path = "layered_wave_heating/spatial_interference.rs"]
mod spatial_interference;

use fs_couple::layered_wave_heating::{
    LayeredWaveDrive, LayeredWaveHeating, LayeredWaveHeatingError,
};
use fs_em::plane_wave::{
    DielectricResponse, Polarization,
    layers::{DielectricLayer, IlluminationBudget, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS as EPS, ELECTRIC_PERMITTIVITY_PROPERTY as EPS_NAME,
        MAGNETIC_PERMEABILITY_DIMS as MU, MAGNETIC_PERMEABILITY_PROPERTY as MU_NAME,
        dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, SPECIFIC_HEAT_CAPACITY_DIMS,
        SPECIFIC_HEAT_CAPACITY_PROPERTY as CP_NAME, ScalarAdmissibility, ScalarPropertyRequirement,
        THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY as K_NAME,
        resolve_material_state_point,
    },
    thermal::{
        conduction::LayerConductionConfig,
        transport::{
            stack::{
                LayerStackSourceQuery,
                transient::{
                    BoundLayeredThermalBody, LayeredThermalIds, PreparedLayeredThermalBody,
                    TransientLayerSourceQuery,
                    implicit::{LayeredThermalSolveControls, LayeredThermalWorkspace},
                },
            },
            storage::{HeatStorageAmount, HeatStorageComponentQuery},
        },
    },
};
use fs_math::c64::C64;
use fs_qty::{
    Dims, QuantitySpec, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{DepositionBudget, DepositionSimulation, Event, StepConfig};

const E: f64 = 1e-11;
const M: f64 = 1e-6;
const OMEGA: f64 = 1e9;
const AREA: f64 = 0.6;
const RHO: f64 = 1000.0;
const CP: f64 = 200.0;
const K: f64 = 2.0;
const COHERENCE: IlluminationBudget = IlluminationBudget {
    maximum_modes: 16,
    maximum_layer_mode_pairs: 1024,
};
fn k() -> f64 {
    OMEGA * E.sqrt() * M.sqrt()
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn at(t: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), pressure)
        .unwrap()
}
fn card(name: &str, epsilon_ratio: f64, loss_ratio: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (EPS_NAME, EPS, E * (epsilon_ratio - loss_ratio)),
        (MU_NAME, MU, M),
        ("delta", EPS, 2.0 * E * loss_ratio),
        ("tau", Time::DIMS, 1.0 / OMEGA),
        (K_NAME, THERMAL_CONDUCTIVITY_DIMS, K),
        (CP_NAME, SPECIFIC_HEAT_CAPACITY_DIMS, CP),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        280.0,
                        360.0,
                    )
                    .with_quantity(
                        "pressure",
                        quantity(QuantityKind::Pressure),
                        90_000.0,
                        200_000.0,
                    ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored coherent wave-to-heat oracle".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "numerical film".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn state(card: &MaterialCard, debye: bool, pressure: f64) -> ResolvedMaterialStatePoint {
    let mut fields: Vec<(&str, Dims)> = vec![(EPS_NAME, EPS), (MU_NAME, MU)];
    if debye {
        fields.extend([("delta", EPS), ("tau", Time::DIMS)]);
    }
    let requirements: Vec<_> = fields
        .into_iter()
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap()
        })
        .collect();
    resolve_material_state_point(
        card,
        &at(300.0, pressure),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn debye(card: &MaterialCard, pressure: f64) -> ResolvedDebyeCapacitor {
    ResolvedDebyeCapacitor::try_from_material(
        state(card, true, pressure),
        EPS_NAME,
        &[DebyePoleProperties {
            permittivity_increment: "delta".into(),
            relaxation_time: "tau".into(),
        }],
        0.013,
        0.004,
    )
    .unwrap()
}
fn body(
    cards: &[&MaterialCard],
    widths: &[f64],
    counts: &[usize],
    area: f64,
) -> PreparedLayeredThermalBody {
    let lower = at(280.0, 100_000.0);
    let upper = at(360.0, 100_000.0);
    let conductivity = PropertyKey::new(K_NAME, THERMAL_CONDUCTIVITY_DIMS);
    let capacity = PropertyKey::new(CP_NAME, SPECIFIC_HEAT_CAPACITY_DIMS);
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let queries: Vec<_> = cards
        .iter()
        .zip(widths)
        .zip(counts)
        .enumerate()
        .map(|(i, ((&card, &width), &count))| {
            let mass = RHO * area * width;
            TransientLayerSourceQuery {
                conduction: LayerStackSourceQuery {
                    id: 100 + i as u64,
                    card,
                    conductivity: &conductivity,
                    thickness_m: width,
                    temperature_axis: "T",
                    lower: &lower,
                    upper: &upper,
                    selection,
                },
                capacity: HeatStorageComponentQuery {
                    id: 200 + i as u64,
                    card,
                    capacity: &capacity,
                    amount: HeatStorageAmount::MassKg(mass),
                    temperature_axis: "T",
                    lower: &lower,
                    upper: &upper,
                    selection,
                },
                mass_kg: mass,
                cell_count: count,
            }
        })
        .collect();
    PreparedLayeredThermalBody::try_new(&queries, &vec![0.0; cards.len() - 1], area, 300.0).unwrap()
}
fn ids() -> LayeredThermalIds {
    LayeredThermalIds {
        first_cell: 10,
        first_material: 30,
        first_link: 50,
        body_id: 70,
    }
}
fn conduction() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 96,
        absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-12,
    }
}
fn initial<'a>(
    bound: &'a BoundLayeredThermalBody<'_>,
    sources: &[f64],
) -> DepositionSimulation<'a> {
    let n = bound.prepared_body().cell_count();
    let mut events = bound
        .birth_events(0.0, &vec![300.0; n], conduction(), || Ok(()))
        .unwrap();
    let mut i = 0;
    for event in &mut events {
        if let Event::Birth(cell) = &mut event.event {
            cell.boundary.source_w = sources[i];
            i += 1;
        }
    }
    let materials: Vec<_> = bound.materials().into_iter().collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget {
            max_cells: n,
            max_links: n - 1,
            max_events: 2 * n - 1,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e} (tol {tolerance:e})"
    );
}

/// G3: native sublayer loss integrates the original film, in authored order.
/// Partial coherence is applied to complete fields before cell heat is formed.
#[test]
fn exact_sublayers_recover_original_films_and_full_coherence() {
    let exterior = state(&card("exterior", 1.0, 0.0), false, 100_000.0);
    let cards = [card("first", 3.0, 0.7), card("second", 1.8, 0.25)];
    let models = cards.each_ref().map(|c| debye(c, 100_000.0));
    let responses = models.each_ref().map(DielectricResponse::Debye);
    let widths = [0.31 / k(), 0.42 / k()];
    let native = body(&[&cards[0], &cards[1]], &widths, &[2, 3], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let owner = LayeredWaveHeating::try_new(
        &bound,
        &exterior,
        DielectricResponse::Lossless(&exterior),
        OMEGA,
        &responses,
        || false,
    )
    .unwrap();
    assert_eq!(owner.cell_ids(), bound.cell_ids());
    assert_eq!(owner.optical_layers().layer_count(), 5);
    for i in 0..5 {
        let original = if i < 2 {
            models[0].material()
        } else {
            models[1].material()
        };
        assert!(std::ptr::eq(
            owner.optical_layers().layer(i).unwrap().response.material(),
            original
        ));
    }
    let unsplit = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        OMEGA,
        &std::array::from_fn::<_, 2, _>(|i| DielectricLayer {
            response: responses[i],
            thickness_m: widths[i],
        }),
    )
    .unwrap();
    let factors = [
        [C64::new(0.8, 0.4), C64::new(-0.3, 0.7)],
        [C64::new(0.2, -0.6), C64::new(0.9, 0.1)],
    ];
    for polarization in [Polarization::Te, Polarization::Tm] {
        let prepared = owner.prepare(0.73, polarization, || false).unwrap();
        let reference = unsplit.prepare_illumination(0.73, polarization).unwrap();
        let mut expected = [0.0; 5];
        for input in factors {
            let mut actual = [0.0; 5];
            let mut whole = [0.0; 2];
            let power = prepared
                .cell_power_into(LayeredWaveDrive::Coherent(input), &mut actual, || false)
                .unwrap();
            let original = reference.evaluate(input, &mut whole).unwrap();
            near(actual[..2].iter().sum(), AREA * whole[0], 2e-17);
            near(actual[2..].iter().sum(), AREA * whole[1], 2e-17);
            near(power.absorbed_w_m2, original.power.absorbed_w_m2, 3e-17);
            for i in 0..5 {
                expected[i] += actual[i];
            }
        }
        let mut partial = [0.0; 5];
        prepared
            .cell_power_into(
                LayeredWaveDrive::Coherence {
                    factors: &factors,
                    budget: COHERENCE,
                },
                &mut partial,
                || false,
            )
            .unwrap();
        for i in 0..5 {
            near(partial[i], expected[i], 2e-18);
        }
        let inputs = [C64::new(0.8, 0.4), C64::new(-0.3, 0.7)];
        let diagonal = [[inputs[0], C64::ZERO], [C64::ZERO, inputs[1]]];
        let mut independent = [0.0; 5];
        prepared
            .cell_power_into(
                LayeredWaveDrive::Independent(inputs),
                &mut independent,
                || false,
            )
            .unwrap();
        prepared
            .cell_power_into(
                LayeredWaveDrive::Coherence {
                    factors: &diagonal,
                    budget: COHERENCE,
                },
                &mut partial,
                || false,
            )
            .unwrap();
        assert_eq!(independent, partial);
    }
}

/// Independent constant-capacity Fourier control-volume equations for one
/// homogeneous insulated slab. These constants are not recovered from the
/// production Jacobian or its temperature result.
fn thermal_reference(powers: &[f64], width: f64, area: f64, dt: f64) -> Vec<f64> {
    let n = powers.len();
    let dx = width / n as f64;
    let capacity = RHO * area * dx * CP;
    let gdt = area * K / dx * dt;
    let mut diagonal: Vec<_> = (0..n)
        .map(|i| capacity + gdt * (usize::from(i > 0) + usize::from(i + 1 < n)) as f64)
        .collect();
    let mut rhs: Vec<_> = powers.iter().map(|p| p * dt).collect();
    for i in 1..n {
        let factor = -gdt / diagonal[i - 1];
        diagonal[i] += factor * gdt;
        rhs[i] -= factor * rhs[i - 1];
    }
    let mut change = vec![0.0; n];
    change[n - 1] = rhs[n - 1] / diagonal[n - 1];
    for i in (0..n - 1).rev() {
        change[i] = (rhs[i] + gdt * change[i + 1]) / diagonal[i];
    }
    change.into_iter().map(|v| 300.0 + v).collect()
}

/// G1: analytical weak-loss standing-wave integrals and the actual native
/// enthalpy endpoint. Source phase changes WHICH cells heat, not only total W.
#[test]
fn standing_wave_cell_integrals_drive_native_temperature_and_conservation() {
    let exterior = state(&card("exterior", 1.0, 0.0), false, 100_000.0);
    let loss = 1e-7;
    let card = card("weak absorber", 1.0, loss);
    let model = debye(&card, 100_000.0);
    let h = core::f64::consts::PI / k();
    let native = body(&[&card], &[h], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let old_cells = initial.cells().clone();
    let owner = LayeredWaveHeating::try_new(
        &bound,
        &exterior,
        DielectricResponse::Lossless(&exterior),
        OMEGA,
        &[DielectricResponse::Debye(&model)],
        || false,
    )
    .unwrap();
    let prepared = owner.prepare(1.0, Polarization::Te, || false).unwrap();
    let amplitude = 1e6;
    let dt = 0.02;
    let mut profiles = Vec::new();
    for sign in [1.0, -1.0] {
        let input = [C64::from_re(amplitude), C64::from_re(sign * amplitude)];
        let mut powers = [0.0; 4];
        prepared
            .cell_power_into(LayeredWaveDrive::Coherent(input), &mut powers, || false)
            .unwrap();
        for (i, &power) in powers.iter().enumerate() {
            let dx = h / 4.0;
            let center = (i as f64 + 0.5) * dx;
            // E=a exp(-ikz)+b exp(-ik(h-z)), in the zero-loss matched limit.
            let integral = 2.0
                * amplitude
                * amplitude
                * (dx + sign * (k() * (h - 2.0 * center)).cos() * (k() * dx).sin() / k());
            let reference = AREA * 0.5 * OMEGA * E * loss * integral;
            near(power, reference, 2e-6 * reference);
        }
        let mut workspace = LayeredThermalWorkspace::new(4).unwrap();
        let result = prepared
            .advance_to(
                &initial,
                dt,
                LayeredWaveDrive::Coherent(input),
                conduction(),
                LayeredThermalSolveControls::default(),
                &mut workspace,
                || false,
            )
            .unwrap();
        let combined: Vec<_> = powers.iter().zip(authored).map(|(&p, s)| p + s).collect();
        let temperatures = thermal_reference(&combined, h, AREA, dt);
        for (i, &id) in owner.cell_ids().iter().enumerate() {
            near(
                result.thermal.temperature_k(id).unwrap(),
                temperatures[i],
                2e-10,
            );
            assert_eq!(
                result.thermal.cells()[&id].boundary,
                initial.cells()[&id].boundary
            );
            near(result.cell_energy_j[i].1, powers[i] * dt, 1e-13);
            assert_eq!(result.cell_energy_j[i].0, id);
        }
        near(
            result.optical.absorbed_j,
            powers.iter().sum::<f64>() * dt,
            1e-11,
        );
        near(
            result.heat_report.energy.external_heat_j,
            result.optical.absorbed_j + authored.iter().sum::<f64>() * dt,
            1e-8,
        );
        near(
            result.heat_report.energy.stored_energy_change_j,
            result.heat_report.energy.external_heat_j,
            2e-8,
        );
        assert!(result.optical.residual_j.abs() < 1e-6);
        assert_eq!(result.thermal.time_s(), dt);
        profiles.push(powers);
    }
    assert!(profiles[0][1] > 2.0 * profiles[0][0]);
    assert!(profiles[1][0] > 2.0 * profiles[1][1]);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &old_cells);
}

/// G3: reversed material order and source faces reverse the complete heat map.
/// Area scales W and mass together, without changing normal-flux conventions.
#[test]
fn reversed_cells_and_area_scaling_preserve_physical_source_orientation() {
    let exterior = state(&card("exterior", 1.0, 0.0), false, 100_000.0);
    let cards = [card("a", 2.7, 0.5), card("b", 1.4, 0.2)];
    let models = cards.each_ref().map(|c| debye(c, 100_000.0));
    let h = [0.45 / k(), 0.21 / k()];
    let input = [C64::new(1.0, 0.2), C64::new(-0.4, 0.7)];
    let mut results = Vec::new();
    for (reverse, scale) in [(false, 1.0), (true, 1.0), (false, 3.0)] {
        let order = if reverse { [1, 0] } else { [0, 1] };
        let widths = order.map(|i| h[i]);
        let counts = order.map(|i| [2, 3][i]);
        let native = body(&order.map(|i| &cards[i]), &widths, &counts, AREA * scale);
        let laws = native.thermal_laws(false).unwrap();
        let bound = native.bind(&laws, ids()).unwrap();
        let owner = LayeredWaveHeating::try_new(
            &bound,
            &exterior,
            DielectricResponse::Lossless(&exterior),
            OMEGA,
            &order.map(|i| DielectricResponse::Debye(&models[i])),
            || false,
        )
        .unwrap();
        let prepared = owner.prepare(0.8, Polarization::Tm, || false).unwrap();
        let mut heat = [0.0; 5];
        prepared
            .cell_power_into(
                LayeredWaveDrive::Coherent(if reverse { [input[1], input[0]] } else { input }),
                &mut heat,
                || false,
            )
            .unwrap();
        results.push(heat);
    }
    for i in 0..5 {
        near(results[1][4 - i], results[0][i], 2e-17);
        near(results[2][i], 3.0 * results[0][i], 2e-17);
    }
}

/// G1/G3: partial source covariance reaches actual cell energy, and an opaque
/// body receives its two face loads without inventing uniform interior heat.
#[test]
fn partial_coherence_and_opaque_limits_reach_native_heat() {
    let exterior = state(&card("exterior", 1.0, 0.0), false, 100_000.0);
    let card = card("opaque", 2.0, 1.0);
    let model = debye(&card, 100_000.0);
    let native = body(&[&card], &[1e5 / k()], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial(&bound, &[0.0; 4]);
    let owner = LayeredWaveHeating::try_new(
        &bound,
        &exterior,
        DielectricResponse::Lossless(&exterior),
        OMEGA,
        &[DielectricResponse::Debye(&model)],
        || false,
    )
    .unwrap();
    let prepared = owner.prepare(0.8, Polarization::Tm, || false).unwrap();
    let factors = [
        [C64::new(300.0, 100.0), C64::new(-200.0, 400.0)],
        [C64::new(100.0, -200.0), C64::new(300.0, 200.0)],
    ];
    let mut expected = [0.0; 4];
    for factor in factors {
        let mut heat = [0.0; 4];
        prepared
            .cell_power_into(LayeredWaveDrive::Coherent(factor), &mut heat, || false)
            .unwrap();
        for i in 0..4 {
            expected[i] += heat[i];
        }
    }
    assert!(expected[0] > 0.0 && expected[3] > 0.0);
    assert_eq!(expected[1], 0.0);
    assert_eq!(expected[2], 0.0);
    let dt = 0.1;
    let mut workspace = LayeredThermalWorkspace::new(4).unwrap();
    let result = prepared
        .advance_to(
            &initial,
            dt,
            LayeredWaveDrive::Coherence {
                factors: &factors,
                budget: COHERENCE,
            },
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    for (i, (_, q)) in result.cell_energy_j.iter().enumerate() {
        near(*q, expected[i] * dt, 2e-14);
    }
    near(
        result.heat_report.energy.external_heat_j,
        expected.iter().sum::<f64>() * dt,
        1e-8,
    );
    near(
        result.heat_report.energy.stored_energy_change_j,
        result.optical.absorbed_j,
        2e-8,
    );
}

/// G0/G4: source/geometry ownership, bounded work and late thermal cancellation
/// are refused transactionally; the same accepted initial state remains usable.
#[test]
fn source_work_and_timeline_refusals_preserve_outputs_and_accepted_history() {
    let exterior = state(&card("exterior", 1.0, 0.0), false, 100_000.0);
    let card = card("absorber", 2.0, 0.5);
    let model = debye(&card, 100_000.0);
    let wrong_pressure = debye(&card, 150_000.0);
    let changed_card = self::card("different source with same numbers", 2.0, 0.5);
    let changed = debye(&changed_card, 100_000.0);
    let native = body(&[&card], &[0.5 / k()], &[2], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    for response in [
        DielectricResponse::Debye(&wrong_pressure),
        DielectricResponse::Debye(&changed),
    ] {
        assert!(
            LayeredWaveHeating::try_new(
                &bound,
                &exterior,
                DielectricResponse::Lossless(&exterior),
                OMEGA,
                &[response],
                || false
            )
            .is_err()
        );
    }
    let large = body(&[&card], &[0.5 / k()], &[65], AREA);
    let large_laws = large.thermal_laws(false).unwrap();
    let large_bound = large.bind(&large_laws, ids()).unwrap();
    assert!(
        LayeredWaveHeating::try_new(
            &large_bound,
            &exterior,
            DielectricResponse::Lossless(&exterior),
            OMEGA,
            &[DielectricResponse::Debye(&model)],
            || false
        )
        .is_err()
    );
    let initial = initial(&bound, &[0.0; 2]);
    let before = initial.cells().clone();
    let owner = LayeredWaveHeating::try_new(
        &bound,
        &exterior,
        DielectricResponse::Lossless(&exterior),
        OMEGA,
        &[DielectricResponse::Debye(&model)],
        || false,
    )
    .unwrap();
    let prepared = owner.prepare(1.0, Polarization::Te, || false).unwrap();
    let input = [C64::new(300.0, 100.0), C64::new(-200.0, 400.0)];
    let mut sentinel = [123.0; 2];
    assert!(
        prepared
            .cell_power_into(
                LayeredWaveDrive::Coherence {
                    factors: &[input],
                    budget: IlluminationBudget {
                        maximum_modes: 0,
                        maximum_layer_mode_pairs: 0
                    },
                },
                &mut sentinel,
                || false
            )
            .is_err()
    );
    assert_eq!(sentinel, [123.0; 2]);
    let mut polls = 0;
    assert!(matches!(
        prepared.cell_power_into(
            LayeredWaveDrive::Coherence {
                factors: &[input],
                budget: COHERENCE,
            },
            &mut sentinel,
            || {
                polls += 1;
                polls == 2
            }
        ),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    assert_eq!(sentinel, [123.0; 2]);
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let mut total_polls = 0;
    let accepted = prepared
        .advance_to(
            &initial,
            0.01,
            LayeredWaveDrive::Coherent(input),
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
    let mut polls = 0;
    let refused = prepared.advance_to(
        &initial,
        0.01,
        LayeredWaveDrive::Coherent(input),
        conduction(),
        LayeredThermalSolveControls::default(),
        &mut workspace,
        || {
            polls += 1;
            polls == total_polls
        },
    );
    assert!(matches!(refused, Err(LayeredWaveHeatingError::Cancelled)));
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
    let retry = prepared
        .advance_to(
            &initial,
            0.01,
            LayeredWaveDrive::Coherent(input),
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(retry.thermal.cells(), accepted.thermal.cells());
    assert_eq!(retry.cell_energy_j, accepted.cell_energy_j);
    // Numerically identical law objects are not the original bound timeline.
    let other_laws = native.thermal_laws(false).unwrap();
    let other_bound = native.bind(&other_laws, ids()).unwrap();
    let other_initial = self::initial(&other_bound, &[0.0; 2]);
    assert!(
        prepared
            .advance_to(
                &other_initial,
                0.01,
                LayeredWaveDrive::Coherent(input),
                conduction(),
                LayeredThermalSolveControls::default(),
                &mut workspace,
                || false
            )
            .is_err()
    );
}
