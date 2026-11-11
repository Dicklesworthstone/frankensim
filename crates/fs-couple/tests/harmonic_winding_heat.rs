#![cfg(feature = "harmonic-winding-heating")]
//! G1/G3/G4: actual coherent AC fields, spatial heat and native caloric receivers.
//! The manufactured two-square modes have an independent scalar-field oracle.

#[path = "../../fs-em/tests/nonlinear_winding_circuit/common.rs"]
mod field;

use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve,
        magnetic::{MagneticRegion, MagneticThermalNetwork},
        material::{ConductionLink, MaterialCell},
    },
    material_circuit::MaterialCircuitSolve,
};
use fs_em::{
    planar::{
        eddy::{
            EddyConductor,
            harmonic::{
                HarmonicWindingControls, PlanarHarmonicWindings,
                exposure::{HarmonicWindingExposureLimits, HarmonicWindingTone},
            },
        },
        nonlinear::PlanarNonlinearMagnetostatics,
    },
    winding::{
        WindingConfig,
        incremental::{IncrementalWindingConfig, PlanarIncrementalWindingBasis},
    },
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, PropertyKey, PropertyValue, QueryPoint,
    SelectionPolicy,
};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_math::c64::C64;
use fs_qty::Density;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const DEPTH: f64 = 0.25;
const SIGMA: f64 = 18.0;
const CAPACITY: f64 = 6.0 * DEPTH;
const WINDING: CellId = CellId(7);
const EDDY: CellId = CellId(9);
const MATERIAL: MaterialId = MaterialId(2);
const OMEGA: [f64; 2] = [13.0, 17.0];
const CURRENT: [[C64; 2]; 2] = [
    [C64::new(0.2, -0.07), C64::new(-0.03, 0.09)],
    [C64::new(-0.13, 0.08), C64::new(0.04, 0.06)],
];

fn card(sigma: f64, varying: bool) -> MaterialCard {
    let original = field::source(false, sigma);
    if varying {
        return original;
    }
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut claim = source.clone();
        for basis in [
            ElectricalPropertyBasis::Resistivity,
            ElectricalPropertyBasis::Conductivity,
        ] {
            if claim.key.name() == basis.property_name() {
                claim.value = PropertyValue::Scalar {
                    value: if basis == ElectricalPropertyBasis::Resistivity {
                        sigma.recip()
                    } else {
                        sigma
                    },
                    dims: basis.dimensions(),
                };
                claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &field::material::point(200.0),
            upper: &field::material::point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 1000,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 1e-7,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-8,
            power_tolerance_w: 1e-8,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-13,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn exposure_limits() -> HarmonicWindingExposureLimits {
    HarmonicWindingExposureLimits {
        maximum_lines: 8,
        maximum_pair_products: 100_000,
        maximum_workspace_bytes: 1_000_000,
        power_tolerance_w: 1e-8,
    }
}
fn basis(problem: &PlanarNonlinearMagnetostatics) -> PlanarIncrementalWindingBasis<'_> {
    PlanarIncrementalWindingBasis::new(
        problem,
        vec![
            vec![24.0, 0.0, 0.0, 0.0, -24.0, 0.0, 0.0, 0.0],
            vec![0.0, -12.0, 0.0, 0.0, 0.0, 12.0, 0.0, 0.0],
        ],
        &field::boundary(),
        vec![0.0; 2],
        vec![0.02, 0.03],
        &[0.0; 8],
        IncrementalWindingConfig {
            winding: WindingConfig {
                max_windings: 4,
                max_basis_values: 10_000,
                pairing_tolerance_h: 1e-10,
                current_balance_tolerance_turns: 1e-10,
                inverse_residual_tolerance: 1e-10,
            },
            minimum_smooth_field_margin_t: 1e-8,
        },
        field::material::controls().nonlinear,
        &mut || false,
    )
    .unwrap()
}
fn model<'b, 'm>(
    basis: &'b PlanarIncrementalWindingBasis<'m>,
    card: &MaterialCard,
) -> PlanarHarmonicWindings<'b, 'm> {
    let conductor = field::material::conductor(card, 300.0, ElectricalPropertyBasis::Resistivity);
    let eddies = [2, 3, 6, 7]
        .into_iter()
        .map(|triangle| EddyConductor {
            triangle,
            region: if triangle < 4 { 7 } else { 9 },
            conductor: conductor.clone(),
        })
        .collect();
    let copper = field::conductors(card, 300.0)
        .into_iter()
        .filter(|s| [0, 1, 4, 5].contains(&s.triangle))
        .collect();
    PlanarHarmonicWindings::new(
        basis,
        eddies,
        &[0, 1, 4, 5],
        copper,
        HarmonicWindingControls {
            relative_tolerance: 1e-11,
            absolute_residual_a_m: 1e-12,
            restart: 8,
            max_cycles: 50,
            max_krylov_scalars: 10_000,
            max_field_scalars: 10_000,
            impedance_tolerance_ohm: 1e-8,
            power_tolerance_w: 1e-8,
            current_balance_tolerance_a: 1e-10,
            max_cell_skin_ratio: 1000.0,
        },
        &mut || false,
    )
    .unwrap()
}
fn boundary() -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    }
}
fn initial(
    curve: &ResolvedHeatCapacityCurve,
    volume: f64,
    linked: bool,
    event: Option<f64>,
) -> DepositionSimulation<'_> {
    let mut events: Vec<_> = [WINDING, EDDY]
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MATERIAL,
                mass_kg: 2.0 * volume,
                specific_h_j_kg: 300.0,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 1.0,
                },
                boundary: boundary(),
            }),
        })
        .collect();
    if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 4,
                a: WINDING,
                b: EDDY,
                area_m2: 1.0,
                conductance_w_k: 123.0,
            }),
        });
    }
    if let Some(time_s) = event {
        events.push(TimedEvent {
            time_s,
            event: Event::Boundary(WINDING, boundary()),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MATERIAL, curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 5,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn cells<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
    key: &'a PropertyKey,
) -> [MaterialCell<'a>; 2] {
    [WINDING, EDDY].map(|cell| MaterialCell {
        cell,
        material_id: MATERIAL,
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((key, MaterialPropertySelection::SingleClaimOnly)),
    })
}
fn routes() -> Vec<MagneticRegion> {
    (0..8)
        .map(|triangle| MagneticRegion::Cell {
            triangle,
            cell: if [0, 1, 4, 5].contains(&triangle) {
                WINDING
            } else {
                EDDY
            },
        })
        .collect()
}
fn signal(amplitude: C64, omega: f64, time: f64) -> f64 {
    let (sin, cos) = (omega * time).sin_cos();
    amplitude.re * cos + amplitude.im * sin
}
fn integral(start: f64, end: f64, mut value: impl FnMut(f64) -> f64) -> f64 {
    let n = 4096;
    let dt = (end - start) / f64::from(n);
    let mut sum = 0.0;
    for i in 0..=n {
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * value(start + f64::from(i) * dt);
    }
    sum * dt / 3.0
}
// Independent one-hat elimination: K=4*depth, C=sigma*depth/36,
// F I=2*depth*(I0-I1/2), with the opposite hat in the return square.
fn reference_heat(start: f64, end: f64) -> [f64; 2] {
    let a = [0, 1].map(|j| {
        (CURRENT[j][0] - CURRENT[j][1].scale(0.5)) / C64::new(2.0, -OMEGA[j] * SIGMA / 72.0)
    });
    let copper = integral(start, end, |time| {
        [288.0, 72.0]
            .into_iter()
            .enumerate()
            .map(|(k, r)| {
                let current: f64 = (0..2).map(|j| signal(CURRENT[j][k], OMEGA[j], time)).sum();
                r * DEPTH / SIGMA * current * current
            })
            .sum()
    });
    let eddy = integral(start, end, |time| {
        let derivative: f64 = (0..2)
            .map(|j| {
                signal(
                    C64::new(OMEGA[j] * a[j].im, -OMEGA[j] * a[j].re),
                    OMEGA[j],
                    time,
                )
            })
            .sum();
        SIGMA * DEPTH / 18.0 * derivative * derivative
    });
    [copper, eddy]
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-8 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn coherent_spatial_heat_drives_actual_temperatures_and_native_conduction() {
    let card = card(SIGMA, false);
    let curve = caloric(&card);
    let context = field::material::point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&card, &curve, &context, &key);
    let regions = routes();
    let problem = field::problem(&card, DEPTH);
    let basis = basis(&problem);
    let model = model(&basis, &card);
    let frequencies = OMEGA.map(|w| model.at_frequency(w, &mut || false).unwrap());
    let tones = [0, 1].map(|j| HarmonicWindingTone {
        frequency: &frequencies[j],
        currents_a: &CURRENT[j],
    });
    let exposure = model
        .prepare_exposure(&tones, exposure_limits(), &mut || false)
        .unwrap();
    for linked in [false, true] {
        let links: Vec<_> = if linked {
            vec![ConductionLink {
                link: 4,
                distance_a_m: 0.4,
                distance_b_m: 0.6,
            }]
        } else {
            vec![]
        };
        let layout = MagneticThermalNetwork {
            cells: &cells,
            links: &links,
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        for duration in [1e-4, 0.081, 0.31] {
            let initial = initial(&curve, DEPTH, linked, None);
            let candidate = layout
                .step_harmonic_winding_candidate(
                    &initial,
                    duration,
                    &exposure,
                    controls(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            let heat = reference_heat(0.0, duration);
            let average = (heat[0] + heat[1]) / (2.0 * CAPACITY);
            let difference =
                (heat[0] - heat[1]) / (CAPACITY + if linked { 2.0 * duration } else { 0.0 });
            near(
                candidate.thermal().temperature_k(WINDING).unwrap() - 300.0,
                average + 0.5 * difference,
            );
            near(
                candidate.thermal().temperature_k(EDDY).unwrap() - 300.0,
                average - 0.5 * difference,
            );
            let report = candidate.report();
            near(report.deposited_joule_j, heat[0] + heat[1]);
            near(report.thermal.stored_energy_change_j, heat[0] + heat[1]);
            near(report.electrical.dissipated_heat_j, heat[0] + heat[1]);
            near(report.non_electrical_heat_j, 0.0);
            near(
                report.electrical.supplied_work_j,
                report.thermal.stored_energy_change_j + report.electrical.magnetic_energy_change_j,
            );
            assert!(report.total_energy_residual_j.abs() < 1e-7);
            assert_eq!(report.source.magnetic.len(), 8);
            near(initial.temperature_k(WINDING).unwrap(), 300.0);
        }
    }
}

#[test]
fn consecutive_intervals_retain_source_phase_and_exact_adiabatic_heat() {
    let card = card(SIGMA, false);
    let curve = caloric(&card);
    let context = field::material::point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&card, &curve, &context, &key);
    let regions = routes();
    let layout = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let problem = field::problem(&card, DEPTH);
    let basis = basis(&problem);
    let model = model(&basis, &card);
    let frequencies = OMEGA.map(|w| model.at_frequency(w, &mut || false).unwrap());
    let tones = [0, 1].map(|j| HarmonicWindingTone {
        frequency: &frequencies[j],
        currents_a: &CURRENT[j],
    });
    let exposure = model
        .prepare_exposure(&tones, exposure_limits(), &mut || false)
        .unwrap();
    let initial = initial(&curve, DEPTH, false, None);
    let whole = layout
        .step_harmonic_winding_candidate(&initial, 0.31, &exposure, controls(), &[0.0; 2], || false)
        .unwrap();
    let first = layout
        .step_harmonic_winding_candidate(&initial, 0.07, &exposure, controls(), &[0.0; 2], || false)
        .unwrap();
    let second = layout
        .step_harmonic_winding_candidate(
            first.thermal(),
            0.31,
            &exposure,
            controls(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    for cell in [WINDING, EDDY] {
        near(
            second.thermal().temperature_k(cell).unwrap(),
            whole.thermal().temperature_k(cell).unwrap(),
        );
    }
    near(
        first.report().deposited_joule_j + second.report().deposited_joule_j,
        whole.report().deposited_joule_j,
    );
    near(
        first.report().electrical.supplied_work_j + second.report().electrical.supplied_work_j,
        whole.report().electrical.supplied_work_j,
    );
    assert_eq!(second.report().electrical.window.start_s(), 0.07);
    near(
        second.report().deposited_joule_j,
        reference_heat(0.07, 0.31).into_iter().sum(),
    );
}

#[test]
fn equal_frequency_destructive_interference_deposits_zero_heat() {
    let card = card(SIGMA, false);
    let curve = caloric(&card);
    let context = field::material::point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&card, &curve, &context, &key);
    let regions = routes();
    let layout = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let problem = field::problem(&card, DEPTH);
    let basis = basis(&problem);
    let model = model(&basis, &card);
    let frequency = model.at_frequency(OMEGA[0], &mut || false).unwrap();
    let reversed = CURRENT[0].map(|i| i.scale(-1.0));
    let tones = [
        HarmonicWindingTone {
            frequency: &frequency,
            currents_a: &CURRENT[0],
        },
        HarmonicWindingTone {
            frequency: &frequency,
            currents_a: &reversed,
        },
    ];
    let exposure = model
        .prepare_exposure(&tones, exposure_limits(), &mut || false)
        .unwrap();
    let initial = initial(&curve, DEPTH, false, None);
    let candidate = layout
        .step_harmonic_winding_candidate(&initial, 0.1, &exposure, controls(), &[0.0; 2], || false)
        .unwrap();
    assert_eq!(candidate.report().deposited_joule_j, 0.0);
    assert_eq!(candidate.report().electrical.supplied_work_j, 0.0);
    for cell in [WINDING, EDDY] {
        assert_eq!(
            candidate.thermal().cells()[&cell].specific_h_j_kg.to_bits(),
            initial.cells()[&cell].specific_h_j_kg.to_bits()
        );
    }
}

#[test]
fn changing_laws_foreign_receivers_geometry_events_and_cancellation_refuse_transactionally() {
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let context = field::material::point(300.0);
    for varying in [false, true] {
        let card = card(SIGMA, varying);
        let curve = caloric(&card);
        let cells = cells(&card, &curve, &context, &key);
        let regions = routes();
        let layout = MagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let problem = field::problem(&card, DEPTH);
        let basis = basis(&problem);
        let model = model(&basis, &card);
        let frequency = model.at_frequency(OMEGA[0], &mut || false).unwrap();
        let exposure = model
            .prepare_exposure(
                &[HarmonicWindingTone {
                    frequency: &frequency,
                    currents_a: &CURRENT[0],
                }],
                exposure_limits(),
                &mut || false,
            )
            .unwrap();
        let initial = initial(&curve, DEPTH, false, None);
        if varying {
            assert!(matches!(
                layout.step_harmonic_winding_candidate(
                    &initial,
                    0.1,
                    &exposure,
                    controls(),
                    &[0.0; 2],
                    || false
                ),
                Err(CircuitEnthalpyError::Invalid(_))
            ));
            near(initial.temperature_k(WINDING).unwrap(), 300.0);
            continue;
        }
        let mut calls = 0;
        layout
            .step_harmonic_winding_candidate(
                &initial,
                0.1,
                &exposure,
                controls(),
                &[0.0; 2],
                || {
                    calls += 1;
                    false
                },
            )
            .unwrap();
        for stop in [1, calls / 2] {
            let mut call = 0;
            assert!(matches!(
                layout.step_harmonic_winding_candidate(
                    &initial,
                    0.1,
                    &exposure,
                    controls(),
                    &[0.0; 2],
                    || {
                        call += 1;
                        call >= stop
                    }
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
        }
        let invalid_volume = self::initial(&curve, DEPTH * 1.01, false, None);
        assert!(
            layout
                .step_harmonic_winding_candidate(
                    &invalid_volume,
                    0.1,
                    &exposure,
                    controls(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
        let event = self::initial(&curve, DEPTH, false, Some(0.05));
        assert!(
            layout
                .step_harmonic_winding_candidate(
                    &event,
                    0.1,
                    &exposure,
                    controls(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
        let mut budget = controls();
        budget.maximum_circuit_trials = 1;
        assert!(
            layout
                .step_harmonic_winding_candidate(
                    &initial,
                    0.1,
                    &exposure,
                    budget,
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
        let foreign = self::card(SIGMA * 2.0, false);
        let foreign_curve = caloric(&foreign);
        let foreign_cells = self::cells(&foreign, &foreign_curve, &context, &key);
        let foreign_thermal = self::initial(&foreign_curve, DEPTH, false, None);
        let foreign_layout = MagneticThermalNetwork {
            cells: &foreign_cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        assert!(
            foreign_layout
                .step_harmonic_winding_candidate(
                    &foreign_thermal,
                    0.1,
                    &exposure,
                    controls(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
        assert_eq!(initial.time_s(), 0.0);
        near(initial.temperature_k(WINDING).unwrap(), 300.0);
        near(initial.temperature_k(EDDY).unwrap(), 300.0);
    }
}
