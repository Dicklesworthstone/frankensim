//! Manufactured coupled-frame/caloric regressions, not measured material data.
//! References use physical consistent inertia, section integrals, continuous
//! damped motion and independently solved finite-cell heat equations.

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    material::{ConductionLink, MaterialCell, MaterialNetwork},
};
use fs_couple::frame_vibration::{
    FrameVibrationBasis, FrameVibrationBudget, FrameVibrationComponent,
    thermal::{FrameLayerCell, FrameThermalSolve},
};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig, PreparedModalBudget};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_ports::PortedModule;
use fs_qty::{Density, DynViscosity, Pressure, QuantitySpec, Temperature};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, FrameViscousSection, Geometry};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const RATE: u32 = 256;
const LENGTH: f64 = 1.0;
const RADIUS: [f64; 2] = [0.06, 0.1];
const RHO: f64 = 3.0;
const YOUNG: f64 = 1200.0;
const CELLS: [CellId; 2] = [CellId(11), CellId(27)];
const DESTINATIONS: [FrameLayerCell; 2] = [
    FrameLayerCell {
        element: 0,
        layer: 0,
        cell: CELLS[0],
    },
    FrameLayerCell {
        element: 0,
        layer: 1,
        cell: CELLS[1],
    },
];

#[derive(Clone, Copy)]
struct Spec {
    axial: f64,
    bending: f64,
    cp: f64,
    mechanical_upper: f64,
    curved_eta: bool,
}
fn specimens() -> [Spec; 2] {
    [
        Spec {
            axial: 2.0,
            bending: 10.0,
            cp: 0.5,
            mechanical_upper: 600.0,
            curved_eta: false,
        },
        Spec {
            axial: 8.0,
            bending: 40.0,
            cp: 1.0,
            mechanical_upper: 600.0,
            curved_eta: false,
        },
    ]
}
fn areas() -> [f64; 2] {
    let a = core::f64::consts::PI * RADIUS[0].powi(2);
    [a, core::f64::consts::PI * RADIUS[1].powi(2) - a]
}
fn moments() -> [f64; 2] {
    let a = core::f64::consts::PI * RADIUS[0].powi(4) / 4.0;
    [a, core::f64::consts::PI * RADIUS[1].powi(4) / 4.0 - a]
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with("omega", 2048.0)
        .unwrap()
}
fn material(index: usize, spec: Spec) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", RHO, Density::DIMS),
        ("young_modulus", YOUNG, Pressure::DIMS),
        ("axial_eta", spec.axial, DynViscosity::DIMS),
        ("bending_eta", spec.bending, DynViscosity::DIMS),
        (
            "specific_heat_capacity",
            spec.cp,
            SPECIFIC_HEAT_CAPACITY_DIMS,
        ),
        (
            "thermal_conductivity",
            0.1 * (index + 1) as f64,
            THERMAL_CONDUCTIVITY_DIMS,
        ),
    ] {
        let viscosity = name.ends_with("_eta");
        let mut validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            200.0,
            if viscosity {
                spec.mechanical_upper
            } else {
                600.0
            },
        );
        if viscosity {
            validity = validity.with("omega", 0.0, 4096.0);
        }
        let curved = viscosity && spec.curved_eta;
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if curved {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(200.0, 0.5 * value), (600.0, 2.5 * value)],
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value, dims }
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if curved {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured layered frame heat fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("manufactured layer {index}"),
            phase: "solid".into(),
            process: "linear frame calorimetry fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(600.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-7,
    )
    .unwrap();
    (card, curve)
}
fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 2,
        max_connections: 3,
        max_setup_terms: 1024,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: 100.0,
        maximum_abs_connection_force_n: 1000.0,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}
fn prepared() -> PreparedModalBudget {
    PreparedModalBudget {
        max_modes: 2,
        max_dense_bytes: 8 * 1024 * 1024,
        max_multiply_terms: 1_000_000_000,
    }
}
fn solve() -> FrameThermalSolve {
    FrameThermalSolve {
        maximum_cells: 2,
        maximum_material_trials: 128,
        maximum_newton_iterations: 12,
        difference_step_j: 1e-7,
        cell_energy_tolerance_j: 1e-10,
        thermal_energy_tolerance_j: 1e-10,
        total_energy_tolerance_j: 2e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn fixture<R>(
    specs: [Spec; 2],
    bending: bool,
    conduction: bool,
    authored_zeta: f64,
    run: impl FnOnce(FrameVibrationComponent, DepositionSimulation<'_>, &MaterialNetwork<'_>) -> R,
) -> R {
    let sources = [material(0, specs[0]), material(1, specs[1])];
    let context = point(300.0);
    let requirement = [
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        ("axial_eta", DynViscosity::DIMS),
        ("bending_eta", DynViscosity::DIMS),
    ]
    .map(|(key, dims)| {
        ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap()
    });
    let layers = sources
        .iter()
        .enumerate()
        .map(|(i, (card, _))| ConcentricMaterialLayer {
            outer_radius_m: RADIUS[i],
            material: resolve_material_state_point(
                card,
                &context,
                &requirement,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        })
        .collect();
    let section = FrameViscousSection::try_new(
        ConcentricElasticSection::try_new(0.0, layers).unwrap(),
        "axial_eta",
        "bending_eta",
    )
    .unwrap();
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(LENGTH, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(
        1,
        if bending {
            [true, true, false]
        } else {
            [false, true, true]
        },
    );
    frame
        .add_material_viscous_element(0, 1, 5, section)
        .unwrap();
    let basis = FrameVibrationBasis::from_frame(
        &frame,
        1,
        FrameVibrationBudget {
            max_frame_dofs: 6,
            max_dense_entries: 36,
        },
        || false,
    )
    .unwrap();
    let mut component = basis
        .coupled_component(
            RATE,
            &[authored_zeta],
            ModalAcousticTimeBudget::audible_reference(),
            18,
            || false,
        )
        .unwrap();
    let dof = if bending { 5 } else { 3 };
    component
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.01 / basis.modes()[0].displacement[dof],
            velocity_m_sqrt_kg_per_s: 0.03 / basis.modes()[0].displacement[dof],
        }])
        .unwrap();
    let a = areas();
    let mut events: Vec<_> = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELLS[i],
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(i as u64),
                mass_kg: RHO * a[i] * LENGTH,
                specific_h_j_kg: 100.0 * specs[i].cp,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [LENGTH, 0.0, 0.0],
                    volume_m3: a[i] * LENGTH,
                    surface_area_m2: 2.0
                        * core::f64::consts::PI
                        * (RADIUS[i] + if i == 0 { 0.0 } else { RADIUS[0] })
                        * LENGTH
                        + 2.0 * a[i],
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    if conduction {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 1,
                a: CELLS[0],
                b: CELLS[1],
                area_m2: 2.0 * core::f64::consts::PI * RADIUS[0] * LENGTH,
                conductance_w_k: 0.0,
            }),
        });
    }
    let laws = [
        (MaterialId(0), sources[0].1.curve() as &dyn ThermalLaw),
        (MaterialId(1), sources[1].1.curve() as &dyn ThermalLaw),
    ];
    let thermal = DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1024,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells: Vec<_> = sources
        .iter()
        .enumerate()
        .map(|(i, (card, curve))| MaterialCell {
            cell: CELLS[i],
            material_id: MaterialId(i as u64),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: conduction
                .then_some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect();
    let links = if conduction {
        vec![ConductionLink {
            link: 1,
            distance_a_m: 0.03,
            distance_b_m: 0.02,
        }]
    } else {
        vec![]
    };
    let network = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &links,
    };
    run(component, thermal, &network)
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance={tolerance:e}"
    );
}
fn oscillator(mass: f64, stiffness: f64, damping: f64, time: f64) -> (f64, f64, f64) {
    let a = damping / (2.0 * mass);
    let w = (stiffness / mass - a * a).sqrt();
    let (s, c) = (w * time).sin_cos();
    let b = (0.03 + a * 0.01) / w;
    let q = (-a * time).exp() * (0.01 * c + b * s);
    let v = (-a * time).exp() * (-0.01 * w * s + b * w * c) - a * q;
    (q, v, 0.5 * mass * v * v + 0.5 * stiffness * q * q)
}

#[test]
fn layered_axial_heat_matches_continuous_work_and_finite_caloric_storage() {
    let specs = specimens();
    let area = areas();
    let mass = RHO * area.iter().sum::<f64>() * LENGTH / 3.0;
    let stiffness = YOUNG * area.iter().sum::<f64>() / LENGTH;
    let viscous: [f64; 2] = std::array::from_fn(|i| specs[i].axial * area[i] / LENGTH);
    let damping: f64 = viscous.iter().sum();
    let initial = oscillator(mass, stiffness, damping, 0.0).2;
    for preparation in [None, Some(prepared())] {
        fixture(
            specs,
            false,
            false,
            0.0,
            |component, mut thermal, network| {
                let gate = CancelGate::new_clock_free();
                let mut owner = component
                    .into_thermal_system(
                        &thermal,
                        network,
                        &DESTINATIONS,
                        config(),
                        preparation,
                        &gate,
                    )
                    .unwrap();
                let initial_thermal = thermal.energy_j();
                let (mut heat, mut entropy) = ([0.0; 2], 0.0);
                for sample in 1..=64 {
                    let report = owner
                        .step(&mut thermal, &[0.0], solve(), &[0.0; 2], &gate)
                        .unwrap();
                    let expected =
                        oscillator(mass, stiffness, damping, sample as f64 / f64::from(RATE));
                    let shape = owner.basis().modes()[0].displacement[3];
                    let state = owner.mechanical().components()[0].states()[0];
                    near(shape * state.displacement_m_sqrt_kg, expected.0, 2e-10);
                    near(shape * state.velocity_m_sqrt_kg_per_s, expected.1, 2e-9);
                    for i in 0..2 {
                        assert_eq!(report.layer_heat[i].destination, DESTINATIONS[i]);
                        assert_eq!(report.layer_heat[i].bending_heat_j, 0.0);
                        heat[i] += report.layer_heat[i].total_j();
                        near(
                            heat[i],
                            (initial - expected.2) * viscous[i] / damping,
                            1e-10,
                        );
                        near(
                            thermal.temperature_k(CELLS[i]).unwrap(),
                            300.0 + heat[i] / (RHO * area[i] * LENGTH * specs[i].cp),
                            1e-8,
                        );
                        let point = report.sections[0].1.layers()[i].material.query_point();
                        assert_eq!(
                            point.iter().find(|(key, _)| key == "T").unwrap().1,
                            thermal.temperature_k(CELLS[i]).unwrap()
                        );
                    }
                    entropy += report.entropy.thermal_storage_change_j_k;
                    assert!(report.entropy.material_loss_generation_j_k >= 0.0);
                    assert!(!report.entropy.has_unclassified_sources);
                    assert_eq!(report.external_mechanical_loss_j, 0.0);
                    near(report.combined_energy_residual_j, 0.0, 2e-10);
                }
                near(
                    owner.mechanical().total_energy_j().unwrap() + thermal.energy_j()
                        - initial_thermal,
                    initial,
                    2e-10,
                );
                let expected_entropy: f64 = (0..2)
                    .map(|i| {
                        RHO * area[i]
                            * LENGTH
                            * specs[i].cp
                            * (thermal.temperature_k(CELLS[i]).unwrap() / 300.0).ln()
                    })
                    .sum();
                near(entropy, expected_entropy, 2e-14);
            },
        );
    }
}

#[test]
fn bending_heat_uses_original_eta_i_and_both_curvature_ports() {
    let specs = specimens();
    fixture(
        specs,
        true,
        false,
        0.0,
        |component, mut thermal, network| {
            let gate = CancelGate::new_clock_free();
            let mut owner = component
                .into_thermal_system(
                    &thermal,
                    network,
                    &DESTINATIONS,
                    config(),
                    Some(prepared()),
                    &gate,
                )
                .unwrap();
            let area = areas();
            let inertia = moments();
            // Only the end rotation is free: exact Hermite consistent inertia,
            // bending stiffness and viscous moment coefficient, all in physical SI.
            let mass = RHO * area.iter().sum::<f64>() * LENGTH.powi(3) / 105.0;
            let stiffness = 4.0 * YOUNG * inertia.iter().sum::<f64>() / LENGTH;
            let local: [f64; 2] =
                std::array::from_fn(|i| 4.0 * specs[i].bending * inertia[i] / LENGTH);
            let damping: f64 = local.iter().sum();
            let initial = oscillator(mass, stiffness, damping, 0.0).2;
            let mut heat = [0.0; 2];
            for sample in 1..=64 {
                let report = owner
                    .step(&mut thermal, &[0.0], solve(), &[0.0; 2], &gate)
                    .unwrap();
                let expected =
                    oscillator(mass, stiffness, damping, sample as f64 / f64::from(RATE));
                let shape = owner.basis().modes()[0].displacement[5];
                let state = owner.mechanical().components()[0].states()[0];
                near(shape * state.displacement_m_sqrt_kg, expected.0, 2e-10);
                near(shape * state.velocity_m_sqrt_kg_per_s, expected.1, 2e-9);
                let ports = owner.mechanical().last_window_transfers().unwrap();
                near(
                    ports[1].energy_into().value(),
                    3.0 * ports[2].energy_into().value(),
                    1e-13,
                );
                for i in 0..2 {
                    assert_eq!(report.layer_heat[i].axial_heat_j, 0.0);
                    heat[i] += report.layer_heat[i].bending_heat_j;
                    near(heat[i], (initial - expected.2) * local[i] / damping, 1e-10);
                    near(
                        thermal.temperature_k(CELLS[i]).unwrap(),
                        300.0 + heat[i] / (RHO * area[i] * LENGTH * specs[i].cp),
                        1e-8,
                    );
                }
            }
            // eta*A would be a materially different, physically incorrect split.
            let wrong_share = specs[0].bending * area[0]
                / (specs[0].bending * area[0] + specs[1].bending * area[1]);
            assert!((heat[0] / heat.iter().sum::<f64>() - wrong_share).abs() > 0.05);
        },
    );
}

#[test]
fn data_only_heat_capacity_swap_and_sourced_conduction_share_the_original_trajectory() {
    let area = areas();
    let mut endpoints = Vec::new();
    for cp_scale in [1.0, 2.0] {
        let specs = specimens().map(|mut s| {
            s.cp *= cp_scale;
            s
        });
        fixture(
            specs,
            false,
            true,
            0.0,
            |component, mut thermal, network| {
                let gate = CancelGate::new_clock_free();
                let mut owner = component
                    .into_thermal_system(
                        &thermal,
                        network,
                        &DESTINATIONS,
                        config(),
                        Some(prepared()),
                        &gate,
                    )
                    .unwrap();
                let capacity: [f64; 2] =
                    std::array::from_fn(|i| RHO * area[i] * LENGTH * specs[i].cp);
                let local: [f64; 2] = std::array::from_fn(|i| specs[i].axial * area[i] / LENGTH);
                let damping: f64 = local.iter().sum();
                let mass = RHO * area.iter().sum::<f64>() * LENGTH / 3.0;
                let stiffness = YOUNG * area.iter().sum::<f64>() / LENGTH;
                let conductance =
                    2.0 * core::f64::consts::PI * RADIUS[0] * LENGTH / (0.03 / 0.1 + 0.02 / 0.2);
                let dt = 1.0 / f64::from(RATE);
                let g = dt * conductance;
                let determinant = (capacity[0] + g) * (capacity[1] + g) - g * g;
                let mut excess = [0.0; 2];
                let mut energy = oscillator(mass, stiffness, damping, 0.0).2;
                let mut entropy = 0.0;
                for sample in 1..=64 {
                    let next_energy = oscillator(mass, stiffness, damping, sample as f64 * dt).2;
                    let heat = energy - next_energy;
                    energy = next_energy;
                    let rhs: [f64; 2] = std::array::from_fn(|i| {
                        capacity[i] * excess[i] + heat * local[i] / damping
                    });
                    excess = [
                        ((capacity[1] + g) * rhs[0] + g * rhs[1]) / determinant,
                        (g * rhs[0] + (capacity[0] + g) * rhs[1]) / determinant,
                    ];
                    let report = owner
                        .step(&mut thermal, &[0.0], solve(), &[0.0; 2], &gate)
                        .unwrap();
                    for i in 0..2 {
                        near(
                            thermal.temperature_k(CELLS[i]).unwrap(),
                            300.0 + excess[i],
                            1e-8,
                        );
                    }
                    near(report.source.links[0].conductance_w_k, conductance, 1e-13);
                    assert!(report.entropy.conduction_generation_j_k >= 0.0);
                    let expected_generation: f64 = report
                        .layer_heat
                        .iter()
                        .map(|h| h.total_j() / thermal.temperature_k(h.destination.cell).unwrap())
                        .sum();
                    near(
                        report.entropy.material_loss_generation_j_k,
                        expected_generation,
                        1e-18,
                    );
                    entropy += report.entropy.thermal_storage_change_j_k;
                }
                let expected_entropy: f64 = (0..2)
                    .map(|i| capacity[i] * (thermal.temperature_k(CELLS[i]).unwrap() / 300.0).ln())
                    .sum();
                near(entropy, expected_entropy, 2e-14);
                endpoints.push((
                    owner.mechanical().components()[0].states().to_vec(),
                    [
                        thermal.temperature_k(CELLS[0]).unwrap(),
                        thermal.temperature_k(CELLS[1]).unwrap(),
                    ],
                ));
            },
        );
    }
    assert_eq!(endpoints[0].0, endpoints[1].0);
    assert!(endpoints[1].1[0] < endpoints[0].1[0] && endpoints[1].1[1] < endpoints[0].1[1]);
}

#[test]
fn refused_heat_preserves_nonzero_continuation_and_deterministic_retry() {
    fixture(
        specimens(),
        false,
        false,
        0.02,
        |component, mut thermal, network| {
            let gate = CancelGate::new_clock_free();
            let mut duplicate = component
                .basis()
                .coupled_component(
                    RATE,
                    &[0.02],
                    ModalAcousticTimeBudget::audible_reference(),
                    18,
                    || false,
                )
                .unwrap();
            duplicate.restore_states(component.states()).unwrap();
            let mut expected_thermal = thermal.clone();
            let mut expected = duplicate
                .into_thermal_system(
                    &expected_thermal,
                    network,
                    &DESTINATIONS,
                    config(),
                    Some(prepared()),
                    &gate,
                )
                .unwrap();
            let mut actual = component
                .into_thermal_system(
                    &thermal,
                    network,
                    &DESTINATIONS,
                    config(),
                    Some(prepared()),
                    &gate,
                )
                .unwrap();
            for _ in 0..3 {
                let report = actual
                    .step(&mut thermal, &[0.001], solve(), &[0.0; 2], &gate)
                    .unwrap();
                assert!(report.external_mechanical_loss_j > 0.0);
                let comparison = expected
                    .step(&mut expected_thermal, &[0.001], solve(), &[0.0; 2], &gate)
                    .unwrap();
                assert_eq!(report.mechanical, comparison.mechanical);
            }
            let saved = actual.mechanical().last_frame().unwrap().clone();
            let states = actual.mechanical().components()[0].states().to_vec();
            let cells = thermal.cells().clone();
            let time = thermal.time_s();
            let mut exhausted = solve();
            exhausted.maximum_material_trials = 1;
            assert!(
                actual
                    .step(&mut thermal, &[0.001], exhausted, &[0.0; 2], &gate)
                    .is_err()
            );
            assert_eq!(actual.mechanical().last_frame(), Some(&saved));
            assert_eq!(
                actual.mechanical().components()[0].states(),
                states.as_slice()
            );
            assert_eq!(thermal.cells(), &cells);
            assert_eq!(thermal.time_s(), time);
            let cancel = CancelGate::new_clock_free();
            cancel.request();
            assert!(matches!(
                actual.step(&mut thermal, &[0.001], solve(), &[0.0; 2], &cancel),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            let retried = actual
                .step(&mut thermal, &[0.001], solve(), &[0.0; 2], &gate)
                .unwrap();
            let comparison = expected
                .step(&mut expected_thermal, &[0.001], solve(), &[0.0; 2], &gate)
                .unwrap();
            assert_eq!(retried.mechanical, comparison.mechanical);
            assert_eq!(retried.layer_heat, comparison.layer_heat);
            assert_eq!(retried.entropy, comparison.entropy);
            assert_eq!(thermal.cells(), expected_thermal.cells());
        },
    );
    let mut narrow = specimens();
    narrow[0].mechanical_upper = 300.000001;
    fixture(
        narrow,
        false,
        false,
        0.0,
        |component, mut thermal, network| {
            let gate = CancelGate::new_clock_free();
            let mut owner = component
                .into_thermal_system(&thermal, network, &DESTINATIONS, config(), None, &gate)
                .unwrap();
            let cells = thermal.cells().clone();
            let states = owner.mechanical().components()[0].states().to_vec();
            assert!(
                owner
                    .step(&mut thermal, &[0.0], solve(), &[0.0; 2], &gate)
                    .is_err()
            );
            assert_eq!(owner.mechanical().samples_rendered(), 0);
            assert_eq!(
                owner.mechanical().components()[0].states(),
                states.as_slice()
            );
            assert_eq!(thermal.cells(), &cells);
            assert_eq!(thermal.time_s(), 0.0);
        },
    );
    let mut curved = specimens();
    curved[0].curved_eta = true;
    fixture(curved, false, false, 0.0, |component, thermal, network| {
        assert!(
            component
                .into_thermal_system(
                    &thermal,
                    network,
                    &DESTINATIONS,
                    config(),
                    None,
                    &CancelGate::new_clock_free()
                )
                .is_err()
        );
    });
    fixture(
        specimens(),
        false,
        false,
        0.0,
        |component, thermal, network| {
            let mut wrong = DESTINATIONS;
            wrong[1].cell = wrong[0].cell;
            assert!(
                component
                    .into_thermal_system(
                        &thermal,
                        network,
                        &wrong,
                        config(),
                        None,
                        &CancelGate::new_clock_free()
                    )
                    .is_err()
            );
        },
    );
}
