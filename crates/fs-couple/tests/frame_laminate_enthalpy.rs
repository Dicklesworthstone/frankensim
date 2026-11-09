//! Actual finite caloric cells drive native laminate-frame shape and prestress.
//! Manufactured sources test the physical equations, not a measured alloy fit.

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve,
    material::{ConductionLink, MaterialCell, MaterialNetwork},
};
use fs_couple::frame_enthalpy::laminate::{
    LaminateCaloricReference, LaminateFrameCell, ReciprocalLaminateFrame,
};
use fs_couple::material_circuit::MaterialCircuitSolve;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::laminate::{
    ReferenceSurface,
    strength::{LaminaStrengthProperties, LaminaStrengthSelection},
    thermal::source::{
        LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput, LaminateThermalSource,
    },
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, MaterialPropertySelection,
    SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{
    Density, Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::{Frame2d, Geometry};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;
use std::sync::Arc;

const WIDTH: f64 = 0.01;
const RHO: f64 = 1000.0;
const CONDUCTIVITY: f64 = 0.02;
const REFERENCE: f64 = 300.0;
const RANGE: [f64; 2] = [280.0, 330.0];

#[derive(Clone, Copy)]
struct Ply {
    e: f64,
    slope: f64,
    alpha: f64,
    nu: f64,
    cp: f64,
    h: f64,
}
fn ply(alpha: f64) -> Ply {
    Ply {
        e: 1e8,
        slope: 0.0,
        alpha,
        nu: 0.2,
        cp: 100.0,
        h: 0.002,
    }
}
fn kelvin() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", kelvin(), t).unwrap()
}
fn card(i: usize, p: Ply) -> Arc<MaterialCard> {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        (DENSITY_PROPERTY, Density::DIMS, RHO, 0.0),
        ("young_modulus", Pressure::DIMS, p.e, p.slope),
        ("poisson_ratio", Dims::NONE, p.nu, 0.0),
        (ALPHA, INVERSE_TEMPERATURE_DIMS, p.alpha, 0.0),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            p.cp,
            0.0,
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            CONDUCTIVITY,
            0.0,
        ),
        ("tension", Pressure::DIMS, 20e6, -1e5),
        ("compression", Pressure::DIMS, 20e6, -1e5),
        ("shear", Pressure::DIMS, 20e6, -1e5),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if slope == 0.0 {
                    PropertyValue::Scalar { value, dims }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: RANGE.map(|t| (t, value + slope * (t - REFERENCE))).to_vec(),
                        dims,
                    }
                },
                interpolation: if slope == 0.0 {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    kelvin(),
                    RANGE[0],
                    RANGE[1],
                ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "analytical finite-heat laminate fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("finite heated layer {i}"),
                phase: "solid".into(),
                process: "manufactured source".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}
fn with_cx<R>(gate: &CancelGate, run: impl FnOnce(&Cx<'_>) -> R) -> R {
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 70,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        maximum_circuit_trials: 128,
        maximum_newton_iterations: 12,
        difference_step_j: 1e-7,
        cell_energy_tolerance_j: 2e-10,
        thermal_energy_tolerance_j: 3e-10,
        total_energy_tolerance_j: 5e-10,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-10,
            power_tolerance_w: 1e-10,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-13,
            linear_restart: 2,
            forcing_minimum: 1e-5,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn mapping(layers: usize) -> Vec<LaminateFrameCell> {
    (0..3)
        .flat_map(|station| {
            (0..layers).map(move |layer| LaminateFrameCell {
                element: 0,
                station,
                layer,
                cell: CellId(layer as u64),
            })
        })
        .collect()
}
fn fixture<R>(
    plies: &[Ply],
    length: f64,
    temperatures: &[f64],
    powers: &[f64],
    linked: bool,
    cx: &Cx<'_>,
    run: impl FnOnce(Frame2d, DepositionSimulation<'_>, &MaterialNetwork<'_>) -> R,
) -> R {
    fixture_with_mass(
        plies,
        length,
        temperatures,
        powers,
        linked,
        FrameMassFormulation::Lumped,
        cx,
        run,
    )
}
#[allow(clippy::too_many_arguments)]
fn fixture_with_mass<R>(
    plies: &[Ply],
    length: f64,
    temperatures: &[f64],
    powers: &[f64],
    linked: bool,
    mass: FrameMassFormulation,
    cx: &Cx<'_>,
    run: impl FnOnce(Frame2d, DepositionSimulation<'_>, &MaterialNetwork<'_>) -> R,
) -> R {
    let cards: Vec<_> = plies.iter().enumerate().map(|(i, &p)| card(i, p)).collect();
    let curves: Vec<_> = cards
        .iter()
        .map(|card| {
            ResolvedHeatCapacityCurve::try_new_exact(
                card,
                HeatCapacitySourceQuery {
                    heat_capacity: &PropertyKey::new(
                        "specific_heat_capacity",
                        SPECIFIC_HEAT_CAPACITY_DIMS,
                    ),
                    density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                    temperature_axis: "T",
                    lower: &point(RANGE[0]),
                    upper: &point(RANGE[1]),
                    selection: SelectionPolicy::SingleClaimOnly,
                },
                SolidLiquidPhase::Solid,
                0.0,
            )
            .unwrap()
        })
        .collect();
    let source = Arc::new(
        LaminateThermalSource::try_new(
            cards
                .iter()
                .zip(plies)
                .map(|(card, p)| LaminaSourceInput {
                    card: Arc::clone(card),
                    reference: point(REFERENCE),
                    thickness_m: p.h,
                    orientation_rad: 0.0,
                    elasticity: LaminaElasticSource::Isotropic(
                        IsotropicElasticParameters::YoungPoisson {
                            young_modulus: ScalarPropertyRequirement::try_new(
                                "young_modulus",
                                Pressure::DIMS,
                                ScalarAdmissibility::StrictlyPositive,
                            )
                            .unwrap(),
                            poisson_ratio: ScalarPropertyRequirement::try_new(
                                "poisson_ratio",
                                Dims::NONE,
                                ScalarAdmissibility::Finite,
                            )
                            .unwrap(),
                        },
                    ),
                    expansion: LaminaExpansionSource::Isotropic,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                })
                .collect(),
            ReferenceSurface::Offset(0.007),
            || false,
        )
        .unwrap(),
    );
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, mass);
    frame.add_node(0.0, 0.0);
    frame.add_node(length, 0.0);
    frame.fix(0, [true; 3]);
    frame
        .add_laminate_thermal_element(
            0,
            1,
            source,
            WIDTH,
            &vec![temperatures.to_vec(); 3],
            0.05,
            cx,
        )
        .unwrap();
    frame = frame
        .with_laminate_thermoelasticity(0, &vec![RANGE; plies.len()], cx)
        .unwrap();
    let strength = LaminaStrengthSelection {
        properties: LaminaStrengthProperties {
            tensile: ["tension".into(), "tension".into()],
            compressive: ["compression".into(), "compression".into()],
            shear: "shear".into(),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    frame = frame
        .with_laminate_strength(0, &vec![strength; plies.len()], cx)
        .unwrap();
    let mut events: Vec<_> = plies
        .iter()
        .enumerate()
        .map(|(i, p)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(i as u64),
                mass_kg: RHO * WIDTH * length * p.h,
                specific_h_j_kg: p.cp * (temperatures[i] - RANGE[0]),
                geometry: CellGeometry {
                    start_m: [0.0, i as f64 * p.h, 0.0],
                    end_m: [length, i as f64 * p.h, 0.0],
                    volume_m3: WIDTH * length * p.h,
                    surface_area_m2: 2.0 * (WIDTH * length + WIDTH * p.h + length * p.h),
                },
                boundary: Boundary {
                    ambient_k: REFERENCE,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: powers[i],
                },
            }),
        })
        .collect();
    if linked {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 8,
                a: CellId(0),
                b: CellId(1),
                area_m2: WIDTH * length,
                conductance_w_k: 0.0,
            }),
        });
    }
    let laws: Vec<_> = curves
        .iter()
        .enumerate()
        .map(|(i, c)| (MaterialId(i as u64), c.curve() as &dyn ThermalLaw))
        .collect();
    let thermal = DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 128,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let context = point(REFERENCE);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells: Vec<_> = cards
        .iter()
        .zip(&curves)
        .enumerate()
        .map(|(i, (card, curve))| MaterialCell {
            cell: CellId(i as u64),
            material_id: MaterialId(i as u64),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: linked.then_some((&key, MaterialPropertySelection::SingleClaimOnly)),
        })
        .collect();
    let links = if linked {
        vec![ConductionLink {
            link: 8,
            distance_a_m: plies[0].h / 2.0,
            distance_b_m: plies[1].h / 2.0,
        }]
    } else {
        vec![]
    };
    run(
        frame,
        thermal,
        &MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &links,
        },
    )
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}; allowance {tolerance:e}"
    );
}

#[path = "frame_laminate_enthalpy/dynamics.rs"]
mod dynamics;

#[test]
fn finite_source_heat_distinguishes_free_expansion_and_restrained_prestress() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        for alpha in [1e-3, -1e-3] {
            let p = ply(alpha);
            let length = 0.2;
            let volume = WIDTH * p.h * length;
            let heat = 2.0;
            fixture(
                &[p],
                length,
                &[REFERENCE],
                &[heat],
                false,
                cx,
                |frame, thermal, source| {
                    for restrained in [false, true] {
                        let mut raw = frame.clone();
                        raw.fix(1, [restrained, true, true]);
                        let model = ReciprocalLaminateFrame::try_new(
                            raw,
                            &thermal,
                            source,
                            &mapping(1),
                            &[0.0; 6],
                            LaminateCaloricReference::StressFreeSolid,
                            3,
                            cx,
                        )
                        .unwrap();
                        let before_mass = model
                            .frame()
                            .laminate_thermal_section(0, 1)
                            .unwrap()
                            .prepared()
                            .reference_mass_moments();
                        let candidate = model
                            .step_candidate(&thermal, 1.0, controls(), &[heat], cx)
                            .unwrap();
                        let cp = RHO * p.cp * volume;
                        let b = p.e * alpha * alpha * volume;
                        let a = cp - REFERENCE * b;
                        let delta = if restrained {
                            2.0 * heat / (a + (a * a - 2.0 * b * heat).sqrt())
                        } else {
                            heat / cp
                        };
                        near(
                            candidate.thermal().temperature_k(CellId(0)).unwrap(),
                            REFERENCE + delta,
                            2e-8,
                        );
                        near(
                            candidate.frame().frame().displacements()[3],
                            if restrained {
                                0.0
                            } else {
                                alpha * delta * length
                            },
                            2e-11,
                        );
                        let station = candidate
                            .frame()
                            .frame()
                            .laminate_thermal_response(0, 1)
                            .unwrap();
                        near(
                            station.axial_force_n,
                            if restrained {
                                -p.e * WIDTH * p.h * alpha * delta
                            } else {
                                0.0
                            },
                            2e-8,
                        );
                        let report = candidate.report();
                        near(report.external_thermal_heat_j, heat, 1e-11);
                        near(report.mechanical.external_work_j, 0.0, 1e-14);
                        near(report.total_energy_residual_j, 0.0, 5e-10);
                        near(
                            report.mechanical.reversible_heat_j,
                            if restrained {
                                b * (REFERENCE * delta + 0.5 * delta * delta)
                            } else {
                                0.0
                            },
                            1e-9,
                        );
                        assert_eq!(
                            before_mass,
                            candidate
                                .frame()
                                .frame()
                                .laminate_thermal_section(0, 1)
                                .unwrap()
                                .prepared()
                                .reference_mass_moments()
                        );
                        assert_eq!(model.time_s(), 0.0);
                        assert_eq!(thermal.time_s(), 0.0);
                        assert_eq!(thermal.temperature_k(CellId(0)).unwrap(), REFERENCE);
                    }
                },
            );
        }
    });
}

#[test]
fn held_load_softening_closes_native_work_and_refreshes_original_strength() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let mut p = ply(1e-3);
        p.slope = -1e6;
        p.cp = 500.0;
        let length = 0.2;
        let force = 10.0;
        let heat = 10.0;
        fixture(
            &[p],
            length,
            &[REFERENCE],
            &[heat],
            false,
            cx,
            |mut frame, thermal, source| {
                frame.fix(1, [false, true, true]);
                let load = [0.0, 0.0, 0.0, force, 0.0, 0.0];
                let model = ReciprocalLaminateFrame::try_new(
                    frame,
                    &thermal,
                    source,
                    &mapping(1),
                    &load,
                    LaminateCaloricReference::StressFreeSolid,
                    3,
                    cx,
                )
                .unwrap();
                let initial_strength = model
                    .frame()
                    .assess_laminate_frame_strength(cx)
                    .unwrap()
                    .utilization();
                let result = model
                    .step_candidate(&thermal, 1.0, controls(), &[heat], cx)
                    .unwrap();
                let area = WIDTH * p.h;
                let volume = area * length;
                let sigma = force / area;
                let strain =
                    |t: f64| sigma / (p.e + p.slope * (t - REFERENCE)) + p.alpha * (t - REFERENCE);
                let energy = |t: f64| {
                    let e = p.e + p.slope * (t - REFERENCE);
                    volume
                        * (0.5 * sigma * sigma / e
                            + t * (sigma * p.alpha - 0.5 * p.slope * (sigma / e).powi(2)))
                };
                let residual = |t: f64| {
                    RHO * p.cp * volume * (t - REFERENCE) + energy(t)
                        - energy(REFERENCE)
                        - force * length * (strain(t) - strain(REFERENCE))
                        - heat
                };
                let (mut lo, mut hi) = (REFERENCE, 320.0);
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if residual(mid) > 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                let t = 0.5 * (lo + hi);
                let expected_work = force * length * (strain(t) - strain(REFERENCE));
                near(result.thermal().temperature_k(CellId(0)).unwrap(), t, 2e-8);
                near(
                    result.frame().frame().displacements()[3],
                    length * strain(t),
                    2e-11,
                );
                near(
                    result.report().mechanical.external_work_j,
                    expected_work,
                    1e-9,
                );
                near(
                    result.report().mechanical.stress_work_j,
                    expected_work,
                    1e-9,
                );
                near(
                    result.report().mechanical.internal_energy_change_j,
                    energy(t) - energy(REFERENCE),
                    1e-9,
                );
                near(result.report().total_energy_residual_j, 0.0, 5e-10);
                let hot_strength = result
                    .frame()
                    .frame()
                    .assess_laminate_frame_strength(cx)
                    .unwrap()
                    .utilization();
                assert!(hot_strength > initial_strength);
                near(hot_strength, sigma / (20e6 - 1e5 * (t - REFERENCE)), 1e-10);
            },
        );
    });
}

// Independent continuous two-layer entropy equations after free equilibrium.
// Psi_eq=E*b*L*h*(alpha1*dT1-alpha0*dT0)^2/16. This small test-only
// reference uses physical temperatures and 2x2 scalar algebra, not native
// frame strains, thermal residuals, source exchange or numerical reports.
fn continuum(plies: [Ply; 2], length: f64, mut t: [f64; 2], duration: f64) -> [f64; 2] {
    let volume = WIDTH * length * plies[0].h;
    let capacity = RHO * plies[0].cp * volume;
    let k = plies[0].e * volume / 8.0;
    let a = plies.map(|p| p.alpha);
    let g = CONDUCTIVITY * WIDTH * length / plies[0].h;
    let rate = |t: [f64; 2]| {
        let m = [
            [capacity - t[0] * k * a[0] * a[0], t[0] * k * a[0] * a[1]],
            [t[1] * k * a[0] * a[1], capacity - t[1] * k * a[1] * a[1]],
        ];
        let q = g * (t[1] - t[0]);
        let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        [
            (m[1][1] + m[0][1]) * q / det,
            -(m[1][0] + m[0][0]) * q / det,
        ]
    };
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let k1 = rate(t);
        let k2 = rate(std::array::from_fn(|i| t[i] + 0.5 * dt * k1[i]));
        let k3 = rate(std::array::from_fn(|i| t[i] + 0.5 * dt * k2[i]));
        let k4 = rate(std::array::from_fn(|i| t[i] + dt * k3[i]));
        t = std::array::from_fn(|i| t[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    t
}

#[test]
fn closed_bimetal_heat_exchange_changes_curvature_and_converges_in_physical_time() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let plies = [
            Ply {
                e: 1e10,
                slope: 0.0,
                alpha: 1e-4,
                nu: 0.13,
                cp: 500.0,
                h: 0.001,
            },
            Ply {
                e: 1e10,
                slope: 0.0,
                alpha: 1.5e-4,
                nu: 0.36,
                cp: 500.0,
                h: 0.001,
            },
        ];
        let length = 0.05;
        let initial = [295.0, 305.0];
        let duration = 5.0;
        let exact = continuum(plies, length, initial, duration);
        fixture(
            &plies,
            length,
            &initial,
            &[0.0, 0.0],
            true,
            cx,
            |frame, thermal, source| {
                let start = ReciprocalLaminateFrame::try_new(
                    frame,
                    &thermal,
                    source,
                    &mapping(2),
                    &[0.0; 6],
                    LaminateCaloricReference::StressFreeSolid,
                    6,
                    cx,
                )
                .unwrap();
                let initial_energy = thermal.energy_j()
                    + start
                        .frame()
                        .laminate_thermoelastic_energy(0, cx)
                        .unwrap()
                        .internal_energy_j;
                let initial_tip = start.frame().displacements()[4];
                let mut errors = Vec::new();
                for steps in [8, 16, 32] {
                    let (mut native, mut coupled) = (thermal.clone(), start.clone());
                    let mut total_heat = 0.0;
                    for step in 1..=steps {
                        let result = coupled
                            .step_candidate(
                                &native,
                                duration * step as f64 / steps as f64,
                                controls(),
                                &[0.0; 2],
                                cx,
                            )
                            .unwrap();
                        near(result.report().external_thermal_heat_j, 0.0, 1e-12);
                        near(result.report().total_energy_residual_j, 0.0, 5e-10);
                        assert_eq!(result.report().mechanical.external_work_j, 0.0);
                        total_heat += result.report().mechanical.reversible_heat_j;
                        (native, coupled, _) = result.into_parts();
                    }
                    let t = [
                        native.temperature_k(CellId(0)).unwrap(),
                        native.temperature_k(CellId(1)).unwrap(),
                    ];
                    let free =
                        std::array::from_fn::<_, 2, _>(|i| plies[i].alpha * (t[i] - REFERENCE));
                    let curvature = 3.0 * (free[0] - free[1]) / (4.0 * plies[0].h);
                    let response = coupled.frame().laminate_thermal_response(0, 1).unwrap();
                    near(response.curvature_per_m, curvature, 2e-8);
                    near(response.axial_force_n, 0.0, 2e-7);
                    near(response.bending_moment_n_m, 0.0, 2e-9);
                    near(
                        coupled.frame().displacements()[4],
                        0.5 * length * length * curvature,
                        2e-11,
                    );
                    assert!((coupled.frame().displacements()[4] - initial_tip).abs() > 1e-4);
                    assert!(response.strain_energy_j_m > 0.0);
                    assert!(total_heat.abs() > 1e-5);
                    near(
                        native.energy_j()
                            + coupled
                                .frame()
                                .laminate_thermoelastic_energy(0, cx)
                                .unwrap()
                                .internal_energy_j,
                        initial_energy,
                        2e-8,
                    );
                    errors.push((t[0] - exact[0]).abs().max((t[1] - exact[1]).abs()));
                    assert_eq!(
                        native.cells()[&CellId(0)].mass_kg,
                        thermal.cells()[&CellId(0)].mass_kg
                    );
                    assert_eq!(
                        native.cells()[&CellId(1)].geometry,
                        thermal.cells()[&CellId(1)].geometry
                    );
                }
                assert!(
                    errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
                    "time errors {errors:?}"
                );
                assert!(errors[2] < 0.03, "time errors {errors:?}");
            },
        );
    });
}

fn changed_mass<'a>(
    thermal: &DepositionSimulation<'_>,
    source: &MaterialNetwork<'a>,
    factor: f64,
) -> DepositionSimulation<'a> {
    let events = thermal
        .cells()
        .values()
        .map(|cell| {
            let mut cell = cell.clone();
            cell.mass_kg *= factor;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    let laws: Vec<_> = source
        .cells
        .iter()
        .map(|c| (c.material_id, c.curve as &dyn ThermalLaw))
        .collect();
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 128,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn original_inventory_chart_capacity_and_transaction_refusals_are_not_fallbacks() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        fixture(
            &[ply(1e-3)],
            0.2,
            &[REFERENCE],
            &[2.0],
            false,
            cx,
            |mut frame, thermal, source| {
                frame.fix(1, [true; 3]);
                let locations = mapping(1);
                let make = |frame: Frame2d,
                            state: &DepositionSimulation<'_>,
                            locations: &[LaminateFrameCell]| {
                    ReciprocalLaminateFrame::try_new(
                        frame,
                        state,
                        source,
                        locations,
                        &[0.0; 6],
                        LaminateCaloricReference::StressFreeSolid,
                        3,
                        cx,
                    )
                };
                assert!(make(frame.clone(), &thermal, &locations[..2]).is_err());
                let wrong = changed_mass(&thermal, source, 1.001);
                assert!(make(frame.clone(), &wrong, &locations).is_err());
                let clone = (*source.cells[0].curve).clone();
                let original = &source.cells[0];
                let cloned_cells = [MaterialCell {
                    cell: original.cell,
                    material_id: original.material_id,
                    curve: &clone,
                    point: original.point,
                    conductivity: None,
                }];
                let cloned_network = MaterialNetwork {
                    cells: &cloned_cells,
                    resistors: &[],
                    links: &[],
                };
                assert!(
                    ReciprocalLaminateFrame::try_new(
                        frame.clone(),
                        &thermal,
                        &cloned_network,
                        &locations,
                        &[0.0; 6],
                        LaminateCaloricReference::StressFreeSolid,
                        3,
                        cx
                    )
                    .is_err()
                );
                let model = make(frame, &thermal, &locations).unwrap();
                assert!(
                    model
                        .step_candidate(&wrong, 1.0, controls(), &[2.0], cx)
                        .is_err()
                );
                assert!(
                    model
                        .step_candidate(&thermal, 100.0, controls(), &[200.0], cx)
                        .is_err()
                );
                let gate = CancelGate::new_clock_free();
                gate.request();
                with_cx(&gate, |cancelled| {
                    assert!(matches!(
                        model.step_candidate(&thermal, 1.0, controls(), &[2.0], cancelled),
                        Err(CircuitEnthalpyError::Cancelled)
                    ))
                });
                let first = model
                    .step_candidate(&thermal, 1.0, controls(), &[2.0], cx)
                    .unwrap();
                let retry = model
                    .step_candidate(&thermal, 1.0, controls(), &[2.0], cx)
                    .unwrap();
                assert_eq!(first.thermal().cells(), retry.thermal().cells());
                assert_eq!(
                    first.frame().frame().displacements(),
                    retry.frame().frame().displacements()
                );
                assert_eq!(first.report().mechanical, retry.report().mechanical);
                assert_eq!(thermal.time_s(), 0.0);
                assert_eq!(model.time_s(), 0.0);
                let (hot, next, _) = first.into_parts();
                assert!(
                    model
                        .step_candidate(&hot, 2.0, controls(), &[2.0], cx)
                        .is_err()
                );
                assert!(
                    next.step_candidate(&hot, 2.0, controls(), &[2.0], cx)
                        .is_ok()
                );
            },
        );
        let mut insufficient = ply(1e-3);
        insufficient.cp = 1.0;
        fixture(
            &[insufficient],
            0.2,
            &[REFERENCE],
            &[0.0],
            false,
            cx,
            |frame, thermal, source| {
                assert!(
                    ReciprocalLaminateFrame::try_new(
                        frame,
                        &thermal,
                        source,
                        &mapping(1),
                        &[0.0; 6],
                        LaminateCaloricReference::StressFreeSolid,
                        3,
                        cx
                    )
                    .is_err()
                );
            },
        );
    });
}

#[test]
fn whole_frame_temperature_refresh_matches_native_single_member_updates_atomically() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let length = 0.2;
        fixture(
            &[ply(1e-4)],
            length,
            &[REFERENCE],
            &[0.0],
            false,
            cx,
            |mut frame, _thermal, _source| {
                let source = Arc::new(
                    frame
                        .laminate_thermal_section(0, 0)
                        .unwrap()
                        .source()
                        .source()
                        .clone(),
                );
                frame.add_node(2.0 * length, 0.0);
                frame
                    .add_laminate_thermal_element(
                        1,
                        2,
                        source,
                        WIDTH,
                        &vec![vec![REFERENCE]; 3],
                        0.05,
                        cx,
                    )
                    .unwrap();
                frame = frame
                    .with_laminate_thermoelasticity(1, &[RANGE], cx)
                    .unwrap();
                let initial: Vec<_> = (0..2)
                    .map(|element| {
                        frame
                            .laminate_thermoelastic_response(element, 1, cx)
                            .unwrap()
                    })
                    .collect();
                let mass: Vec<_> = (0..2)
                    .map(|element| {
                        frame
                            .laminate_thermal_section(element, 1)
                            .unwrap()
                            .prepared()
                            .reference_mass_moments()
                    })
                    .collect();
                let field = vec![
                    vec![vec![310.0], vec![315.0], vec![320.0]],
                    vec![vec![305.0], vec![312.0], vec![325.0]],
                ];
                let mut sequential = frame
                    .with_laminate_temperatures(0, &field[0], cx)
                    .unwrap()
                    .with_laminate_temperatures(1, &field[1], cx)
                    .unwrap();
                let mut batch = frame.with_laminate_temperature_field(&field, cx).unwrap();
                for element in 0..2 {
                    assert_eq!(
                        batch.section_deformations(element),
                        sequential.section_deformations(element)
                    );
                    assert_eq!(
                        batch.element_forces(element),
                        sequential.element_forces(element)
                    );
                    assert_eq!(
                        batch
                            .laminate_thermoelastic_response(element, 1, cx)
                            .unwrap(),
                        sequential
                            .laminate_thermoelastic_response(element, 1, cx)
                            .unwrap()
                    );
                    assert_eq!(
                        batch
                            .laminate_thermal_section(element, 1)
                            .unwrap()
                            .prepared()
                            .reference_mass_moments(),
                        mass[element]
                    );
                }
                batch.static_load_with_cx(&[0.0; 9], 1, cx).unwrap();
                sequential.static_load_with_cx(&[0.0; 9], 1, cx).unwrap();
                assert_eq!(batch.displacements(), sequential.displacements());
                assert_eq!(
                    batch.assess_laminate_element_strength(0, cx).unwrap(),
                    sequential.assess_laminate_element_strength(0, cx).unwrap()
                );
                let mut unsupported = field.clone();
                unsupported[1][1][0] = RANGE[1] + 1.0;
                assert!(
                    frame
                        .with_laminate_temperature_field(&unsupported, cx)
                        .is_err()
                );
                assert!(
                    frame
                        .with_laminate_temperature_field(&field[..1], cx)
                        .is_err()
                );
                for element in 0..2 {
                    assert_eq!(
                        frame
                            .laminate_thermoelastic_response(element, 1, cx)
                            .unwrap(),
                        initial[element]
                    );
                    assert_eq!(
                        frame
                            .laminate_thermal_section(element, 1)
                            .unwrap()
                            .prepared()
                            .reference_mass_moments(),
                        mass[element]
                    );
                }
                assert!(frame.displacements().iter().all(|&value| value == 0.0));
            },
        );
    });
}
