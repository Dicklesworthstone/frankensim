//! G1/G3/G4: original two-axis shear, pinned material and finite caloric state.
//! Synthetic claims define independent references, not identified film data.

use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::tangential::planar::material_layer::{
    LayerThermalSolve, PlanarLayerThermalBinding, PlanarMaterialLayerSystem,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedShearLayer, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, resolve_material_state_point,
};
use fs_qty::{Density, DynViscosity, Pressure, Temperature};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use std::cell::Cell;
use std::collections::BTreeMap;
use std::convert::Infallible;

const AREA: f64 = 2e-6;
const THICKNESS: f64 = 1e-6;
const RHO: f64 = 1000.0;
const CP: f64 = 10.0;
const FILM: CellId = CellId(10);
const DT: f64 = 1.0 / 48_000.0;
const K_SLOPE: f64 = -2e4 * AREA / THICKNESS;
const K_INTERCEPT: f64 = 14e6 * AREA / THICKNESS;

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with("omega", 0.0)
        .unwrap()
}

fn affine() -> Option<Vec<(f64, f64)>> {
    Some(vec![(280.0, 8.4e6), (340.0, 7.2e6)])
}

fn coating(
    modulus: Option<Vec<(f64, f64)>>,
    yield_pa: f64,
    variable_loss: bool,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    coating_with_memory(modulus, yield_pa, variable_loss, None)
}

fn coating_with_memory(
    modulus: Option<Vec<(f64, f64)>>,
    yield_pa: f64,
    variable_loss: bool,
    memory: Option<(Vec<(f64, f64)>, f64)>,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, cold, hot, varies) in [
        ("density", Density::DIMS, RHO, RHO, false),
        ("shear_modulus", Pressure::DIMS, 8e6, 8e6, modulus.is_some()),
        (
            "shear_yield_stress",
            Pressure::DIMS,
            if variable_loss { 120.0 } else { yield_pa },
            if variable_loss { 60.0 } else { yield_pa },
            variable_loss,
        ),
        (
            "plastic_shear_viscosity",
            DynViscosity::DIMS,
            if variable_loss { 7.0 } else { 6.0 },
            if variable_loss { 4.0 } else { 6.0 },
            variable_loss,
        ),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            CP,
            CP,
            false,
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            0.2,
            0.2,
            false,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if varies {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: if name == "shear_modulus" {
                            modulus.clone().unwrap()
                        } else {
                            vec![(280.0, cold), (340.0, hot)]
                        },
                        dims,
                    }
                } else {
                    PropertyValue::Scalar { value: cold, dims }
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 340.0)
                    .with("omega", 0.0, 1e6),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if varies {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured planar thermal shear layer".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    if let Some((knots, tau)) = memory {
        for (name, dims, value, interpolation) in [
            ("memory_modulus", Pressure::DIMS, PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS, knots,
                dims: Pressure::DIMS }, InterpolationPolicy::LinearInside),
            ("memory_time", fs_qty::Time::DIMS, PropertyValue::Scalar {
                value: tau, dims: fs_qty::Time::DIMS }, InterpolationPolicy::ConstantWithinValidity),
        ] {
            claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims), value,
                validity: ValidityDomain::unconstrained().with("T", 280.0, 340.0)
                    .with("omega", 0.0, 1e6), uncertainty: UncertaintyModel::Unstated,
                interpolation, observations: vec![], provenance: Provenance {
                    source: "manufactured affine shear memory".into(), license: "CC0-1.0".into(),
                    artifact: None } }).unwrap();
        }
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured-film".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(340.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP,
        1e-8,
    )
    .unwrap();
    (card, caloric)
}

fn source(card: &MaterialCard) -> ResolvedShearLayer {
    source_with_memory(card, false)
}

fn source_with_memory(card: &MaterialCard, memory: bool) -> ResolvedShearLayer {
    let mut requirements: Vec<_> = [
        ("density", Density::DIMS),
        ("shear_modulus", Pressure::DIMS),
        ("shear_yield_stress", Pressure::DIMS),
        ("plastic_shear_viscosity", DynViscosity::DIMS),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
    ]
    .into_iter()
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    })
    .collect();
    if memory {
        for (name, dims) in [("memory_modulus", Pressure::DIMS), ("memory_time", fs_qty::Time::DIMS)] {
            requirements.push(ScalarPropertyRequirement::try_new(name, dims,
                ScalarAdmissibility::StrictlyPositive).unwrap());
        }
    }
    let resolved = resolve_material_state_point(
        card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let layer = ResolvedShearLayer::new(resolved, AREA, THICKNESS, 0.05).unwrap();
    if memory { layer.with_relaxation(vec![("memory_modulus".into(), "memory_time".into())]).unwrap() }
        else { layer }
}

fn thermal(curve: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let volume = AREA * THICKNESS;
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: FILM,
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: RHO * volume,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.002, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 2.0 * AREA + 0.006 * THICKNESS,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(0), curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-14,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn binding<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
) -> PlanarLayerThermalBinding<'a> {
    PlanarLayerThermalBinding {
        cell: FILM,
        material: card,
        caloric: curve,
        temperature_axis: "T",
    }
}

fn owner<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    heat: &DepositionSimulation<'_>,
    connected: bool,
) -> PlanarMaterialLayerSystem<'a> {
    PlanarMaterialLayerSystem::new(
        network(connected),
        attachments(0, IDENTITY),
        Some(attachments(1, IDENTITY)),
        source(card),
        binding(card, curve),
        heat,
        config(),
        matrix_config(),
        &CancelGate::new(),
    )
    .unwrap()
}

fn controls() -> LayerThermalSolve {
    LayerThermalSolve {
        maximum_cells: 1,
        maximum_mechanical_trials: 256,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-11,
        cell_energy_tolerance_j: 1e-13,
        thermal_energy_tolerance_j: 1e-13,
        total_energy_tolerance_j: 2e-12,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-17,
            relative_tolerance: 1e-11,
            linear_restart: 1,
            ..Default::default()
        },
    }
}

fn other(endpoint: &DepositionSimulation<'_>) -> Result<StepCoefficients<f64>, Infallible> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: endpoint.temperature_k(FILM).unwrap(),
    })
}

fn input(speed: [f64; 2]) -> PlanarTangentialModalDrive {
    PlanarTangentialModalDrive {
        contact_temperature_k: None,
        ..drive(speed)
    }
}

fn heat_energy(heat: &DepositionSimulation<'_>) -> f64 {
    heat.cells()
        .values()
        .map(|c| c.mass_kg * c.specific_h_j_kg)
        .sum()
}

fn norm2(x: [f64; 2]) -> f64 {
    x[0] * x[0] + x[1] * x[1]
}

// Independent exact undamped oscillators, without asking the production
// response operator for either its compliance or its free displacement.
fn free_response(
    owner: &PlanarMaterialLayerSystem<'_>,
    speed: [f64; 2],
) -> ([[f64; 2]; 2], [f64; 2]) {
    let mut c = [[0.0; 2]; 2];
    let mut u = speed.map(|v| v * DT);
    for body in 0..2 {
        for mode in 0..2 {
            let w = OMEGAS[body][mode];
            let state = owner.mechanics().network().components()[body].states()[mode];
            let q_free = state.displacement_m_sqrt_kg * ((w * DT).cos() - 1.0)
                + state.velocity_m_sqrt_kg_per_s * (w * DT).sin() / w;
            for i in 0..2 {
                u[i] -= if body == 0 { 1.0 } else { -1.0 } * SHAPES[body][i][mode] * q_free;
                for j in 0..2 {
                    c[i][j] +=
                        SHAPES[body][i][mode] * SHAPES[body][j][mode] * (1.0 - (w * DT).cos())
                            / (w * w);
                }
            }
        }
    }
    (c, u)
}

fn elastic(c: [[f64; 2]; 2], free: [f64; 2], old: [f64; 2], t: f64) -> ([f64; 2], [f64; 2]) {
    let half_k = 0.5 * (K_INTERCEPT + K_SLOPE * t);
    let a = [
        [1.0 + half_k * c[0][0], half_k * c[0][1]],
        [half_k * c[1][0], 1.0 + half_k * c[1][1]],
    ];
    let b = [
        half_k * (2.0 * old[0] + free[0]),
        half_k * (2.0 * old[1] + free[1]),
    ];
    let determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let f = [
        (a[1][1] * b[0] - a[0][1] * b[1]) / determinant,
        (a[0][0] * b[1] - a[1][0] * b[0]) / determinant,
    ];
    let e = std::array::from_fn(|i| old[i] + free[i] - c[i][0] * f[0] - c[i][1] * f[1]);
    (f, e)
}

#[test]
fn temperature_dependent_memory_closes_full_two_axis_thermal_root_and_preserves_refused_history() {
    const K: f64 = 16e6;
    const MEMORY_INTERCEPT: f64 = 14.4e6;
    const MEMORY_SLOPE: f64 = -2e4;
    const TAU: f64 = 0.0005;
    let (card, caloric) = coating_with_memory(None, 1e8, false,
        Some((vec![(280.0, 4.4e6), (340.0, 3.8e6)], TAU)));
    let mut heat = thermal(caloric.curve());
    let gate = CancelGate::new();
    let mut owner = PlanarMaterialLayerSystem::new(network(false), attachments(0, IDENTITY),
        Some(attachments(1, IDENTITY)), source_with_memory(&card, true),
        binding(&card, caloric.curve()), &heat, config(), matrix_config(), &gate).unwrap();
    assert_eq!(owner.stiffness_slope_n_m_k(), 0.0,
        "thermal memory must remain active even when equilibrium stiffness is constant");
    close(owner.relaxation_stiffness_slopes_n_m_k()[0], MEMORY_SLOPE, 1e-10);
    for speed in [[0.0012, -0.0007], [0.0, 0.0], [-0.0010, 0.0006]] {
        let old = owner.mechanics().state().elastic_slip_m();
        let old_r = owner.mechanics().state().relaxation_slip_m()[0];
        let old_t = heat.temperature_k(FILM).unwrap();
        let old_h = heat_energy(&heat);
        let old_body = owner.mechanics().network().total_energy_j().unwrap();
        let (c, free) = free_response(&owner, speed);
        // Independent two-axis midpoint Maxwell equations, followed by a
        // scalar caloric root; no production source/return-map calls below.
        let candidate = |t: f64| {
            let memory = MEMORY_INTERCEPT + MEMORY_SLOPE*t;
            let w = 1.0/(1.0+DT/(2.0*TAU));
            let half = 0.5*(K+memory*w);
            let a = [[1.0+half*c[0][0], half*c[0][1]],
                [half*c[1][0], 1.0+half*c[1][1]]];
            let b: [f64; 2] = std::array::from_fn(|i| K*old[i]+memory*w*old_r[i]+half*free[i]);
            let det = a[0][0]*a[1][1]-a[0][1]*a[1][0];
            let force = [(a[1][1]*b[0]-a[0][1]*b[1])/det,
                (a[0][0]*b[1]-a[1][0]*b[0])/det];
            let de: [f64; 2] = std::array::from_fn(|i| free[i]-c[i][0]*force[0]-c[i][1]*force[1]);
            let e: [f64; 2] = std::array::from_fn(|i| old[i]+de[i]);
            let r: [f64; 2] = std::array::from_fn(|i| 2.0*w*(old_r[i]+0.5*de[i])-old_r[i]);
            let mean: [f64; 2] = std::array::from_fn(|i| 0.5*(r[i]+old_r[i]));
            let loss = DT*memory/TAU*norm2(mean);
            let coefficient = 0.5*MEMORY_SLOPE*(t-old_t)*norm2(old_r);
            let entropy_energy = -0.5*MEMORY_SLOPE*(t*norm2(r)-old_t*norm2(old_r));
            let residual = RHO*AREA*THICKNESS*CP*(t-old_t)-loss+coefficient+entropy_energy;
            (force, e, r, loss, coefficient, entropy_energy, residual)
        };
        let (mut lo, mut hi) = (280.0, 340.0);
        assert!(candidate(lo).6 < 0.0 && candidate(hi).6 > 0.0);
        for _ in 0..64 {
            let mid = 0.5*(lo+hi);
            if candidate(mid).6 < 0.0 { lo=mid; } else { hi=mid; }
        }
        let expected_t = 0.5*(lo+hi);
        let (force, e, r, loss, coefficient, entropy_energy, _) = candidate(expected_t);
        let frame = owner.step_with_layer_heat(&mut heat, &[0.0; 4], input(speed),
            controls(), &[0.0], other, &gate).unwrap();
        close(heat.temperature_k(FILM).unwrap(), expected_t, 2e-7);
        for axis in 0..2 {
            close(frame.mechanical.traction_n[axis], force[axis], 2e-10);
            close(owner.mechanics().state().elastic_slip_m()[axis], e[axis], 2e-14);
            close(owner.mechanics().state().relaxation_slip_m()[0][axis], r[axis], 2e-14);
        }
        close(frame.relaxation_heat_j, loss, 2e-16);
        close(frame.coefficient_free_energy_change_j, coefficient, 2e-16);
        close(frame.mechanical_entropy_energy_change_j, entropy_energy, 2e-16);
        close(frame.reversible_heat_j, -coefficient-entropy_energy, 2e-16);
        close(frame.layer_internal_energy_change_j,
            0.5*K*(norm2(e)-norm2(old))+0.5*MEMORY_INTERCEPT*(norm2(r)-norm2(old_r)), 2e-16);
        assert_eq!(frame.plastic_heat_j, 0.0);
        assert!(frame.relaxation_heat_j > 0.0);
        close(heat_energy(&heat)-old_h, frame.relaxation_heat_j+frame.reversible_heat_j, 1e-13);
        close(frame.mechanical.network_energy_j-old_body+frame.layer_internal_energy_change_j
            +heat_energy(&heat)-old_h, frame.mechanical.drive_work_j, 2e-12);
        assert!(frame.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
    }
    let state = owner.mechanics().state().clone();
    let source = owner.source().clone();
    let cells = heat.cells().clone();
    let sample = owner.mechanics().network().samples_rendered();
    assert!(owner.step_with_layer_heat(&mut heat, &[0.0; 4], input([0.0002, 0.0001]),
        controls(), &[0.0], |_| Err::<StepCoefficients<()>, _>("manufactured late source refusal"),
        &gate).is_err());
    assert_eq!(owner.mechanics().state(), &state);
    assert_eq!(owner.source(), &source);
    assert_eq!(heat.cells(), &cells);
    assert_eq!(owner.mechanics().network().samples_rendered(), sample);
}

#[test]
fn loaded_two_axis_thermoelastic_heating_and_cooling_match_an_independent_root() {
    let (card, caloric) = coating(affine(), 1e8, false);
    let mut heat = thermal(caloric.curve());
    let mut owner = owner(&card, caloric.curve(), &heat, false);
    let gate = CancelGate::new();
    let mut signs = [false; 2];
    for (step, speed) in [[0.0012, -0.0007], [0.0013, -0.0007], [-0.0010, 0.0006]]
        .into_iter()
        .enumerate()
    {
        let old = owner.mechanics().state().elastic_slip_m();
        let old_t = heat.temperature_k(FILM).unwrap();
        let old_k = owner.source().stiffness_n_m();
        let old_h = heat_energy(&heat);
        let old_body = owner.mechanics().network().total_energy_j().unwrap();
        let (c, free) = free_response(&owner, speed);
        assert!(c[0][1].abs() > 0.3 * c[0][0]);
        let residual = |t| {
            let (_, e) = elastic(c, free, old, t);
            RHO * AREA * THICKNESS * CP * (t - old_t) - 0.5 * K_SLOPE * t * (norm2(e) - norm2(old))
        };
        let (mut lo, mut hi) = (280.0, 340.0);
        assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if residual(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let expected_t = 0.5 * (lo + hi);
        let (expected_f, expected_e) = elastic(c, free, old, expected_t);
        let frame = owner
            .step_with_layer_heat(
                &mut heat,
                &[0.0; 4],
                input(speed),
                controls(),
                &[0.0],
                other,
                &gate,
            )
            .unwrap();
        let t = heat.temperature_k(FILM).unwrap();
        close(t, expected_t, 2e-7);
        for i in 0..2 {
            close(frame.mechanical.traction_n[i], expected_f[i], 2e-10);
            close(
                owner.mechanics().state().elastic_slip_m()[i],
                expected_e[i],
                2e-14,
            );
        }
        close(
            frame.coefficient_free_energy_change_j,
            0.5 * (owner.source().stiffness_n_m() - old_k) * norm2(old),
            2e-17,
        );
        close(
            frame.mechanical_entropy_energy_change_j,
            -0.5 * K_SLOPE * (t * norm2(expected_e) - old_t * norm2(old)),
            2e-16,
        );
        close(
            frame.reversible_heat_j,
            0.5 * K_SLOPE * t * (norm2(expected_e) - norm2(old)),
            2e-16,
        );
        close(
            frame.layer_internal_energy_change_j,
            0.5 * K_INTERCEPT * (norm2(expected_e) - norm2(old)),
            2e-16,
        );
        assert_eq!(frame.plastic_heat_j, 0.0);
        assert_eq!(frame.source, t);
        close(heat_energy(&heat) - old_h, frame.reversible_heat_j, 1e-13);
        close(
            frame.mechanical.network_energy_j - old_body
                + frame.layer_internal_energy_change_j
                + heat_energy(&heat)
                - old_h,
            frame.mechanical.drive_work_j,
            2e-12,
        );
        assert!(frame.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
        if step > 0 {
            assert!(norm2(old) > 0.0);
            signs[usize::from(frame.reversible_heat_j > 0.0)] = true;
        }
    }
    assert!(
        signs.into_iter().all(|s| s),
        "loaded strain increase must cool and release must heat"
    );
}

#[test]
fn nonproportional_plastic_history_heats_the_same_source_and_closes_internal_energy() {
    let (card, caloric) = coating(affine(), 100.0, true);
    let mut heat = thermal(caloric.curve());
    let mut owner = owner(&card, caloric.curve(), &heat, true);
    let original = owner.source().clone();
    let gate = CancelGate::new();
    let (mut work, mut network_loss, mut plastic, mut reversible) = (0.0, 0.0, 0.0, 0.0);
    for n in 0..32 {
        let speed = [[0.001, 0.0], [0.0, 0.001], [-0.001, 0.0], [0.0, -0.001]][n / 8];
        let old_e = owner.mechanics().state().elastic_slip_m();
        let old_p = owner.mechanics().state().plastic_slip_m();
        let frame = owner
            .step_with_layer_heat(
                &mut heat,
                &[1e-4, -2e-4, 1e-4, -1e-4],
                input(speed),
                controls(),
                &[0.0],
                other,
                &gate,
            )
            .unwrap();
        let state = owner.mechanics().state();
        let e = state.elastic_slip_m();
        let p = state.plastic_slip_m();
        let dp = [p[0] - old_p[0], p[1] - old_p[1]];
        let distance = dp[0].hypot(dp[1]);
        let expected_heat = owner.source().yield_force_n() * distance
            + owner.source().plastic_resistance_n_s_m() * distance * distance / DT;
        close(frame.plastic_heat_j, expected_heat, 1e-16);
        assert_eq!(
            frame.plastic_heat_j,
            frame.mechanical.friction_dissipation_j
        );
        for i in 0..2 {
            // History consumes the constitutive displacement. The original
            // mechanical owner separately reports its actual endpoint mismatch.
            let compatibility = frame.mechanical.compatibility_residual_m[i];
            assert!(compatibility.abs() <= matrix_config().displacement_tolerance_m);
            close(
                e[i],
                old_e[i] + frame.mechanical.slip_distance_m[i] - compatibility - dp[i],
                1e-15,
            );
            close(
                frame.mechanical.traction_n[i],
                0.5 * owner.source().stiffness_n_m() * (old_e[i] + e[i]),
                1e-10,
            );
        }
        close(
            frame.reversible_heat_j
                + frame.coefficient_free_energy_change_j
                + frame.mechanical_entropy_energy_change_j,
            0.0,
            1e-17,
        );
        close(
            frame.deposited_heat_j,
            frame.plastic_heat_j + frame.reversible_heat_j,
            1e-17,
        );
        assert_eq!(
            owner.source().temperature_kelvin(),
            heat.temperature_k(FILM).unwrap()
        );
        work += frame.mechanical.external_work_j
            + frame.mechanical.drive_work_j
            + frame.external_heat_j;
        network_loss += frame.mechanical.network_dissipation_j;
        plastic += frame.plastic_heat_j;
        reversible += frame.reversible_heat_j;
        // U_layer = (k - T*k') |e|²/2 = constant intercept * |e|²/2.
        close(
            frame.mechanical.network_energy_j
                + 0.5 * K_INTERCEPT * norm2(e)
                + heat_energy(&heat)
                + network_loss,
            work,
            2e-11,
        );
        assert!(frame.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j);
    }
    assert!(plastic > 0.0);
    assert!(
        owner
            .mechanics()
            .state()
            .plastic_slip_m()
            .iter()
            .all(|v| v.abs() > 1e-9)
    );
    assert!(heat.temperature_k(FILM).unwrap() > 300.0);
    close(heat_energy(&heat), plastic + reversible, 2e-12);
    for (a, b) in owner
        .source()
        .resolved()
        .properties()
        .iter()
        .zip(original.resolved().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    assert_eq!(owner.mechanics().network().samples_rendered(), 32);
    assert_eq!(heat.time_s(), 32.0 * DT);
}

#[test]
fn constant_modulus_and_loss_recover_original_planar_motion_without_aging_trials() {
    let (card, caloric) = coating(None, 100.0, false);
    let mut heat = thermal(caloric.curve());
    let mut heated = owner(&card, caloric.curve(), &heat, false);
    let s = source(&card);
    let law = ThirdBodyLayerContact::new(
        s.stiffness_n_m(),
        s.yield_force_n(),
        s.plastic_resistance_n_s_m(),
        300.0,
        s.maximum_elastic_slip_m(),
    )
    .unwrap();
    let mut ordinary = system(law, false, IDENTITY);
    let gate = CancelGate::new();
    let mut generated = 0.0;
    for n in 0..24 {
        let speed = if n < 12 {
            [0.001, -0.0006]
        } else {
            [-0.0005, 0.0008]
        };
        let frame = heated
            .step_with_layer_heat(
                &mut heat,
                &[0.0; 4],
                input(speed),
                controls(),
                &[0.0],
                other,
                &gate,
            )
            .unwrap();
        let reference = *ordinary.step(&[0.0; 4], drive(speed)).unwrap();
        assert_eq!(frame.mechanical, reference);
        assert_eq!(states(heated.mechanics()), states(&ordinary));
        assert_eq!(
            heated.mechanics().state().elastic_slip_m(),
            ordinary.state().elastic_slip_m()
        );
        assert_eq!(
            heated.mechanics().state().plastic_slip_m(),
            ordinary.state().plastic_slip_m()
        );
        assert_eq!(frame.coefficient_free_energy_change_j, 0.0);
        assert_eq!(frame.mechanical_entropy_energy_change_j, 0.0);
        assert_eq!(frame.reversible_heat_j, 0.0);
        generated += frame.plastic_heat_j;
    }
    assert!(generated > 0.0);
    close(heat_energy(&heat), generated, 2e-12);
    assert_eq!(
        heated.mechanics().network().samples_rendered(),
        ordinary.network().samples_rendered()
    );
    assert_eq!(heat.time_s(), 24.0 * DT);
}

#[test]
fn source_affinity_geometry_and_explicit_load_temperature_boundaries_refuse() {
    let (card, caloric) = coating(affine(), 100.0, false);
    let heat = thermal(caloric.curve());
    let gate = CancelGate::new();
    let make = |source, binding, thermal: &DepositionSimulation<'_>| {
        PlanarMaterialLayerSystem::new(
            network(false),
            attachments(0, IDENTITY),
            Some(attachments(1, IDENTITY)),
            source,
            binding,
            thermal,
            config(),
            matrix_config(),
            &gate,
        )
    };
    let (kink, kink_caloric) = coating(
        Some(vec![(280.0, 8.4e6), (300.0, 8e6), (340.0, 7.6e6)]),
        100.0,
        false,
    );
    let kink_heat = thermal(kink_caloric.curve());
    assert!(
        make(
            source(&kink),
            binding(&kink, kink_caloric.curve()),
            &kink_heat
        )
        .is_err()
    );
    let (foreign, _) = coating(None, 200.0, false);
    assert!(make(source(&card), binding(&foreign, caloric.curve()), &heat).is_err());
    let wider = ResolvedShearLayer::new(
        source(&card).resolved().clone(),
        2.0 * AREA,
        THICKNESS,
        0.05,
    )
    .unwrap();
    assert!(make(wider, binding(&card, caloric.curve()), &heat).is_err());
    let mut owner = owner(&card, caloric.curve(), &heat, false);
    let mut heat = heat;
    let initial = states(owner.mechanics());
    for drive in [
        PlanarTangentialModalDrive {
            contact_temperature_k: Some(300.0),
            ..input([0.001, 0.0])
        },
        PlanarTangentialModalDrive {
            normal_force_n: 0.0,
            ..input([0.001, 0.0])
        },
        PlanarTangentialModalDrive {
            nominal_pressure_pa: Some(2.0 / AREA),
            ..input([0.001, 0.0])
        },
    ] {
        assert!(
            owner
                .step_with_layer_heat(
                    &mut heat,
                    &[0.0; 4],
                    drive,
                    controls(),
                    &[0.0],
                    other,
                    &gate
                )
                .is_err()
        );
        assert_eq!(states(owner.mechanics()), initial);
        assert_eq!(owner.mechanics().network().samples_rendered(), 0);
        assert_eq!(heat.time_s(), 0.0);
    }
    let mut outside = thermal(caloric.curve());
    assert!(
        owner
            .step_with_layer_heat(
                &mut outside,
                &[0.0; 4],
                input([0.001, 0.0]),
                controls(),
                &[60.0 * RHO * AREA * THICKNESS * CP],
                other,
                &gate
            )
            .is_err()
    );
    let accepted = owner
        .step_with_layer_heat(
            &mut heat,
            &[0.0; 4],
            PlanarTangentialModalDrive {
                nominal_pressure_pa: Some(1.0 / AREA),
                ..input([0.001, 0.0])
            },
            controls(),
            &[0.0],
            other,
            &gate,
        )
        .unwrap();
    assert!(accepted.plastic_heat_j > 0.0);
}

#[test]
fn loaded_late_source_budget_and_cancel_failures_preserve_modal_film_and_heat_retry() {
    let (card, caloric) = coating(affine(), 100.0, true);
    let (mut heat, mut reference_heat) = (thermal(caloric.curve()), thermal(caloric.curve()));
    let mut owner = owner(&card, caloric.curve(), &heat, true);
    let mut reference = self::owner(&card, caloric.curve(), &reference_heat, true);
    let gate = CancelGate::new();
    let drive = input([0.001, -0.0007]);
    for _ in 0..4 {
        owner
            .step_with_layer_heat(
                &mut heat,
                &[0.0; 4],
                drive,
                controls(),
                &[0.0],
                other,
                &gate,
            )
            .unwrap();
        reference
            .step_with_layer_heat(
                &mut reference_heat,
                &[0.0; 4],
                drive,
                controls(),
                &[0.0],
                other,
                &gate,
            )
            .unwrap();
    }
    let saved_states = states(owner.mechanics());
    let saved_source = owner.source().clone();
    let saved_history = owner.mechanics().state().clone();
    let saved_frame = owner.mechanics().last_frame().copied();
    let saved_cells = heat.cells().clone();
    let saved_time = heat.time_s();
    assert!(norm2(saved_history.plastic_slip_m()) > 0.0);
    let unchanged = |owner: &PlanarMaterialLayerSystem<'_>, heat: &DepositionSimulation<'_>| {
        assert_eq!(states(owner.mechanics()), saved_states);
        assert_eq!(owner.source(), &saved_source);
        assert_eq!(owner.mechanics().state(), &saved_history);
        assert_eq!(owner.mechanics().last_frame().copied(), saved_frame);
        assert_eq!(owner.mechanics().network().samples_rendered(), 4);
        assert_eq!(heat.cells(), &saved_cells);
        assert_eq!(heat.time_s(), saved_time);
    };
    let calls = Cell::new(0);
    assert!(
        owner
            .step_with_layer_heat(
                &mut heat,
                &[0.0; 4],
                drive,
                controls(),
                &[0.0],
                |endpoint| {
                    calls.set(calls.get() + 1);
                    if calls.get() == 2 {
                        Err("late transport refusal")
                    } else {
                        Ok(other(endpoint).unwrap())
                    }
                },
                &gate
            )
            .is_err()
    );
    assert_eq!(calls.get(), 2);
    unchanged(&owner, &heat);
    let mut tiny = controls();
    tiny.maximum_mechanical_trials = 1;
    assert!(
        owner
            .step_with_layer_heat(&mut heat, &[0.0; 4], drive, tiny, &[0.0], other, &gate)
            .is_err()
    );
    unchanged(&owner, &heat);
    let cancelled = CancelGate::new();
    assert!(
        owner
            .step_with_layer_heat(
                &mut heat,
                &[0.0; 4],
                drive,
                controls(),
                &[0.0],
                |endpoint| {
                    cancelled.request();
                    other(endpoint)
                },
                &cancelled
            )
            .is_err()
    );
    unchanged(&owner, &heat);
    let retry = owner
        .step_with_layer_heat(
            &mut heat,
            &[0.0; 4],
            drive,
            controls(),
            &[0.0],
            other,
            &gate,
        )
        .unwrap();
    let expected = reference
        .step_with_layer_heat(
            &mut reference_heat,
            &[0.0; 4],
            drive,
            controls(),
            &[0.0],
            other,
            &gate,
        )
        .unwrap();
    assert_eq!(retry.mechanical, expected.mechanical);
    assert_eq!(retry.deposited_heat_j, expected.deposited_heat_j);
    assert_eq!(states(owner.mechanics()), states(reference.mechanics()));
    assert_eq!(owner.mechanics().state(), reference.mechanics().state());
    assert_eq!(owner.source(), reference.source());
    assert_eq!(heat.cells(), reference_heat.cells());
    assert_eq!(heat.time_s(), reference_heat.time_s());
}
