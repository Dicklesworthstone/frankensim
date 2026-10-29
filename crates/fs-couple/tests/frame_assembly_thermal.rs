//! Independent physical two-frame references on one mechanical/caloric clock.
//! Manufactured source cards test coupling and material substitution, not an
//! empirical metal fit. Native port blocks deliberately reverse component order.

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    material::{ConductionLink, MaterialCell, MaterialNetwork},
};
use fs_couple::frame_vibration::{
    FrameVibrationBasis, FrameVibrationBudget,
    thermal::{FrameLayerCell, FrameThermalSolve},
    thermoelastic::{FrameThermoelasticBinding, FrameThermoelasticSystem, FrameConnectionHeatDestination},
};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, PreparedModalBudget,
};
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
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, resolve_material_state_point,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Density, DynViscosity, Pressure, QuantitySpec, Temperature};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, FrameViscousSection, Geometry};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

const RATE: u32 = 256;
const LENGTH: f64 = 1.0;
const RADIUS: f64 = 0.1;
const CELLS: [CellId; 2] = [CellId(11), CellId(27)];
const STARTS: [usize; 2] = [7, 0];
const INITIAL: [f64; 4] = [0.01, 0.03, -0.006, -0.02];
const BRIDGE_K: f64 = 3.0;
const BRIDGE_C: f64 = 0.02;
const LINK_G: f64 = 0.015; // .003 / (.03/.2 + .02/.4) [W/K]

#[derive(Clone, Copy)]
struct Spec {
    rho: f64,
    young: f64,
    eta: f64,
    cp: f64,
    e_slope: f64,
    eta_slope: f64,
    alpha: f64,
}
fn specs(live: bool) -> [Spec; 2] {
    [0, 1].map(|i| Spec {
        rho: [3.0, 5.0][i],
        young: [1000.0, 1300.0][i],
        eta: [1.2, 4.5][i],
        cp: [0.01, 0.02][i],
        e_slope: if live { [-1.0, 2.0][i] } else { 0.0 },
        eta_slope: if live { 0.003 } else { 0.0 },
        alpha: if live { [0.00002, -0.00003][i] } else { 0.0 },
    })
}
fn area() -> f64 { core::f64::consts::PI * RADIUS * RADIUS }
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t).unwrap()
        .with("omega", 2048.0).unwrap()
}
fn material(i: usize, s: Spec) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, value, dims, slope) in [
        ("density", s.rho, Density::DIMS, 0.0),
        ("young_modulus", s.young, Pressure::DIMS, s.e_slope),
        ("eta", s.eta, DynViscosity::DIMS, s.eta * s.eta_slope),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, s.alpha, INVERSE_TEMPERATURE_DIMS, 0.0),
        ("specific_heat_capacity", s.cp, SPECIFIC_HEAT_CAPACITY_DIMS, 0.0),
        ("thermal_conductivity", 0.2 * (i + 1) as f64, THERMAL_CONDUCTIVITY_DIMS, 0.0),
    ] {
        let curved = slope != 0.0;
        let mut validity = ValidityDomain::unconstrained().with_quantity(
            "T", QuantitySpec::dimensional(Temperature::DIMS), 200.0, 600.0,
        );
        if name == "eta" { validity = validity.with("omega", 0.0, 4096.0); }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: if curved {
                PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value - 100.0 * slope), (600.0, value + 300.0 * slope)],
                    dims,
                }
            } else { PropertyValue::Scalar { value, dims } },
            validity, uncertainty: UncertaintyModel::Unstated,
            interpolation: if curved { InterpolationPolicy::LinearInside }
                else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![],
            provenance: Provenance {
                source: "manufactured shared thermoelastic frame reference".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId { chemistry: format!("frame {i}"), phase: "solid".into(),
            process: "analytical fixture".into(), revision: 0 }, claims, vec![],
    ).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
            lower: &point(200.0), upper: &point(600.0), selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid, 0.0, 1e-7,
    ).unwrap();
    (card, curve)
}
fn basis(card: &MaterialCard) -> FrameVibrationBasis {
    let requirements = [
        ("density", Density::DIMS), ("young_modulus", Pressure::DIMS),
        ("eta", DynViscosity::DIMS),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS),
    ].map(|(key, dims)| ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap());
    let section = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
        outer_radius_m: RADIUS,
        material: resolve_material_state_point(card, &point(300.0), &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap(),
    }]).unwrap();
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0); frame.add_node(LENGTH, 0.0);
    frame.fix(0, [true; 3]); frame.fix(1, [false, true, true]);
    frame.add_material_viscous_element(0, 1, 5, FrameViscousSection::try_new(section, "eta", "eta").unwrap()).unwrap();
    FrameVibrationBasis::from_frame(&frame, 1,
        FrameVibrationBudget { max_frame_dofs: 6, max_dense_entries: 36 }, || false).unwrap()
}
fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 2, max_connections: 13, max_setup_terms: 32768,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: 100.0, maximum_abs_connection_force_n: 1000.0,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}
fn controls() -> FrameThermalSolve {
    FrameThermalSolve {
        maximum_cells: 2, maximum_material_trials: 128, maximum_newton_iterations: 12,
        difference_step_j: 1e-7, cell_energy_tolerance_j: 1e-10,
        thermal_energy_tolerance_j: 1e-10, total_energy_tolerance_j: 2e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13, relative_tolerance: 1e-14,
            forcing_minimum: 1e-4, linear_restart: 2, ..NewtonKrylovConfig::default()
        },
    }
}
fn fixture<R>(s: [Spec; 2], linked: bool,
    run: impl FnOnce(&[FrameVibrationBasis; 2], DepositionSimulation<'_>, &MaterialNetwork<'_>) -> R,
) -> R {
    let sources = [material(0, s[0]), material(1, s[1])];
    let bases = sources.each_ref().map(|(card, _)| basis(card));
    let context = point(300.0);
    let mut events: Vec<_> = (0..2).map(|i| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CELLS[i], segment_id: i as u64, body_id: i as u64, material_id: MaterialId(i as u64),
            mass_kg: s[i].rho * area() * LENGTH, specific_h_j_kg: 100.0 * s[i].cp,
            geometry: CellGeometry {
                start_m: [0.0, i as f64, 0.0], end_m: [LENGTH, i as f64, 0.0],
                volume_m3: area() * LENGTH,
                surface_area_m2: 2.0 * core::f64::consts::PI * RADIUS * LENGTH + 2.0 * area(),
            },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }),
    }).collect();
    if linked {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: 3, a: CELLS[0], b: CELLS[1], area_m2: 0.003, conductance_w_k: 0.0,
        }) });
    }
    let laws = [(MaterialId(0), sources[0].1.curve() as &dyn ThermalLaw),
        (MaterialId(1), sources[1].1.curve() as &dyn ThermalLaw)];
    let thermal = DepositionSimulation::new(events, &laws,
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 1024 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 0.01, energy_tolerance_j: 1e-10 }, || Ok(())).unwrap().0;
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells: Vec<_> = sources.iter().enumerate().map(|(i, (card, curve))| MaterialCell {
        cell: CELLS[i], material_id: MaterialId(i as u64), curve: curve.curve(),
        point: ThermalMaterialPoint { material: card, context: &context, temperature_axis: "T", liquid_fraction_axis: None },
        conductivity: linked.then_some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
    }).collect();
    let links = if linked { vec![ConductionLink { link: 3, distance_a_m: 0.03, distance_b_m: 0.02 }] } else { vec![] };
    run(&bases, thermal, &MaterialNetwork { cells: &cells, resistors: &[], links: &links })
}
fn parts<'n, 'a>(bases: &[FrameVibrationBasis; 2], thermal: &DepositionSimulation<'_>,
    source: &'n MaterialNetwork<'a>, coupled: bool, cached: bool,
) -> (Vec<FrameThermoelasticBinding<'n, 'a>>, CoupledModalSystem) {
    parts_on_cells(bases, thermal, source, coupled, cached, CELLS)
}
fn parts_on_cells<'n, 'a>(bases: &[FrameVibrationBasis; 2], thermal: &DepositionSimulation<'_>,
    source: &'n MaterialNetwork<'a>, coupled: bool, cached: bool, destinations: [CellId; 2],
) -> (Vec<FrameThermoelasticBinding<'n, 'a>>, CoupledModalSystem) {
    let gate = CancelGate::new_clock_free();
    let shapes = bases.each_ref().map(|b| b.modes()[0].displacement[3]);
    let mut bindings = Vec::new(); let mut models = Vec::new(); let mut blocks = Vec::new();
    for i in 0..2 {
        let mut component = bases[i].coupled_component(RATE, &[0.0],
            ModalAcousticTimeBudget::audible_reference(), 18, || false).unwrap();
        component.restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: INITIAL[2 * i] / shapes[i],
            velocity_m_sqrt_kg_per_s: INITIAL[2 * i + 1] / shapes[i],
        }]).unwrap();
        let (binding, model, ports) = component.into_thermoelastic_parts(
            thermal, source, &[FrameLayerCell { element: 0, layer: 0, cell: destinations[i] }],
            0.05, i, STARTS[i], config(), &gate,
        ).unwrap();
        assert_eq!(ports.len(), 6);
        bindings.push(binding); models.push(model); blocks.push(ports);
    }
    let first = blocks.remove(0);
    let mut connections = blocks.remove(0);
    connections.push(ModalConnection {
        left: ModalAttachment { component: 0, shapes: vec![shapes[0]] },
        right: ModalAttachment { component: 1, shapes: vec![shapes[1]] },
        stiffness_n_m: if coupled { BRIDGE_K } else { 0.0 },
        damping_n_s_m: if coupled { BRIDGE_C } else { 0.0 }, rest_extension_m: 0.0,
    });
    connections.extend(first);
    let mechanical = if cached {
        CoupledModalSystem::new_prepared(models, connections, config(), PreparedModalBudget {
            max_modes: 2, max_dense_bytes: 8 * 1024 * 1024, max_multiply_terms: 1_000_000_000,
        }, &gate)
    } else { CoupledModalSystem::new(models, connections, config(), &gate) }.unwrap();
    (bindings, mechanical)
}
fn assembled<'n, 'a>(bases: &[FrameVibrationBasis; 2], thermal: &DepositionSimulation<'_>,
    source: &'n MaterialNetwork<'a>, coupled: bool, cached: bool, reverse: bool,
) -> FrameThermoelasticSystem<'n, 'a> {
    let (mut bindings, mechanical) = parts(bases, thermal, source, coupled, cached);
    if reverse { bindings.reverse(); }
    FrameThermoelasticBinding::into_shared_system(bindings, mechanical, thermal, &CancelGate::new_clock_free()).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e} (tolerance {tolerance:e})");
}

// Independent continuous SI equations, including three distinct heat integrals.
// Coefficients are frozen only over this reference interval, as in the public
// first-order split. No modal coordinate, native port or native work is reused.
fn reference_interval(s: [Spec; 2], temperatures: [f64; 2], coupled: bool,
    mut y: [f64; 7], dt: f64,
) -> [f64; 7] {
    let mass = s.map(|s| s.rho * area() * LENGTH / 3.0);
    let stiffness: [f64; 2] = std::array::from_fn(|i| (s[i].young + s[i].e_slope * (temperatures[i] - 300.0)) * area() / LENGTH);
    let damping: [f64; 2] = std::array::from_fn(|i| s[i].eta * (1.0 + s[i].eta_slope * (temperatures[i] - 300.0)) * area() / LENGTH);
    let rest: [f64; 2] = std::array::from_fn(|i| LENGTH * s[i].alpha * (temperatures[i] - 300.0));
    let rhs = |y: [f64; 7]| {
        let bridge = if coupled { BRIDGE_K * (y[0] - y[2]) + BRIDGE_C * (y[1] - y[3]) } else { 0.0 };
        [y[1], (-stiffness[0] * (y[0] - rest[0]) - damping[0] * y[1] - bridge) / mass[0],
         y[3], (-stiffness[1] * (y[2] - rest[1]) - damping[1] * y[3] + bridge) / mass[1],
         damping[0] * y[1] * y[1], damping[1] * y[3] * y[3],
         if coupled { BRIDGE_C * (y[1] - y[3]).powi(2) } else { 0.0 }]
    };
    let h = dt / 256.0;
    let shifted = |y: [f64; 7], k: [f64; 7], a: f64| std::array::from_fn(|i| y[i] + a * k[i]);
    for _ in 0..256 {
        let k1 = rhs(y); let k2 = rhs(shifted(y, k1, 0.5 * h));
        let k3 = rhs(shifted(y, k2, 0.5 * h)); let k4 = rhs(shifted(y, k3, h));
        for i in 0..7 { y[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0; }
    }
    y
}
fn physical_state(owner: &FrameThermoelasticSystem<'_, '_>) -> [f64; 4] {
    let mut out = [0.0; 4];
    for frame in owner.frames() {
        let i = frame.component_index();
        let shape = frame.basis().modes()[0].displacement[3];
        let state = owner.mechanical().components()[i].states()[0];
        out[2 * i] = shape * state.displacement_m_sqrt_kg;
        out[2 * i + 1] = shape * state.velocity_m_sqrt_kg_per_s;
    }
    out
}

#[test]
fn two_frames_share_native_motion_finite_heat_and_cross_frame_conduction() {
    let s = specs(false);
    for cached in [false, true] {
        for coupled in [false, true] {
            fixture(s, true, |bases, mut thermal, source| {
                let gate = CancelGate::new_clock_free();
                let mut owner = assembled(bases, &thermal, source, coupled, cached, true);
                assert_eq!(owner.frames().len(), 2);
                assert_eq!(owner.frames()[0].connection_range(), 7..13);
                assert_eq!(owner.frames()[1].connection_range(), 0..6);
                assert_eq!(owner.frames()[0].member_range(), 0..1);
                assert_eq!(owner.frames()[1].member_range(), 1..2);
                let initial_energy = owner.mechanical_internal_energy_j().unwrap() + thermal.energy_j();
                let capacity = s.map(|s| s.rho * area() * LENGTH * s.cp);
                let dt = 1.0 / f64::from(RATE);
                let g = dt * LINK_G;
                let determinant = capacity[0] * capacity[1] + g * (capacity[0] + capacity[1]);
                let mut excess = [0.0; 2];
                let mut y = [INITIAL[0], INITIAL[1], INITIAL[2], INITIAL[3], 0.0, 0.0, 0.0];
                for sample in 1..=32 {
                    let next = reference_interval(s, [300.0; 2], coupled, y, dt);
                    let heat = [next[4] - y[4], next[5] - y[5]];
                    let rhs = [capacity[0] * excess[0] + heat[0], capacity[1] * excess[1] + heat[1]];
                    excess = [((capacity[1] + g) * rhs[0] + g * rhs[1]) / determinant,
                        (g * rhs[0] + (capacity[0] + g) * rhs[1]) / determinant];
                    let report = owner.step(&mut thermal, &[0.0; 2], controls(), &[0.0; 2], &gate).unwrap();
                    for i in 0..2 {
                        assert_eq!(report.members[i].component, i);
                        assert_eq!(report.members[i].element, 0);
                        near(report.layers[i].material_heat_j, heat[i], 2e-11);
                        near(thermal.temperature_k(CELLS[i]).unwrap(), 300.0 + excess[i], 2e-7);
                    }
                    for (a, b) in physical_state(&owner).iter().zip(&next[..4]) { near(*a, *b, 2e-9); }
                    near(report.external_mechanical_loss_j, next[6] - y[6], 2e-11);
                    near(report.source.links[0].conductance_w_k, LINK_G, 1e-14);
                    near(report.total_energy_residual_j, 0.0, 2e-10);
                    near(owner.mechanical_internal_energy_j().unwrap() + thermal.energy_j() + next[6], initial_energy, 2e-9);
                    assert_eq!(owner.mechanical().samples_rendered(), sample);
                    near(thermal.time_s(), sample as f64 * dt, 0.0);
                    y = next;
                }
                assert!(excess[0] > 0.0 && excess[1] > 0.0);
            });
        }
    }
}

#[test]
fn live_material_updates_reach_both_frames_and_data_only_heat_capacity_changes_the_assembly() {
    let mut endpoints = Vec::new();
    for cp_scale in [1.0, 2.0] {
        let mut s = specs(true); s[0].cp *= cp_scale;
        fixture(s, true, |bases, mut thermal, source| {
            let gate = CancelGate::new_clock_free();
            let mut owner = assembled(bases, &thermal, source, true, true, false);
            let initial_energy = owner.mechanical_internal_energy_j().unwrap() + thermal.energy_j();
            let mut outside = 0.0;
            for _ in 0..32 {
                let temperatures = CELLS.map(|c| thermal.temperature_k(c).unwrap());
                let old = physical_state(&owner);
                let reference = reference_interval(s, temperatures, true,
                    [old[0], old[1], old[2], old[3], 0.0, 0.0, 0.0], 1.0 / f64::from(RATE));
                let report = owner.step(&mut thermal, &[0.0; 2], controls(), &[0.0; 2], &gate).unwrap();
                for (a, b) in physical_state(&owner).iter().zip(&reference[..4]) { near(*a, *b, 2e-9); }
                for i in 0..2 {
                    let t = thermal.temperature_k(CELLS[i]).unwrap();
                    let law = owner.mechanical().connection_coefficients(STARTS[i] + 3).unwrap();
                    near(law.stiffness_n_m, (s[i].young + s[i].e_slope * (t - 300.0)) * area() / LENGTH, 2e-10);
                    near(law.rest_extension_m, LENGTH * s[i].alpha * (t - 300.0), 2e-13);
                    near(owner.mechanical().connection_damping_n_s_m(STARTS[i]).unwrap(),
                        s[i].eta * (1.0 + s[i].eta_slope * (t - 300.0)) * area() / LENGTH, 2e-12);
                    near(report.layers[i].material_heat_j, reference[4 + i], 2e-11);
                    // Independent affine-E, constant-alpha internal energy at
                    // actual strain, not a sum of the report's own fields.
                    let x = physical_state(&owner)[2 * i] / LENGTH;
                    let e = s[i].young + s[i].e_slope * (t - 300.0);
                    let elastic = x - s[i].alpha * (t - 300.0);
                    let psi = area() * LENGTH * 0.5 * e * elastic * elastic;
                    let entropy = area() * LENGTH * (e * s[i].alpha * elastic - 0.5 * s[i].e_slope * elastic * elastic);
                    near(report.members[i].layers[0].internal_energy_j, psi + t * entropy, 2e-10);
                }
                outside += report.external_mechanical_loss_j;
                near(owner.mechanical_internal_energy_j().unwrap() + thermal.energy_j() + outside, initial_energy, 3e-9);
            }
            endpoints.push((physical_state(&owner), CELLS.map(|c| thermal.temperature_k(c).unwrap())));
        });
    }
    assert!((endpoints[0].1[0] - endpoints[1].1[0]).abs() > 1e-4);
    assert!((endpoints[0].1[1] - endpoints[1].1[1]).abs() > 1e-5);
    assert!((endpoints[0].0[2] - endpoints[1].0[2]).abs() > 1e-9);
}

#[test]
fn shared_failure_preserves_every_frame_heat_clock_and_retry_independent_of_binding_order() {
    fixture(specs(true), true, |bases, mut thermal, source| {
        let gate = CancelGate::new_clock_free();
        let mut replay_thermal = thermal.clone();
        let mut actual = assembled(bases, &thermal, source, true, true, false)
            .with_connection_heat(&thermal, heat_routes(), 2, &gate).unwrap();
        let mut reversed_routes = heat_routes(); reversed_routes.reverse();
        let mut replay = assembled(bases, &replay_thermal, source, true, true, true)
            .with_connection_heat(&replay_thermal, reversed_routes, 2, &gate).unwrap();
        for _ in 0..3 {
            let a = actual.step(&mut thermal, &[0.0; 2], controls(), &[0.0; 2], &gate).unwrap();
            let b = replay.step(&mut replay_thermal, &[0.0; 2], controls(), &[0.0; 2], &gate).unwrap();
            assert_eq!(a.mechanical, b.mechanical); assert_eq!(a.layers, b.layers);
            assert_eq!(a.connection_heat, b.connection_heat);
        }
        let state = physical_state(&actual);
        let cells = thermal.cells().clone(); let time = thermal.time_s();
        let sample = actual.mechanical().samples_rendered();
        let last = actual.mechanical().last_frame().unwrap().clone();
        let coefficients: Vec<_> = (0..13).map(|i| actual.mechanical().connection_coefficients(i).unwrap()).collect();
        let mut exhausted = controls(); exhausted.maximum_material_trials = 1;
        assert!(actual.step(&mut thermal, &[0.0; 2], exhausted, &[0.0; 2], &gate).is_err());
        let cancelled = CancelGate::new_clock_free(); cancelled.request();
        assert!(matches!(actual.step(&mut thermal, &[0.0; 2], controls(), &[0.0; 2], &cancelled), Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(physical_state(&actual), state); assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), time); assert_eq!(actual.mechanical().samples_rendered(), sample);
        assert_eq!(actual.mechanical().last_frame(), Some(&last));
        for (i, before) in coefficients.iter().enumerate() {
            let after = actual.mechanical().connection_coefficients(i).unwrap();
            assert_eq!(before.stiffness_n_m, after.stiffness_n_m);
            assert_eq!(before.damping_n_s_m, after.damping_n_s_m);
            assert_eq!(before.rest_extension_m, after.rest_extension_m);
        }
        let a = actual.step(&mut thermal, &[0.0; 2], controls(), &[0.0; 2], &gate).unwrap();
        let b = replay.step(&mut replay_thermal, &[0.0; 2], controls(), &[0.0; 2], &gate).unwrap();
        assert_eq!(a.mechanical, b.mechanical); assert_eq!(a.layers, b.layers);
            assert_eq!(a.connection_heat, b.connection_heat);
        assert_eq!(thermal.cells(), replay_thermal.cells());
    });
}

#[test]
fn duplicate_frames_and_foreign_thermal_networks_are_not_silently_coowned() {
    fixture(specs(false), true, |bases, thermal, source| {
        let gate = CancelGate::new_clock_free();
        let (mut bindings, mechanical) = parts(bases, &thermal, source, true, false);
        let (mut duplicates, _) = parts(bases, &thermal, source, true, false);
        bindings[1] = duplicates.remove(0);
        assert!(FrameThermoelasticBinding::into_shared_system(bindings, mechanical, &thermal, &gate).is_err());
        let duplicate_bases = [bases[0].clone(), bases[0].clone()];
        let (bindings, mechanical) = parts_on_cells(&duplicate_bases, &thermal, source,
            true, false, [CELLS[0]; 2]);
        assert!(FrameThermoelasticBinding::into_shared_system(bindings, mechanical, &thermal, &gate).is_err());
        let (mut bindings, mechanical) = parts(bases, &thermal, source, true, false);
        let other = MaterialNetwork { cells: source.cells, resistors: source.resistors, links: source.links };
        let (mut foreign, _) = parts(bases, &thermal, &other, true, false);
        bindings[1] = foreign.remove(1);
        assert!(FrameThermoelasticBinding::into_shared_system(bindings, mechanical, &thermal, &gate).is_err());
    });
}

fn heat_routes() -> Vec<FrameConnectionHeatDestination> {
    // Deliberately reverse cell order; setup canonicalizes the authored map.
    vec![
        FrameConnectionHeatDestination { connection: 6, cell: CELLS[1], fraction: 0.75 },
        FrameConnectionHeatDestination { connection: 6, cell: CELLS[0], fraction: 0.25 },
    ]
}

#[test]
fn connector_loss_heats_declared_material_cells_without_depositing_elastic_work_twice() {
    let s = specs(false);
    for cached in [false, true] {
        fixture(s, true, |bases, mut thermal, source| {
            let gate = CancelGate::new_clock_free();
            let mut owner = assembled(bases, &thermal, source, true, cached, true)
                .with_connection_heat(&thermal, heat_routes(), 2, &gate).unwrap();
            assert_eq!(owner.connection_heat_destinations()[0].cell, CELLS[0]);
            let initial_energy = owner.mechanical_internal_energy_j().unwrap() + thermal.energy_j();
            let capacity = s.map(|s| s.rho * area() * LENGTH * s.cp);
            let dt = 1.0 / f64::from(RATE);
            let g = dt * LINK_G;
            let determinant = capacity[0] * capacity[1] + g * (capacity[0] + capacity[1]);
            let mut excess = [0.0; 2];
            let mut y = [INITIAL[0], INITIAL[1], INITIAL[2], INITIAL[3], 0.0, 0.0, 0.0];
            for _ in 0..32 {
                let next = reference_interval(s, [300.0; 2], true, y, dt);
                let bridge_heat = next[6] - y[6];
                let shares = [0.25, 0.75];
                let heat: [f64; 2] = std::array::from_fn(|i| next[4 + i] - y[4 + i] + shares[i] * bridge_heat);
                let rhs = [capacity[0] * excess[0] + heat[0], capacity[1] * excess[1] + heat[1]];
                excess = [((capacity[1] + g) * rhs[0] + g * rhs[1]) / determinant,
                    (g * rhs[0] + (capacity[0] + g) * rhs[1]) / determinant];
                let report = owner.step(&mut thermal, &[0.0; 2], controls(), &[0.0; 2], &gate).unwrap();
                assert_eq!(report.connection_heat.len(), 2);
                for i in 0..2 {
                    assert_eq!(report.connection_heat[i].cell, CELLS[i]);
                    assert_eq!(report.connection_heat[i].connection, 6);
                    near(report.connection_heat[i].heat_j, shares[i] * bridge_heat, 2e-11);
                    near(thermal.temperature_k(CELLS[i]).unwrap(), 300.0 + excess[i], 2e-7);
                    near(report.layers[i].material_heat_j, next[4 + i] - y[4 + i], 2e-11);
                }
                for (a, b) in physical_state(&owner).iter().zip(&next[..4]) { near(*a, *b, 2e-9); }
                near(report.deposited_material_heat_j, heat.iter().sum(), 2e-11);
                near(report.external_mechanical_loss_j, 0.0, 2e-12);
                near(report.entropy.material_loss_generation_j_k,
                    heat[0] / (300.0 + excess[0]) + heat[1] / (300.0 + excess[1]), 2e-13);
                near(owner.mechanical_internal_energy_j().unwrap() + thermal.energy_j(), initial_energy, 2e-9);
                y = next;
            }
            assert!(y[6] > 0.0);
        });
    }
}

#[test]
fn connector_heating_changes_live_materials_and_the_next_native_motion() {
    fixture(specs(true), true, |bases, mut hot_thermal, source| {
        let gate = CancelGate::new_clock_free();
        let mut sink_thermal = hot_thermal.clone();
        let mut hot = assembled(bases, &hot_thermal, source, true, true, false)
            .with_connection_heat(&hot_thermal, heat_routes(), 2, &gate).unwrap();
        let mut sink = assembled(bases, &sink_thermal, source, true, true, false);
        let initial_energy = hot.mechanical_internal_energy_j().unwrap() + hot_thermal.energy_j();
        for sample in 0..32 {
            let report = hot.step(&mut hot_thermal, &[0.0; 2], controls(), &[0.0; 2], &gate).unwrap();
            sink.step(&mut sink_thermal, &[0.0; 2], controls(), &[0.0; 2], &gate).unwrap();
            near(report.external_mechanical_loss_j, 0.0, 2e-12);
            near(hot.mechanical_internal_energy_j().unwrap() + hot_thermal.energy_j(), initial_energy, 3e-9);
            if sample == 0 {
                // The first frozen-T motion is unchanged. Only its accepted
                // heat destinations differ; later material feedback is real.
                assert_eq!(physical_state(&hot), physical_state(&sink));
                for cell in CELLS {
                    assert!(hot_thermal.temperature_k(cell).unwrap() > sink_thermal.temperature_k(cell).unwrap());
                }
            }
        }
        assert!((physical_state(&hot)[2] - physical_state(&sink)[2]).abs() > 1e-9);
        for i in 0..2 {
            assert_ne!(hot.mechanical().connection_damping_n_s_m(STARTS[i]),
                sink.mechanical().connection_damping_n_s_m(STARTS[i]));
        }
    });
}

#[test]
fn incomplete_duplicate_or_already_owned_connector_heat_is_refused() {
    fixture(specs(false), true, |bases, thermal, source| {
        let gate = CancelGate::new_clock_free();
        let route = |connection, cell, fraction| FrameConnectionHeatDestination { connection, cell, fraction };
        let cases = vec![
            (vec![], 2),
            (vec![route(6, CELLS[0], 0.9)], 2),
            (vec![route(6, CELLS[0], 0.5), route(6, CELLS[0], 0.5)], 2),
            (vec![route(6, CellId(999), 1.0)], 2),
            (vec![route(0, CELLS[0], 1.0)], 2),
            (vec![route(13, CELLS[0], 1.0)], 2),
            (vec![route(6, CELLS[0], f64::NAN)], 2),
            (vec![route(6, CELLS[0], -1.0)], 2),
            (heat_routes(), 1),
        ];
        let original = thermal.cells().clone();
        for (routes, bound) in cases {
            let owner = assembled(bases, &thermal, source, true, false, false);
            assert!(owner.with_connection_heat(&thermal, routes, bound, &gate).is_err());
        }
        let owner = assembled(bases, &thermal, source, true, false, false)
            .with_connection_heat(&thermal, heat_routes(), 2, &gate).unwrap();
        assert!(owner.with_connection_heat(&thermal, heat_routes(), 2, &gate).is_err());
        assert_eq!(thermal.cells(), &original);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
