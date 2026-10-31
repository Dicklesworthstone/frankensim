//! Native capillary body coupling using the original manufactured liquid/wetting fixture.
#[path = "capillary_modal/axisymmetric.rs"]
mod axisymmetric;

use core::f64::consts::PI;
use fs_contact::capillary_bridge::*;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::liquid::surface::adhesion::{
    EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY, EquilibriumAngleOrigin, resolve_young_dupre_wetting,
};
use fs_material::liquid::surface::*;
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, Dims};

const VOLUME: f64 = 1e-10;
const WIDTH: f64 = 1e-3;
fn claim(name: &str, value: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 330.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: format!("manufactured planar bridge: {name}"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn sources(
    theta: f64,
    sigma: f64,
    temperature: f64,
    young: bool,
) -> (ResolvedLiquidSurfaceState, [ResolvedLiquidWettingState; 2]) {
    let (material, interfaces) = source_cards(theta, sigma, young);
    resolved_sources(&material, &interfaces, temperature, young)
}
fn source_cards(theta: f64, sigma: f64, young: bool) -> (MaterialCard, [InterfaceSystemCard; 2]) {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        (DYNAMIC_VISCOSITY_PROPERTY, 1e-3, DYNAMIC_VISCOSITY_DIMS),
        (SURFACE_TENSION_PROPERTY, sigma, SURFACE_TENSION_DIMS),
    ] {
        claims.insert_claim(claim(name, value, dims)).unwrap();
    }
    let id = |name: &str, phase: &str| MaterialStateId {
        chemistry: name.into(),
        phase: phase.into(),
        process: "manufactured verification".into(),
        revision: 0,
    };
    let card = MaterialCard::assemble(id("bridge liquid", "liquid"), claims, vec![]).unwrap();
    let interfaces = std::array::from_fn(|i| {
        let mut claims = ClaimSet::new();
        let (name, value, dims) = if young {
            (
                EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
                sigma * (1.0 + theta.cos()),
                SURFACE_TENSION_DIMS,
            )
        } else {
            (EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, theta, Dims::NONE)
        };
        claims.insert_claim(claim(name, value, dims)).unwrap();
        InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: id(&format!("plate {i}"), "solid"),
                texture_frame: format!("plate {i}"),
            },
            SurfaceSpec {
                material: card.id().clone(),
                texture_frame: "liquid".into(),
            },
            SystemContext {
                medium: "declared vapor".into(),
                third_body: None,
                environment: "manufactured isothermal system".into(),
                history: "clean".into(),
            },
            claims,
            vec![],
        )
        .unwrap()
    });
    (card, interfaces)
}
fn resolved_sources(
    card: &MaterialCard,
    interfaces: &[InterfaceSystemCard; 2],
    temperature: f64,
    young: bool,
) -> (ResolvedLiquidSurfaceState, [ResolvedLiquidWettingState; 2]) {
    let point = QueryPoint::new().with("T", temperature).unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let liquid = resolve_liquid_surface_state(
        card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        selection.clone(),
    )
    .unwrap();
    let wetting = std::array::from_fn(|i| {
        if young {
            resolve_young_dupre_wetting(&interfaces[i], &liquid, selection.clone()).unwrap()
        } else {
            resolve_liquid_wetting_state(&interfaces[i], &point, selection.clone()).unwrap()
        }
    });
    (liquid, wetting)
}
fn limits() -> BridgeLimits {
    BridgeLimits {
        minimum_gap_m: 1e-5,
        maximum_gap_m: 0.02,
        maximum_wetted_span_m: 0.02,
        maximum_abs_force_n: 0.1,
        maximum_abs_surface_energy_j: 1e-3,
        maximum_abs_step_work_j: 1e-3,
        maximum_transition_release_j: 1e-3,
        energy_absolute_tolerance_j: 1e-20,
        energy_relative_tolerance: 1e-12,
        gravity_m_s2: 0.0,
        maximum_bond_number: 0.01,
        maximum_capillary_number: 0.001,
        maximum_weber_number: 0.001,
    }
}
fn model(theta: f64, sigma: f64, young: bool) -> SourcePlanarCapillaryBridge {
    let (liquid, wetting) = sources(theta, sigma, 300.0, young);
    SourcePlanarCapillaryBridge::new(liquid, wetting, VOLUME, WIDTH, limits()).unwrap()
}

use fs_contact::capillary_bridge::load_line::BridgeLoadLineConfig;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::capillary::{
    CapillaryModalConfig, CapillaryModalSystem,
};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingConfig};
use fs_math::{c64::C64, det};
const MA: f64 = 2.0;
const MB: f64 = 3.0;
const RATE: u32 = 1000;
fn network(moving: bool) -> CoupledModalSystem {
    let models = [(MA, -0.0006), (MB, 0.0004)]
        .into_iter()
        .map(|(mass, velocity)| {
            let mut m = ModalAcousticTimeModel::try_new_with_free_coordinates(
                RATE,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            m.restore_states(&[ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: if moving {
                    det::sqrt(mass) * velocity
                } else {
                    0.0
                },
            }])
            .unwrap();
            m
        })
        .collect();
    CoupledModalSystem::new(
        models,
        vec![],
        ModalCouplingConfig {
            max_modes: 2,
            max_connections: 0,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 0.01,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-17,
            energy_relative_tolerance: 1e-10,
        },
        &CancelGate::new(),
    )
    .unwrap()
}
fn config() -> CapillaryModalConfig {
    CapillaryModalConfig {
        contact: BridgeLoadLineConfig {
            max_iterations: 32,
            max_evaluations: 161,
            max_backtracks: 16,
            gap_absolute_tolerance_m: 2e-16,
            gap_relative_tolerance: 1e-11,
            maximum_gap_increment_m: 0.01,
            minimum_compatibility_slope: 1e-4,
        },
        force_absolute_tolerance_n: 1e-12,
        force_relative_tolerance: 1e-10,
    }
}
fn assigned(source: &SourcePlanarCapillaryBridge) -> [&MaterialStateId; 2] {
    [
        &source.wetting()[0].surfaces()[0].material,
        &source.wetting()[1].surfaces()[0].material,
    ]
}
fn column() -> Vec<f64> {
    vec![-1.0 / det::sqrt(MA), 1.0 / det::sqrt(MB)]
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}");
}
#[test]
fn actual_native_bodies_exchange_surface_work_and_retain_liquid_through_rupture() {
    let source = model(PI / 3.0, 0.05, true);
    let gap = 0.5 * (source.formation_gap_m() + source.rupture_gap_m());
    let state = source.initialize(gap, true).unwrap();
    let mut owner = CapillaryModalSystem::new(
        network(true),
        &source,
        assigned(&source),
        state,
        column(),
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let initial = owner.total_energy_j().unwrap();
    let mut release = 0.0;
    let mut ruptures = 0;
    for _ in 0..700 {
        let a = owner.components()[0].states()[0].velocity_m_sqrt_kg_per_s;
        let b = owner.components()[1].states()[0].velocity_m_sqrt_kg_per_s;
        let frame = *owner.step(&[0.0; 2], 0.0).unwrap();
        let va = owner.components()[0].states()[0].velocity_m_sqrt_kg_per_s;
        let vb = owner.components()[1].states()[0].velocity_m_sqrt_kg_per_s;
        close(det::sqrt(MA) * va + det::sqrt(MB) * vb, 0.0, 2e-14);
        close(
            det::sqrt(MA) * (va - a),
            -owner.sample_period_s() * frame.applied_mean_force_n,
            2e-14,
        );
        close(
            det::sqrt(MB) * (vb - b),
            owner.sample_period_s() * frame.applied_mean_force_n,
            2e-14,
        );
        release += frame.topology_release_j + frame.network_dissipation_j;
        close(owner.total_energy_j().unwrap() + release, initial, 2e-15);
        ruptures += usize::from(frame.bridge.event == BridgeEvent::Ruptured);
        assert!(frame.evaluations <= config().contact.max_evaluations);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
    assert_eq!(ruptures, 1);
    assert!(!owner.state().attached());
    assert!(std::ptr::eq(owner.source(), &source));
    assert_eq!(
        owner
            .last_frame()
            .unwrap()
            .bridge
            .profile
            .detached_cap_volumes_m3,
        Some([0.5 * VOLUME; 2])
    );
}
#[test]
fn prescribed_coalescence_keeps_work_release_and_source_assignment_distinct() {
    let source = model(PI / 3.0, 0.05, true);
    let old = source
        .initialize(source.formation_gap_m() + 5e-9, false)
        .unwrap();
    let mut owner = CapillaryModalSystem::new(
        network(false),
        &source,
        assigned(&source),
        old,
        vec![0.0; 2],
        config(),
        &CancelGate::new(),
    )
    .unwrap();
    let initial = owner.total_energy_j().unwrap();
    let rate = -1e-8 / owner.sample_period_s();
    let frame = *owner.step(&[0.0; 2], rate).unwrap();
    assert_eq!(frame.bridge.event, BridgeEvent::Reformed);
    assert!(frame.topology_release_j > 0.0);
    close(
        frame.prescribed_work_j,
        frame.bridge.work_on_liquid_j,
        1e-19,
    );
    close(
        owner.total_energy_j().unwrap() - initial + frame.topology_release_j,
        frame.prescribed_work_j,
        1e-18,
    );
    close(frame.network_energy_j, 0.0, 1e-20);
    let materials = assigned(&source);
    assert!(
        CapillaryModalSystem::new(
            network(false),
            &source,
            [materials[1], materials[0]],
            old,
            column(),
            config(),
            &CancelGate::new()
        )
        .is_err()
    );
}
#[test]
fn native_refusal_keeps_body_bridge_and_observation_for_identical_retry() {
    let source = model(PI / 3.0, 0.05, true);
    let old = source
        .initialize(
            0.5 * (source.formation_gap_m() + source.rupture_gap_m()),
            true,
        )
        .unwrap();
    let build = || {
        CapillaryModalSystem::new(
            network(true),
            &source,
            assigned(&source),
            old,
            column(),
            config(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let mut owner = build();
    let mut reference = build();
    owner.step(&[0.0; 2], 0.0).unwrap();
    reference.step(&[0.0; 2], 0.0).unwrap();
    let state = owner.state();
    let frame = *owner.last_frame().unwrap();
    let bodies: Vec<_> = owner
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect();
    assert!(owner.step(&[0.0; 2], -1e6).is_err());
    let gate = CancelGate::new();
    gate.request();
    assert!(owner.step_under_gate(&[0.0; 2], 0.0, &gate).is_err());
    assert_eq!(owner.state(), state);
    assert_eq!(*owner.last_frame().unwrap(), frame);
    assert_eq!(
        owner
            .components()
            .iter()
            .flat_map(|m| m.states().iter().copied())
            .collect::<Vec<_>>(),
        bodies
    );
    assert_eq!(
        owner.step(&[0.0; 2], 0.0).unwrap(),
        reference.step(&[0.0; 2], 0.0).unwrap()
    );
}

#[test]
fn material_source_refresh_preserves_moving_native_history_and_changes_subsequent_force() {
    let (original, surfaces) = source_cards(PI / 3.0, 0.05, true);
    let curves = |original: &ClaimSet| {
        let mut claims = ClaimSet::new();
        for (_, original) in original.claims_ordered() {
            let mut changed = original.clone();
            if matches!(
                original.key.name(),
                SURFACE_TENSION_PROPERTY | EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY
            ) {
                let PropertyValue::Scalar { value, dims } = original.value else {
                    unreachable!()
                };
                changed.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(300.0, value), (320.0, 0.8 * value)],
                    dims,
                };
                changed.interpolation = InterpolationPolicy::LinearInside;
            }
            claims.insert_claim(changed).unwrap();
        }
        claims
    };
    let material =
        MaterialCard::assemble(original.id().clone(), curves(original.claims()), vec![]).unwrap();
    let surfaces = surfaces.map(|card| {
        InterfaceSystemCard::assemble(
            card.surface_a().clone(),
            card.surface_b().clone(),
            card.context().clone(),
            curves(card.claims()),
            vec![],
        )
        .unwrap()
    });
    let (liquid, wetting) = resolved_sources(&material, &surfaces, 300.0, true);
    let source =
        SourcePlanarCapillaryBridge::new(liquid, wetting, VOLUME, WIDTH, limits()).unwrap();
    let gap = 0.5 * (source.formation_gap_m() + source.rupture_gap_m());
    let initial_state = source.initialize(gap, true).unwrap();
    let gate = CancelGate::new_clock_free();
    let build = || {
        CapillaryModalSystem::new(
            network(true),
            &source,
            assigned(&source),
            initial_state,
            column(),
            config(),
            &gate,
        )
        .unwrap()
    };
    let mut owner = build();
    let mut control = build();
    owner.step(&[0.0; 2], 0.0).unwrap();
    control.step(&[0.0; 2], 0.0).unwrap();
    let state = owner.state();
    let gap = owner.gap_m();
    let before_energy = owner.total_energy_j().unwrap();
    let frame = *owner.last_frame().unwrap();
    let bodies = owner
        .components()
        .iter()
        .map(|m| m.states().to_vec())
        .collect::<Vec<_>>();
    let query = QueryPoint::new().with("T", 320.0).unwrap();
    let stale = source
        .refresh_at_query_point(initial_state, &material, surfaces.each_ref(), &query, &gate)
        .unwrap();
    let refresh = source
        .refresh_at_query_point(state, &material, surfaces.each_ref(), &query, &gate)
        .unwrap();
    assert!(owner.apply_source_refresh(&stale, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(owner.apply_source_refresh(&refresh, &cancelled).is_err());
    assert_eq!(owner.state(), state);
    assert_eq!(owner.total_energy_j().unwrap(), before_energy);
    assert_eq!(owner.samples_rendered(), 1);
    assert_eq!(owner.last_frame(), Some(&frame));
    assert_eq!(
        owner
            .components()
            .iter()
            .map(|m| m.states().to_vec())
            .collect::<Vec<_>>(),
        bodies
    );

    let parameter_work = owner.apply_source_refresh(&refresh, &gate).unwrap();
    assert_eq!(parameter_work, refresh.parameter_work_j());
    close(
        owner.total_energy_j().unwrap() - before_energy,
        parameter_work,
        1e-21,
    );
    assert_eq!(owner.samples_rendered(), 1);
    assert_eq!(owner.gap_m(), gap);
    assert_eq!(owner.state(), refresh.state());
    assert_eq!(owner.last_frame(), Some(&frame));
    assert_eq!(
        owner
            .components()
            .iter()
            .map(|m| m.states().to_vec())
            .collect::<Vec<_>>(),
        bodies
    );
    assert!(owner.apply_source_refresh(&refresh, &gate).is_err());
    let changed = *owner.step(&[0.0; 2], 0.0).unwrap();
    let frozen = *control.step(&[0.0; 2], 0.0).unwrap();
    assert!(changed.applied_mean_force_n.abs() < frozen.applied_mean_force_n.abs());
    assert_ne!(
        owner.components()[0].states(),
        control.components()[0].states()
    );
    let mut exported = changed.network_dissipation_j + changed.topology_release_j;
    for _ in 0..20 {
        let next = *owner.step(&[0.0; 2], 0.0).unwrap();
        exported += next.network_dissipation_j + next.topology_release_j;
    }
    close(
        owner.total_energy_j().unwrap() + exported,
        before_energy + parameter_work,
        1e-15,
    );
}
