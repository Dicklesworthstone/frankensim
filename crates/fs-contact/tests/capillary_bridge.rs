//! G0/G1/G3/G4: manufactured source cards and independent planar geometry/work.
//! This is verification of the declared 2-D model, not a measured rupture law.
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

#[path = "capillary_bridge/axisymmetric.rs"]
mod axisymmetric;

const VOLUME: f64 = 1e-10;
const WIDTH: f64 = 1e-3;

#[path = "capillary_bridge/load_line.rs"]
mod load_line;

#[path = "capillary_bridge/axisymmetric_load_line.rs"]
mod axisymmetric_load_line;
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

fn temperature_cards(
    young: bool,
    density_factor: f64,
    angle_factor: f64,
) -> (MaterialCard, [InterfaceSystemCard; 2]) {
    let (card, interfaces) = source_cards(0.7, 0.05, young);
    let curves = |original: &ClaimSet| {
        let mut claims = ClaimSet::new();
        for (_, original) in original.claims_ordered() {
            let factor = match original.key.name() {
                "density" => density_factor,
                DYNAMIC_VISCOSITY_PROPERTY => 2.0,
                SURFACE_TENSION_PROPERTY => 0.8,
                EQUILIBRIUM_CONTACT_ANGLE_PROPERTY => angle_factor,
                EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY => {
                    0.8 * (1.0 + (0.7 * angle_factor).cos()) / (1.0 + 0.7_f64.cos())
                }
                _ => unreachable!(),
            };
            let mut changed = original.clone();
            let PropertyValue::Scalar { value, dims } = original.value else {
                unreachable!()
            };
            changed.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(300.0, value), (320.0, factor * value)],
                dims,
            };
            changed.interpolation = InterpolationPolicy::LinearInside;
            claims.insert_claim(changed).unwrap();
        }
        claims
    };
    let material =
        MaterialCard::assemble(card.id().clone(), curves(card.claims()), vec![]).unwrap();
    let interfaces = interfaces.map(|card| {
        InterfaceSystemCard::assemble(
            card.surface_a().clone(),
            card.surface_b().clone(),
            card.context().clone(),
            curves(card.claims()),
            vec![],
        )
        .unwrap()
    });
    (material, interfaces)
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
fn near(a: f64, b: f64, rel: f64, scale: f64) {
    assert!(
        (a - b).abs() <= rel * scale.max(a.abs()).max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let n = 4096;
    let h = (b - a) / f64::from(n);
    let mut sum = f(a) + f(b);
    for i in 1..n {
        sum += f(a + f64::from(i) * h) * if i % 2 == 0 { 2.0 } else { 4.0 };
    }
    sum * h / 3.0
}

#[test]
fn g1_original_circle_geometry_young_laplace_force_and_energy_derivatives() {
    for young in [false, true] {
        for theta in [0.2, 0.7, 1.3] {
            let m = model(theta, 0.05, young);
            let h = 0.5 * (m.formation_gap_m() + m.rupture_gap_m());
            let p = m.evaluate(m.initialize(h, true).unwrap()).unwrap();
            let theta = m.surface_law().contact_angle_radians();
            let (s, c) = theta.sin_cos();
            let radius = h / (2.0 * c);
            // Independent direct circular profile integration, with no use of
            // the owner's area/neck/energy coefficient formulas.
            let x = |y: f64| p.wetted_span_m / 2.0 + radius * s - (radius * radius - y * y).sqrt();
            near(
                WIDTH * simpson(|y| 2.0 * x(y), -h / 2.0, h / 2.0),
                VOLUME,
                2e-9,
                0.0,
            );
            near(
                WIDTH
                    * simpson(
                        |y| 2.0 * radius / (radius * radius - y * y).sqrt(),
                        -h / 2.0,
                        h / 2.0,
                    ),
                p.liquid_gas_area_m2,
                2e-8,
                0.0,
            );
            near(2.0 * x(0.0), p.neck_width_m.unwrap(), 2e-14, 1e-6);
            near(p.liquid_minus_gas_pressure_pa, -0.05 / radius, 2e-14, 0.0);
            let force =
                -p.liquid_minus_gas_pressure_pa * p.wetted_span_m * WIDTH + 2.0 * 0.05 * WIDTH * s;
            near(p.normal_force_n, -force, 2e-14, 0.0);
            near(
                p.surface_free_energy_j,
                m.surface_law()
                    .surface_energy_j(p.liquid_gas_area_m2, p.wetted_area_m2)
                    .unwrap(),
                2e-14,
                1e-8,
            );
            let dh = h * 1e-5;
            let lo = m.evaluate(m.initialize(h - dh, true).unwrap()).unwrap();
            let hi = m.evaluate(m.initialize(h + dh, true).unwrap()).unwrap();
            near(
                -(hi.surface_free_energy_j - lo.surface_free_energy_j) / (2.0 * dh),
                p.normal_force_n,
                2e-9,
                0.0,
            );
            near(
                (hi.normal_force_n - lo.normal_force_n) / (2.0 * dh),
                p.force_gap_derivative_n_per_m,
                2e-9,
                0.0,
            );
            let port = p.point_port([0.0, 0.0, 1.0], [0.01, 0.0, 0.0]).unwrap();
            assert_eq!(port.action_force_n, [0.0, 0.0, p.normal_force_n]);
            assert_eq!(port.reaction_force_n, port.action_force_n.map(|v| -v));
            near(
                port.action_moment_n_m[1],
                -0.01 * p.normal_force_n,
                1e-15,
                0.0,
            );
            assert_eq!(port.tangent_n_per_m, -p.force_gap_derivative_n_per_m);
        }
    }
}

#[test]
fn g3_rupture_and_coalescence_keep_all_liquid_and_close_the_hysteresis_cycle() {
    let m = model(PI / 3.0, 0.05, true);
    let gate = CancelGate::new_clock_free();
    let far = 1.2 * m.rupture_gap_m();
    let close = 0.8 * m.formation_gap_m();
    let start = m.initialize(far, false).unwrap();
    let initial = m.evaluate(start).unwrap();
    let theta = m.surface_law().contact_angle_radians();
    let segment = theta - theta.sin() * theta.cos();
    let cap_radius = ((VOLUME / WIDTH) / (2.0 * segment)).sqrt();
    near(
        initial.surface_free_energy_j,
        4.0 * 0.05 * WIDTH * cap_radius * segment,
        2e-14,
        0.0,
    );
    near(
        m.formation_gap_m(),
        2.0 * cap_radius * (1.0 - theta.cos()),
        2e-14,
        0.0,
    );
    let joined = m.advance(start, close, 100.0, &gate).unwrap();
    assert_eq!(joined.event, BridgeEvent::Reformed);
    assert!(joined.irreversible_release_j > 0.0);
    assert!(joined.profile.attached);
    assert_eq!(joined.profile.liquid_volume_m3, VOLUME);
    let split = m.advance(joined.state(), far, 100.0, &gate).unwrap();
    assert_eq!(split.event, BridgeEvent::Ruptured);
    assert!(split.irreversible_release_j > 0.0);
    assert_eq!(
        split.profile.detached_cap_volumes_m3,
        Some([VOLUME / 2.0; 2])
    );
    assert_eq!(
        split
            .profile
            .detached_cap_volumes_m3
            .unwrap()
            .iter()
            .sum::<f64>(),
        VOLUME
    );
    assert_eq!(split.profile.normal_force_n, 0.0);
    assert!(
        split.mean_normal_force_n < 0.0,
        "endpoint force would incorrectly discard rupture-crossing work"
    );
    for (step, delta) in [(joined, close - far), (split, far - close)] {
        near(
            step.mean_normal_force_n * delta,
            -step.work_on_liquid_j,
            2e-15,
            1e-8,
        );
        near(
            step.surface_energy_change_j + step.irreversible_release_j,
            step.work_on_liquid_j,
            2e-14,
            1e-8,
        );
        assert!(step.energy_balance_residual_j.abs() < 1e-20);
    }
    near(
        joined.work_on_liquid_j + split.work_on_liquid_j,
        joined.irreversible_release_j + split.irreversible_release_j,
        2e-14,
        1e-8,
    );
    println!(
        "planar bridge: V={VOLUME:e} m3, touching gap={:e} m, geometric zero-neck gap={:e} m, coalescence release={:e} J, rupture release={:e} J, cycle work={:e} J",
        m.formation_gap_m(),
        m.rupture_gap_m(),
        joined.irreversible_release_j,
        split.irreversible_release_j,
        joined.work_on_liquid_j + split.work_on_liquid_j
    );
    assert_eq!(
        m.evaluate(split.state()).unwrap().surface_free_energy_j,
        initial.surface_free_energy_j
    );
    let free = m.advance(split.state(), 1.1 * far, 100.0, &gate).unwrap();
    assert_eq!(free.event, BridgeEvent::None);
    assert_eq!(free.work_on_liquid_j, 0.0);
    assert_eq!(free.irreversible_release_j, 0.0);
    assert_eq!(free.mean_normal_force_n, 0.0);
}

#[test]
fn g3_event_subdivision_and_adjacent_float_work_do_not_lose_surface_work() {
    let m = model(0.7, 0.05, false);
    let gate = CancelGate::new_clock_free();
    let near_gap = 0.8 * m.formation_gap_m();
    let far = 1.1 * m.rupture_gap_m();
    let start = m.initialize(near_gap, true).unwrap();
    let direct = m.advance(start, far, 100.0, &gate).unwrap();
    let a = m
        .advance(start, 0.9 * m.rupture_gap_m(), 100.0, &gate)
        .unwrap();
    let b = m
        .advance(a.state(), m.rupture_gap_m(), 100.0, &gate)
        .unwrap();
    let c = m.advance(b.state(), far, 100.0, &gate).unwrap();
    assert_eq!(direct.state(), c.state());
    assert_eq!(b.event, BridgeEvent::Ruptured);
    near(
        direct.work_on_liquid_j,
        a.work_on_liquid_j + b.work_on_liquid_j + c.work_on_liquid_j,
        2e-14,
        1e-8,
    );
    assert_eq!(direct.irreversible_release_j, b.irreversible_release_j);
    let direct = m.advance(c.state(), near_gap, 100.0, &gate).unwrap();
    let a = m
        .advance(c.state(), m.formation_gap_m(), 100.0, &gate)
        .unwrap();
    let b = m.advance(a.state(), near_gap, 100.0, &gate).unwrap();
    assert_eq!(a.event, BridgeEvent::Reformed);
    assert_eq!(direct.state(), b.state());
    near(
        direct.work_on_liquid_j,
        a.work_on_liquid_j + b.work_on_liquid_j,
        2e-14,
        1e-8,
    );
    let next = f64::from_bits(near_gap.to_bits() + 1);
    let tiny = m.advance(start, next, 1.0, &gate).unwrap();
    assert!(tiny.work_on_liquid_j > 0.0);
    near(
        tiny.mean_normal_force_n,
        m.evaluate(start).unwrap().normal_force_n,
        2e-14,
        0.0,
    );
}

#[test]
fn g3_geometric_and_surface_tension_scaling_are_physical() {
    let m = model(0.7, 0.05, true);
    let h = 0.5 * (m.formation_gap_m() + m.rupture_gap_m());
    let p = m.evaluate(m.initialize(h, true).unwrap()).unwrap();
    for scale in [0.5_f64, 3.0] {
        let (liquid, wetting) = sources(0.7, 0.05, 300.0, true);
        let scaled = SourcePlanarCapillaryBridge::new(
            liquid,
            wetting,
            VOLUME * scale.powi(3),
            WIDTH * scale,
            limits(),
        )
        .unwrap();
        let q = scaled
            .evaluate(scaled.initialize(h * scale, true).unwrap())
            .unwrap();
        near(q.normal_force_n, p.normal_force_n * scale, 2e-14, 0.0);
        near(
            q.surface_free_energy_j,
            p.surface_free_energy_j * scale * scale,
            2e-14,
            1e-8,
        );
        near(
            q.liquid_minus_gas_pressure_pa,
            p.liquid_minus_gas_pressure_pa / scale,
            2e-14,
            0.0,
        );
        near(
            scaled.rupture_gap_m(),
            m.rupture_gap_m() * scale,
            2e-14,
            0.0,
        );
    }
    let stronger = model(0.7, 0.15, true);
    let q = stronger
        .evaluate(stronger.initialize(h, true).unwrap())
        .unwrap();
    near(q.normal_force_n, 3.0 * p.normal_force_n, 3e-14, 0.0);
    near(
        q.surface_free_energy_j,
        3.0 * p.surface_free_energy_j,
        3e-14,
        1e-8,
    );
}

#[test]
fn g0_g4_limits_source_mismatch_and_cancellation_refuse_without_mutating_history() {
    let m = model(0.7, 0.05, true);
    let gate = CancelGate::new_clock_free();
    let h = 0.8 * m.formation_gap_m();
    let state = m.initialize(h, true).unwrap();
    let snapshot = m.evaluate(state).unwrap();
    assert!(m.initialize(m.rupture_gap_m(), true).is_err());
    assert!(m.initialize(0.9 * m.formation_gap_m(), false).is_err());
    assert!(m.advance(state, 2.0 * h, 1e-10, &gate).is_err());
    assert!(m.advance(state, 0.0, 1.0, &gate).is_err());
    assert!(
        m.advance(state, 2.0 * limits().maximum_gap_m, 1.0, &gate)
            .is_err()
    );
    assert!(m.advance(state, h, f64::NAN, &gate).is_err());
    let (liquid, wetting) = sources(0.7, 0.05, 301.0, true);
    let other = SourcePlanarCapillaryBridge::new(liquid, wetting, VOLUME, WIDTH, limits()).unwrap();
    assert!(matches!(
        other.evaluate(state),
        Err(BridgeError::StateMismatch)
    ));
    let (liquid, _) = sources(0.7, 0.05, 300.0, true);
    assert!(
        SourcePlanarCapillaryBridge::new(liquid, other.wetting().clone(), VOLUME, WIDTH, limits())
            .is_err()
    );
    let (liquid, wetting) = sources(0.7, 0.05, 300.0, true);
    let low = SourcePlanarCapillaryBridge::new(
        liquid,
        wetting,
        VOLUME,
        WIDTH,
        BridgeLimits {
            maximum_abs_step_work_j: 1e-12,
            ..limits()
        },
    )
    .unwrap();
    assert!(matches!(
        low.advance(
            low.initialize(h, true).unwrap(),
            0.9 * low.rupture_gap_m(),
            100.0,
            &gate
        ),
        Err(BridgeError::Limit {
            field: "absolute bridge step work",
            ..
        })
    ));
    let (liquid, wetting) = sources(0.7, 0.05, 300.0, true);
    let gravity = SourcePlanarCapillaryBridge::new(
        liquid,
        wetting,
        VOLUME,
        WIDTH,
        BridgeLimits {
            gravity_m_s2: 9.81,
            ..limits()
        },
    )
    .unwrap();
    assert!(gravity.initialize(h, true).is_err());
    for angle in [PI / 2.0, 2.0] {
        let (liquid, wetting) = sources(angle, 0.05, 300.0, false);
        assert!(
            SourcePlanarCapillaryBridge::new(liquid, wetting, VOLUME, WIDTH, limits()).is_err()
        );
    }
    let cancel = CancelGate::new_clock_free();
    cancel.request();
    assert!(matches!(
        m.advance(state, 1.1 * m.rupture_gap_m(), 100.0, &cancel),
        Err(BridgeError::Cancelled)
    ));
    assert_eq!(m.evaluate(state).unwrap(), snapshot);
    assert_eq!(
        m.advance(state, 1.1 * m.rupture_gap_m(), 100.0, &gate)
            .unwrap(),
        m.advance(state, 1.1 * m.rupture_gap_m(), 100.0, &gate)
            .unwrap()
    );
}

#[test]
fn g1_original_adhesion_and_liquid_sources_are_retained_without_sample_copies() {
    let m = model(0.7, 0.05, true);
    let liquid_id = m.liquid().material().identity();
    let liquid_address = m.liquid().material().properties().as_ptr();
    let surface_addresses = m
        .wetting()
        .each_ref()
        .map(|w| w.resolved().properties().as_ptr());
    assert_ne!(
        m.wetting()[0].resolved().identity(),
        m.wetting()[1].resolved().identity()
    );
    for wall in m.wetting() {
        assert_eq!(wall.angle_origin(), EquilibriumAngleOrigin::YoungDupre);
        assert!(
            wall.resolved()
                .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                .is_some()
        );
        assert!(
            wall.resolved()
                .property(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY)
                .is_none()
        );
        assert_eq!(
            wall.resolved().query_point(),
            m.liquid().material().query_point()
        );
    }
    let gate = CancelGate::new_clock_free();
    let mut state = m.initialize(1.1 * m.rupture_gap_m(), false).unwrap();
    for h in [
        0.8 * m.formation_gap_m(),
        1.1 * m.rupture_gap_m(),
        0.8 * m.formation_gap_m(),
    ] {
        state = m.advance(state, h, 100.0, &gate).unwrap().state();
        assert_eq!(m.liquid().material().identity(), liquid_id);
        assert_eq!(m.liquid().material().properties().as_ptr(), liquid_address);
        assert_eq!(
            m.wetting()
                .each_ref()
                .map(|w| w.resolved().properties().as_ptr()),
            surface_addresses
        );
    }
}

#[test]
fn g1_source_refresh_changes_surface_force_and_accounts_parameter_work_on_both_branches() {
    let gate = CancelGate::new_clock_free();
    for young in [false, true] {
        let (material, interfaces) = temperature_cards(young, 1.0, 1.0);
        let (liquid, wetting) = resolved_sources(&material, &interfaces, 300.0, young);
        let source =
            SourcePlanarCapillaryBridge::new(liquid, wetting, VOLUME, WIDTH, limits()).unwrap();
        for attached in [true, false] {
            let gap = if attached {
                0.8 * source.formation_gap_m()
            } else {
                1.2 * source.rupture_gap_m()
            };
            let state = source.initialize(gap, attached).unwrap();
            let before = source.evaluate(state).unwrap();
            let refresh = source
                .refresh_at_query_point(
                    state,
                    &material,
                    interfaces.each_ref(),
                    &QueryPoint::new().with("T", 320.0).unwrap(),
                    &gate,
                )
                .unwrap();
            let after = refresh.source().evaluate(refresh.state()).unwrap();
            assert_eq!(refresh.previous_state(), state);
            assert_eq!(refresh.state().gap_m(), gap);
            assert_eq!(refresh.state().attached(), attached);
            assert_eq!(after.liquid_volume_m3, before.liquid_volume_m3);
            assert_eq!(
                refresh.source().liquid().density_kg_m3(),
                source.liquid().density_kg_m3()
            );
            near(
                after.normal_force_n,
                0.8 * before.normal_force_n,
                2e-14,
                1e-8,
            );
            near(
                after.surface_free_energy_j,
                0.8 * before.surface_free_energy_j,
                2e-14,
                1e-10,
            );
            near(
                refresh.parameter_work_j(),
                -0.2 * before.surface_free_energy_j,
                2e-14,
                1e-10,
            );
            assert_eq!(refresh.source().liquid().dynamic_viscosity_pa_s(), 2e-3);
            assert_ne!(refresh.source().identity(), source.identity());
            for (old, new) in source
                .liquid()
                .material()
                .properties()
                .iter()
                .zip(refresh.source().liquid().material().properties())
            {
                assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
            }
            for wall in refresh.source().wetting() {
                assert_eq!(
                    wall.resolved().query_point(),
                    refresh.source().liquid().material().query_point()
                );
            }
            let held = refresh
                .source()
                .advance(refresh.state(), gap, 1.0, &gate)
                .unwrap();
            assert_eq!(held.work_on_liquid_j, 0.0);
            assert_eq!(held.irreversible_release_j, 0.0);
            let back = refresh
                .source()
                .refresh_at_query_point(
                    refresh.state(),
                    &material,
                    interfaces.each_ref(),
                    &QueryPoint::new().with("T", 300.0).unwrap(),
                    &gate,
                )
                .unwrap();
            assert_eq!(back.state().gap_m(), state.gap_m());
            assert_eq!(back.state().attached(), state.attached());
            assert_eq!(back.source().evaluate(back.state()).unwrap(), before);
            assert_eq!(back.parameter_work_j(), -refresh.parameter_work_j());
            assert_eq!(source.evaluate(state).unwrap(), before);
        }
    }
}

#[test]
fn g3_source_refresh_resolves_wetting_and_refuses_mass_topology_source_or_domain_substitution() {
    let gate = CancelGate::new_clock_free();
    for young in [false, true] {
        let (material, interfaces) = temperature_cards(young, 1.0, 2.0);
        let build = |temperature| {
            let (liquid, wetting) = resolved_sources(&material, &interfaces, temperature, young);
            SourcePlanarCapillaryBridge::new(liquid, wetting, VOLUME, WIDTH, limits()).unwrap()
        };
        let cold = build(300.0);
        let hot = build(320.0);
        let gap = 0.5 * cold.formation_gap_m().min(hot.formation_gap_m());
        let state = cold.initialize(gap, true).unwrap();
        let hot_point = QueryPoint::new().with("T", 320.0).unwrap();
        let refresh = cold
            .refresh_at_query_point(state, &material, interfaces.each_ref(), &hot_point, &gate)
            .unwrap();
        assert_eq!(
            refresh.source().evaluate(refresh.state()).unwrap(),
            hot.evaluate(hot.initialize(gap, true).unwrap()).unwrap()
        );
        assert!(
            refresh.source().wetting()[0].contact_angle_radians()
                > cold.wetting()[0].contact_angle_radians()
        );

        let (large, small, to) = if cold.rupture_gap_m() > hot.rupture_gap_m() {
            (&cold, &hot, 320.0)
        } else {
            (&hot, &cold, 300.0)
        };
        let state = large
            .initialize(0.5 * (large.rupture_gap_m() + small.rupture_gap_m()), true)
            .unwrap();
        assert!(
            large
                .refresh_at_query_point(
                    state,
                    &material,
                    interfaces.each_ref(),
                    &QueryPoint::new().with("T", to).unwrap(),
                    &gate
                )
                .is_err()
        );
        let (small, large, to) = if cold.formation_gap_m() < hot.formation_gap_m() {
            (&cold, &hot, 320.0)
        } else {
            (&hot, &cold, 300.0)
        };
        let state = small
            .initialize(
                0.5 * (small.formation_gap_m() + large.formation_gap_m()),
                false,
            )
            .unwrap();
        assert!(
            small
                .refresh_at_query_point(
                    state,
                    &material,
                    interfaces.each_ref(),
                    &QueryPoint::new().with("T", to).unwrap(),
                    &gate
                )
                .is_err()
        );
    }

    let (material, interfaces) = temperature_cards(false, 0.99, 1.0);
    let (liquid, wetting) = resolved_sources(&material, &interfaces, 300.0, false);
    let source =
        SourcePlanarCapillaryBridge::new(liquid, wetting, VOLUME, WIDTH, limits()).unwrap();
    let state = source
        .initialize(0.8 * source.formation_gap_m(), true)
        .unwrap();
    let before = source.evaluate(state).unwrap();
    assert!(matches!(
        source.refresh_at_query_point(
            state,
            &material,
            interfaces.each_ref(),
            &QueryPoint::new().with("T", 320.0).unwrap(),
            &gate
        ),
        Err(BridgeError::Invalid { .. })
    ));
    assert!(
        source
            .refresh_at_query_point(
                state,
                &material,
                interfaces.each_ref(),
                &QueryPoint::new().with("T", 350.0).unwrap(),
                &gate
            )
            .is_err()
    );
    let (other, _) = source_cards(0.7, 0.05, false);
    assert!(
        source
            .refresh_at_query_point(
                state,
                &other,
                interfaces.each_ref(),
                &QueryPoint::new().with("T", 300.0).unwrap(),
                &gate
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        source.refresh_at_query_point(
            state,
            &material,
            interfaces.each_ref(),
            &QueryPoint::new().with("T", 300.0).unwrap(),
            &cancelled
        ),
        Err(BridgeError::Cancelled)
    ));
    assert_eq!(source.evaluate(state).unwrap(), before);
}
