//! Manufactured reference inventories exercise the native grid and couplings.
use super::*;
use crate::{Edge, Instrument, Nonlinear};
use std::cell::Cell;

fn near(actual: f64, expected: f64, atol: f64, rtol: f64) {
    assert!(
        (actual - expected).abs() <= atol + rtol * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}

fn spec(thickness: &dyn Fn(f64) -> f64) -> PlateSpec<'_> {
    PlateSpec {
        radius: 0.04,
        hole: 0.0,
        grid: 0.01,
        edge: Edge::Clamped,
        thickness,
        rest_height: &|_| 0.0,
        young: 2e6,
        poisson: 0.3,
        density: 1000.0,
        tension: 1.0,
        bending: true,
        nonlinear: Nonlinear::None,
        t60: [(100.0, 0.2), (1000.0, 0.1)],
        support: None,
        muffle: None,
    }
}

fn half_velocity(body: &Body, dt: f64) -> Vec<f64> {
    body.w
        .iter()
        .zip(&body.wp)
        .map(|(w, wp)| (w - wp) / dt)
        .collect()
}

#[test]
fn field_is_sampled_once_and_constant_fields_preserve_original_operators() {
    let thickness = |r: f64| 0.001 * (1.0 + 2.0 * r);
    let s = spec(&thickness);
    let original = Body::plate_without_surface_loss("scalar", &s).unwrap();
    let calls = Cell::new(0);
    let density = |_: [f64; 2]| {
        calls.set(calls.get() + 1);
        1000.0
    };
    let mut field = Body::plate_with_reference_density(
        "field",
        &s,
        PlateReferenceDensity::Volumetric(&density),
    )
    .unwrap();
    let areal = |[x, y]: [f64; 2]| 1000.0 * thickness((x * x + y * y).sqrt());
    let integrated =
        Body::plate_with_reference_density("areal", &s, PlateReferenceDensity::Areal(&areal))
            .unwrap();
    assert_eq!(calls.get(), field.n);
    assert_eq!(original.mass, field.mass);
    assert_eq!(field.mass, integrated.mass);
    assert_eq!(original.xy, field.xy);
    assert_eq!(original.stiff.val, field.stiff.val);
    assert_eq!(original.lap.val, field.lap.val);
    assert!(!field.surface_t60);
    assert_eq!(field.s1, 0.0);
    assert!(field.s0.iter().all(|s| *s == 0.0));
    for (i, &[x, y]) in field.xy.iter().enumerate() {
        let old_mass = s.density * thickness((x * x + y * y).sqrt()) * s.grid * s.grid;
        assert_eq!(field.mass[i].to_bits(), old_mass.to_bits());
        field.w[i] = 1e-7 * (0.2 * i as f64).sin();
        field.wp[i] = field.w[i];
    }
    field.sync();
    let mut instrument = Instrument::new(1e6);
    instrument.add(field);
    assert!(instrument.stability().0 < 1.0);
    let mut replay = instrument.clone();
    for _ in 0..3 {
        instrument.try_step(&mut []).unwrap();
        replay.try_step(&mut []).unwrap();
    }
    assert_eq!(instrument.bodies[0].w, replay.bodies[0].w);
    assert_eq!(calls.get(), original.n);
    assert_eq!(instrument.bodies[0].mass, original.mass);
}

#[test]
fn light_and_redistributed_density_reach_the_actual_stiffness_mass_clock() {
    let thickness = |_| 1e-5;
    let mut s = spec(&thickness);
    s.bending = false;
    s.density = 1e-5;
    let original = Body::plate_without_surface_loss("light scalar", &s).unwrap();
    assert!(original.mass.iter().all(|m| *m > 0.0 && *m < 1e-12));
    for mass in &original.mass {
        assert_eq!(*mass, s.density * 1e-5 * s.grid * s.grid);
    }
    let rho = |[x, _]: [f64; 2]| {
        1e-5 * if x < 0.0 {
            0.6
        } else if x > 0.0 {
            1.4
        } else {
            1.0
        }
    };
    let denser = |p: [f64; 2]| 4.0 * rho(p);
    let graded = Body::plate_with_reference_density(
        "light graded",
        &s,
        PlateReferenceDensity::Volumetric(&rho),
    )
    .unwrap();
    let heavy = Body::plate_with_reference_density(
        "fourfold graded",
        &s,
        PlateReferenceDensity::Volumetric(&denser),
    )
    .unwrap();
    let inventory: f64 = original.mass.iter().sum();
    near(graded.mass.iter().sum(), inventory, 0.0, 2e-14);
    let first_moment: f64 = graded
        .xy
        .iter()
        .zip(&graded.mass)
        .map(|(p, m)| p[0] * m)
        .sum();
    assert!(first_moment > 0.0);
    assert_eq!(graded.stiff.val, original.stiff.val);
    assert_eq!(graded.lap.val, original.lap.val);
    let admission = |body: Body, rate| {
        let mut ins = Instrument::new(rate);
        ins.add(body);
        ins.stability().0
    };
    let scalar_fraction = admission(original, 1e4);
    let graded_fraction = admission(graded.clone(), 1e4);
    assert!(graded_fraction.is_finite() && graded_fraction > 1.0);
    assert!((graded_fraction / scalar_fraction - 1.0).abs() > 0.01);
    near(admission(heavy, 1e4), 0.25 * graded_fraction, 0.0, 3e-12);
    let resolved_fraction = admission(graded, 1e9);
    near(resolved_fraction, graded_fraction * 1e-10, 0.0, 3e-12);
    assert!(resolved_fraction < 1.0);
}

#[test]
fn tiny_physical_mass_cannot_collapse_a_large_clock_rate_to_zero() {
    let thickness = |_| 0.001;
    let mut s = spec(&thickness);
    s.bending = false;
    s.tension = 1e-10;
    let ordinary_density = |_: [f64; 2]| 1.0;
    let tiny_density = |_: [f64; 2]| 1e-190;
    let mut ordinary = Body::plate_with_reference_density(
        "ordinary metric",
        &s,
        PlateReferenceDensity::Areal(&ordinary_density),
    )
    .unwrap();
    let mut tiny = Body::plate_with_reference_density(
        "tiny metric",
        &s,
        PlateReferenceDensity::Areal(&tiny_density),
    )
    .unwrap();
    let ordinary_rate = ordinary.max_rate2();
    let original_position = tiny.w.clone();
    let original_gradient = tiny.g.clone();
    let tiny_rate = tiny.max_rate2();
    assert!(tiny_rate.is_finite() && tiny_rate > 1e180);
    near(tiny_rate / ordinary_rate, 1e190, 0.0, 5e-12);
    near(tiny.max_rate2_with_shell(), tiny_rate, 0.0, 5e-12);
    assert_eq!(tiny.w, original_position);
    assert_eq!(tiny.g, original_gradient);
    let mut clock = Instrument::new(1e90);
    clock.add(tiny.clone());
    assert!(clock.stability().0 > 1.0);
    clock.rate = 1e94;
    assert!(clock.stability().0 < 1.0);

    // The next case has valid positive masses but an unrepresentable physical
    // inverse-mass operator. Both estimates must retain a refusal-sized rate.
    s.tension = 1.0;
    let unresolved = |_: [f64; 2]| 1e-310;
    let mut overflow = Body::plate_with_reference_density(
        "overflowed operator",
        &s,
        PlateReferenceDensity::Areal(&unresolved),
    )
    .unwrap();
    let saved = overflow.w.clone();
    assert_eq!(overflow.max_rate2(), f64::INFINITY);
    assert_eq!(overflow.max_rate2_with_shell(), f64::INFINITY);
    assert_eq!(overflow.w, saved);
    let mut clock = Instrument::new(1e100);
    clock.add(overflow);
    assert_eq!(clock.stability().0, f64::INFINITY);
    clock.rate = 1e200;
    assert_eq!(clock.stability().0, f64::INFINITY);
    near(spectral_norm(&[3e-200, 4e-200]), 5e-200, 0.0, 2e-15);
    near(spectral_norm(&[3e180, 4e180]), 5e180, 0.0, 2e-15);
}

#[test]
fn graded_mass_has_the_original_acceleration_impulse_and_kinetic_balance() {
    let thickness = |_| 0.001;
    let mut s = spec(&thickness);
    s.edge = Edge::Free;
    s.bending = false;
    s.tension = 0.0;
    // An areal inventory is authoritative even when the unused scalar density
    // has no value; neither thickness nor the old scalar is applied twice.
    s.density = f64::NAN;
    let sigma = |[x, y]: [f64; 2]| 0.8 + 4.0 * x + 2.0 * y;
    let mut body =
        Body::plate_with_reference_density("graded free", &s, PlateReferenceDensity::Areal(&sigma))
            .unwrap();
    let (dt, acceleration, initial_velocity) = (1e-4, 0.7, 0.03);
    for i in 0..body.n {
        body.wp[i] = -dt * initial_velocity + 0.5 * dt * dt * acceleration;
        body.bias[i] = body.mass[i] * acceleration;
    }
    body.sync();
    let total_mass: f64 = body.mass.iter().sum();
    let initial_momentum = total_mass * (initial_velocity - 0.5 * dt * acceleration);
    let mut instrument = Instrument::new(1.0 / dt);
    instrument.add(body);
    assert_eq!(instrument.stability().0, 0.0);
    for step in 1..=4 {
        instrument.try_step(&mut []).unwrap();
        let body = &instrument.bodies[0];
        let time = step as f64 * dt;
        let v = initial_velocity + (time - 0.5 * dt) * acceleration;
        let velocities = half_velocity(body, dt);
        let momentum: f64 = body.mass.iter().zip(&velocities).map(|(m, v)| m * v).sum();
        near(
            momentum - initial_momentum,
            time * total_mass * acceleration,
            1e-17,
            3e-12,
        );
        near(
            body.discrete_energy(dt),
            0.5 * total_mass * v * v,
            1e-18,
            3e-12,
        );
        for (u, actual_v) in body.w.iter().zip(&velocities) {
            near(
                *u,
                time * initial_velocity + 0.5 * time * time * acceleration,
                1e-17,
                3e-12,
            );
            near(*actual_v, v, 1e-14, 3e-12);
        }
    }
}

#[test]
fn contact_force_response_and_attached_air_use_the_same_graded_inverse() {
    let thickness = |_| 0.001;
    let mut s = spec(&thickness);
    s.edge = Edge::Free;
    s.bending = false;
    s.tension = 0.0;
    let sigma = |[x, y]: [f64; 2]| 1.0 + 5.0 * x - 3.0 * y;
    let mut body = Body::plate_with_reference_density(
        "air-loaded graded plate",
        &s,
        PlateReferenceDensity::Areal(&sigma),
    )
    .unwrap();
    body.air_mass = 0.003;
    let dt = 1e-5;
    let force = 0.7;
    let node = body.node_at(0.0, 0.0).unwrap();
    body.linear_forces();
    body.f[node] += force;
    body.predict(dt);
    let predicted = body.wn.clone();
    // The exact original rank-one air inertia is gamma * 1*1^T. Its inverse
    // is independent of the native response implementation used by contacts.
    let gamma = body.air_mass * body.h.powi(4) / body.area.powi(2);
    let reciprocal_sum: f64 = body.mass.iter().map(|m| 1.0 / m).sum();
    let coefficient = gamma * (force / body.mass[node]) / (1.0 + gamma * reciprocal_sum);
    for i in 0..body.n {
        let original_force = if i == node { force } else { 0.0 };
        let expected = dt * dt * (original_force - coefficient) / body.mass[i];
        near(predicted[i], expected, 1e-20, 3e-12);
    }
    body.footprint_response(0, &[(0, node, 1.0)], force, &[], 0.0, dt);
    for (contact, predicted) in body.tmp.iter().zip(&predicted) {
        near(*contact, *predicted, 1e-20, 3e-12);
    }
    body.finish(dt);
    body.rotate();
    let velocities = half_velocity(&body, dt);
    let sum: f64 = velocities.iter().sum();
    let momentum: f64 = body
        .mass
        .iter()
        .zip(&velocities)
        .map(|(m, v)| m * v)
        .sum::<f64>()
        + body.n as f64 * gamma * sum;
    near(momentum, dt * force, 1e-18, 3e-12);
    near(
        body.discrete_energy(dt),
        0.5 * dt * force * velocities[node],
        1e-20,
        3e-12,
    );
}

fn viscosity() -> PlaneStressViscosity {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::plane_stress_viscosity::{MaterialPlane, PLANE_STRESS_VISCOSITY_PROPERTIES};
    use fs_material::state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    };
    use fs_qty::DynViscosity;
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value) in PLANE_STRESS_VISCOSITY_PROPERTIES
        .into_iter()
        .zip([8e4, 6e4, 3e4, 1e4, 4e3, -3e3])
    {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, DynViscosity::DIMS),
                value: PropertyValue::Scalar {
                    value,
                    dims: DynViscosity::DIMS,
                },
                validity: ValidityDomain::unconstrained().with("omega", 0.0, 1e9),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic graded plate viscosity".into(),
                    license: "test".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(
                name,
                DynViscosity::DIMS,
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic plate".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let source = resolve_material_state_point(
        &card,
        &QueryPoint::new().with("omega", 1000.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PlaneStressViscosity::from_material(
        source,
        MaterialPlane::Axes12,
        PLANE_STRESS_VISCOSITY_PROPERTIES,
    )
    .unwrap()
}

#[test]
fn sourced_and_support_losses_keep_physical_work_with_spatial_mass() {
    let thickness = |_| 0.001;
    let s = spec(&thickness);
    let source = viscosity();
    let rho = |[x, y]: [f64; 2]| 1000.0 * (1.0 + 5.0 * x + 2.0 * y);
    let mut body = Body::plate_with_reference_density_and_viscosity(
        "viscous graded",
        &s,
        PlateReferenceDensity::Volumetric(&rho),
        &source,
        0.23,
    )
    .unwrap();
    assert_eq!(body.bending_viscosity(), Some((&source, 0.23)));
    let mut clock = Instrument::new(1e6);
    clock.add(body.clone());
    clock.validate_bending_loss_clock().unwrap();
    assert!(clock.stability().0 < 1.0);
    let dt = clock.k();
    for i in 0..body.n {
        body.w[i] = 1e-6 * (0.31 * i as f64).sin();
        body.wp[i] = body.w[i] - dt * 0.001 * (0.29 * i as f64).cos();
    }
    body.sync();
    let initial = body.discrete_energy(dt);
    for _ in 0..8 {
        body.linear_forces();
        body.predict(dt);
        body.finish(dt);
        body.rotate();
        near(
            body.discrete_energy(dt) + body.bending_loss_dissipated_j(),
            initial,
            1e-19,
            2e-11,
        );
    }
    assert!(body.bending_loss_dissipated_j() > 0.0);

    let mut supported = spec(&thickness);
    supported.edge = Edge::Free;
    supported.bending = false;
    supported.tension = 0.0;
    supported.support = Some((0.0, 0.1, 0.0, 20.0));
    let mut body = Body::plate_with_reference_density(
        "graded dashpots",
        &supported,
        PlateReferenceDensity::Volumetric(&rho),
    )
    .unwrap();
    for i in 0..body.n {
        body.wp[i] = -dt * 0.01;
    }
    body.sync();
    let before = body.discrete_energy(dt);
    body.linear_forces();
    body.predict(dt);
    let mut work = 0.0;
    for i in 0..body.n {
        let drag = 20.0 * body.h * body.h;
        near(2.0 * body.mass[i] * body.s0[i], drag, 0.0, 2e-14);
        let centered = (body.wn[i] - body.wp[i]) / (2.0 * dt);
        work += dt * drag * centered * centered;
    }
    body.finish(dt);
    body.rotate();
    near(body.discrete_energy(dt) + work, before, 1e-20, 2e-11);
}

#[test]
fn invalid_or_unrepresentable_inventory_refuses_without_a_mass_floor() {
    let thickness = |_| 0.001;
    let mut s = spec(&thickness);
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        s.density = bad;
        assert!(Body::plate_without_surface_loss("invalid scalar", &s).is_err());
        let sample = |_: [f64; 2]| bad;
        for field in [
            PlateReferenceDensity::Volumetric(&sample),
            PlateReferenceDensity::Areal(&sample),
        ] {
            assert!(Body::plate_with_reference_density("invalid field", &s, field).is_err());
        }
    }
    s.density = f64::from_bits(1);
    assert!(
        Body::plate_without_surface_loss("underflow", &s)
            .err()
            .unwrap()
            .contains("reference mass")
    );
    let tiny = |_: [f64; 2]| f64::from_bits(1);
    assert!(
        Body::plate_with_reference_density(
            "underflow field",
            &s,
            PlateReferenceDensity::Areal(&tiny)
        )
        .err()
        .unwrap()
        .contains("reference mass")
    );
    s.radius = 7.0;
    s.grid = 2.0;
    let overflow = |_: [f64; 2]| f64::MAX;
    assert!(
        Body::plate_with_reference_density(
            "overflow field",
            &s,
            PlateReferenceDensity::Areal(&overflow)
        )
        .err()
        .unwrap()
        .contains("reference mass")
    );
}
