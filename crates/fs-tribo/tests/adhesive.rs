//! G0/G1/G3: junction shear, detached limit, passivity and area scaling.
use fs_tribo::adhesive::{AdhesiveContactState, AdhesiveFrictionLaw};
use fs_tribo::{
    ContactFrame, FrictionLaw, FrictionRegime, InputAuthority, InterfaceMedium, InterfaceSystemRef,
    TangentialSlip,
};

fn interface() -> InterfaceSystemRef {
    InterfaceSystemRef::new(
        "adhesive-A/B",
        "bond-1",
        "synthetic-junction",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap()
}
fn pressure() -> FrictionLaw {
    FrictionLaw::Coulomb {
        static_mu: 0.6,
        kinetic_mu: 0.4,
    }
}
fn slip(v: [f64; 3]) -> TangentialSlip {
    TangentialSlip::new(&ContactFrame::new([0.0, 0.0, 1.0]).unwrap(), v).unwrap()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-14 * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

#[test]
fn attached_junction_transmits_shear_without_positive_net_load() {
    // A normal owner may supply this finite area while its NET normal force
    // is zero or tensile. No net-load clipping belongs in the shear law.
    let law = AdhesiveFrictionLaw::new(
        FrictionLaw::Coulomb {
            static_mu: 0.0,
            kinetic_mu: 0.0,
        },
        3e6,
        2e6,
    )
    .unwrap();
    let contact = AdhesiveContactState::new(2e-6, 0.0).unwrap();
    let identity = interface();
    let at_rest = law.evaluate(&identity, contact, slip([0.0; 3])).unwrap();
    assert_eq!(at_rest.regime, FrictionRegime::Sticking);
    close(at_rest.static_capacity_n, 6.0);
    assert_eq!(at_rest.traction_n, [0.0; 3]);
    assert_eq!(at_rest.dissipated_power_w, 0.0);
    let moving = law
        .evaluate(&identity, contact, slip([0.03, 0.04, 0.0]))
        .unwrap();
    close(moving.traction_n[0], -2.4);
    close(moving.traction_n[1], -3.2);
    close(moving.dissipated_power_w, 0.2);
    assert!(std::ptr::eq(moving.provenance, identity.provenance()));
}

#[test]
fn detached_contact_has_no_capacity_or_frictional_heat() {
    let law = AdhesiveFrictionLaw::new(pressure(), 3e6, 2e6).unwrap();
    let identity = interface();
    let response = law
        .evaluate(
            &identity,
            AdhesiveContactState::new(0.0, 0.0).unwrap(),
            slip([2.0, -1.0, 0.0]),
        )
        .unwrap();
    assert_eq!(response.static_capacity_n, 0.0);
    assert_eq!(response.traction_n, [0.0; 3]);
    assert_eq!(response.dissipated_power_w, 0.0);
}

#[test]
fn zero_junction_strength_recovers_existing_pressure_law() {
    let identity = interface();
    for speed in [1e-8, 0.02, 2.0] {
        for pressure in [
            pressure(),
            FrictionLaw::Hyperbolic {
                static_mu: 0.8,
                kinetic_mu: 0.3,
                characteristic_speed: 0.01,
            },
        ] {
            let law = AdhesiveFrictionLaw::new(pressure.clone(), 0.0, 0.0).unwrap();
            let velocity = slip([speed, 0.0, 0.0]);
            let baseline = pressure.evaluate(&identity, 8.0, velocity).unwrap();
            let actual = law
                .evaluate(
                    &identity,
                    AdhesiveContactState::new(1e-6, 8.0).unwrap(),
                    velocity,
                )
                .unwrap();
            close(actual.static_capacity_n, baseline.static_limit);
            close(actual.traction_n[0], baseline.traction_n()[0]);
            close(actual.dissipated_power_w, baseline.dissipated_power_w());
        }
    }
}

#[test]
fn area_and_load_scale_extensively_and_power_equals_negative_slip_work() {
    let law = AdhesiveFrictionLaw::new(pressure(), 3e6, 2e6).unwrap();
    let identity = interface();
    let c = AdhesiveContactState::new(1e-6, 5.0).unwrap();
    for v in [[0.3, 0.4, 0.0], [-0.4, 0.3, 0.0], [0.0, -0.5, 0.0]] {
        let a = law.evaluate(&identity, c, slip(v)).unwrap();
        let b = law
            .evaluate(
                &identity,
                AdhesiveContactState::new(7e-6, 35.0).unwrap(),
                slip(v),
            )
            .unwrap();
        close(
            a.dissipated_power_w,
            -a.traction_n.iter().zip(v).map(|(f, x)| f * x).sum::<f64>(),
        );
        assert!(a.dissipated_power_w >= 0.0);
        close(b.static_capacity_n, 7.0 * a.static_capacity_n);
        for i in 0..3 {
            close(b.traction_n[i], 7.0 * a.traction_n[i]);
        }
    }
}

#[test]
fn scalar_regularization_is_continuous_odd_and_dissipative() {
    let law = AdhesiveFrictionLaw::new(pressure(), 3e6, 2e6).unwrap();
    let c = AdhesiveContactState::new(2e-6, 5.0).unwrap();
    let scale = 0.01;
    let limit = law.sliding_capacity_n(c, scale).unwrap();
    for v in [0.0, 1e-12, 0.005, scale, 0.02] {
        let force = law.regularized_traction_1d(c, v, scale).unwrap();
        close(force, -law.regularized_traction_1d(c, -v, scale).unwrap());
        assert!(-force * v >= 0.0);
        close(force, -limit * (v / scale).min(1.0));
    }
    let a = law
        .regularized_traction_1d(c, scale * (1.0 - 1e-10), scale)
        .unwrap();
    let b = law
        .regularized_traction_1d(c, scale * (1.0 + 1e-10), scale)
        .unwrap();
    assert!((a - b).abs() < 2e-10 * limit);
}

#[test]
fn invalid_inputs_refuse_and_subnormal_slip_does_not_overflow_force() {
    assert!(AdhesiveContactState::new(0.0, 1.0).is_err());
    assert!(AdhesiveContactState::new(1.0, -1.0).is_err());
    assert!(AdhesiveContactState::new(f64::NAN, 0.0).is_err());
    assert!(AdhesiveFrictionLaw::new(pressure(), 1.0, 2.0).is_err());
    let law = AdhesiveFrictionLaw::new(pressure(), 3.0, 2.0).unwrap();
    let identity = interface();
    let r = law
        .evaluate(
            &identity,
            AdhesiveContactState::new(1.0, 0.0).unwrap(),
            slip([f64::MIN_POSITIVE * 0.5, 0.0, 0.0]),
        )
        .unwrap();
    close(r.traction_n[0], -2.0);
    assert!(r.dissipated_power_w > 0.0);
    assert!(
        law.static_capacity_n(AdhesiveContactState::new(f64::MAX, 0.0).unwrap())
            .is_err()
    );
}
