//! The original source bundle feeds every parameter of the actual speed law.

use super::*;
use fs_contact::interface_binding::thermal_friction::rate_dependent::{
    SourceThermalFriction, ThermalFrictionCurve,
};
use fs_tribo::FrictionLaw;
use fs_tribo::sliding_patch::rate_dependent::{HertzFrictionPatch, HertzFrictionQuadrature};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const SLOPE: Dims = Dims([-1, 0, 1, 0, 0, 0]);

fn parameter(name: &str, dims: Dims, values: [f64; 2]) -> PropertyClaim {
    let mut claim = curve(name, values);
    claim.key = PropertyKey::new(name, dims);
    claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, values[0]), (320.0, values[1])],
        dims,
    };
    claim
}

fn rate_card() -> InterfaceSystemCard {
    card(vec![
        curve(STATIC, [0.9, 0.6]),
        curve(KINETIC, [0.4, 0.2]),
        parameter("speed", SPEED, [0.08, 0.12]),
        parameter("viscous", SLOPE, [0.1, 0.3]),
        parameter("slope", SLOPE, [-0.3, -0.1]),
        parameter("exponent", Dims::NONE, [0.6, 0.9]),
    ])
}

fn typed(name: &str, dims: Dims) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
}

fn curves() -> [ThermalFrictionCurve; 5] {
    [
        ThermalFrictionCurve::Coulomb {
            static_mu: requirement(STATIC),
            kinetic_mu: requirement(KINETIC),
        },
        ThermalFrictionCurve::VelocityDependent {
            static_mu: requirement(STATIC),
            mu_zero: requirement(KINETIC),
            slope_per_speed: typed("slope", SLOPE),
        },
        ThermalFrictionCurve::Stribeck {
            static_mu: requirement(STATIC),
            kinetic_mu: requirement(KINETIC),
            characteristic_speed: typed("speed", SPEED),
            viscous_per_speed: typed("viscous", SLOPE),
        },
        ThermalFrictionCurve::Hyperbolic {
            static_mu: requirement(STATIC),
            kinetic_mu: requirement(KINETIC),
            characteristic_speed: typed("speed", SPEED),
        },
        ThermalFrictionCurve::ShiftedPower {
            static_mu: requirement(STATIC),
            kinetic_mu: requirement(KINETIC),
            characteristic_speed: typed("speed", SPEED),
            exponent: typed("exponent", Dims::NONE),
        },
    ]
}

fn source(card: &InterfaceSystemCard, curve: ThermalFrictionCurve) -> SourceThermalFriction {
    SourceThermalFriction::resolve(
        bind_dry_interface_system_card(card, InputAuthority::SyntheticFixture).unwrap(),
        query(280.0),
        curve,
        MaterialPropertySelection::SingleClaimOnly,
        axes(),
    )
    .unwrap()
}

#[test]
fn all_original_curves_refresh_complete_parameters_and_actual_mixed_hertz_traction() {
    let card = rate_card();
    let quadrature = HertzFrictionQuadrature {
        radial_panels: 8,
        angular_panels: 8,
    };
    for (family, curve) in curves().into_iter().enumerate() {
        let source = source(&card, curve);
        let saved = source.clone();
        for (temperature, s, k, speed, viscous, slope, exponent) in [
            (280.0, 0.9, 0.4, 0.08, 0.1, -0.3, 0.6),
            (320.0, 0.6, 0.2, 0.12, 0.3, -0.1, 0.9),
        ] {
            let expected = match family {
                0 => FrictionLaw::Coulomb {
                    static_mu: s,
                    kinetic_mu: k,
                },
                1 => FrictionLaw::VelocityDependent {
                    static_mu: s,
                    mu_zero: k,
                    slope_per_speed: slope,
                },
                2 => FrictionLaw::Stribeck {
                    static_mu: s,
                    kinetic_mu: k,
                    characteristic_speed: speed,
                    viscous_per_speed: viscous,
                },
                3 => FrictionLaw::Hyperbolic {
                    static_mu: s,
                    kinetic_mu: k,
                    characteristic_speed: speed,
                },
                _ => FrictionLaw::ShiftedPower {
                    static_mu: s,
                    kinetic_mu: k,
                    characteristic_speed: speed,
                    exponent,
                },
            };
            let binding = source.hertz_at(temperature, None, 10.0, 0.02).unwrap();
            assert_eq!(binding.law().friction_law(), &expected);
            assert_eq!(binding.owner_identity(), source.identity());
            assert_eq!(
                binding.friction_source().state_coordinate("T"),
                Some(temperature)
            );
            for property in binding.friction_source().properties() {
                let initial = source
                    .initial_source()
                    .property(property.requirement().name())
                    .unwrap();
                assert_eq!(
                    property.answer().receipt.selected,
                    initial.answer().receipt.selected
                );
                card.claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
            let actual = binding
                .law()
                .evaluate([0.07, -0.03], 4.0, quadrature)
                .unwrap();
            let expected = HertzFrictionPatch::new(10.0, 0.02, expected)
                .unwrap()
                .evaluate([0.07, -0.03], 4.0, quadrature)
                .unwrap();
            assert_eq!(actual, expected);
        }
        assert_eq!(source, saved);
        assert_eq!(
            source.interface().interface().provenance().authority(),
            InputAuthority::SyntheticFixture
        );
    }
}

#[test]
fn coulomb_limit_preserves_original_pressure_resolved_force_moment_and_heat() {
    let card = rate_card();
    let old = resolve_at(&card, query(280.0)).unwrap();
    let source = source(&card, curves()[0].clone());
    let cold = source.hertz_at(280.0, None, 10.0, 0.02).unwrap();
    let old = old.hertz_at(280.0, None, 10.0, 0.02).unwrap();
    let actual = cold
        .law()
        .evaluate(
            [0.003, -0.004],
            0.25,
            HertzFrictionQuadrature {
                radial_panels: 1,
                angular_panels: 1,
            },
        )
        .unwrap();
    let parameters = old.law().parameters();
    let capacity = parameters.normal_force_n * parameters.kinetic_mu;
    // Independent circular-Hertz pressure moments at |u| = a*|omega|.
    let force = capacity * 9.0 * core::f64::consts::PI / 32.0;
    let moment = capacity * parameters.contact_radius_m * 9.0 * core::f64::consts::PI / 128.0;
    close(actual.force_n[0], -0.6 * force);
    close(actual.force_n[1], 0.8 * force);
    close(actual.spin_moment_nm, -moment);
    close(actual.dissipated_power_w, force * 0.005 + moment * 0.25);
}

#[test]
fn typed_roles_source_support_and_native_physical_domains_cannot_be_bypassed() {
    let dataset = rate_card();
    let bound = bind_dry_interface_system_card(&dataset, InputAuthority::SyntheticFixture).unwrap();
    let attempt = |curve, point| {
        SourceThermalFriction::resolve(
            bound.clone(),
            point,
            curve,
            MaterialPropertySelection::SingleClaimOnly,
            axes(),
        )
    };
    let wrong_unit = ThermalFrictionCurve::Stribeck {
        static_mu: requirement(STATIC),
        kinetic_mu: requirement(KINETIC),
        characteristic_speed: requirement("speed"),
        viscous_per_speed: typed("viscous", SLOPE),
    };
    assert!(attempt(wrong_unit, query(280.0)).is_err());
    assert!(
        attempt(
            ThermalFrictionCurve::Coulomb {
                static_mu: requirement(STATIC),
                kinetic_mu: typed(STATIC, Dims::NONE),
            },
            query(280.0)
        )
        .is_err()
    );
    // The original curve owns local speed. A single source speed cannot stand
    // in for all quadrature points across this translating and spinning patch.
    let speed_query = query(280.0)
        .with_quantity("sliding-speed", QuantitySpec::dimensional(SPEED), 0.1)
        .unwrap();
    assert!(attempt(curves()[2].clone(), speed_query).is_err());
    let bound_source = source(&dataset, curves()[2].clone());
    let saved = bound_source.clone();
    assert!(bound_source.hertz_at(321.0, None, 10.0, 0.02).is_err());
    assert!(
        bound_source
            .hertz_at(280.0, Some(100.0), 10.0, 0.02)
            .is_err()
    );
    assert!(bound_source.hertz_at(280.0, None, -1.0, 0.02).is_err());
    assert!(bound_source.hertz_at(280.0, None, 10.0, 0.0).is_err());
    assert_eq!(bound_source, saved);
    let invalid = card(vec![curve(STATIC, [0.9, 0.1]), curve(KINETIC, [0.4, 0.2])]);
    let source = source(&invalid, curves()[0].clone());
    assert!(source.coefficients_at(320.0, None).is_err());
    assert_eq!(
        source.coefficients_at(280.0, None).unwrap().1,
        FrictionLaw::Coulomb {
            static_mu: 0.9,
            kinetic_mu: 0.4
        }
    );
}

#[test]
fn family_roles_authority_and_ordered_interface_remain_bound_when_values_coincide() {
    let card = card(vec![curve(STATIC, [0.5, 0.5]), curve(KINETIC, [0.5, 0.5])]);
    let forward = source(&card, curves()[0].clone());
    let swapped = source(
        &card,
        ThermalFrictionCurve::Coulomb {
            static_mu: requirement(KINETIC),
            kinetic_mu: requirement(STATIC),
        },
    );
    assert_ne!(forward.identity(), swapped.identity());
    let shared = source(
        &card,
        ThermalFrictionCurve::Coulomb {
            static_mu: requirement(STATIC),
            kinetic_mu: requirement(STATIC),
        },
    );
    assert_eq!(shared.initial_source().properties().len(), 1);
    assert_eq!(
        forward.coefficients_at(300.0, None).unwrap().1,
        shared.coefficients_at(300.0, None).unwrap().1
    );
    let caller = SourceThermalFriction::resolve(
        bind_dry_interface_system_card(&card, InputAuthority::CallerDeclared).unwrap(),
        query(280.0),
        curves()[0].clone(),
        MaterialPropertySelection::SingleClaimOnly,
        axes(),
    )
    .unwrap();
    assert_ne!(caller.identity(), forward.identity());
    let reversed = InterfaceSystemCard::assemble(
        card.surface_b().clone(),
        card.surface_a().clone(),
        card.context().clone(),
        card.claims().clone(),
        vec![],
    )
    .unwrap();
    assert_ne!(
        source(&reversed, curves()[0].clone()).identity(),
        forward.identity()
    );
}
