//! Calibrated kinetics on genuine source-resolved reduced tensors.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::plane_stress_viscosity::{
    MaterialPlane, PLANE_STRESS_VISCOSITY_PROPERTIES as KEYS, PlaneStressViscosity,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_qty::semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, DynViscosity, QuantitySpec, Temperature};

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn moisture() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass),
        ValueForm::Static,
    ))
}
fn shift() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 1.3,
            c2_k: 60.0,
        },
        300.0,
        [280.0, 330.0],
    )
    .unwrap()
}

// variation: 0 = constants, 1 = already temperature-dependent coefficients,
// 2 = independent moisture curves held at their original fixed moisture.
fn source(
    values: [f64; 6],
    axis: &str,
    typed: bool,
    variation: usize,
    second_temperature: bool,
    band: (f64, f64),
) -> PlaneStressViscosity {
    let mut validity = ValidityDomain::unconstrained().with("omega", band.0, band.1);
    let mut point = QueryPoint::new().with("omega", band.0).unwrap();
    if typed {
        validity = validity.with_quantity(axis, absolute(), 280.0, 320.0);
        point = point.with_quantity(axis, absolute(), 300.0).unwrap();
    } else {
        validity = validity.with(axis, 280.0, 320.0);
        point = point.with(axis, 300.0).unwrap();
    }
    if variation == 2 {
        validity = validity.with_quantity("moisture", moisture(), 0.0, 0.2);
        point = point.with_quantity("moisture", moisture(), 0.1).unwrap();
    }
    if second_temperature {
        validity = validity.with_quantity("other_temperature", absolute(), 280.0, 320.0);
        point = point
            .with_quantity("other_temperature", absolute(), 300.0)
            .unwrap();
    }
    // Extra query coordinates carry no source dependence. They still cannot
    // impersonate the declared absolute-temperature owner of a shift.
    point = point
        .with_quantity(
            "temperature_interval",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap();
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, value) in KEYS.into_iter().zip(values) {
        let property = match variation {
            1 => PropertyValue::Curve {
                abscissa: axis.into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, 0.5 * value), (320.0, 1.5 * value)],
                dims: DynViscosity::DIMS,
            },
            2 => PropertyValue::Curve {
                abscissa: "moisture".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, value), (0.2, 2.0 * value)],
                dims: DynViscosity::DIMS,
            },
            _ => PropertyValue::Scalar {
                value,
                dims: DynViscosity::DIMS,
            },
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, DynViscosity::DIMS),
                value: property,
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: if variation == 0 {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic fixed-temperature reduced tensor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(
                key,
                DynViscosity::DIMS,
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "calibrated viscosity snapshot".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let state = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    PlaneStressViscosity::from_material(state, MaterialPlane::Axes12, KEYS).unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 3e-12 * (a.abs() + b.abs()).max(1e-300),
        "{a:e} != {b:e}"
    );
}

#[test]
fn temperature_scales_original_tensor_stress_heat_and_support_without_changing_source() {
    for (axis, typed) in [("T", false), ("T", true), ("film_temperature", true)] {
        let reference = source(
            [4.0, 9.0, 16.0, 1.0, 0.5, -0.75],
            axis,
            typed,
            0,
            false,
            (10.0, 10000.0),
        );
        let before = reference.clone();
        assert_eq!(
            shift()
                .shift_plane_stress_viscosity(&reference, axis, 300.0)
                .unwrap()
                .law(),
            &reference
        );
        for temperature in [280.0, 330.0] {
            let shifted = shift()
                .shift_plane_stress_viscosity(&reference, axis, temperature)
                .unwrap();
            let a = shifted.state().factor();
            assert_eq!(shifted.state().temperature_k(), temperature);
            assert_eq!(shifted.law().source(), reference.source());
            assert_eq!(shifted.law().properties(), reference.properties());
            assert_eq!(shifted.law().plane(), reference.plane());
            assert_eq!(shifted.law().omega_band_rad_s(), (10.0 / a, 10000.0 / a));
            for (x, v) in shifted
                .law()
                .matrix_pa_s()
                .iter()
                .zip(reference.matrix_pa_s())
            {
                assert_eq!(*x, a * v);
            }
            let (stress, power) = reference.stress_and_power([0.3, -0.7, 1.1]).unwrap();
            let (next, heat) = shifted.law().stress_and_power([0.3, -0.7, 1.1]).unwrap();
            for i in 0..3 {
                near(next[i], a * stress[i]);
            }
            near(heat, a * power);
            shifted.law().validate_frequency(10000.0 / a).unwrap();
            assert!(
                shifted
                    .law()
                    .validate_frequency((10000.0 / a).next_up())
                    .is_err()
            );
            assert!(shifted.law().matrix_pa_s()[5] < 0.0);
        }
        assert_eq!(reference, before);
    }
}

#[test]
fn cached_weak_channels_and_held_moisture_survive_the_same_positive_shift() {
    let a = 3.7;
    let b = a + a * 2.0_f64.powi(-48);
    let weak = source([a, b, 0.0, a, 0.0, 0.0], "T", false, 0, false, (0.0, 1e6));
    for temperature in [280.0, 330.0] {
        let shifted = shift()
            .shift_plane_stress_viscosity(&weak, "T", temperature)
            .unwrap();
        let expected = shifted.state().factor() * (b - a);
        let actual = shifted.law().stress_and_power([-1.0, 1.0, 0.0]).unwrap().1;
        near(actual, expected);
        assert_eq!(
            shifted.law().stress_and_power([0.0, 0.0, 1.0]).unwrap().1,
            0.0
        );
    }
    let wet = source(
        [4.0, 9.0, 16.0, 1.0, 0.5, -0.75],
        "film_temperature",
        true,
        2,
        false,
        (10.0, 10000.0),
    );
    let shifted = shift()
        .shift_plane_stress_viscosity(&wet, "film_temperature", 315.0)
        .unwrap();
    assert_eq!(
        shifted.law().source().state_coordinate("moisture"),
        Some(0.1)
    );
    assert_eq!(shifted.law().source(), wet.source());
    near(
        shifted.law().stress_and_power([0.3, -0.7, 1.1]).unwrap().1,
        shifted.state().factor() * wet.stress_and_power([0.3, -0.7, 1.1]).unwrap().1,
    );
}

#[test]
fn wrong_reference_axes_existing_temperature_dependence_and_second_owners_refuse() {
    let values = [4.0, 9.0, 16.0, 1.0, 0.5, -0.75];
    let reference = source(values, "T", false, 0, false, (10.0, 10000.0));
    for axis in ["missing_temperature", "temperature_interval", "omega"] {
        assert!(
            shift()
                .shift_plane_stress_viscosity(&reference, axis, 315.0)
                .is_err()
        );
    }
    let wrong_reference = TemperatureShift::new(
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 10000.0,
        },
        301.0,
        [280.0, 330.0],
    )
    .unwrap();
    assert!(
        wrong_reference
            .shift_plane_stress_viscosity(&reference, "T", 315.0)
            .is_err()
    );
    for (variation, second) in [(1, false), (0, true)] {
        let original = source(values, "T", true, variation, second, (10.0, 10000.0));
        assert!(
            shift()
                .shift_plane_stress_viscosity(&original, "T", 315.0)
                .is_err()
        );
    }
    assert!(
        shift()
            .shift_plane_stress_viscosity(&reference, "T", 331.0)
            .is_err()
    );
}

#[test]
fn nonzero_overflow_underflow_and_unrepresentable_frequency_support_refuse() {
    for (value, temperature) in [(8e307, 280.0), (f64::from_bits(1), 330.0)] {
        let reference = source(
            [value, 0.0, 0.0, 0.0, 0.0, 0.0],
            "T",
            false,
            0,
            false,
            (0.0, 1e6),
        );
        let before = reference.clone();
        assert!(
            shift()
                .shift_plane_stress_viscosity(&reference, "T", temperature)
                .is_err()
        );
        assert_eq!(reference, before);
    }
    for (band, temperature) in [((f64::from_bits(1), 1.0), 280.0), ((1e308, 1.1e308), 330.0)] {
        let reference = source([0.0; 6], "T", false, 0, false, band);
        assert!(
            shift()
                .shift_plane_stress_viscosity(&reference, "T", temperature)
                .is_err()
        );
    }
    let zero = source([0.0; 6], "T", false, 0, false, (0.0, f64::MAX));
    let shifted = shift()
        .shift_plane_stress_viscosity(&zero, "T", 330.0)
        .unwrap();
    assert_eq!(shifted.law().matrix_pa_s(), &[0.0; 9]);
    assert_eq!(shifted.law().omega_band_rad_s(), (0.0, f64::MAX));
}
