//! Source-owned finite thermoelastic waves: independent moduli and refusal oracles.
use super::*;
use fs_couple::finite_solid::thermal::conservative::ConservativeThermalField;
use fs_couple::finite_solid::thermal::thermodynamics::{
    CaloricReference, FiniteThermoelasticSourceError,
};
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::thermoelastic::{FiniteThermoelasticState, ThermalReferenceJet};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};

const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn caloric_key() -> PropertyKey {
    PropertyKey::new("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0]))
}

// All quantities share one immutable card. Separate slopes at 300 K let the
// refusal fixture distinguish a true source kink from a redundant source knot.
fn thermodynamic_card(
    stretches: [f64; 3],
    right_slope_scale: f64,
    capacities: [f64; 3],
) -> MaterialCard {
    let base = source_card(6.0, 1.0, false, false);
    let mut claims = ClaimSet::new();
    let validity = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        290.0,
        310.0,
    );
    // The chart's density metadata varies with T; original mechanics still
    // owns rho0=6 at its 300 K datum and actual rho=rho0/det(F).
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("density", Density::DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(290.0, 6.5), (300.0, 6.0), (310.0, 5.5)],
                dims: Density::DIMS,
            },
            validity: validity.clone(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: base.models()[0].provenance.clone(),
        })
        .unwrap();
    for (axis, datum) in [0.2, -0.1, 0.05].into_iter().enumerate() {
        let change = (1.0 + datum) * (stretches[axis] - 1.0);
        claims
            .insert_claim(PropertyClaim {
                key: key(axis),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (290.0, datum - change),
                        (300.0, datum),
                        (310.0, datum + right_slope_scale * change),
                    ],
                    dims: Dims::NONE,
                },
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: base.models()[0].provenance.clone(),
            })
            .unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: caloric_key(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: [290.0, 300.0, 310.0].into_iter().zip(capacities).collect(),
                dims: Dims([2, 0, -2, -1, 0, 0]),
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: base.models()[0].provenance.clone(),
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &caloric_key(),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(290.0),
            upper: &point(310.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        0.01,
    )
    .unwrap()
}

fn h_at_300(source: &ResolvedHeatCapacityCurve) -> f64 {
    source
        .curve()
        .knots()
        .iter()
        .find(|knot| knot.temperature_k == 300.0)
        .unwrap()
        .specific_enthalpy_j_kg
}

#[test]
fn same_source_entropy_constraint_has_independent_bulk_increment_and_original_mass() {
    let card = thermodynamic_card([1.01; 3], 1.0, [1.0; 3]);
    let material = bind_expansion(&card, Q, 300.0);
    let chart = caloric(&card);
    let state = material
        .adiabatic_acoustic_state(
            &I,
            ReferenceDeformation::identity(),
            &chart,
            h_at_300(&chart),
            CaloricReference::StressFreeZeroPressure,
        )
        .unwrap();
    assert!(std::ptr::eq(state.material(), &material));
    assert!(std::ptr::eq(state.caloric_source(), &chart));
    assert_eq!(state.caloric_state().temperature_k(), 300.0);
    let response = state.thermodynamics();
    near(response.free_energy_j_m3(), 0.0, 1e-12);
    near(response.entropy_j_m3_k(), 0.0, 1e-12);
    near(response.internal_energy_j_m3(), 0.0, 1e-12);
    // At the stress-free datum, alpha=.001/K, K=40 Pa and mu=15 Pa.
    // beta=-3*K*alpha*I, C_F=rho0*Cp-9*T*K*alpha^2, and only the
    // longitudinal modulus increases by T*(3*K*alpha)^2/C_F.
    let total_capacity = 6.0 - 9.0 * 300.0 * 40.0 * 0.001_f64.powi(2);
    near(
        response.mechanical_heat_capacity_j_m3_k(),
        total_capacity - 6.0,
        1e-11,
    );
    near(state.native().caloric_heat_capacity_j_m3_k(), 6.0, 1e-12);
    near(
        state.native().fixed_deformation_heat_capacity_j_m3_k(),
        total_capacity,
        1e-11,
    );
    for (index, beta) in response
        .stress_temperature_derivative_pa_k()
        .iter()
        .enumerate()
    {
        near(*beta, if index % 4 == 0 { -0.12 } else { 0.0 }, 1e-11);
    }
    let bulk_increment = 300.0 * 0.12_f64.powi(2) / total_capacity;
    let waves = state.directional_waves([0.3, -0.4, 0.5]).unwrap();
    let held = material
        .acoustic_state(&I)
        .unwrap()
        .directional_waves(waves.direction)
        .unwrap();
    for i in 0..3 {
        for j in 0..3 {
            near(
                waves.acoustic_tensor_pa[i][j] - held.acoustic_tensor_pa[i][j],
                bulk_increment * waves.direction[i] * waves.direction[j],
                1e-10,
            );
        }
    }
    for (mode, modulus) in waves
        .modes
        .into_iter()
        .zip([15.0, 15.0, 60.0 + bulk_increment])
    {
        near(mode.directional_modulus_pa, modulus, 1e-10);
        let speed = (modulus / 6.0).sqrt();
        near(mode.phase_speed_m_s.unwrap(), speed, 1e-11);
        for (value, n) in mode
            .energy_velocity_m_s
            .unwrap()
            .into_iter()
            .zip(waves.direction)
        {
            near(value, speed * n, 1e-10);
        }
    }

    // Rotated anisotropic expansion and a noncommuting manufacturing shear:
    // compose the authored linear stretches independently, then compare the
    // source field's finite-strain entropy/temperature response and actual mass.
    let stretches = [1.01, 0.996, 1.006];
    let card = thermodynamic_card(stretches, 1.0, [10.0; 3]);
    let material = bind_expansion(&card, Q, 305.0);
    let chart = caloric(&card);
    let before =
        ReferenceDeformation::try_new([1.05, 0.12, 0.0, 0.0, 0.94, 0.0, 0.0, 0.0, 1.02]).unwrap();
    let f = [1.12, 0.20, 0.02, 0.0, 0.98, 0.03, 0.0, 0.0, 1.01];
    let alpha = stretches.map(|s| (s - 1.0) / 10.0);
    let thermal: [f64; 9] = core::array::from_fn(|ij| {
        (0..3)
            .map(|k| Q[ij / 3][k] * (1.0 + 5.0 * alpha[k]) * Q[ij % 3][k])
            .sum()
    });
    let derivative: [f64; 9] =
        core::array::from_fn(|ij| (0..3).map(|k| Q[ij / 3][k] * alpha[k] * Q[ij % 3][k]).sum());
    let compose = |left: [f64; 9]| {
        core::array::from_fn(|ij| {
            (0..3)
                .map(|k| left[3 * (ij / 3) + k] * before.deformation_gradient()[3 * k + ij % 3])
                .sum()
        })
    };
    let jet = ThermalReferenceJet::try_new(
        ReferenceDeformation::try_new(compose(thermal)).unwrap(),
        compose(derivative),
        [0.0; 9],
    )
    .unwrap();
    let expected = FiniteThermoelasticState::try_new(material.law().law(), &f, 305.0, jet)
        .unwrap()
        .adiabatic_acoustics(6.0, 60.0)
        .unwrap();
    with_cx(|_, cx| {
        let field = ConservativeThermalField::uniform(material, before, 2, cx).unwrap();
        let actual = field
            .adiabatic_acoustic_state(
                1,
                &f,
                &chart,
                150.0,
                CaloricReference::StressFreeZeroPressure,
            )
            .unwrap();
        assert_eq!(actual.reference_before_thermal(), &before);
        near(
            actual.thermodynamics().internal_energy_j_m3(),
            expected.source().internal_energy_j_m3(),
            1e-9,
        );
        near(
            actual.native().fixed_deformation_heat_capacity_j_m3_k(),
            expected.fixed_deformation_heat_capacity_j_m3_k(),
            1e-9,
        );
        near(actual.caloric_state().bulk_density_kg_m3(), 5.75, 1e-12);
        near(
            actual.native().acoustic().reference_density_kg_m3(),
            6.0,
            0.0,
        );
        near(
            actual.native().acoustic().current_density_kg_m3(),
            6.0 / det(&f),
            1e-12,
        );
        for (a, b) in actual
            .thermodynamics()
            .stress_temperature_derivative_pa_k()
            .iter()
            .zip(expected.source().stress_temperature_derivative_pa_k())
        {
            near(*a, *b, 1e-10);
        }
        assert!(
            field
                .adiabatic_acoustic_state(
                    2,
                    &f,
                    &chart,
                    150.0,
                    CaloricReference::StressFreeZeroPressure,
                )
                .is_err()
        );
    });
}

#[test]
fn source_and_caloric_kinks_endpoints_identity_and_negative_total_capacity_are_refused() {
    let reference = ReferenceDeformation::identity();
    let convention = CaloricReference::StressFreeZeroPressure;
    let card = thermodynamic_card([1.01; 3], 1.0, [1.0; 3]);
    let chart = caloric(&card);
    let material = bind_expansion(&card, Q, 300.0);
    let wrong_chart = caloric(&thermodynamic_card([1.01; 3], 1.0, [2.0; 3]));
    let error = material
        .adiabatic_acoustic_state(
            &I,
            reference,
            &wrong_chart,
            h_at_300(&wrong_chart),
            convention,
        )
        .unwrap_err();
    assert!(error.to_string().contains("same solid material card"));
    let error = material
        .adiabatic_acoustic_state(&I, reference, &chart, 11.0, convention)
        .unwrap_err();
    assert!(error.to_string().contains("actual caloric temperature"));
    let endpoint = bind_expansion(&card, Q, 290.0);
    let error = endpoint
        .adiabatic_acoustic_state(&I, reference, &chart, 0.0, convention)
        .unwrap_err();
    assert!(error.to_string().contains("two-sided expansion-source"));

    let kink_card = thermodynamic_card([1.01; 3], 2.0, [1.0; 3]);
    let kink = bind_expansion(&kink_card, Q, 300.0);
    let kink_chart = caloric(&kink_card);
    let error = kink
        .adiabatic_acoustic_state(
            &I,
            reference,
            &kink_chart,
            h_at_300(&kink_chart),
            convention,
        )
        .unwrap_err();
    assert!(error.to_string().contains("slope-changing expansion knot"));

    let kink_card = thermodynamic_card([1.01; 3], 1.0, [1.0, 1.0, 2.0]);
    let kink = bind_expansion(&kink_card, Q, 300.0);
    let kink_chart = caloric(&kink_card);
    let error = kink
        .adiabatic_acoustic_state(
            &I,
            reference,
            &kink_chart,
            h_at_300(&kink_chart),
            convention,
        )
        .unwrap_err();
    assert!(error.to_string().contains("slope-changing caloric knot"));

    let insufficient = thermodynamic_card([1.01; 3], 1.0, [0.001; 3]);
    let cold_capacity = caloric(&insufficient);
    let material = bind_expansion(&insufficient, Q, 300.0);
    let error = material
        .adiabatic_acoustic_state(
            &I,
            reference,
            &cold_capacity,
            h_at_300(&cold_capacity),
            convention,
        )
        .unwrap_err();
    assert!(matches!(error, FiniteThermoelasticSourceError::Material(_)));
}
