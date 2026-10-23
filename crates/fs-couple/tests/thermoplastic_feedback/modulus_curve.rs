//! Affine sourced moduli change native stress, entropy and actual thermal evolution.
use super::*;
use fs_material::SmallStrainLaw;
use fs_material::plastic::J2State;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::thermal::ThermoplasticSource;
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
};

const E0: f64 = 1.2e9;
const H0: f64 = 1e8;
const NU: f64 = 0.2;
const RHO: f64 = 1000.0;

fn affine(cp: f64, de: f64, dh: f64, alpha: Option<f64>) -> Fixture {
    let mut data = properties(cp, false);
    for (index, base, slope) in [(0, E0, de), (3, H0, dh)] {
        curve(
            &mut data[index],
            "T",
            [280.0, 296.0, 300.0, 309.0, 320.0]
                .map(|t| (t, base + slope * (t - 300.0)))
                .to_vec(),
        );
    }
    if let Some(alpha) = alpha {
        data.push(scalar(ALPHA, INVERSE_TEMPERATURE_DIMS, alpha));
    }
    Fixture::from_properties(
        data,
        cp,
        QueryPoint::new().with("T", 300.0).unwrap(),
        "T",
        SolidLiquidPhase::Solid,
    )
}
fn source(fixture: &Fixture, expansion: bool) -> ThermoplasticSource<'_> {
    let mut query = fixture.query();
    query.thermal_expansion = expansion.then_some(MaterialPropertySelection::SingleClaimOnly);
    ThermoplasticSource::try_new(query, 0.0).unwrap()
}
fn tolerances() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-4,
        yield_pa: 1e-3,
    }
}
fn inner(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    (0..6)
        .map(|i| {
            if i < 3 {
                a[i] * b[i]
            } else {
                2.0 * a[i] * b[i]
            }
        })
        .sum()
}
// Independent tensor elasticity, rather than the native stress/energy getter.
fn elastic_stress(t: f64, e: &[f64; 6], de: f64) -> [f64; 6] {
    let youngs = E0 + de * (t - 300.0);
    let mu = youngs / (2.0 * (1.0 + NU));
    let lambda = youngs * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
    let trace = e[0] + e[1] + e[2];
    core::array::from_fn(|i| 2.0 * mu * e[i] + if i < 3 { lambda * trace } else { 0.0 })
}
fn elastic_energy(t: f64, e: &[f64; 6], de: f64) -> f64 {
    0.5 * inner(e, &elastic_stress(t, e, de))
}

#[test]
fn g1_affine_moduli_close_actual_free_energy_entropy_and_heating_cooling_work() {
    let cp = 5000.0;
    let (de, dh, alpha) = (-1e7, -1e6, 1e-5);
    let fixture = affine(cp, de, dh, Some(alpha));
    let source = source(&fixture, true);
    let p = 0.002;
    let state = J2State {
        plastic_strain: [p, -0.5 * p, -0.5 * p, 0.0, 0.0, 0.0],
        alpha: p,
    };
    let e0 = [0.0004, 0.0005, 0.0003, 0.0001, -0.00005, 0.00002];
    let e1 = [0.0007, 0.0006, 0.0005, 0.00015, -0.00002, 0.00008];
    for (t0, t1, old_e, new_e) in [(293.0, 314.0, e0, e1), (314.0, 293.0, e1, e0)] {
        let total = |t: f64, e: &[f64; 6]| {
            core::array::from_fn(|i| {
                e[i] + state.plastic_strain[i] + if i < 3 { alpha * (t - 300.0) } else { 0.0 }
            })
        };
        let old = total(t0, &old_e);
        let new = total(t1, &new_e);
        let increment = source
            .advance_with_energy(
                cp * (t0 - 300.0),
                cp * (t1 - 300.0),
                &old,
                &state,
                &new,
                tolerances(),
            )
            .unwrap();
        assert_eq!(increment.mechanical.state, state);
        let sigma0 = elastic_stress(t0, &old_e, de);
        let sigma1 = elastic_stress(t1, &new_e, de);
        let psi0 = elastic_energy(t0, &old_e, de) + 0.5 * (H0 + dh * (t0 - 300.0)) * p * p;
        let psi1 = elastic_energy(t1, &new_e, de) + 0.5 * (H0 + dh * (t1 - 300.0)) * p * p;
        let entropy = |t: f64, e: &[f64; 6], sigma: &[f64; 6]| {
            alpha * (sigma[0] + sigma[1] + sigma[2])
                - de / (E0 + de * (t - 300.0)) * elastic_energy(t, e, de)
                - 0.5 * dh * p * p
        };
        let storage = t1 * entropy(t1, &new_e, &sigma1) - t0 * entropy(t0, &old_e, &sigma0);
        let coefficient = de * (t1 - t0) / (E0 + de * (t0 - 300.0))
            * elastic_energy(t0, &old_e, de)
            + 0.5 * dh * (t1 - t0) * p * p;
        let heat = (sigma1[0] + sigma1[1] + sigma1[2]) * alpha * (t1 - t0) - coefficient - storage;
        let delta = core::array::from_fn(|i| new[i] - old[i]);
        let delta_e = core::array::from_fn(|i| new_e[i] - old_e[i]);
        close(increment.stored_change_j_m3(), psi1 - psi0, 1e-6);
        close(increment.coefficient_storage_change_j_m3, coefficient, 1e-6);
        close(increment.thermoelastic_storage_change_j_m3, storage, 1e-6);
        close(increment.reversible_heat_j_m3, heat, 1e-6);
        close(increment.mechanical_work_j_m3, inner(&sigma1, &delta), 1e-6);
        close(
            increment.mechanical.energy.numerical_dissipation_j_m3,
            elastic_energy(t1, &delta_e, de),
            1e-6,
        );
        close(
            increment.mechanical.energy.physical_dissipation_j_m3,
            0.0,
            1e-12,
        );
        close(increment.mechanical.energy.residual_j_m3, 0.0, 1e-6);
        close(
            increment.mechanical_work_j_m3
                - increment.stored_change_j_m3()
                - storage
                - increment.caloric_heat_j_m3()
                - increment.mechanical.energy.numerical_dissipation_j_m3,
            0.0,
            1e-6,
        );
    }
}

#[test]
fn g3_affine_surface_moduli_keep_the_typed_axis_source_and_coefficient_work() {
    let cp = 5000.0;
    let mut data = properties(cp, false);
    data[2] = scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6);
    pressure_surface(&mut data[0], [1.6e9, 1.2e9, 1.4e9, 1.0e9]);
    pressure_surface(&mut data[3], [1.3e8, 0.9e8, 1.1e8, 0.7e8]);
    let axis = "specimen_kelvin";
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    for index in [0, 3] {
        if let PropertyValue::Surface2D { axes, .. } = &mut data[index].value {
            axes[0] = axis.into();
        }
        data[index].validity = ValidityDomain::unconstrained()
            .with_quantity(axis, absolute, 280.0, 320.0)
            .with_quantity("pressure", pressure, 1e5, 3e5);
    }
    let context = QueryPoint::new()
        .with_quantity(axis, absolute, 300.0)
        .unwrap()
        .with_quantity("pressure", pressure, 2e5)
        .unwrap();
    let fixture = Fixture::from_properties(data, cp, context, axis, SolidLiquidPhase::Solid);
    let source = source(&fixture, false);
    let (_, hot) = source.resolve_enthalpy(10.0 * cp).unwrap();
    close(hot.law().elastic.youngs, 1.2e9, 1e-6);
    close(hot.law().hardening, 0.9e8, 1e-7);
    assert_eq!(hot.material().axis_quantities().get(axis), Some(&absolute));
    assert!(
        !hot.material()
            .query_point()
            .iter()
            .any(|(name, _)| name == "T")
    );
    for name in [YOUNG_MODULUS_PROPERTY, BASIS.property_name()] {
        assert_eq!(
            hot.material()
                .property(name)
                .unwrap()
                .answer()
                .receipt
                .selected,
            source
                .reference()
                .material()
                .property(name)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    let strain = [0.001, 0.001, 0.001, 0.0, 0.0, 0.0];
    let state = source.reference().initial_state();
    let increment = source
        .advance_with_energy(0.0, 10.0 * cp, &strain, &state, &strain, tolerances())
        .unwrap();
    // Fixed strain under affine E: psi changes by -250 J/m3, Ts by +250.
    // Internal energy is unchanged; no coefficient loss may become heat.
    close(increment.coefficient_storage_change_j_m3, -250.0, 1e-6);
    close(increment.stored_change_j_m3(), -250.0, 1e-6);
    close(increment.thermoelastic_storage_change_j_m3, 250.0, 1e-6);
    close(increment.reversible_heat_j_m3, 0.0, 1e-6);
    close(increment.mechanical_work_j_m3, 0.0, 1e-12);
    close(
        increment.mechanical.energy.numerical_dissipation_j_m3,
        0.0,
        1e-12,
    );
    for stress in &increment.mechanical.stress_pa[..3] {
        close(*stress, 2e6, 1e-5);
    }
}

#[test]
fn g3_modulus_data_changes_actual_adiabatic_temperature_and_warming_stress() {
    context(|gate, cx| {
        let cp = 500.0;
        let strain = [0.001, 0.001, 0.001, 0.0, 0.0, 0.0];
        let mut temperatures = Vec::new();
        for de in [0.0, -1e7] {
            let fixture = affine(cp, de, 0.0, None);
            let original = fixture.point();
            let ambient = environment(300.0, 0.0);
            let loaded = original
                .advance(
                    strain,
                    1.0,
                    config(),
                    |input| ambient.advance(cx, input),
                    || gate.is_requested(),
                )
                .unwrap();
            let entropy = -de * 2.5e-6;
            let expected = 300.0 * RHO * cp / (RHO * cp + entropy);
            close(loaded.next.caloric_state().temperature_k(), expected, 2e-8);
            close(loaded.reversible_heat_j, -VOLUME * expected * entropy, 1e-9);
            close(loaded.coefficient_storage_change_j, 0.0, 1e-12);
            close(loaded.combined_residual_j, 0.0, 1e-8);
            let hot = environment(320.0, 200.0);
            let heated = loaded
                .next
                .advance(
                    strain,
                    1.0,
                    config(),
                    |input| hot.advance(cx, input),
                    || gate.is_requested(),
                )
                .unwrap();
            assert!(heated.next.caloric_state().temperature_k() > expected);
            close(heated.reversible_heat_j, 0.0, 1e-8);
            close(
                heated.coefficient_storage_change_j + heated.thermoelastic_storage_change_j,
                0.0,
                1e-8,
            );
            close(heated.combined_residual_j, 0.0, 1e-8);
            if de < 0.0 {
                assert!(heated.coefficient_storage_change_j < 0.0);
                assert!(heated.next.stress_pa()[0] < loaded.next.stress_pa()[0]);
            }
            assert_eq!(original.time_s(), 0.0);
            assert_eq!(original.history().alpha, 0.0);
            temperatures.push(expected);
        }
        assert!(temperatures[1] < temperatures[0]);
    });
}

#[test]
fn g4_uniaxial_tangent_conversion_must_not_hide_nonaffine_hardening_energy() {
    for (de, et, admitted) in [(-1e7, 2e7, false), (-1e7, 0.0, true), (0.0, 2e7, true)] {
        let mut data = properties(500.0, false);
        curve(
            &mut data[0],
            "T",
            vec![(280.0, E0 - 20.0 * de), (320.0, E0 + 20.0 * de)],
        );
        let basis = HardeningPropertyBasis::UniaxialTangentModulus;
        data[3] = scalar(basis.property_name(), Pressure::DIMS, et);
        let fixture = Fixture::from_properties(
            data,
            500.0,
            QueryPoint::new().with("T", 300.0).unwrap(),
            "T",
            SolidLiquidPhase::Solid,
        );
        let mut query = fixture.query();
        query.hardening_basis = basis;
        let result = ThermoplasticSource::try_new(query, 0.0);
        if admitted {
            let source = result.unwrap();
            assert_eq!(source.reference().hardening_basis(), basis);
        } else {
            assert!(matches!(result,Err(FeedbackError::Invalid(message))
                if message.contains("native plastic-strain modulus")));
        }
    }
}

#[test]
fn g4_modulus_expansion_cross_term_refuses_negative_capacity_without_changing_history() {
    let cp = 100.0;
    let alpha = 1e-4;
    let control = affine(cp, 0.0, 0.0, Some(alpha));
    let mut query = control.query();
    query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
    let accepted = CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE).unwrap();
    let old_state = accepted.caloric_state();
    let old_identity = accepted.resolved().material().identity();
    let bad = affine(cp, -4e7, 0.0, Some(alpha));
    let mut query = bad.query();
    query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
    assert!(matches!(CaloricPoint::try_new(query,0.0,VOLUME,SURFACE),
        Err(FeedbackError::Invalid(message)) if message.contains("effective heat capacity")));
    // A state inside BOTH elastic and total strain bounds has negative
    // fixed-history capacity, despite the old constant-modulus gate passing.
    let t = 310.0;
    let e = 0.03 / 3.0_f64.sqrt();
    let bulk = (E0 - 4e7 * (t - 300.0)) / (3.0 * (1.0 - 2.0 * NU));
    let bulk_slope = -4e7 / (3.0 * (1.0 - 2.0 * NU));
    let capacity =
        RHO * cp + t * (2.0 * alpha * 3.0 * bulk_slope * 3.0 * e - 9.0 * bulk * alpha * alpha);
    assert!(capacity < 0.0);
    assert!(3.0_f64.sqrt() * (e + alpha * (t - 300.0)) < 0.04);
    assert_eq!(accepted.caloric_state(), old_state);
    assert_eq!(accepted.resolved().material().identity(), old_identity);
    assert_eq!(accepted.time_s(), 0.0);
}
