//! G0/G1/G3: original fracture-plane viscosity reaches the native timed law.
use super::*;
use fs_material::brittle::initiation::TensileInitiationSettings;
use fs_material::brittle::source::rate::*;

const NORMAL_ETA: f64 = 1e8;
const TANGENTIAL_ETA: f64 = 2e8;

fn viscous_properties() -> Vec<PropertyClaim> {
    let mut claims = properties();
    claims.push(scalar(
        BRITTLE_NORMAL_VISCOSITY_PROPERTY,
        BRITTLE_VISCOSITY_DIMS,
        NORMAL_ETA,
    ));
    claims.push(scalar(
        BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY,
        BRITTLE_VISCOSITY_DIMS,
        TANGENTIAL_ETA,
    ));
    claims
}

fn resolved(
    source: &MaterialCard,
    width: f64,
) -> Result<ResolvedBrittleCrackBand, BrittleSourceError> {
    resolve_viscous_brittle_crack_band(
        source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        settings(width),
        1.0,
    )
}

#[test]
fn g1_source_viscosity_adds_physical_stress_and_separate_loss_to_original_fracture() {
    let mut claims = viscous_properties();
    claims.push(scalar(DENSITY_PROPERTY, Density::DIMS, 2300.0));
    let source = card(claims);
    let prepared = resolve_viscous_brittle_dynamic_material(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
        1.0,
    )
    .unwrap();
    assert_eq!(prepared.material().card_identity(), source.content_hash());
    assert_eq!(prepared.material().properties().len(), 10);
    near(prepared.reference_density_kg_m3().unwrap(), 2300.0);
    let viscosity = prepared.law().viscosity().unwrap();
    assert_eq!(viscosity.normal_pa_s_m(), NORMAL_ETA);
    assert_eq!(viscosity.tangential_pa_s_m(), TANGENTIAL_ETA);

    let strain = [0.0002, 0.0, 0.0, 0.0, 0.0, 0.0];
    let dt = 0.02;
    let virgin = prepared.law().initial_state();
    let step = prepared
        .law()
        .trial_over_time(&strain, &virgin, dt)
        .unwrap();
    let elastic_modulus = 1e9 * 0.8 / (1.2 * 0.6);
    let jump = 0.01 * strain[0];
    near(
        step.stress_pa[0],
        elastic_modulus * strain[0] + NORMAL_ETA * jump / dt,
    );
    near(
        step.viscous_dissipated_j_m3,
        NORMAL_ETA * jump * jump / (dt * 0.01),
    );
    assert_eq!(step.equilibrium.energy.dissipated, 0.0);
    near(
        step.total_work_j_m3(),
        step.equilibrium.energy.stored_change + step.viscous_dissipated_j_m3,
    );
    assert_eq!(virgin.strain(), &[0.0; 6]);
    assert_eq!(virgin.damage(), 0.0);
    assert!(prepared.law().trial(&strain, &virgin).is_err());
}

#[test]
fn g3_explicit_viscosities_are_not_rescaled_with_the_geometric_band() {
    let source = card(viscous_properties());
    let jump = 1e-6;
    let duration = 0.02;
    for width in [0.01, 0.02, 0.04] {
        let prepared = resolved(&source, width).unwrap();
        let strain = [jump / width, 0.0, 0.0, 0.0, 0.0, 0.0];
        let step = prepared
            .law()
            .trial_over_time(&strain, &prepared.law().initial_state(), duration)
            .unwrap();
        assert_eq!(step.equilibrium.state.damage(), 0.0);
        assert_eq!(
            prepared.law().viscosity().unwrap().normal_pa_s_m(),
            NORMAL_ETA
        );
        // Same undamaged physical jump path: the per-area dashpot law is fixed.
        // No assertion of arbitrary-mesh objectivity after damage is made.
        near(
            width * step.viscous_dissipated_j_m3,
            NORMAL_ETA * jump * jump / duration,
        );
        assert_eq!(
            prepared.law().parameters().mode_i_fracture_energy_j_m2,
            100.0
        );
    }
}

#[test]
fn g0_missing_wrong_dimensional_or_nonphysical_viscosity_refuses() {
    for index in [7, 8] {
        let mut claims = viscous_properties();
        claims.remove(index);
        assert!(resolved(&card(claims), 0.01).is_err());
        for dims in [Dims([-1, 1, -1, 0, 0, 0]), Pressure::DIMS, Dims::NONE] {
            let mut claims = viscous_properties();
            let name = if index == 7 {
                BRITTLE_NORMAL_VISCOSITY_PROPERTY
            } else {
                BRITTLE_TANGENTIAL_VISCOSITY_PROPERTY
            };
            claims[index] = scalar(name, dims, NORMAL_ETA);
            assert!(resolved(&card(claims), 0.01).is_err());
        }
    }
    for value in [0.0, -1.0] {
        let mut claims = viscous_properties();
        claims[7] = scalar(
            BRITTLE_NORMAL_VISCOSITY_PROPERTY,
            BRITTLE_VISCOSITY_DIMS,
            value,
        );
        assert!(resolved(&card(claims), 0.01).is_err());
    }
    let source = card(viscous_properties());
    for maximum in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            resolve_viscous_brittle_crack_band(
                &source,
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                settings(0.01),
                maximum,
            )
            .is_err()
        );
    }
}

#[test]
fn g0_rate_conditioned_measurements_cannot_be_frozen_into_constant_dashpots() {
    for axis in ["strain_rate", "frequency"] {
        let mut claims = viscous_properties();
        let quantity = QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0]));
        claims[7].validity =
            ValidityDomain::unconstrained().with_quantity(axis, quantity, 0.0, 2.0);
        let point = QueryPoint::new()
            .with_quantity(axis, quantity, 1.0)
            .unwrap();
        assert!(
            resolve_viscous_brittle_crack_band(
                &card(claims),
                &point,
                MaterialPropertySelection::SingleClaimOnly,
                settings(0.01),
                1.0,
            )
            .is_err()
        );
    }
}

#[test]
fn g3_temperature_refresh_retains_rate_model_claims_and_damaged_history() {
    let mut claims = viscous_properties();
    for (index, eta) in [(7, NORMAL_ETA), (8, TANGENTIAL_ETA)] {
        claims[index].value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, eta), (400.0, 2.0 * eta)],
            dims: BRITTLE_VISCOSITY_DIMS,
        };
        claims[index].validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
            400.0,
        );
        claims[index].interpolation = InterpolationPolicy::LinearInside;
    }
    let source = card(claims);
    let query = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    let cold = resolve_viscous_brittle_crack_band(
        &source,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        settings(0.01),
        1.0,
    )
    .unwrap();
    let warm = cold.at_temperature(&source, "T", 400.0).unwrap();
    assert_eq!(
        warm.law().viscosity().unwrap().normal_pa_s_m(),
        2.0 * NORMAL_ETA
    );
    assert_eq!(
        warm.law().viscosity().unwrap().tangential_pa_s_m(),
        2.0 * TANGENTIAL_ETA
    );
    assert_eq!(warm.law().viscosity().unwrap().maximum_rate_m_s(), 1.0);
    for (old, new) in cold
        .material()
        .properties()
        .iter()
        .zip(warm.material().properties())
    {
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    let strain = [0.003, 0.0, 0.0, 0.0, 0.0, 0.0];
    let loaded = cold
        .law()
        .trial_over_time(&strain, &cold.law().initial_state(), 0.1)
        .unwrap();
    let state = loaded.equilibrium.state;
    assert!(state.damage() > 0.0 && !state.failed());
    let transferred = cold.transfer_temperature(&state, &warm, "T", 1e-6).unwrap();
    assert_eq!(transferred.state.strain(), state.strain());
    assert_eq!(transferred.state.damage(), state.damage());
    assert_eq!(transferred.state.dissipated_j_m3(), state.dissipated_j_m3());
    assert_eq!(transferred.parameter_work_j_m3, 0.0);
    assert_eq!(transferred.damage_dissipated_j_m3, 0.0);
    assert!(cold.at_temperature(&source, "T", 401.0).is_err());

    let dropped =
        ResolvedBrittleCrackBand::try_from_material(warm.material().clone(), settings(0.01))
            .unwrap();
    assert!(cold.require_temperature_successor(&dropped, "T").is_err());
    let different_bound = ResolvedBrittleCrackBand::try_from_viscous_material(
        warm.material().clone(),
        settings(0.01),
        0.5,
    )
    .unwrap();
    assert!(
        cold.require_temperature_successor(&different_bound, "T")
            .is_err()
    );
    assert!(
        cold.tensile_initiator(TensileInitiationSettings {
            maximum_strain: 0.05,
            relative_stress_tolerance: 1e-8,
            fraction_tolerance: 1e-8,
            maximum_bisections: 64,
        })
        .is_err()
    );
}
