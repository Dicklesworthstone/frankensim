//! A source-selected intrinsic kinetic fit releases heat into an actual finite
//! thermal receiver. All coefficients are manufactured, not material evidence.
use super::*;
use fs_material::field_storage::hysteresis::source::MAGNETIC_VISCOSITY_DIMS;

const ETA: f64 = 0.6;

fn kinetic_card(bounds: [f64; 2]) -> (MaterialCard, PropertyKey) {
    let original = card(true);
    let mut claims = original.claims().clone();
    let key = PropertyKey::new("intrinsic_switching_viscosity", MAGNETIC_VISCOSITY_DIMS);
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Scalar {
                value: ETA,
                dims: MAGNETIC_VISCOSITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                bounds[0],
                bounds[1],
            ),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured intrinsic switching heat fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    (
        MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap(),
        key,
    )
}

// Independent collinear backward-Euler minimizer of quadratic storage,
// dry pinning and quadratic kinetic dissipation inside the original bounds.
fn kinetic_play(old: f64, h: f64, dt: f64) -> f64 {
    let a = CHI.recip() + ETA / dt;
    let drive = h + ETA / dt * old;
    old.clamp((drive - HC) / a, (drive + HC) / a)
        .clamp(-1.0, 1.0)
}

#[test]
fn intrinsic_relaxation_deposits_full_physical_heat_once_and_refreshes_all_sources() {
    let (card, key) = kinetic_card([200.0, 500.0]);
    let curve = caloric(&card);
    let volume = 2.0 * DEPTH;
    let mut thermal = thermal_state(&curve, volume);
    let problem = problem_with_viscosity(&card, true, Some(key.clone()));
    let mut field = state(&card, &problem, false);
    let context = point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes(8);
    let network = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let mut z = 0.0;
    let mut expected_temperature = 300.0;
    let mut saw_hold_relaxation = false;
    let mut previous_h = 0.0;
    for (h, dt) in [(6.0, 0.01), (6.0, 0.1), (6.0, 0.2), (0.0, 0.1), (-6.0, 0.2)] {
        let next_z = kinetic_play(z, h, dt);
        let pinning = volume * HC * (next_z - z).abs();
        let viscous = volume * ETA * (next_z - z).powi(2) / dt;
        let heat = pinning + viscous;
        expected_temperature += heat / (DENSITY * CP * volume);
        let mut currents = [6.0 * h; 8];
        currents[4..].fill(-6.0 * h);
        let candidate = network
            .step_hysteretic_eddy_candidate(
                &thermal,
                field.time_s() + dt,
                &field,
                &currents,
                controls(),
                field_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let report = candidate.report();
        near(report.electrical.magnetic.pinning_heat_j, pinning);
        near(report.electrical.magnetic.viscous_heat_j, viscous);
        near(report.electrical.hysteresis_heat_j, heat);
        assert_eq!(report.electrical.joule_heat_j, 0.0);
        near(candidate.deposited_heat_j(), heat);
        near(report.thermal.stored_energy_change_j, heat);
        near(
            candidate.thermal().temperature_k(CELL).unwrap(),
            expected_temperature,
        );
        near(report.total_energy_residual_j, 0.0);
        assert!(report.electrical.numerical_dissipation_j > 0.0);
        near(
            candidate.field().magnetic_state().cell_states()[0].polarizations_t()[0][0],
            next_z,
        );
        for sample in &report.source.magnetic {
            assert_eq!(sample.material.properties().len(), 6);
            near(
                sample.material.property_by_key(&key).unwrap().value_si(),
                ETA,
            );
            near(sample.material.query().axes()["T"], expected_temperature);
        }
        saw_hold_relaxation |= h == previous_h && viscous > 0.0;
        (thermal, field, _) = candidate.into_parts();
        z = next_z;
        previous_h = h;
    }
    assert!(saw_hold_relaxation);
}

#[test]
fn finite_thermal_owner_refuses_a_kinetic_fit_without_full_caloric_domain_support() {
    let (card, key) = kinetic_card([299.0, 301.0]);
    let curve = caloric(&card);
    let thermal = thermal_state(&curve, 2.0 * DEPTH);
    let problem = problem_with_viscosity(&card, true, Some(key));
    let field = state(&card, &problem, false);
    let context = point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes(8);
    let network = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let before = field.magnetic_state().cell_states().to_vec();
    let mut currents = [36.0; 8];
    currents[4..].fill(-36.0);
    assert!(
        network
            .step_hysteretic_eddy_candidate(
                &thermal,
                0.1,
                &field,
                &currents,
                controls(),
                field_controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    assert_eq!(field.magnetic_state().cell_states(), before);
    assert_eq!(field.time_s(), 0.0);
    near(thermal.temperature_k(CELL).unwrap(), 300.0);
}
