use super::*;
use crate::dynamics::CaloricDomain;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};

fn gas(slope: f64, formation: f64) -> FrozenGas {
    let mut a = [0.0; 9];
    a[2] = 3.5;
    a[3] = slope;
    a[7] = formation;
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 1.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 0.001, 50.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "Manufactured thermally perfect wall reference".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(1.0),
        ElementalReferenceIdV1::new("wall-reference").unwrap(),
    )
    .unwrap();
    let species = Nasa9StandardStateModelV1::new(
        SpeciesId::new("gas").unwrap(),
        MolarMass::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
        convention,
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    let source = FrozenIdealGasMixtureModelV1::new(
        vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    FrozenGas::new(
        source,
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn constant_cp_wall_compression_matches_analytic_shock_from_weak_to_strong() {
    let gas = gas(0.0, 0.0);
    for velocity in [1e-12_f64, 1e-5, 0.2, 2.0, 10.0] {
        let state = gas.conserved(1.0, velocity, 1.0).unwrap();
        let primitive = gas.primitive(state, || false).unwrap();
        let flux = gas.wall_flux(state, primitive, 1.0, &mut || false).unwrap();
        let k = 0.5 * 2.4 * velocity * velocity;
        let exact = 1.0 + 0.5 * k + 0.5 * (k * k + 4.0 * k * (1.0 + 1.0 / 6.0)).sqrt();
        close(flux.momentum, exact, 2e-10);
        assert_eq!(flux.mass, 0.0);
        assert_eq!(flux.energy, 0.0);
        assert!(!flux.used_hll_fallback);
    }
}

#[test]
fn constant_cp_retreat_matches_isentropic_rarefaction_and_never_pulls() {
    let gas = gas(0.0, 0.0);
    for speed in [0.01_f64, 0.5, 2.0, 3.0] {
        let state = gas.conserved(1.0, speed, 0.4).unwrap();
        let p = gas.primitive(state, || false).unwrap();
        let wall = gas.stopped_wall(p, -speed, &mut || false).unwrap();
        let factor = 1.0 - 0.2 * speed / (1.4_f64 * 0.4).sqrt();
        close(wall.pressure, 0.4 * factor.powf(7.0), 2e-10);
        close(wall.temperature, 0.4 * factor * factor, 2e-10);
        assert!(wall.pressure > 0.0 && wall.pressure < p.pressure);
        let f = gas.wall_flux(state, p, -1.0, &mut || false).unwrap();
        assert_eq!(f.mass, 0.0);
        assert_eq!(f.energy, 0.0);
        assert_eq!(f.momentum, wall.pressure);
    }
}

#[test]
fn variable_cp_wall_shock_closes_independent_mass_momentum_energy_and_entropy() {
    let gas = gas(0.4, -1e18);
    let q = 3.0;
    let state = gas.conserved(0.8, q, 1.6).unwrap();
    let before = gas.primitive(state, || false).unwrap();
    let wall = gas.stopped_wall(before, q, &mut || false).unwrap();
    let w = q / (wall.density / before.density - 1.0);
    let j1 = before.density * (q + w);
    let j2 = wall.density * w;
    close(j1, j2, 2e-11);
    close(
        before.pressure + j1 * (q + w),
        wall.pressure + j2 * w,
        2e-10,
    );
    // Cp/R = 3.5+0.4T, so this independently integrates the declared Cv/R.
    let delta_e = 2.5 * (wall.temperature - before.temperature)
        + 0.2 * (wall.temperature.powi(2) - before.temperature.powi(2));
    let hugoniot =
        0.5 * (wall.pressure + before.pressure) * (1.0 / before.density - 1.0 / wall.density);
    close(delta_e, hugoniot, 2e-10);
    let entropy = 3.5 * (wall.temperature / before.temperature).ln()
        + 0.4 * (wall.temperature - before.temperature)
        - (wall.pressure / before.pressure).ln();
    assert!(entropy > 0.0);
    let frozen_gamma = (3.5 + 0.4 * before.temperature) / (2.5 + 0.4 * before.temperature);
    let k = 0.5 * (frozen_gamma + 1.0) * before.density * q * q;
    let b = (frozen_gamma - 1.0) / (frozen_gamma + 1.0) * before.pressure;
    let frozen_pressure =
        before.pressure + 0.5 * k + 0.5 * (k * k + 4.0 * k * (before.pressure + b)).sqrt();
    assert!((wall.pressure / frozen_pressure - 1.0).abs() > 0.005);
}

#[test]
fn variable_cp_rarefaction_matches_independent_simpson_integral_and_entropy() {
    let gas = gas(0.4, 1e18);
    let state = gas.conserved(1.0, -1.0, 4.0).unwrap();
    let before = gas.primitive(state, || false).unwrap();
    let wall = gas.stopped_wall(before, -1.0, &mut || false).unwrap();
    let n = 2_000;
    let h = (before.temperature - wall.temperature) / n as f64;
    let f = |t: f64| ((3.5 + 0.4 * t) * (2.5 + 0.4 * t) / t).sqrt();
    let mut sum = f(wall.temperature) + f(before.temperature);
    for i in 1..n {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(wall.temperature + i as f64 * h);
    }
    close(h * sum / 3.0, 1.0, 2e-10);
    let entropy = 3.5 * (wall.temperature / before.temperature).ln()
        + 0.4 * (wall.temperature - before.temperature)
        - (wall.pressure / before.pressure).ln();
    assert!(entropy.abs() < 1e-10);
}

#[test]
fn wall_is_orientation_and_formation_reference_invariant() {
    let reference = gas(0.4, 0.0);
    for velocity in [-1.0_f64, 0.0, 2.0] {
        let state = reference.conserved(1.0, velocity, 3.0).unwrap();
        let p = reference.primitive(state, || false).unwrap();
        let expected = reference.wall_flux(state, p, 1.0, &mut || false).unwrap();
        let reverse = reference.conserved(1.0, -velocity, 3.0).unwrap();
        let reversed = reference
            .wall_flux(
                reverse,
                reference.primitive(reverse, || false).unwrap(),
                -1.0,
                &mut || false,
            )
            .unwrap();
        assert_eq!(expected, reversed);
        for formation in [-1e18, 1e18] {
            let shifted = gas(0.4, formation);
            let q = shifted.conserved(1.0, velocity, 3.0).unwrap();
            let f = shifted
                .wall_flux(q, shifted.primitive(q, || false).unwrap(), 1.0, &mut || {
                    false
                })
                .unwrap();
            assert_eq!(f, expected);
        }
    }
}

#[test]
fn wall_refuses_unsupported_hot_cold_and_cancelled_paths_without_extrapolation() {
    let gas = gas(0.0, 0.0);
    for (velocity, pressure, normal) in [(30.0, 1.0, 1.0), (4.0, 0.4, -1.0)] {
        let state = gas.conserved(1.0, velocity, pressure).unwrap();
        let primitive = gas.primitive(state, || false).unwrap();
        assert!(matches!(
            gas.wall_flux(state, primitive, normal, &mut || false),
            Err(Error::CaloricPath(_))
        ));
    }
    let state = gas.conserved(1.0, -1.0, 3.0).unwrap();
    let p = gas.primitive(state, || false).unwrap();
    let mut polls = 0;
    assert_eq!(
        gas.wall_flux(state, p, 1.0, &mut || {
            polls += 1;
            polls >= 4
        }),
        Err(Error::Cancelled)
    );
}
