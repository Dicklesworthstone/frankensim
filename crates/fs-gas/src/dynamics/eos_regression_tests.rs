//! Focused source-energy regressions; synthetic polynomials are analytic oracles.
use super::super::{Boundaries, Boundary, Euler1d, StepControl};
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Length, MolarMass, Time};
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};

fn gas(formation: f64) -> FrozenGas {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[3] = 0.4;
    coefficients[7] = formation;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
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
        sources: vec![],
        provenance: Provenance {
            source: "analytic variable-Cp shock".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    };
    let species = Nasa9StandardStateModelV1::new(
        SpeciesId::new("caloric-shock").unwrap(),
        MolarMass::new(R),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(1.0),
            ElementalReferenceIdV1::new("synthetic-frozen").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}

#[test]
fn conservative_construction_refuses_erased_caloric_energy() {
    let gas = gas(0.0);
    // Previously this returned Ok with cached T=1 K, while decoding the
    // returned conservative state yielded Tmin=0.001 K: sensible energy had
    // vanished when added to kinetic energy, with all fields still finite.
    assert!(gas.conserved(1.0, 1e10, 1.0).is_err());
}

#[test]
fn admitted_construction_round_trips_within_declared_caloric_precision() {
    let gas = gas(-1e18);
    for temperature in [0.1, 1.0, 20.0] {
        for velocity in [-30.0, 0.0, 45.0] {
            let state = gas.conserved(0.8, velocity, 0.8 * temperature).unwrap();
            let recovered = gas.primitive(state, || false).unwrap();
            let error =
                gas.sensible(recovered.temperature).unwrap() - gas.sensible(temperature).unwrap();
            let tolerance = gas.domain().absolute_energy_tolerance
                + gas.domain().relative_energy_tolerance * gas.sensible(temperature).unwrap().abs();
            // Construction and inversion each have this declared tolerance.
            assert!(error.abs() <= 2.0 * tolerance);
        }
    }
}

#[test]
fn variable_cp_moving_shock_retains_source_heating_and_reference_invariance() {
    let normal = gas(0.0);
    let shifted = gas(1e18);
    let t1 = 1.0_f64;
    let u1 = 3.0 * ((3.5 + 0.4 * t1) / (2.5 + 0.4 * t1) * t1).sqrt();
    let b = u1 + t1 / u1;
    let mut lower = t1 / u1;
    let mut upper = u1;
    // Independent Rankine-Hugoniot energy root for analytic Cp=3.5+0.4*T.
    // Factoring u-u1 removes the reversible upstream root.
    for _ in 0..100 {
        let u = 0.5 * (lower + upper);
        let t = u * (b - u);
        let secant_cp = 3.5 + 0.2 * (t + t1);
        let residual = secant_cp * (t1 / u1 - u) + 0.5 * (u + u1);
        if residual > 0.0 {
            lower = u;
        } else {
            upper = u;
        }
    }
    let u2 = 0.5 * (lower + upper);
    let t2 = u2 * (b - u2);
    let rho2 = u1 / u2;
    let p2 = rho2 * t2;
    assert!(t2 > t1);
    let ds = normal
        .source()
        .specific_entropy_change(
            Temperature::new(t1),
            Pressure::new(1.0),
            Temperature::new(t2),
            Pressure::new(p2),
        )
        .unwrap()
        .value();
    assert!(ds > 0.0);
    let speed = 0.7;
    let initial_position = 0.35;
    let n = 160;
    let dx = 1.0 / n as f64;
    let make = |gas: FrozenGas| {
        let states = (0..n)
            .map(|i| {
                if (i as f64 + 0.5) * dx < initial_position {
                    gas.conserved(1.0, u1 + speed, 1.0).unwrap()
                } else {
                    gas.conserved(rho2, u2 + speed, p2).unwrap()
                }
            })
            .collect();
        Euler1d::new(
            gas,
            Length::new(dx),
            1.0,
            states,
            Boundaries::Ends {
                left: Boundary::Transmissive,
                right: Boundary::Transmissive,
            },
            || false,
        )
        .unwrap()
    };
    let mut a = make(normal);
    let mut b = make(shifted);
    let end = 0.1;
    let mut steps = 0;
    while a.time().value() < end {
        let dt = Time::new(end - a.time().value());
        let report = a.step(dt, StepControl::default(), || false).unwrap();
        let shifted_report = b.step(dt, StepControl::default(), || false).unwrap();
        assert_eq!(report, shifted_report);
        assert!(report.conservation_residual.energy.abs() < 2e-11);
        steps += 1;
        assert!(steps < 2000);
    }
    assert_eq!(a.states(), b.states());
    let face = a
        .states()
        .windows(2)
        .enumerate()
        .max_by(|(_, l), (_, r)| {
            (l[1].density - l[0].density)
                .abs()
                .total_cmp(&(r[1].density - r[0].density).abs())
        })
        .unwrap()
        .0
        + 1;
    let observed_position = face as f64 * dx;
    let predicted_position = initial_position + speed * end;
    assert!(
        (observed_position - predicted_position).abs() < 3.0 * dx,
        "moving shock: observed {observed_position}, predicted {predicted_position}"
    );
    let downstream = a.gas().primitive(a.states()[n * 3 / 4], || false).unwrap();
    assert!((downstream.temperature - t2).abs() < 1e-7);
}
