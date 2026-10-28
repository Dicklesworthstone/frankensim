use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Length, MolarMass, Pressure, Temperature, Time};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};

pub(super) fn source(
    name: &str,
    mass: f64,
    cp: f64,
    cp_linear: f64,
    formation: f64,
) -> Nasa9StandardStateModelV1 {
    let a = [0.0, 0.0, cp, cp_linear, 0.0, 0.0, 0.0, formation, 0.0];
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
        sources: vec![],
        provenance: Provenance {
            source: "manufactured multicomponent Euler caloric source".into(),
            license: "fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(1.0),
            ElementalReferenceIdV1::new("multicomponent-euler-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
pub(super) fn bank(models: Vec<Nasa9StandardStateModelV1>) -> IdealGasEquilibriumV1 {
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1; ids.len()]],
        )
        .unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap()
}
pub(super) fn gas(models: Vec<Nasa9StandardStateModelV1>) -> MixtureGas {
    MixtureGas::new(
        bank(models),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}
fn standard(formation: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, formation),
        source("B", 2.0 * R, 3.5, 0.0, -2.0 * formation),
    ])
}
fn flow(
    gas: MixtureGas,
    count: usize,
    boundaries: Boundaries,
    initial: impl Fn(f64) -> ([f64; 2], f64, f64),
) -> MixtureEuler1d {
    let states = (0..count)
        .map(|i| {
            let (d, u, p) = initial((i as f64 + 0.5) / count as f64);
            gas.conserved(&d, u, p).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        0.7,
        states,
        boundaries,
        || false,
    )
    .unwrap()
}
fn run(flow: &mut MixtureEuler1d, end: f64) -> usize {
    let mut steps = 0;
    while flow.time().value() < end {
        let r = flow
            .step(
                Time::new(end - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        assert!(r.maximum_cfl <= 0.40000000000001);
        for value in r.conservation_residual.species_masses[..flow.gas().species_count()]
            .iter()
            .chain([
                &r.conservation_residual.momentum,
                &r.conservation_residual.energy,
            ])
        {
            assert!(value.abs() < 2e-12, "residual {value}");
        }
        steps += 1;
        assert!(steps < 10000);
    }
    steps
}

#[test]
fn caloric_mixture_uses_species_masses_variable_cp_and_frozen_sound_speed() {
    let gas = gas(vec![
        source("A", R, 3.5, 0.25, 7.0),
        source("B", 2.0 * R, 4.0, 0.1, -12.0),
    ]);
    let densities = [0.7, 0.3];
    let t = 2.4;
    let p = (0.7 + 0.3 / 2.0) * t;
    let q = gas.conserved(&densities, 0.8, p).unwrap();
    let primitive = gas.primitive(q.as_cell(), || false).unwrap();
    let cv = 0.7 * (2.5 + 0.25 * t) + 0.15 * (3.0 + 0.1 * t);
    let gamma = 1.0 + 0.85 / cv;
    assert!((primitive.temperature - t).abs() < 2e-12);
    assert!((primitive.sound_speed.powi(2) - gamma * p).abs() < 2e-12);
    let thermal = 0.7 * (2.5 * (t - 0.001) + 0.125 * (t * t - 1e-6))
        + 0.15 * (3.0 * (t - 0.001) + 0.05 * (t * t - 1e-6));
    assert!((q.energy - (thermal + 0.5 * 0.8 * 0.8)).abs() < 2e-12);
    assert_eq!(gas.source().models()[0].species().as_str(), "A");
    assert!(gas.conserved(&[-0.1, 1.1], 0.0, 1.0).is_err());
    assert!(gas.conserved(&[0.0, 0.0], 0.0, 1.0).is_err());
    assert!(gas.conserved(&[1.0], 0.0, 1.0).is_err());
    assert!(gas.conserved(&[1.0, 0.0], 1e16, 1.0).is_err());
    let zero = gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    assert_eq!(zero.partial_densities[1], 0.0);
    let mut flux = [0.0; 4];
    gas.hllc_flux(q.as_cell(), q.as_cell(), &mut flux, || false)
        .unwrap();
    assert!((flux[0] - 0.7 * 0.8).abs() < 2e-13);
    assert!((flux[1] - 0.3 * 0.8).abs() < 2e-13);
    assert!((flux[2] - (0.8 * 0.8 + p)).abs() < 2e-13);
    assert!((flux[3] - (q.energy + p) * 0.8).abs() < 2e-13);
}

#[test]
fn stationary_material_contact_and_tapered_rest_state_are_preserved() {
    // A spatially balanced field remains balanced to the caloric inversion
    // tolerance. Use a tighter declared root budget for this roundoff test.
    let mut domain = CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0));
    domain.absolute_energy_tolerance = 1e-13;
    domain.relative_energy_tolerance = 1e-14;
    let gas = MixtureGas::new(
        bank(vec![
            source("A", R, 3.5, 0.2, 0.0),
            source("B", 2.0 * R, 5.0, 0.4, -5.0),
        ]),
        domain,
        || false,
    )
    .unwrap();
    let n = 32;
    let positions: Vec<_> = (0..=n).map(|i| (i as f64 / n as f64).powi(2)).collect();
    let areas: Vec<_> = positions.iter().map(|x| 1.0 + 0.5 * x).collect();
    let volumes: Vec<_> = positions
        .windows(2)
        .map(|x| (x[1] - x[0]) * (1.0 + 0.25 * (x[1] + x[0])))
        .collect();
    let mesh = Mesh1d::new(positions, areas, volumes).unwrap();
    let states = (0..n)
        .map(|i| {
            gas.conserved(if i < n / 2 { &[1.0, 0.0] } else { &[0.0, 2.0] }, 0.0, 1.0)
                .unwrap()
        })
        .collect();
    let mut f = MixtureEuler1d::new_on_mesh(
        gas,
        mesh,
        states,
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    let original = f.states.clone();
    run(&mut f, 0.02);
    for (a, b) in f.states.iter().zip(&original) {
        assert!((a - b).abs() < 2e-11, "rest drift {}", a - b);
    }
    for i in 0..n {
        let p = f.gas().primitive(f.cell(i).unwrap(), || false).unwrap();
        assert!(p.velocity.abs() < 2e-11);
        assert!((p.pressure - 1.0).abs() < 2e-11);
    }
}

#[test]
fn smooth_species_contact_converges_and_conserves_every_species() {
    let mut errors = Vec::new();
    for n in [40, 80, 160] {
        let mut f = flow(standard(0.0), n, Boundaries::Periodic, |x| {
            let a = 0.5 + 0.2 * (2.0 * std::f64::consts::PI * x).sin();
            ([a, 1.0 - a], 0.6, 1.0)
        });
        let initial = f.totals();
        run(&mut f, 0.2);
        let end = f.totals();
        for k in 0..2 {
            assert!((end.species_masses[k] - initial.species_masses[k]).abs() < 2e-13);
        }
        assert!((end.energy - initial.energy).abs() < 2e-13);
        let mut error = 0.0;
        for i in 0..n {
            let x = (i as f64 + 0.5) / n as f64;
            let expected = 0.5 + 0.2 * (2.0 * std::f64::consts::PI * (x - 0.12)).sin();
            let q = f.cell(i).unwrap();
            error += (q.partial_densities[0] - expected).abs() / n as f64;
            assert!((q.density() - 1.0).abs() < 5e-12);
            let p = f.gas().primitive(q, || false).unwrap();
            assert!((p.pressure - 1.0).abs() < 5e-11);
            assert!((p.velocity - 0.6).abs() < 5e-11);
        }
        eprintln!("species contact n={n}: L1={error:.9e}");
        errors.push(error);
    }
    assert!(
        errors[0] / errors[1] > 3.0 && errors[1] / errors[2] > 3.0,
        "{errors:?}"
    );
}

#[test]
fn shock_matches_existing_frozen_euler_for_uniform_composition() {
    use crate::dynamics::{Euler1d, FrozenGas};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{Composition, CompositionBasis};
    let gas = standard(13.0);
    // Y=(0.6,0.4), M=(R,2R) gives mole fractions (0.75,0.25).
    let frozen = FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            gas.source().models().to_vec(),
            Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        )
        .unwrap(),
        gas.domain(),
        || false,
    )
    .unwrap();
    let n = 120;
    let mut mixed = flow(
        gas,
        n,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        |x| {
            let (rho, p) = if x < 0.5 { (1.0, 1.0) } else { (0.125, 0.1) };
            ([0.6 * rho, 0.4 * rho], 0.0, p)
        },
    );
    let states = (0..n)
        .map(|i| {
            let (rho, p) = if i < n / 2 { (1.0, 1.0) } else { (0.125, 0.1) };
            frozen.conserved(rho, 0.0, p).unwrap()
        })
        .collect();
    let mut single = Euler1d::new(
        frozen,
        Length::new(1.0 / n as f64),
        0.7,
        states,
        mixed.boundaries,
        || false,
    )
    .unwrap();
    for _ in 0..200 {
        let a = mixed
            .step(Time::new(0.0005), StepControl::default(), || false)
            .unwrap();
        let b = single
            .step(Time::new(0.0005), StepControl::default(), || false)
            .unwrap();
        assert_eq!(a.elapsed, b.elapsed);
    }
    for i in 0..n {
        let q = mixed.cell(i).unwrap();
        let a = mixed.gas().primitive(q, || false).unwrap();
        let b = single
            .gas()
            .primitive(single.states()[i], || false)
            .unwrap();
        for (a, b) in [
            (a.density, b.density),
            (a.velocity, b.velocity),
            (a.pressure, b.pressure),
            (a.temperature, b.temperature),
        ] {
            assert!((a - b).abs() < 2e-10, "mixture/frozen {a} {b}");
        }
        assert!((q.partial_densities[0] / q.density() - 0.6).abs() < 2e-12);
    }
}

#[test]
fn independent_formation_reference_shifts_do_not_change_transport() {
    let init = |x: f64| {
        if x < 0.5 {
            ([0.8, 0.2], 0.1, 1.0)
        } else {
            ([0.1, 0.2], -0.1, 0.3)
        }
    };
    let boundaries = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    let mut a = flow(standard(0.0), 40, boundaries, init);
    let mut b = flow(standard(1e18), 40, boundaries, init);
    assert_eq!(a.states, b.states);
    run(&mut a, 0.02);
    run(&mut b, 0.02);
    assert_eq!(a.states, b.states);
    assert_eq!(a.time(), b.time());
}

#[test]
fn rarefaction_keeps_structural_species_zeros_and_refusals_are_atomic() {
    let mut f = flow(
        standard(0.0),
        80,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        |x| ([1.0, 0.0], if x < 0.5 { -1.5 } else { 1.5 }, 1.0),
    );
    run(&mut f, 0.025);
    for i in 0..f.cell_count() {
        let q = f.cell(i).unwrap();
        assert_eq!(q.partial_densities[1], 0.0);
        assert!(q.partial_densities[0] > 0.0);
        let p = f.gas().primitive(q, || false).unwrap();
        assert!(p.temperature >= 0.001);
    }
    let states = f.states.clone();
    let time = f.time();
    let mut calls = 0;
    let e = f.step(Time::new(0.01), StepControl::default(), || {
        calls += 1;
        calls > 150
    });
    assert_eq!(e, Err(Error::Cancelled));
    assert_eq!(f.states, states);
    assert_eq!(f.time(), time);
    assert!(
        f.step(
            Time::new(0.01),
            StepControl {
                cfl: 0.46,
                ..StepControl::default()
            },
            || false
        )
        .is_err()
    );
    assert_eq!(f.states, states);
    let mut duplicate = f.duplicate(&mut || false).unwrap();
    duplicate.time = 10.0;
    duplicate.copy_physical_from(&f).unwrap();
    assert_eq!(duplicate.time(), f.time());
    assert_eq!(duplicate.states, f.states);
    duplicate.validate_states(&mut || false).unwrap();
}
