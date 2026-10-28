use super::*;
use crate::transport::{WilkeGasViscosityComponentV1, WilkeGasViscosityMixtureV1};
use crate::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, InitialStatePolicy,
    NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, LawId, LawParameter, Provenance};
use fs_qty::MolarMass;
use fs_qty::semantic::{Composition, CompositionBasis};
use std::collections::BTreeMap;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}

// Manufactured NASA and collision records exercise the physics independently
// of a real-species parameter catalog; these are not measured gas properties.
fn source(name: &str, mass: f64, cp_r: f64, cp_linear: f64) -> Nasa9StandardStateModelV1 {
    let coefficients = [0.0, 0.0, cp_r, cp_linear, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 101325.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 30.0, 10000.0),
        sources: vec![],
        provenance: Provenance {
            source: "manufactured kinetic transport conformance".into(),
            license: "fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(101325.0),
            ElementalReferenceIdV1::new("kinetic-transport-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

fn bank(models: Vec<Nasa9StandardStateModelV1>) -> IdealGasEquilibriumV1 {
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

fn collision(name: &str, diameter_angstrom: f64, epsilon_k: f64) -> LennardJonesSpeciesV1 {
    LennardJonesSpeciesV1::new(
        id(name),
        Length::new(diameter_angstrom * 1e-10),
        Temperature::new(epsilon_k),
    )
    .unwrap()
}

fn domain() -> KineticTransportDomainV1 {
    KineticTransportDomainV1 {
        minimum_temperature_k: 100.0,
        maximum_temperature_k: 1000.0,
        minimum_pressure_pa: 1000.0,
        maximum_pressure_pa: 1e6,
    }
}

fn model(thermo: IdealGasEquilibriumV1, data: Vec<LennardJonesSpeciesV1>) -> KineticGasTransportV1 {
    KineticGasTransportV1::new(thermo, "manufactured-lj", 3, data, domain(), || false).unwrap()
}

fn binary(diameter_scale: f64) -> KineticGasTransportV1 {
    model(
        bank(vec![
            source("A", 0.04, 2.5, 0.0),
            source("B", 0.004, 2.5, 0.0),
        ]),
        vec![
            collision("A", 3.5 * diameter_scale, 100.0),
            collision("B", 2.5 * diameter_scale, 20.0),
        ],
    )
}

fn close(actual: f64, expected: f64, relative_tolerance: f64) {
    assert!(
        (actual - expected).abs() <= relative_tolerance * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_neufeld_table_i_keeps_full_viscosity_correction_and_its_domain() {
    // Independently evaluated Table I reference values, including the small
    // oscillatory viscosity correction (about 0.2% at the high endpoint).
    for (t, expected11, expected22) in [
        (0.3, 2.650176361097787, 2.8455432577357387),
        (1.0, 1.4404663995933595, 1.59314511051151),
        (10.0, 0.741854874843624, 0.8240731640708785),
        (100.0, 0.5167176976723339, 0.5860653933147579),
    ] {
        let (omega11, omega22) = collision_integrals(t).unwrap();
        close(omega11, expected11, 3e-14);
        close(omega22, expected22, 3e-14);
    }
    for t in [0.3_f64.next_down(), 100.0_f64.next_up(), f64::NAN] {
        assert!(matches!(
            collision_integrals(t),
            Err(KineticTransportErrorV1::CollisionDomain { .. })
        ));
    }
}

#[test]
fn g1_si_prefactors_match_independent_engineering_unit_references() {
    let m = binary(1.0);
    let mut w = m.workspace().unwrap();
    let v = m
        .evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &[0.25, 0.75],
            &mut w,
            || false,
        )
        .unwrap();
    // Engineering expressions use M in g/mol, sigma in angstroms and p in atm:
    // mu=2.6693e-6 sqrt(M T)/(sigma^2 Omega22) [Pa s];
    // D=1.8583e-7 T^1.5 sqrt(1/Ma+1/Mb)/(p sigma_ab^2 Omega11) [m^2/s].
    // Their rounded prefactors differ from exact-SI constants by up to 0.026%.
    close(v.pure_viscosities()[0], 2.297292237150358e-5, 3e-4);
    close(v.pure_viscosities()[1], 1.907874283365763e-5, 3e-4);
    close(v.binary_diffusivities()[1], 7.070445273615895e-5, 3e-4);
    assert_eq!(v.binary_diffusivities()[0], 0.0);
    assert_eq!(v.binary_diffusivities()[3], 0.0);
    assert_eq!(v.binary_diffusivities()[1], v.binary_diffusivities()[2]);
    assert_eq!(v.maximum_binary_diffusivity(), v.binary_diffusivities()[1]);
    assert_eq!(v.temperature(), Temperature::new(300.0));
    assert_eq!(v.pressure(), Pressure::new(101325.0));
}

#[test]
fn g1_monatomic_prandtl_schmidt_and_absent_species_limits() {
    let m = model(
        bank(vec![
            source("A", 0.04, 2.5, 0.0),
            source("B", 0.04, 2.5, 0.0),
        ]),
        vec![collision("A", 3.5, 100.0), collision("B", 3.5, 100.0)],
    );
    let mut w = m.workspace().unwrap();
    for x in [[0.25, 0.75], [1.0, 0.0], [1.0 - 1e-14, 1e-14]] {
        let v = m
            .evaluate(
                Temperature::new(300.0),
                Pressure::new(101325.0),
                &x,
                &mut w,
                || false,
            )
            .unwrap();
        let mu = v.dynamic_viscosity().value();
        close(mu, v.pure_viscosities()[0], 3e-15);
        close(
            v.thermal_conductivity().value(),
            v.pure_conductivities()[0],
            3e-15,
        );
        close(
            (2.5 * R / 0.04) * mu / v.thermal_conductivity().value(),
            2.0 / 3.0,
            3e-15,
        );
        let density = 101325.0 * 0.04 / (R * 300.0);
        let (omega11, omega22) = collision_integrals(3.0).unwrap();
        close(
            mu / (density * v.binary_diffusivities()[1]),
            (5.0 / 6.0) * omega11 / omega22,
            3e-14,
        );
    }
    let pure = model(
        bank(vec![source("A", 0.04, 2.5, 0.0)]),
        vec![collision("A", 3.5, 100.0)],
    );
    let mut w = pure.workspace().unwrap();
    let v = pure
        .evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &[1.0],
            &mut w,
            || false,
        )
        .unwrap();
    assert_eq!(v.binary_diffusivities(), &[0.0]);
    assert_eq!(v.maximum_binary_diffusivity(), 0.0);
    assert_eq!(v.dynamic_viscosity().value(), v.pure_viscosities()[0]);
    close(
        v.thermal_conductivity().value(),
        v.pure_conductivities()[0],
        3e-15,
    );
}

#[test]
fn g1_wilke_agrees_with_shared_rule_and_conductivity_uses_declared_mixing() {
    let m = binary(1.0);
    let mut w = m.workspace().unwrap();
    let x = [0.25, 0.75];
    let v = m
        .evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &x,
            &mut w,
            || false,
        )
        .unwrap();
    let components = m
        .thermodynamics()
        .models()
        .iter()
        .zip(v.pure_viscosities())
        .map(|(source, &mu)| {
            WilkeGasViscosityComponentV1::new(
                source.species().clone(),
                source.molar_mass(),
                DynViscosity::new(mu),
            )
            .unwrap()
        })
        .collect();
    let reference = WilkeGasViscosityMixtureV1::new(
        Temperature::new(300.0),
        components,
        Composition::new(CompositionBasis::MoleFraction, x.to_vec()).unwrap(),
    )
    .unwrap()
    .evaluate()
    .unwrap();
    close(
        v.dynamic_viscosity().value(),
        reference.dynamic_viscosity().value(),
        3e-15,
    );
    let lambda = v.pure_conductivities();
    let arithmetic = x[0] * lambda[0] + x[1] * lambda[1];
    let harmonic = 1.0 / (x[0] / lambda[0] + x[1] / lambda[1]);
    close(
        v.thermal_conductivity().value(),
        0.5 * (arithmetic + harmonic),
        3e-15,
    );
    assert!(lambda[0] < v.thermal_conductivity().value());
    assert!(v.thermal_conductivity().value() < lambda[1]);
}

#[test]
fn g1_pressure_collision_size_and_original_cp_drive_transport_changes() {
    let m = binary(1.0);
    let larger = binary(2.0);
    let mut w = m.workspace().unwrap();
    let mut w_large = larger.workspace().unwrap();
    let first = m
        .evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &[0.25, 0.75],
            &mut w,
            || false,
        )
        .unwrap();
    let (mu, lambda, diffusion) = (
        first.dynamic_viscosity().value(),
        first.thermal_conductivity().value(),
        first.binary_diffusivities()[1],
    );
    let dense = m
        .evaluate(
            Temperature::new(300.0),
            Pressure::new(202650.0),
            &[0.25, 0.75],
            &mut w,
            || false,
        )
        .unwrap();
    assert_eq!(dense.dynamic_viscosity().value(), mu);
    assert_eq!(dense.thermal_conductivity().value(), lambda);
    close(dense.binary_diffusivities()[1], 0.5 * diffusion, 3e-15);
    let size = larger
        .evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &[0.25, 0.75],
            &mut w_large,
            || false,
        )
        .unwrap();
    close(size.dynamic_viscosity().value(), 0.25 * mu, 3e-15);
    close(size.thermal_conductivity().value(), 0.25 * lambda, 3e-15);
    close(size.binary_diffusivities()[1], 0.25 * diffusion, 3e-15);
    let hot = m
        .evaluate(
            Temperature::new(600.0),
            Pressure::new(101325.0),
            &[0.25, 0.75],
            &mut w,
            || false,
        )
        .unwrap();
    assert!(hot.dynamic_viscosity().value() > mu);
    assert!(hot.binary_diffusivities()[1] > 2.0 * diffusion);

    let original = bank(vec![source("A", 0.04, 3.0, 0.001)]);
    let varying = model(original.clone(), vec![collision("A", 3.5, 100.0)]);
    assert_eq!(varying.thermodynamics(), &original);
    let mut w = varying.workspace().unwrap();
    for t in [300.0, 600.0, 1000.0] {
        let v = varying
            .evaluate(
                Temperature::new(t),
                Pressure::new(101325.0),
                &[1.0],
                &mut w,
                || false,
            )
            .unwrap();
        close(
            v.thermal_conductivity().value() / v.dynamic_viscosity().value(),
            (3.0 + 0.001 * t + 1.25) * R / 0.04,
            3e-15,
        );
    }
}

#[test]
fn g3_canonical_identity_keeps_exact_source_bank_and_refuses_missing_records() {
    let m = binary(1.0);
    let original = m.thermodynamics().clone();
    let reversed = model(
        original.clone(),
        vec![collision("B", 2.5, 20.0), collision("A", 3.5, 100.0)],
    );
    assert_eq!(m, reversed);
    assert_eq!(m.model_id(), "manufactured-lj");
    assert_eq!(m.model_version(), 3);
    assert_eq!(m.species_count(), 2);
    assert_eq!(m.species_data()[0].species(), &id("A"));
    assert_eq!(m.species_data()[0].diameter(), Length::new(3.5 * 1e-10));
    assert_eq!(
        m.species_data()[0].well_depth_over_k(),
        Temperature::new(100.0)
    );
    for records in [
        vec![collision("A", 3.5, 100.0)],
        vec![collision("A", 3.5, 100.0), collision("A", 2.5, 20.0)],
        vec![collision("A", 3.5, 100.0), collision("C", 2.5, 20.0)],
    ] {
        assert!(matches!(
            KineticGasTransportV1::new(original.clone(), "fixture", 1, records, domain(), || false),
            Err(KineticTransportErrorV1::Invalid(_))
        ));
    }
    let changed = model(
        bank(vec![
            source("A", 0.16, 2.5, 0.0),
            source("B", 0.004, 2.5, 0.0),
        ]),
        m.species_data().to_vec(),
    );
    assert_ne!(changed.thermodynamics(), m.thermodynamics());
    let mut w = m.workspace().unwrap();
    let v = m
        .evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &[0.25, 0.75],
            &mut w,
            || false,
        )
        .unwrap();
    let (mu, lambda, diffusion) = (
        v.pure_viscosities()[0],
        v.pure_conductivities()[0],
        v.binary_diffusivities()[1],
    );
    // A workspace may be reused with a same-sized source bank: there is no
    // temperature-only cache that could return the previous owner's masses.
    let changed_v = changed
        .evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &[0.25, 0.75],
            &mut w,
            || false,
        )
        .unwrap();
    close(changed_v.pure_viscosities()[0], 2.0 * mu, 3e-15);
    close(changed_v.pure_conductivities()[0], 0.5 * lambda, 3e-15);
    assert!(changed_v.binary_diffusivities()[1] < diffusion);
}

#[test]
fn g3_physical_domains_original_calorics_and_cancellation_fail_closed() {
    let m = binary(1.0);
    let mut w = m.workspace().unwrap();
    for (t, p) in [
        (99.0, 101325.0),
        (1001.0, 101325.0),
        (300.0, 999.0),
        (300.0, f64::NAN),
    ] {
        assert_eq!(
            m.evaluate(
                Temperature::new(t),
                Pressure::new(p),
                &[0.25, 0.75],
                &mut w,
                || false
            )
            .unwrap_err(),
            KineticTransportErrorV1::OutsideDomain
        );
    }
    for x in [
        [0.0, 0.0],
        [-0.1, 1.1],
        [f64::MAX, f64::MAX],
        [f64::NAN, 1.0],
    ] {
        assert!(matches!(
            m.evaluate(
                Temperature::new(300.0),
                Pressure::new(101325.0),
                &x,
                &mut w,
                || false
            ),
            Err(KineticTransportErrorV1::Invalid(_))
        ));
    }
    let original = m.thermodynamics().clone();
    let mut collision_domain = domain();
    collision_domain.maximum_temperature_k = 2001.0; // B exceeds T*=100.
    assert!(matches!(
        KineticGasTransportV1::new(
            original.clone(),
            "fixture",
            1,
            m.species_data().to_vec(),
            collision_domain,
            || false
        ),
        Err(KineticTransportErrorV1::CollisionDomain { .. })
    ));
    let mut source_domain = domain();
    source_domain.minimum_temperature_k = 29.0;
    assert!(matches!(
        KineticGasTransportV1::new(
            original,
            "fixture",
            1,
            m.species_data().to_vec(),
            source_domain,
            || false
        ),
        Err(KineticTransportErrorV1::CaloricPath(_))
    ));
    let insufficient_translation = model(
        bank(vec![source("A", 0.04, 2.0, 0.0)]),
        vec![collision("A", 3.5, 100.0)],
    );
    let mut one = insufficient_translation.workspace().unwrap();
    assert!(matches!(
        insufficient_translation.evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &[1.0],
            &mut one,
            || false
        ),
        Err(KineticTransportErrorV1::Invalid(_))
    ));
    assert!(matches!(
        m.evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &[0.25, 0.75],
            &mut one,
            || false
        ),
        Err(KineticTransportErrorV1::Invalid(_))
    ));

    let addresses = (
        w.pure_viscosities.as_ptr(),
        w.pure_conductivities.as_ptr(),
        w.binary_diffusivities.as_ptr(),
    );
    let mut polls = 0;
    assert_eq!(
        m.evaluate(
            Temperature::new(300.0),
            Pressure::new(101325.0),
            &[0.25, 0.75],
            &mut w,
            || {
                polls += 1;
                polls == 3
            }
        )
        .unwrap_err(),
        KineticTransportErrorV1::Cancelled
    );
    let reused = m
        .evaluate(
            Temperature::new(600.0),
            Pressure::new(202650.0),
            &[0.75, 0.25],
            &mut w,
            || false,
        )
        .unwrap();
    let mut fresh = m.workspace().unwrap();
    let reference = m
        .evaluate(
            Temperature::new(600.0),
            Pressure::new(202650.0),
            &[0.75, 0.25],
            &mut fresh,
            || false,
        )
        .unwrap();
    assert_eq!(reused.dynamic_viscosity(), reference.dynamic_viscosity());
    assert_eq!(
        reused.thermal_conductivity(),
        reference.thermal_conductivity()
    );
    assert_eq!(
        reused.binary_diffusivities(),
        reference.binary_diffusivities()
    );
    assert_eq!(
        addresses,
        (
            w.pure_viscosities.as_ptr(),
            w.pure_conductivities.as_ptr(),
            w.binary_diffusivities.as_ptr()
        )
    );
}
