//! G0/G1/G3: actual equilibrium populations drive source-defined conductivity.
#[path = "common/ionized_gas.rs"]
mod fixture;
use fixture::*;
use fs_material::conductor::planar::{PlanarConductivity, PlanarOhmicResponse};
use fs_material::ionized_gas::{IonizedGasConductivity, IonizedGasError};

#[path = "ionized_gas/frozen.rs"]
mod frozen;

#[test]
fn g1_saha_populations_and_particle_masses_drive_current_without_neutral_transport() {
    let m = material(8000.0, 100_000.0);
    let t = m.equilibrium().unwrap().temperature().value();
    let p = m.equilibrium().unwrap().pressure().value();
    let mut expected = 0.0;
    for carrier in m.carriers() {
        let f = m
            .equilibrium()
            .unwrap()
            .species()
            .iter()
            .find(|s| &s.species == carrier.species())
            .unwrap()
            .mole_fraction;
        relative(carrier.number_density_per_m3(), f * p / (KB * t), 1e-14);
        let contribution =
            f * p / (KB * t) * E * E / (carrier.mass_kg() * carrier.relaxation().rate_per_s());
        relative(carrier.conductivity_s_m(), contribution, 1e-14);
        expected += contribution;
    }
    relative(m.conductivity_s_m(), expected, 1e-14);
    assert_eq!(
        m.ohmic_conductivity(),
        PlanarConductivity::Scalar(m.conductivity_s_m())
    );
    assert!(m.has_admitted_ohmic_state());
    let field = m.maximum_electric_field_v_m();
    assert!(m.admits_ohmic_field([0.5 * field, 0.5 * field]));
    assert!(!m.admits_ohmic_field([field, field]));
}

#[test]
fn g3_state_and_collision_sources_have_distinct_physical_effects() {
    let owner = chemistry();
    let state = equilibrium(&owner, 8000.0, 100_000.0);
    let a = IonizedGasConductivity::new(
        &owner,
        state.clone(),
        rates(8000.0, 100_000.0, 1.0),
        limits(),
        || false,
    )
    .unwrap();
    let b = IonizedGasConductivity::new(
        &owner,
        state,
        rates(8000.0, 100_000.0, 2.0),
        limits(),
        || false,
    )
    .unwrap();
    relative(a.conductivity_s_m(), 2.0 * b.conductivity_s_m(), 1e-14);
    relative(
        b.maximum_electric_field_v_m(),
        2.0 * a.maximum_electric_field_v_m(),
        1e-14,
    );
    assert!(
        material(8000.0, 100_000.0).conductivity_s_m()
            > 10.0 * material(6000.0, 100_000.0).conductivity_s_m()
    );
}

#[test]
fn g0_collision_sources_must_cover_exact_present_carriers_and_state() {
    let owner = chemistry();
    let state = equilibrium(&owner, 8000.0, 100_000.0);
    assert!(matches!(
        IonizedGasConductivity::new(&owner, state.clone(), vec![], limits(), || false),
        Err(IonizedGasError::MissingRate(_))
    ));
    let duplicate = vec![source("e-", 8000.0, 100_000.0, 1e10); 2];
    assert!(
        IonizedGasConductivity::new(&owner, state.clone(), duplicate, limits(), || false).is_err()
    );
    assert!(
        IonizedGasConductivity::new(
            &owner,
            state.clone(),
            rates(6000.0, 100_000.0, 1.0),
            limits(),
            || false
        )
        .is_err()
    );
    assert!(
        IonizedGasConductivity::new(
            &owner,
            state.clone(),
            rates(8000.0, 1e6, 1.0),
            limits(),
            || false
        )
        .is_err()
    );
    let extra = vec![
        source("e-", 8000.0, 100_000.0, 1e10),
        source("H+", 8000.0, 100_000.0, 1e9),
        source("H", 8000.0, 100_000.0, 1e9),
    ];
    assert!(IonizedGasConductivity::new(&owner, state, extra, limits(), || false).is_err());
}

#[test]
fn g0_transport_admission_is_independent_from_chemical_equilibrium() {
    let owner = chemistry();
    let state = equilibrium(&owner, 8000.0, 100_000.0);
    let accepted = material(8000.0, 100_000.0);
    let mut domain = limits();
    domain.maximum_charged_mole_fraction = accepted.charged_mole_fraction() * 0.5;
    assert!(
        IonizedGasConductivity::new(
            &owner,
            state.clone(),
            rates(8000.0, 100_000.0, 1.0),
            domain,
            || false
        )
        .is_err()
    );
    domain = limits();
    domain.maximum_coulomb_coupling = accepted.coulomb_coupling() * 0.5;
    assert!(
        IonizedGasConductivity::new(
            &owner,
            state.clone(),
            rates(8000.0, 100_000.0, 1.0),
            domain,
            || false
        )
        .is_err()
    );
    domain = limits();
    domain.maximum_thermal_phase_space_density =
        accepted.maximum_thermal_phase_space_density() * 0.5;
    assert!(
        IonizedGasConductivity::new(
            &owner,
            state.clone(),
            rates(8000.0, 100_000.0, 1.0),
            domain,
            || false
        )
        .is_err()
    );
    assert_eq!(
        IonizedGasConductivity::new(
            &owner,
            state,
            rates(8000.0, 100_000.0, 1.0),
            limits(),
            || true
        ),
        Err(IonizedGasError::Cancelled)
    );
}
