//! Complete reacting entropy uses the original calorics and actual chemical sources.

use super::molecular::molecular;
use super::*;

mod interval;
use crate::multicomponent::dissipation::{ConstantTransport, DissipationOperator};
use crate::multicomponent::{Boundary, MixtureEntropyControl, Reconstruction};
use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1, required_nasa9_parameter_dims};

fn controls() -> ReactiveStepControl {
    ReactiveStepControl {
        transport: StepControl {
            reconstruction: Reconstruction::FirstOrder,
            ..StepControl::default()
        },
        ..ReactiveStepControl::default()
    }
}

// These fixtures have M=R and p0=1 Pa. Integrate their original NASA-9
// Cp/R=a+b*T, s0/R=a*ln(T)+b*T+c independently. In particular, c remains
// physical reaction entropy when a species is created or consumed.
fn source_entropy(flow: &MixtureEuler1d, calorics: &[[f64; 3]]) -> f64 {
    (0..flow.cell_count())
        .map(|i| {
            let q = flow.cell(i).unwrap();
            let t = flow.gas().primitive(q, || false).unwrap().temperature;
            flow.mesh().volumes()[i]
                * q.partial_densities
                    .iter()
                    .zip(calorics)
                    .filter_map(|(&rho, &[a, b, c])| {
                        (rho > 0.0).then(|| rho * ((a - 1.0) * t.ln() + b * t + c - rho.ln()))
                    })
                    .sum::<f64>()
        })
        .sum()
}

#[test]
fn g1_closed_reacting_entropy_admits_product_creation_and_replays_heat_release_duration() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut simulation = uniform(&mechanism);
    let initial_entropy = source_entropy(simulation.flow(), &[[3.5, 0.0, 0.0]; 2]);
    let initial = simulation.flow().totals();
    assert_eq!(initial.species_masses[1], 0.0);
    let report = simulation
        .step_entropy_consistent(
            Time::new(0.2),
            controls(),
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    let step = report.coupled;
    assert!(
        step.coupling_retries > 0,
        "heat release must restart the entire step"
    );
    assert_eq!(step.elapsed, step.transport.elapsed);
    assert_eq!(simulation.flow().time(), step.elapsed);
    assert_eq!(
        step.chemistry.cell_solves,
        2 * simulation.flow().cell_count()
    );
    assert_eq!(step.transport.hll_fallback_faces, 0);
    assert!(step.dissipation.is_none() && report.dissipation_entropy.is_none());

    // Reaction-created B is a real source, not a species continuity failure.
    assert!(step.after.species_masses[1] > 0.05);
    for k in 0..2 {
        let residual = step.after.species_masses[k]
            - initial.species_masses[k]
            - step.chemistry.species_change_kg[k];
        assert!(residual.abs() < 2e-12);
    }
    assert!(report.entropy.maximum_relative_species_residual < 1e-12);
    assert!(step.absolute_energy_residual_j.abs() < 3e-12);
    assert!(step.mass_residual_kg.abs() < 3e-13);

    let measured = source_entropy(simulation.flow(), &[[3.5, 0.0, 0.0]; 2]) - initial_entropy;
    assert!(measured > 0.3);
    assert!((report.entropy.source_entropy_change_j_k - measured).abs() < 3e-10);
    assert!((report.entropy.entropy_change_j_k - measured).abs() < 3e-10);
    assert!(report.transport_entropy.entropy_change_j_k.abs() < 3e-12);
    assert!((step.chemistry.entropy_change_j_k - measured).abs() < 3e-9);

    let volume = simulation.flow().mesh().volumes()[0];
    let reactor = advance_adiabatic(
        &mechanism,
        &[(SpeciesId::new("A").unwrap(), Amount::new(volume / R))],
        Volume::new(volume),
        Temperature::new(1.0),
        step.elapsed,
        chemistry(),
        || false,
    )
    .unwrap();
    for i in 0..simulation.flow().cell_count() {
        let q = simulation.flow().cell(i).unwrap();
        let p = simulation.flow().gas().primitive(q, || false).unwrap();
        // U=rho*Cv*T-40*rho_B is conserved at fixed volume and velocity.
        assert!((p.temperature - (1.0 + 16.0 * q.partial_densities[1])).abs() < 3e-9);
        assert!((p.temperature - reactor.temperature().value()).abs() < 2e-4);
        assert!((p.velocity - 0.7).abs() < 2e-12);
    }
}

fn entropy_datum(model: Nasa9StandardStateModelV1, value: f64) -> Nasa9StandardStateModelV1 {
    let mut card = model.regions()[0].card().clone();
    let parameter = required_nasa9_parameter_dims()
        .into_iter()
        .nth(8)
        .unwrap()
        .0
        .to_string();
    card.parameters.get_mut(&parameter).unwrap().value = value;
    Nasa9StandardStateModelV1::new(
        model.species().clone(),
        model.molar_mass(),
        model.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

fn variable_caloric() -> (ReactingEuler1d, DissipationOperator<ConstantTransport>) {
    let models = bank(vec![
        source("A", R, 3.5, 0.2, 0.0),
        entropy_datum(source("B", R, 4.0, 0.1, -4.0), 1.5),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("entropy-source-isomers").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(20.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        models,
        "variable-caloric-reacting-entropy",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let count = 16;
    let states = (0..count)
        .map(|i| {
            let angle = 2.0 * core::f64::consts::PI * (i as f64 + 0.5) / count as f64;
            let density = 1.0 + 0.1 * angle.cos();
            let fraction = 0.65 + 0.1 * angle.sin();
            let temperature = 1.0 + 0.15 * angle.sin();
            gas.conserved(
                &[density * fraction, density * (1.0 - fraction)],
                0.4 + 0.1 * angle.cos(),
                density * temperature,
            )
            .unwrap()
        })
        .collect();
    let flow = MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        0.1,
        states,
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let operator = molecular(&flow, 0.04, 0.01, 0.04);
    (
        ReactingEuler1d::new(flow, mechanism, chemistry(), || false).unwrap(),
        operator,
    )
}

#[test]
fn g1_variable_caloric_reacting_molecular_flow_retains_physical_species_entropy_datums() {
    let (mut simulation, mut operator) = variable_caloric();
    let calorics = [[3.5, 0.2, 0.0], [4.0, 0.1, 1.5]];
    let initial_entropy = source_entropy(simulation.flow(), &calorics);
    let initial = simulation.flow().totals();
    let original = simulation.accepted.states.clone();
    let report = simulation
        .step_entropy_consistent_with_dissipation(
            &mut operator,
            Time::new(0.002),
            controls(),
            DissipationControl::default(),
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    let step = report.coupled;
    let halves = step.dissipation.unwrap();
    let entropy_halves = report.dissipation_entropy.unwrap();
    assert_eq!(step.elapsed.value(), 0.002);
    assert_eq!(simulation.flow().time(), step.elapsed);
    assert_eq!(halves.first.elapsed.value(), 0.5 * step.elapsed.value());
    assert_eq!(halves.second.elapsed.value(), 0.5 * step.elapsed.value());
    assert_eq!(step.transport.hll_fallback_faces, 0);
    assert_ne!(simulation.accepted.states, original);
    assert!(step.chemistry.species_change_kg[1] > 1e-4);
    for k in 0..2 {
        assert!(
            (step.after.species_masses[k]
                - initial.species_masses[k]
                - step.chemistry.species_change_kg[k])
                .abs()
                < 3e-12
        );
    }
    assert!(report.entropy.maximum_relative_species_residual < 1e-12);
    assert!(step.absolute_energy_residual_j.abs() < 3e-11);
    assert!(step.mass_residual_kg.abs() < 3e-12);
    assert!(step.conservation_residual.momentum.abs() < 3e-12);

    // The B species carries a real 1.5*R molar entropy constant. Removing the
    // entire species reference contribution would erase its reaction entropy.
    let measured = source_entropy(simulation.flow(), &calorics) - initial_entropy;
    let datum_change = 1.5 * (step.after.species_masses[1] - initial.species_masses[1]);
    assert!(datum_change > 1e-4 && measured > datum_change);
    assert!((report.entropy.source_entropy_change_j_k - measured).abs() < 3e-10);
    assert!((report.entropy.entropy_change_j_k - measured).abs() < 3e-10);
    assert!(
        (report.entropy.source_entropy_change_j_k
            - report.entropy.reference_entropy_residual_j_k
            - report.entropy.entropy_change_j_k)
            .abs()
            < 3e-12
    );
    for change in [
        report.transport_entropy.entropy_change_j_k,
        entropy_halves.first.entropy_change_j_k,
        entropy_halves.second.entropy_change_j_k,
    ] {
        assert!(change >= -MixtureEntropyControl::default().absolute_tolerance_j_k);
    }
    let parts = step.chemistry.entropy_change_j_k
        + report.transport_entropy.source_entropy_change_j_k
        + entropy_halves.first.source_entropy_change_j_k
        + entropy_halves.second.source_entropy_change_j_k;
    assert!((parts - measured).abs() < 2e-8);
}

#[test]
fn g4_late_reacting_entropy_cancellation_is_atomic_and_original_owners_are_reusable() {
    let (mut reference, mut reference_operator) = variable_caloric();
    let mut polls = 0;
    let expected = reference
        .step_entropy_consistent_with_dissipation(
            &mut reference_operator,
            Time::new(0.002),
            controls(),
            DissipationControl::default(),
            MixtureEntropyControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 100);
    let (mut simulation, mut operator) = variable_caloric();
    let original = simulation.accepted.states.clone();
    let mut visits = 0;
    let error = simulation
        .step_entropy_consistent_with_dissipation(
            &mut operator,
            Time::new(0.002),
            controls(),
            DissipationControl::default(),
            MixtureEntropyControl::default(),
            || {
                visits += 1;
                visits >= polls
            },
        )
        .unwrap_err();
    assert_eq!(error, ReactiveError::Cancelled);
    assert_eq!(simulation.accepted.states, original);
    assert_eq!(simulation.flow().time().value(), 0.0);
    let actual = simulation
        .step_entropy_consistent_with_dissipation(
            &mut operator,
            Time::new(0.002),
            controls(),
            DissipationControl::default(),
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(actual.coupled.elapsed, expected.coupled.elapsed);
    for (a, b) in simulation
        .accepted
        .states
        .iter()
        .zip(&reference.accepted.states)
    {
        assert!((a - b).abs() < 3e-10);
    }
    assert!(
        (actual.entropy.entropy_change_j_k - expected.entropy.entropy_change_j_k).abs() < 3e-10
    );
}

#[test]
fn g4_reacting_entropy_refuses_unowned_reconstruction_and_advective_boundary_entropy() {
    let mechanism = mechanism(4.0, 20.0, 0.0);
    for (boundaries, reconstruction, message) in [
        (
            Boundaries::Periodic,
            Reconstruction::Muscl,
            "entropy-consistent reacting flow requires first-order reconstruction",
        ),
        (
            Boundaries::Ends {
                left: Boundary::Transmissive,
                right: Boundary::Transmissive,
            },
            Reconstruction::FirstOrder,
            "entropy-consistent reacting flow requires periodic or stationary reflecting ends",
        ),
    ] {
        let mut flow = field(&mechanism, 4, 0.4, |_| [0.7, 0.3]);
        flow.set_boundaries(boundaries, || false).unwrap();
        let mut simulation =
            ReactingEuler1d::new(flow, mechanism.clone(), chemistry(), || false).unwrap();
        let original = simulation.accepted.states.clone();
        let mut control = controls();
        control.transport.reconstruction = reconstruction;
        let error = simulation
            .step_entropy_consistent(
                Time::new(0.002),
                control,
                MixtureEntropyControl::default(),
                || false,
            )
            .unwrap_err();
        assert_eq!(
            error,
            ReactiveError::Flow(FlowError::EntropyDomain(message))
        );
        assert_eq!(simulation.accepted.states, original);
        assert_eq!(simulation.flow().time().value(), 0.0);
    }
}

mod exterior;

mod exterior_molecular;
