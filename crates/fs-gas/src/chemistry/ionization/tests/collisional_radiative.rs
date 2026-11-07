use super::{EV, KB, ME, NA, id, initial};
mod thermal;
use crate::chemistry::collisional_radiative::*;
use fs_matdb::Provenance;
use fs_qty::{Area, Energy, Temperature, Time, Volume};
use fs_thermochem::ReactionId;
use fs_thermochem::electron_collision::{
    ElectronCollisionKindV1, ElectronCrossSectionPointV1, ElectronCrossSectionV1,
    MaxwellianRateControlsV1,
};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::ElectronImpactRateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1, PhotonNumberFluxV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;
use fs_thermochem::kinetics::atomic_ionization::{AtomicIonizationV1, AtomicStatisticalWeightsV1};

const SIGMA: f64 = 1e-22;

pub(super) fn fixture(
    x: f64,
    te: f64,
    impact: f64,
    elemental: f64,
    charge: f64,
) -> (FrozenIonizedGasStateV1, AtomicRadiativeRecombinationV1) {
    let (original, atom) = initial(x, te, 6000.0, elemental, charge);
    let gap = atom.ionization_energy_j_mol() / NA;
    let provenance = Provenance {
        source: "analytic ramp impact and constant bound-free cross-section fixture, not measured H data".into(),
        license: "analytic fixture".into(), artifact: None,
    };
    let collision = ElectronCrossSectionV1::tabulated(
        id("H"),
        ElectronCollisionKindV1::Ionization,
        provenance.clone(),
        vec![
            ElectronCrossSectionPointV1 {
                energy: Energy::new(gap),
                cross_section: Area::new(0.0),
            },
            ElectronCrossSectionPointV1 {
                energy: Energy::new(gap + 1000.0 * EV),
                cross_section: Area::new(impact * 1e-18),
            },
        ],
        Area::new(impact * 1e-18),
    )
    .unwrap();
    let channel = AtomicIonizationV1::new(
        original.thermodynamics().clone(),
        ReactionId::new("analytic-coupled-impact").unwrap(),
        id("H+"),
        ElectronImpactRateV1::new(
            id("e"),
            collision,
            MaxwellianRateControlsV1 {
                relative_tolerance: 1e-11,
                maximum_subdivisions: 4096,
            },
            0.1,
        )
        .unwrap(),
        AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
        1e-12,
    )
    .unwrap();
    let photon = AtomicPhotoionizationV1::new(
        channel,
        provenance.clone(),
        vec![
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap),
                cross_section: Area::new(SIGMA),
            },
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap + 1000.0 * EV),
                cross_section: Area::new(SIGMA),
            },
        ],
    )
    .unwrap();
    let source = AtomicRadiativeRecombinationV1::new(photon, Area::new(SIGMA), provenance).unwrap();
    let state = FrozenIonizedGasStateV1::new(
        original.thermodynamics(),
        original.amounts(),
        Volume::new(0.2),
        original.heavy_temperature(),
        original.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    (state, source)
}
fn controls(tolerance: f64) -> CollisionalRadiativeSolve {
    let mut controls =
        CollisionalRadiativeSolve::new(Temperature::new(1000.0), Temperature::new(50000.0));
    controls.relative_tolerance = tolerance;
    controls
}
fn field(energy: f64, rate: f64) -> [PhotonNumberFluxV1; 1] {
    [PhotonNumberFluxV1 {
        photon_energy: Energy::new(energy),
        number_flux_per_m2_s: rate / SIGMA,
    }]
}

#[test]
fn g1_simultaneous_four_channel_evolution_matches_independent_ode() {
    let (initial, source) = fixture(0.1, 10000.0, 1.0, 0.0, 0.0);
    let photons = field(20.0 * EV, 1000.0);
    // Independent DOP853 integration of x, Ue and all photon/event budgets,
    // rtol=2.5e-13. No production quadrature or time-step code is reused:
    // sigma_impact=C(E-chi) gives k=sqrt(8/pi/me)*C*sqrt(theta)
    // *(chi+2theta)*exp(-chi/theta), and constant sigma_PI gives the
    // complete Milne moments proportional to chi²+2chi theta+2theta²
    // and theta*(chi²+4chi theta+6theta²). The source cutoffs exceed 200 kT.
    for (time, x, te, photo, collision, capture, absorbed, emitted) in [
        (
            2e-5,
            0.10443275846131123,
            13659.196390910463,
            0.01798098860984913,
            0.002768513752706761,
            0.016316743901244695,
            34698.033158534985,
            23563.26188349829,
        ),
        (
            4e-5,
            0.10969305850629041,
            13767.799300512135,
            0.03583668453830872,
            0.008510329833801886,
            0.03465395586582019,
            69154.2881975401,
            50073.05832461097,
        ),
    ] {
        let result = advance_collisional_radiative(
            &initial,
            &source,
            &photons,
            Time::new(time),
            controls(1e-7),
            || false,
        )
        .unwrap();
        let state = result.state();
        let balance = result.balance();
        assert!(
            (state.amounts()[1].value() - x).abs() < 2e-7,
            "{} != {x}",
            state.amounts()[1].value()
        );
        assert!((state.electron_temperature().value() - te).abs() < 0.03);
        for (actual, expected) in [
            (balance.photoionization_extent_mol, photo),
            (balance.collisional_extent_mol, collision),
            (balance.radiative_recombination_extent_mol, capture),
        ] {
            assert!((actual - expected).abs() < 3e-7, "{actual} != {expected}");
        }
        assert!((balance.absorbed_radiation_energy_j - absorbed).abs() < 0.5);
        assert!((balance.emitted_radiation_energy_j - emitted).abs() < 0.5);
        let actual_extent = state.amounts()[1].value() - initial.amounts()[1].value();
        assert!((result.reaction_extent().value() - actual_extent).abs() < 2e-13);
        assert!(
            (balance.chemical_storage_change_j
                - source.source().channel().ionization_energy_j_mol() * actual_extent)
                .abs()
                < 1e-9
        );
        assert_eq!(state.amounts()[1], state.amounts()[2]);
        assert!((state.amounts()[0].value() + state.amounts()[1].value() - 1.0).abs() < 2e-15);
        assert_eq!(state.heavy_temperature(), initial.heavy_temperature());
        assert_eq!(state.thermodynamics(), initial.thermodynamics());
        assert_eq!(result.photon_field(), &photons);
        assert!(balance.maximum_normalized_local_error <= 1.0);
        assert!(balance.energy_residual_j.abs() < 1e-7);
        // This uses the independently published thermodynamic states, not the
        // solver's net-event counter. Both gross radiation channels are large.
        let actual_energy = state.internal_energy().value() - initial.internal_energy().value();
        assert!(
            (actual_energy + balance.emitted_radiation_energy_j
                - balance.absorbed_radiation_energy_j)
                .abs()
                < 3e-7
        );
        assert!(balance.photoionization_extent_mol > 3.0 * actual_extent);
        assert!(balance.radiative_recombination_extent_mol > 3.0 * actual_extent);
    }
}

#[test]
fn g1_neutral_start_and_stiff_radiation_turnover_match_exact_fixed_temperature_limit() {
    // For a constant bound-free source, mean captured energy is 3 kT/2
    // when chi/kT=1+sqrt(7). Giving every photoelectron that same energy
    // holds T exactly while BOTH creation and capture change the population.
    // With zero impact, x'=gamma*(1-x)-beta*x² has an independent Riccati solution.
    let gap = 13.598_434_599_702 * EV;
    let theta = gap / (1.0 + 7.0_f64.sqrt());
    let temperature = theta / KB;
    let photon_energy = gap + 1.5 * theta;
    let gamma = 1000.0;
    let c = 299_792_458.0;
    let alpha = (8.0 / (core::f64::consts::PI * ME)).sqrt() * SIGMA / (ME * c * c)
        * (gap * gap + 2.0 * gap * theta + 2.0 * theta * theta)
        / theta.sqrt();
    let beta = alpha * NA / 0.2;
    let decay = (gamma * gamma + 4.0 * gamma * beta).sqrt();
    let positive_root = 2.0 * gamma / (gamma + decay);
    let negative_root = -(gamma + decay) / (2.0 * beta);
    let (initial, source) = fixture(0.0, 10000.0, 0.0, 0.0, 0.0);
    let photons = field(photon_energy, gamma);
    for time in [2e-5, 5e-3] {
        let attenuation = (-decay * time).exp();
        let expected = positive_root * (1.0 - attenuation)
            / (1.0 - (positive_root / negative_root) * attenuation);
        let mut numerical = controls(1e-6);
        numerical.initial_step_s = time;
        numerical.maximum_step_s = time;
        numerical.energy_tolerance_j = 1e-6;
        let result = advance_collisional_radiative(
            &initial,
            &source,
            &photons,
            Time::new(time),
            numerical,
            || false,
        )
        .unwrap();
        let state = result.state();
        assert!((state.amounts()[1].value() - expected).abs() < 2e-6);
        assert!((state.electron_temperature().value() - temperature).abs() < 2e-5);
        assert_eq!(state.amounts()[1], state.amounts()[2]);
        assert!(state.amounts().iter().all(|n| n.value() >= 0.0));
        assert!(result.balance().collisional_extent_mol.abs() < 1e-12);
        assert!(result.balance().energy_residual_j.abs() < numerical.energy_tolerance_j);
        assert!(
            (state.internal_energy().value() - initial.internal_energy().value()
                + result.balance().emitted_radiation_energy_j
                - result.balance().absorbed_radiation_energy_j)
                .abs()
                < 2e-6
        );
        if time > 1e-3 {
            assert!(decay * time > 20.0);
            assert!(result.balance().photoionization_extent_mol > 5.0 * state.amounts()[1].value());
            assert!(result.balance().radiative_recombination_extent_mol > 1.0);
            assert!(result.balance().rejected_steps > 0);
        }
    }
}

#[test]
fn g3_refinement_and_formation_reference_shifts_preserve_the_coupled_physics() {
    let photons = field(20.0 * EV, 1000.0);
    let (initial, source) = fixture(0.1, 10000.0, 1.0, 0.0, 0.0);
    let coarse = advance_collisional_radiative(
        &initial,
        &source,
        &photons,
        Time::new(2e-5),
        controls(1e-4),
        || false,
    )
    .unwrap();
    let fine = advance_collisional_radiative(
        &initial,
        &source,
        &photons,
        Time::new(2e-5),
        controls(1e-6),
        || false,
    )
    .unwrap();
    let exact = 13659.196390910463;
    assert!(
        (fine.state().electron_temperature().value() - exact).abs()
            < 0.3 * (coarse.state().electron_temperature().value() - exact).abs()
    );
    let (shifted, shifted_source) = fixture(0.1, 10000.0, 1.0, 1e10, 1e6);
    let changed = advance_collisional_radiative(
        &shifted,
        &shifted_source,
        &photons,
        Time::new(2e-5),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert!(
        (changed.state().electron_temperature().value()
            - fine.state().electron_temperature().value())
        .abs()
            < 1e-4
    );
    assert!(
        (changed.state().amounts()[1].value() - fine.state().amounts()[1].value()).abs() < 1e-9
    );
    assert!(
        (changed.balance().emitted_radiation_energy_j - fine.balance().emitted_radiation_energy_j)
            .abs()
            < 1e-3
    );
    assert!(changed.balance().energy_residual_j.abs() < 1e-7);
}

#[test]
fn g0_g4_identities_and_refusals_leave_the_original_state_unchanged() {
    let (neutral, source) = fixture(0.0, 10000.0, 1.0, 0.0, 0.0);
    let dark = advance_collisional_radiative(
        &neutral,
        &source,
        &[],
        Time::new(1.0),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert_eq!(dark.state(), &neutral);
    assert_eq!(dark.balance().absorbed_photons, 0.0);
    assert_eq!(dark.balance().emitted_photons, 0.0);
    let (initial, source) = fixture(0.1, 10000.0, 1.0, 0.0, 0.0);
    let copy = initial.clone();
    let photons = field(20.0 * EV, 1000.0);
    let zero = advance_collisional_radiative(
        &initial,
        &source,
        &photons,
        Time::new(0.0),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert_eq!(zero.state(), &initial);
    for seconds in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            advance_collisional_radiative(
                &initial,
                &source,
                &photons,
                Time::new(seconds),
                controls(1e-6),
                || false
            )
            .is_err()
        );
    }
    let mut limited = controls(1e-6);
    limited.maximum_steps = 1;
    assert_eq!(
        advance_collisional_radiative(
            &initial,
            &source,
            &photons,
            Time::new(4e-5),
            limited,
            || false
        ),
        Err(CollisionalRadiativeError::StepBudget)
    );
    assert_eq!(
        advance_collisional_radiative(
            &initial,
            &source,
            &photons,
            Time::new(4e-5),
            controls(1e-6),
            || true
        ),
        Err(CollisionalRadiativeError::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        advance_collisional_radiative(
            &initial,
            &source,
            &photons,
            Time::new(4e-5),
            controls(1e-6),
            || {
                polls += 1;
                polls > 1000
            }
        ),
        Err(CollisionalRadiativeError::Cancelled)
    );
    let (different, _) = fixture(0.1, 10000.0, 1.0, 1e6, 0.0);
    assert!(matches!(
        advance_collisional_radiative(
            &different,
            &source,
            &photons,
            Time::new(1e-6),
            controls(1e-6),
            || false
        ),
        Err(CollisionalRadiativeError::Invalid(_))
    ));
    assert_eq!(initial, copy);
}
