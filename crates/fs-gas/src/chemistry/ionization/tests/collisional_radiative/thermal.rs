//! Independent thermal limits and coupled source-equation references.

use super::super::{MP, species};
use super::*;
use crate::chemistry::energy::{ElectronHeavyCollisionSources, ElectronHeavyEnergyError};
use fs_evidence::ValidityDomain;
use fs_qty::{Amount, MolarMass, Power};
use fs_thermochem::coulomb::{CoulombLimitsV1, CoulombModelV1, CoulombScreeningV1};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, Nasa9RegionV1, Nasa9StandardStateModelV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

const CV: f64 = 1.5 * R;

fn photon_source(channel: AtomicIonizationV1, sigma: f64) -> AtomicRadiativeRecombinationV1 {
    let gap = channel.ionization_energy_j_mol() / NA;
    let provenance = Provenance {
        source: "analytic constant bound-free source, not measured atom data".into(),
        license: "analytic test fixture".into(),
        artifact: None,
    };
    let photon = AtomicPhotoionizationV1::new(
        channel,
        provenance.clone(),
        vec![
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap),
                cross_section: Area::new(sigma),
            },
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap + 1000.0 * EV),
                cross_section: Area::new(sigma),
            },
        ],
    )
    .unwrap();
    AtomicRadiativeRecombinationV1::new(photon, Area::new(sigma), provenance).unwrap()
}

fn parcel(
    x: f64,
    te: f64,
    th: f64,
    impact: f64,
    sigma: f64,
) -> (FrozenIonizedGasStateV1, AtomicRadiativeRecombinationV1) {
    let (original, source) = fixture(x, te, impact, 0.0, 0.0);
    let state = FrozenIonizedGasStateV1::new(
        original.thermodynamics(),
        original.amounts(),
        Volume::new(1000.0),
        Temperature::new(th),
        original.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    (
        state,
        photon_source(source.source().channel().clone(), sigma),
    )
}

fn collisions(state: &FrozenIonizedGasStateV1, coulomb: bool) -> ElectronHeavyCollisionSources {
    let momentum = ElectronCrossSectionV1::tabulated(
        id("H"),
        ElectronCollisionKindV1::ElasticMomentumTransfer,
        Provenance {
            source: "analytic constant neutral momentum cross section".into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
        vec![
            ElectronCrossSectionPointV1 {
                energy: Energy::new(0.0),
                cross_section: Area::new(1e-20),
            },
            ElectronCrossSectionPointV1 {
                energy: Energy::new(500.0 * EV),
                cross_section: Area::new(1e-20),
            },
        ],
        Area::new(1e-20),
    )
    .unwrap();
    ElectronHeavyCollisionSources {
        momentum_sources: vec![momentum],
        coulomb: coulomb.then(|| {
            CoulombModelV1::new(
                state,
                CoulombScreeningV1::StaticAllChargedDebye,
                CoulombLimitsV1 {
                    minimum_coulomb_logarithm: 4.0,
                    minimum_particles_per_debye_sphere: 10.0,
                    maximum_coulomb_coupling: 0.1,
                    maximum_thermal_phase_space_density: 0.1,
                    maximum_thermal_energy_rest_energy_ratio: 0.001,
                    maximum_quantum_to_classical_ratio: 0.5,
                    maximum_relative_drift_thermal_ratio: 0.1,
                },
                || false,
            )
            .unwrap()
        }),
        maximum_charged_mole_fraction: if coulomb { 1.0 } else { 0.5 },
        maximum_target_thermal_speed_ratio: 0.1,
        rate_controls: MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
    }
}

fn inputs() -> CollisionalRadiativeThermal {
    CollisionalRadiativeThermal::new(Temperature::new(1000.0), Temperature::new(50000.0))
}

fn assert_accounts(result: &CollisionalRadiativeResult) {
    let start = result.initial_state();
    let end = result.state();
    let channel = result.source().source().channel();
    let b = result.balance();
    let ne0 = start.amounts()[channel.electron_index()].value();
    let ne = end.amounts()[channel.electron_index()].value();
    let sensible = CV
        * (ne0 * (end.electron_temperature().value() - start.electron_temperature().value())
            + (ne - ne0) * end.electron_temperature().value());
    let chemical = channel.ionization_energy_j_mol()
        * (end.amounts()[channel.ion_index()].value()
            - start.amounts()[channel.ion_index()].value());
    // Independently read all published original-source heavy caloric changes.
    // This includes spectators, and does not use the integrated transfer target.
    let mut heavy = 0.0;
    for (index, (amount, model)) in start
        .amounts()
        .iter()
        .zip(start.thermodynamics().models())
        .enumerate()
    {
        if Some(index) != start.electron_index() {
            heavy += amount.value()
                * model
                    .molar_internal_energy_change(
                        start.heavy_temperature(),
                        end.heavy_temperature(),
                    )
                    .unwrap()
                    .value();
        }
    }
    assert!(
        (sensible + chemical + b.emitted_radiation_energy_j - b.absorbed_radiation_energy_j
            + b.electron_to_heavy_energy_j
            - b.electron_deposited_energy_j)
            .abs()
            < 3e-7,
        "{b:?}"
    );
    assert!(
        (heavy - b.electron_to_heavy_energy_j - b.heavy_deposited_energy_j).abs() < 2e-7,
        "{b:?}"
    );
    assert!((heavy - b.heavy_sensible_energy_change_j).abs() < 2e-8);
    assert!(
        (b.electron_to_heavy_neutral_energy_j + b.electron_to_heavy_coulomb_energy_j
            - b.electron_to_heavy_energy_j)
            .abs()
            < 1e-10
    );
    for residual in [
        b.energy_residual_j,
        b.electron_sensible_energy_residual_j,
        b.heavy_sensible_energy_residual_j,
    ] {
        assert!(
            residual.abs() <= result.controls().energy_tolerance_j,
            "{b:?}"
        );
    }
    assert_eq!(end.thermodynamics(), start.thermodynamics());
    assert_eq!(end.volume(), start.volume());
    assert!(end.relative_charge_residual() <= result.controls().relative_charge_tolerance);
    assert!(b.maximum_normalized_local_error <= 1.0);
}

#[test]
fn thermal_exchange_only_matches_nonlinear_neutral_relaxation_in_both_directions() {
    // Constant sigma gives the exact Maxwellian momentum moment
    // km=(4/3)*sigma*sqrt(8*kT/(pi*m)), so G=g0*sqrt(Te).
    // At fixed populations and Cv, y=sqrt(Te) satisfies
    // y'=-(g0/2)*(1/Ce+1/Ch)*(y²-Teq), with an exact rational exponential.
    // Both reaction sources are explicitly zero; no chemistry can advance time.
    for (te, th, time) in [(30000.0, 4000.0, 3e-3), (2000.0, 10000.0, 1e-3)] {
        let x = 0.05;
        let (state, source) = parcel(x, te, th, 0.0, 0.0);
        let mut thermal = inputs();
        thermal.collisions = Some(collisions(&state, false));
        let mut c = controls(1e-7);
        c.initial_step_s = time;
        c.maximum_step_s = time;
        let result = advance_collisional_radiative_with_thermal(
            &state,
            &source,
            &[],
            &thermal,
            Time::new(time),
            c,
            || false,
        )
        .unwrap();
        let equilibrium = (x * te + th) / (1.0 + x);
        let g0 = 3.0 * R * x * (1.0 - x) * NA / 1000.0
            * (ME / (MP + ME))
            * (4.0 / 3.0)
            * 1e-20
            * (8.0 * KB / (core::f64::consts::PI * ME)).sqrt();
        let attenuation = (-g0 * (1.0 / (CV * x) + 1.0 / CV) * equilibrium.sqrt() * time).exp();
        let ratio =
            (te.sqrt() - equilibrium.sqrt()) / (te.sqrt() + equilibrium.sqrt()) * attenuation;
        let expected_te = (equilibrium.sqrt() * (1.0 + ratio) / (1.0 - ratio)).powi(2);
        let expected_th = th + x * (te - expected_te);
        assert!((result.state().electron_temperature().value() - expected_te).abs() < 0.01);
        assert!((result.state().heavy_temperature().value() - expected_th).abs() < 0.001);
        assert_eq!(result.state().amounts(), state.amounts());
        assert_eq!(result.reaction_extent(), Amount::new(0.0));
        assert_eq!(result.balance().emitted_photons, 0.0);
        assert_eq!(result.balance().absorbed_photons, 0.0);
        assert_eq!(result.balance().electron_to_heavy_coulomb_energy_j, 0.0);
        assert!(result.balance().electron_to_heavy_energy_j * (te - th) > 0.0);
        assert!(result.balance().accepted_steps > 0);
        assert!(result.balance().accepted_steps < 2000);
        assert_accounts(&result);
    }
}

#[test]
fn thermal_reactions_radiation_both_collision_channels_and_powers_match_independent_ode() {
    let (state, source) = parcel(0.05, 30000.0, 4000.0, 1.0, 1e-20);
    let photons = [PhotonNumberFluxV1 {
        photon_energy: Energy::new(20.0 * EV),
        number_flux_per_m2_s: 3000.0 / 1e-20,
    }];
    let mut thermal = inputs();
    thermal.collisions = Some(collisions(&state, true));
    thermal.electron_power = Power::new(3e8);
    thermal.heavy_power = Power::new(2e7);
    // Independent DOP853 (rtol2.5e-13) and Radau (rtol1e-11) agree.
    // Integrate x, Ue, Uh and every account using the analytic impact/Milne
    // moments above, km's constant-sigma moment and finite-mass Landau G:
    // 4*sqrt(2pi)*kB*ne*ni*C²*ln(lambdaD/b90)/(me*mi*S^(3/2)).
    // S=kB*(Te/me+Th/mi), b90=C/(mu*S), static screening includes both pools.
    for (time, expected) in [
        (
            2e-6,
            [
                0.05740271548326162,
                19420.164844764524,
                4534.898400546476,
                0.005675473750024579,
                0.0017353300472077533,
                8.088313970718891e-6,
                10951.999394564988,
                12.649461831054222,
                208.12842521567413,
                6422.960708564786,
            ],
        ),
        (
            2e-5,
            [
                0.10678908995315124,
                8284.416293149965,
                7671.218564730493,
                0.05522696474967247,
                0.001812754642883156,
                0.0002506294394044328,
                106571.8407206897,
                349.44746943106156,
                420.7962249878713,
                44965.51805429495,
            ],
        ),
    ] {
        let result = advance_collisional_radiative_with_thermal(
            &state,
            &source,
            &photons,
            &thermal,
            Time::new(time),
            controls(1e-7),
            || false,
        )
        .unwrap();
        let end = result.state();
        let b = result.balance();
        let actual = [
            end.amounts()[1].value(),
            end.electron_temperature().value(),
            end.heavy_temperature().value(),
            b.photoionization_extent_mol,
            b.collisional_extent_mol,
            b.radiative_recombination_extent_mol,
            b.absorbed_radiation_energy_j,
            b.emitted_radiation_energy_j,
            b.electron_to_heavy_neutral_energy_j,
            b.electron_to_heavy_coulomb_energy_j,
        ];
        let allowances = [2e-7, 0.08, 0.02, 2e-7, 2e-7, 1e-8, 0.5, 0.03, 0.03, 0.3];
        for ((value, exact), allowance) in actual.into_iter().zip(expected).zip(allowances) {
            assert!(
                (value - exact).abs() < allowance,
                "{actual:?} != {expected:?}"
            );
        }
        assert!((b.electron_deposited_energy_j - 3e8 * time).abs() < 1e-8);
        assert!((b.heavy_deposited_energy_j - 2e7 * time).abs() < 1e-9);
        assert_eq!(result.thermal_inputs(), Some(&thermal));
        assert_eq!(end.amounts()[1], end.amounts()[2]);
        assert!((end.amounts()[0].value() + end.amounts()[1].value() - 1.0).abs() < 2e-15);
        assert!(
            (end.internal_energy().value() - state.internal_energy().value()
                + b.emitted_radiation_energy_j
                - b.absorbed_radiation_energy_j
                - b.electron_deposited_energy_j
                - b.heavy_deposited_energy_j)
                .abs()
                < 3e-7
        );
        assert_accounts(&result);
    }
}

#[test]
fn thermal_held_powers_recover_variable_spectator_calorics_and_reference_shifts() {
    for (elemental, charge) in [(0.0, 0.0), (1e12, 1e8)] {
        let (original, source) = fixture(0.05, 30000.0, 0.0, elemental, charge);
        let background = species("N2", 0.028 / NA, 1.0, elemental * 2.0);
        let mut card = background.regions()[0].card().clone();
        card.parameters.get_mut("a2").unwrap().value = 3.5;
        card.parameters.get_mut("a3").unwrap().value = 1e-4;
        card.validity = ValidityDomain::unconstrained().with("T", 1000.0, 8000.0);
        let nitrogen = Nasa9StandardStateModelV1::new(
            id("N2"),
            MolarMass::new(0.028),
            original.thermodynamics().convention().clone(),
            vec![Nasa9RegionV1::from_card(card).unwrap()],
        )
        .unwrap();
        let ids = vec![id("H"), id("H+"), id("N2"), id("e")];
        let bank = IdealGasEquilibriumV1::new(
            ElementalMatrix::new(
                vec![ElementId::new("H").unwrap(), ElementId::new("N").unwrap()],
                ids.clone(),
                vec![vec![1, 1, 0, 0], vec![0, 0, 2, 0]],
            )
            .unwrap(),
            ChargeVector::new(ids, vec![0, 1, 0, -1]).unwrap(),
            vec![
                original.thermodynamics().models()[0].clone(),
                original.thermodynamics().models()[1].clone(),
                nitrogen,
                original.thermodynamics().models()[2].clone(),
            ],
        )
        .unwrap();
        let state = FrozenIonizedGasStateV1::new(
            &bank,
            &[
                Amount::new(0.95),
                Amount::new(0.05),
                Amount::new(2.0),
                Amount::new(0.05),
            ],
            original.volume(),
            Temperature::new(4000.0),
            original.electron_temperature(),
            1e-12,
            || false,
        )
        .unwrap();
        let channel = AtomicIonizationV1::new(
            bank,
            ReactionId::new("power-only-H-in-N2").unwrap(),
            id("H+"),
            source.source().channel().forward_rate().clone(),
            AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
            1e-12,
        )
        .unwrap();
        let source = photon_source(channel, 0.0);
        for (th, te) in [(6000.0_f64, 35000.0), (3000.0, 25000.0)] {
            let time = 1e-4;
            // Two moles of spectator have Cv/R=2.5+1e-4*T. Thus the
            // independent full heavy primitive is R*(6.5*T+1e-4*T²).
            let heavy_energy = R * (6.5 * (th - 4000.0) + 1e-4 * (th * th - 4000.0_f64.powi(2)));
            let electron_energy = CV * 0.05 * (te - 30000.0);
            let mut thermal = inputs();
            thermal.maximum_heavy_temperature_k = 8000.0;
            thermal.heavy_power = Power::new(heavy_energy / time);
            thermal.electron_power = Power::new(electron_energy / time);
            let mut c = controls(1e-7);
            c.initial_step_s = time;
            let result = advance_collisional_radiative_with_thermal(
                &state,
                &source,
                &[],
                &thermal,
                Time::new(time),
                c,
                || false,
            )
            .unwrap();
            assert!((result.state().heavy_temperature().value() - th).abs() < 2e-9);
            assert!((result.state().electron_temperature().value() - te).abs() < 2e-9);
            assert_eq!(result.state().amounts(), state.amounts());
            assert_eq!(result.balance().accepted_steps, 1);
            assert_eq!(result.balance().electron_to_heavy_energy_j, 0.0);
            assert_eq!(result.balance().emitted_photons, 0.0);
            assert_accounts(&result);
        }
        let mut unsupported = inputs();
        unsupported.maximum_heavy_temperature_k = 9000.0;
        assert!(matches!(
            advance_collisional_radiative_with_thermal(
                &state,
                &source,
                &[],
                &unsupported,
                Time::new(0.0),
                controls(1e-6),
                || false
            ),
            Err(CollisionalRadiativeError::Admission(_))
        ));
    }
}

#[test]
fn thermal_neutral_birth_identity_and_source_refusals_preserve_inputs() {
    let (neutral, source) = parcel(0.0, 10000.0, 6000.0, 0.0, SIGMA);
    let photons = field(20.0 * EV, 1000.0);
    let mut thermal = inputs();
    thermal.collisions = Some(collisions(&neutral, false));
    let created = advance_collisional_radiative_with_thermal(
        &neutral,
        &source,
        &photons,
        &thermal,
        Time::new(1e-6),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert!(created.state().amounts()[2].value() > 0.0);
    assert!(created.state().heavy_temperature().value() > neutral.heavy_temperature().value());
    assert!(created.balance().electron_to_heavy_energy_j > 0.0);
    assert_accounts(&created);
    let mut missing = thermal.clone();
    missing
        .collisions
        .as_mut()
        .unwrap()
        .momentum_sources
        .clear();
    assert!(matches!(
        advance_collisional_radiative_with_thermal(
            &neutral,
            &source,
            &photons,
            &missing,
            Time::new(1e-6),
            controls(1e-6),
            || false
        ),
        Err(CollisionalRadiativeError::Thermal(
            ElectronHeavyEnergyError::MissingSource(_)
        ))
    ));
    for power in [1.0, -1.0] {
        let mut unsupported = thermal.clone();
        unsupported.electron_power = Power::new(power);
        assert!(matches!(
            advance_collisional_radiative_with_thermal(
                &neutral,
                &source,
                &photons,
                &unsupported,
                Time::new(1e-6),
                controls(1e-6),
                || false
            ),
            Err(CollisionalRadiativeError::Invalid(_))
        ));
    }
    // Heavy work remains meaningful at exact zero electron heat capacity.
    let (_, dark_source) = parcel(0.0, 10000.0, 6000.0, 0.0, 0.0);
    let mut heavy_only = inputs();
    heavy_only.heavy_power = Power::new(CV * 1000.0);
    let heated = advance_collisional_radiative_with_thermal(
        &neutral,
        &dark_source,
        &[],
        &heavy_only,
        Time::new(0.01),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert_eq!(
        heated.state().electron_temperature(),
        neutral.electron_temperature()
    );
    assert_eq!(heated.state().amounts(), neutral.amounts());
    assert!((heated.state().heavy_temperature().value() - 6010.0).abs() < 1e-9);
    assert_accounts(&heated);
    let (state, source) = parcel(0.05, 30000.0, 4000.0, 1.0, SIGMA);
    let copy = state.clone();
    let mut thermal = inputs();
    thermal.collisions = Some(collisions(&state, true));
    let identity = advance_collisional_radiative_with_thermal(
        &state,
        &source,
        &photons,
        &thermal,
        Time::new(0.0),
        controls(1e-6),
        || false,
    )
    .unwrap();
    assert_eq!(identity.state(), &state);
    assert_eq!(identity.balance().electron_to_heavy_energy_j, 0.0);
    let stale = parcel(0.06, 30000.0, 4000.0, 1.0, SIGMA).0;
    assert!(matches!(
        advance_collisional_radiative_with_thermal(
            &stale,
            &source,
            &photons,
            &thermal,
            Time::new(1e-6),
            controls(1e-6),
            || false
        ),
        Err(CollisionalRadiativeError::Thermal(
            ElectronHeavyEnergyError::Invalid(_)
        ))
    ));
    assert_eq!(
        advance_collisional_radiative_with_thermal(
            &state,
            &source,
            &photons,
            &thermal,
            Time::new(1e-6),
            controls(1e-6),
            || true
        ),
        Err(CollisionalRadiativeError::Cancelled)
    );
    let mut bounded = controls(1e-8);
    bounded.maximum_steps = 1;
    assert_eq!(
        advance_collisional_radiative_with_thermal(
            &state,
            &source,
            &photons,
            &thermal,
            Time::new(1e-5),
            bounded,
            || false
        ),
        Err(CollisionalRadiativeError::StepBudget)
    );
    assert_eq!(state, copy);
}
