//! Original-source heat/work/entropy conservation for a moving reacting parcel.
use super::*;
use crate::thermochemical_gas::energy::equilibrium::tests::{fixture, isomers};
use crate::thermochemical_gas::tests::{
    point, species_with_energy_reference, transport, transport_in_domain,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::equilibrium::composition::{
    EquilibriumSolveV1, volume::EquilibriumVolumeSolveV1,
};

fn controls() -> DrivenEquilibriumGasSolve {
    let mut chemistry = EquilibriumSolveV1::default();
    chemistry.relative_balance_tolerance = 1e-13;
    chemistry.dimensionless_potential_tolerance = 1e-12;
    DrivenEquilibriumGasSolve {
        equilibrium: EquilibriumInternalEnergySolveV1 {
            volume: EquilibriumVolumeSolveV1 {
                equilibrium: chemistry,
                relative_volume_tolerance: 1e-11,
                maximum_pressure_iterations: 96,
            },
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1500.0,
            energy_tolerance_j: 1e-7,
            maximum_temperature_iterations: 128,
        },
        mass_tolerance_kg: 1e-12,
        energy_tolerance_j: 1e-5,
        entropy_tolerance_j_k: 1e-8,
    }
}
fn feed(amount: f64) -> Vec<(SpeciesId, Amount)> {
    vec![(SpeciesId::new("A").unwrap(), Amount::new(amount))]
}
fn parcel<'a>(
    bank: &'a ThermochemicalGasModel,
    chemistry: &'a IdealGasEquilibriumV1,
) -> DrivenEquilibriumGasParcel<'a> {
    DrivenEquilibriumGasParcel::new(
        bank,
        chemistry,
        &point(300.0, 1e5).with("sample_id", 7.0).unwrap(),
        &feed(2.0),
        Temperature::new(900.0),
        Volume::new(0.02),
        controls(),
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn fraction(t: f64) -> f64 {
    1.0 / (1.0 + (-24000.0 / (R * t)).exp())
}
fn mixing(t: f64) -> f64 {
    let y = fraction(t);
    -y * y.ln() - (1.0 - y) * (1.0 - y).ln()
}
fn energy(t: f64) -> f64 {
    2.0 * (3.0 * R * (t - 900.0) - 24000.0 * (fraction(t) - fraction(900.0)))
}
fn temperature(target: f64) -> f64 {
    let (mut lo, mut hi) = (200.0, 1500.0);
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        if energy(t) < target {
            lo = t;
        } else {
            hi = t;
        }
    }
    0.5 * (lo + hi)
}
fn entropy(t: f64, v: f64) -> f64 {
    2.0 * R * (3.0 * (t / 900.0).ln() + (v / 0.02).ln() + mixing(t) - mixing(900.0))
}

#[test]
fn g1_coupled_heat_work_and_volume_match_independent_reacting_balance() {
    let (bank, chemistry) = isomers();
    let mut gas = parcel(&bank, &chemistry);
    let original = gas.chamber().clone();
    let mass = gas.mass_kg();
    let mut target = 0.0;
    let mut previous_s = 0.0;
    for (q, w, v) in [
        (6000.0, 0.0, 0.02),
        (0.0, 2000.0, 0.03),
        (-3000.0, -1000.0, 0.015),
        (-2000.0, 0.0, 0.02),
    ] {
        target += q - w;
        let step = gas
            .advance(Volume::new(v), Energy::new(q), Energy::new(w), || false)
            .unwrap();
        let t = temperature(target);
        close(gas.sample().temperature.value(), t, 2e-6);
        close(
            gas.chamber().parcel().state().state().pressure,
            2.0 * R * t / v,
            2e-4,
        );
        close(
            gas.chamber().equilibrium().equilibrium().species()[1].mole_fraction,
            fraction(t),
            2e-10,
        );
        close(
            step.stored_energy_change_j,
            q - w,
            controls().energy_tolerance_j,
        );
        assert!(step.first_law_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(step.reference_energy_residual_j.abs() <= controls().energy_tolerance_j);
        close(
            step.internal_energy_change_since_reference_j,
            target,
            controls().energy_tolerance_j,
        );
        close(step.entropy_change_since_reference_j_k, entropy(t, v), 2e-8);
        close(step.entropy_change_j_k, entropy(t, v) - previous_s, 3e-8);
        assert_eq!(gas.energy_input_since_reference_j(), target);
        assert_eq!(gas.mass_kg(), mass);
        assert_eq!(gas.volume_m3(), v);
        assert_eq!(gas.feed(), feed(2.0));
        assert_eq!(gas.chamber().parcel().point().axes()["sample_id"], 7.0);
        assert_eq!(
            gas.chamber().parcel().model().transport_materials(),
            bank.transport_materials()
        );
        previous_s = entropy(t, v);
    }
    close(gas.sample().temperature.value(), 900.0, 2e-6);
    close(
        gas.chamber().balance().parcel_internal_energy_j,
        original.balance().parcel_internal_energy_j,
        controls().energy_tolerance_j,
    );
}

#[test]
fn g1_reservoir_entropy_accepts_heat_and_refuses_impossible_transfer_atomically() {
    let (bank, chemistry) = isomers();
    let mut gas = parcel(&bank, &chemistry);
    let old = gas.chamber().clone();
    let old_sample = gas.sample();
    let refused = gas.advance_with_reservoir(
        Volume::new(0.02),
        Energy::new(1000.0),
        Energy::new(0.0),
        Temperature::new(300.0),
        1e-8,
        || false,
    );
    assert!(matches!(refused, Err(Error::Entropy(_))));
    assert_eq!(gas.chamber(), &old);
    assert_eq!(gas.sample(), old_sample);
    assert_eq!(gas.energy_input_since_reference_j(), 0.0);
    let heating = gas
        .advance_with_reservoir(
            Volume::new(0.02),
            Energy::new(5000.0),
            Energy::new(0.0),
            Temperature::new(1400.0),
            1e-8,
            || false,
        )
        .unwrap();
    close(
        heating.entropy_production_j_k,
        entropy(temperature(5000.0), 0.02) - 5000.0 / 1400.0,
        2e-8,
    );
    assert!(heating.entropy_production_j_k > 0.0);
    let cooling = gas
        .advance_with_reservoir(
            Volume::new(0.02),
            Energy::new(-2000.0),
            Energy::new(0.0),
            Temperature::new(300.0),
            1e-8,
            || false,
        )
        .unwrap();
    assert!(cooling.entropy_production_j_k > 0.0);
    let before = gas.sample().temperature.value();
    let expansion = gas
        .advance_with_reservoir(
            Volume::new(0.04),
            Energy::new(0.0),
            Energy::new(0.0),
            Temperature::new(300.0),
            1e-8,
            || false,
        )
        .unwrap();
    close(gas.sample().temperature.value(), before, 2e-6);
    close(
        expansion.entropy_production_j_k,
        2.0 * R * 2.0_f64.ln(),
        3e-8,
    );
    assert!(expansion.entropy_production_j_k > 0.0);
}

#[test]
fn g0_g4_noop_cancelled_and_invalid_updates_leave_accepted_state_unchanged() {
    let (bank, chemistry) = isomers();
    let mut gas = parcel(&bank, &chemistry);
    let initial = gas.chamber().clone();
    let sample = gas.sample();
    let noop = gas.add_heat(Energy::new(0.0), || false).unwrap();
    assert_eq!(
        (noop.pressure_iterations, noop.composition_iterations),
        (0, 0)
    );
    assert_eq!(gas.chamber(), &initial);
    assert_eq!(gas.sample(), sample);
    for (v, q, w) in [
        (0.0, 0.0, 0.0),
        (f64::NAN, 0.0, 0.0),
        (0.02, f64::INFINITY, 0.0),
        (0.02, 1.0, 1e20),
        (0.02, 1e12, 0.0),
    ] {
        assert!(
            gas.advance(Volume::new(v), Energy::new(q), Energy::new(w), || false)
                .is_err()
        );
        assert_eq!(gas.chamber(), &initial);
        assert_eq!(gas.sample(), sample);
        assert_eq!(gas.energy_input_since_reference_j(), 0.0);
    }
    let mut polls = 0;
    let mut counted = parcel(&bank, &chemistry);
    counted
        .add_heat(Energy::new(5000.0), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut trial = parcel(&bank, &chemistry);
        let before = trial.chamber().clone();
        let mut count = 0;
        assert!(matches!(
            trial.add_heat(Energy::new(5000.0), || {
                count += 1;
                count >= stop
            }),
            Err(Error::Cancelled)
        ));
        assert_eq!(trial.chamber(), &before);
        assert_eq!(trial.energy_input_since_reference_j(), 0.0);
        trial.add_heat(Energy::new(5000.0), || false).unwrap();
        close(
            trial.sample().temperature.value(),
            temperature(5000.0),
            2e-6,
        );
    }
}

#[test]
fn g3_transport_refusal_does_not_commit_chemistry_or_heat() {
    let (bank, chemistry) = isomers();
    let restricted = ThermochemicalGasModel::new(
        bank.thermodynamics().clone(),
        vec![
            transport_in_domain("A", 1.7e-5, 0.025, [200.0, 1000.0], [1000.0, 1e7]),
            transport_in_domain("B", 2e-5, 0.15, [200.0, 1000.0], [1000.0, 1e7]),
        ],
        None,
    )
    .unwrap();
    let mut gas = parcel(&restricted, &chemistry);
    let before = gas.chamber().clone();
    assert!(gas.add_heat(Energy::new(18000.0), || false).is_err());
    assert_eq!(gas.chamber(), &before);
    assert_eq!(gas.energy_input_since_reference_j(), 0.0);
    gas.add_heat(Energy::new(1500.0), || false).unwrap();
    close(gas.sample().temperature.value(), temperature(1500.0), 2e-6);
}

#[test]
fn g3_extensive_scaling_preserves_reacting_temperature_pressure_and_entropy_per_mass() {
    let (bank, chemistry) = isomers();
    for scale in [0.25, 1.0, 4.0] {
        let mut c = controls();
        c.mass_tolerance_kg *= scale;
        c.energy_tolerance_j *= scale;
        c.entropy_tolerance_j_k *= scale;
        c.equilibrium.energy_tolerance_j *= scale;
        let mut gas = DrivenEquilibriumGasParcel::new(
            &bank,
            &chemistry,
            &point(300.0, 1e5),
            &feed(2.0 * scale),
            Temperature::new(900.0),
            Volume::new(0.02 * scale),
            c,
            || false,
        )
        .unwrap();
        let step = gas
            .advance(
                Volume::new(0.03 * scale),
                Energy::new(7000.0 * scale),
                Energy::new(3000.0 * scale),
                || false,
            )
            .unwrap();
        close(gas.sample().temperature.value(), temperature(4000.0), 2e-6);
        close(
            gas.sample().pressure.value(),
            2.0 * R * temperature(4000.0) / 0.03,
            2e-4,
        );
        close(step.stored_energy_change_j / scale, 4000.0, 1e-5);
        close(
            step.entropy_change_j_k / scale,
            entropy(temperature(4000.0), 0.03),
            3e-8,
        );
        assert_eq!(gas.mass_kg(), 0.02 * scale);
    }
}

#[test]
fn g3_small_heat_uses_source_difference_at_a_large_formation_reference() {
    let q = 0.001;
    for offset in [0.0, 1e11] {
        let (bank, chemistry) = fixture(
            vec![species_with_energy_reference("A", 0.01, 4.0, 0.0, offset)],
            vec![transport("A", 1.7e-5, 0.025)],
            &["H"],
            vec![vec![1]],
        );
        let mut c = controls();
        c.equilibrium.energy_tolerance_j = 1e-9;
        let mut gas = DrivenEquilibriumGasParcel::new(
            &bank,
            &chemistry,
            &point(300.0, 1e5),
            &feed(2.0),
            Temperature::new(900.0),
            Volume::new(0.02),
            c,
            || false,
        )
        .unwrap();
        let step = gas.add_heat(Energy::new(q), || false).unwrap();
        close(
            gas.sample().temperature.value(),
            900.0 + q / (6.0 * R),
            3e-10,
        );
        close(step.stored_energy_change_j, q, 2e-8);
        assert!(step.first_law_residual_j.abs() < 2e-8);
        gas.add_heat(Energy::new(-q), || false).unwrap();
        close(gas.sample().temperature.value(), 900.0, 3e-10);
        assert_eq!(gas.energy_input_since_reference_j(), 0.0);
    }
}

#[test]
fn g1_adiabatic_update_refuses_entropy_decrease_before_any_publication() {
    let (bank, chemistry) = isomers();
    let mut gas = parcel(&bank, &chemistry);
    let before = gas.chamber().clone();
    let sample = gas.sample();
    // Work extraction from a rigid, closed, zero-heat equilibrium gas would
    // lower its temperature and entropy. Likewise, cost-free compression
    // cannot retain this energy and halve the volume.
    for (volume, work) in [(0.02, 1000.0), (0.01, 0.0)] {
        let refused = gas.advance(
            Volume::new(volume), Energy::new(0.0), Energy::new(work), || false,
        );
        let Err(Error::AdiabaticEntropy(report)) = refused else {
            panic!("entropy-decreasing adiabatic update was not refused");
        };
        assert!(report.entropy_change_j_k < -controls().entropy_tolerance_j_k);
        assert_eq!(gas.chamber(), &before);
        assert_eq!(gas.sample(), sample);
        assert_eq!(gas.energy_input_since_reference_j(), 0.0);
    }
    // An expansion delivering the same extracted work has a positive
    // entropy increment and must still succeed after the rejected trials.
    let accepted = gas.advance(
        Volume::new(0.03), Energy::new(0.0), Energy::new(1000.0), || false,
    ).unwrap();
    assert!(accepted.entropy_change_j_k > 0.0);
    close(accepted.stored_energy_change_j, -1000.0, controls().energy_tolerance_j);
    close(gas.sample().temperature.value(), temperature(-1000.0), 2e-6);
}
