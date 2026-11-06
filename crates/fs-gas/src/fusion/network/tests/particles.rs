//! Shared open-particle maps, independent continuum calorics and pressure work.

use super::*;
use crate::fusion::volume::FusionVolumeKnot;
use crate::fusion::{FusionEnergyConfinement, FusionParticleExchange};
use fs_qty::{Frequency, Volume};

const CHARGE: [f64; 5] = [1.0, 1.0, 1.0, 2.0, 2.0];
const SUPPLY: [f64; 5] = [1e-6, 2e-5, 1e-5, 4e-6, 1e-6];
const REMOVAL: [f64; 5] = [0.5, 2.0, 3.0, 4.0, 15.0];

fn particle_source() -> FusionParticleExchange {
    FusionParticleExchange::new(
        SUPPLY,
        REMOVAL.map(Frequency::new),
        Temperature::new(2.0 * KELVIN_PER_KEV),
        Temperature::new(KELVIN_PER_KEV),
    )
    .unwrap()
}

fn exact_amounts(before: [f64; 5], time: f64) -> [f64; 5] {
    core::array::from_fn(|i| {
        before[i] * (-REMOVAL[i] * time).exp()
            + SUPPLY[i] / REMOVAL[i] * -(-REMOVAL[i] * time).exp_m1()
    })
}

// Direct continuum ODE, not the production split maps. Analytic independent
// stocks enter U'_s = a U_s + Cv T_in sum(Z S) - U_s sum(Z lambda n)/sum(Z n).
// The final four entries independently integrate exhaust energy and PdV work.
fn reference(
    before: [f64; 5],
    energies: [f64; 2],
    duration: f64,
    volume_ratio: f64,
    steps: u32,
) -> [f64; 6] {
    let compression = -(2.0 / 3.0) * volume_ratio.ln() / duration;
    let rhs = |t: f64, y: [f64; 6]| {
        let stocks = exact_amounts(before, t);
        let mut derivative = [0.0; 6];
        for pool in 0..2 {
            let weights = if pool == 0 { [1.0; 5] } else { CHARGE };
            let inventory: f64 = (0..5).map(|i| weights[i] * stocks[i]).sum();
            let outgoing: f64 = (0..5).map(|i| weights[i] * REMOVAL[i] * stocks[i]).sum();
            let incoming: f64 = (0..5).map(|i| weights[i] * SUPPLY[i]).sum();
            let carried_out = y[pool] * outgoing / inventory;
            let work = compression * y[pool];
            derivative[pool] =
                work + CV_MOLAR * [2.0, 1.0][pool] * KELVIN_PER_KEV * incoming - carried_out;
            derivative[pool + 2] = carried_out;
            derivative[pool + 4] = work;
        }
        derivative
    };
    let mut y = [energies[0], energies[1], 0.0, 0.0, 0.0, 0.0];
    let dt = duration / f64::from(steps);
    for i in 0..steps {
        let t = f64::from(i) * dt;
        let k1 = rhs(t, y);
        let k2 = rhs(
            t + 0.5 * dt,
            core::array::from_fn(|j| y[j] + 0.5 * dt * k1[j]),
        );
        let k3 = rhs(
            t + 0.5 * dt,
            core::array::from_fn(|j| y[j] + 0.5 * dt * k2[j]),
        );
        let k4 = rhs(t + dt, core::array::from_fn(|j| y[j] + dt * k3[j]));
        for j in 0..6 {
            y[j] += dt / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        }
    }
    y
}

#[test]
fn g1_particles_five_ion_sources_match_independent_amount_and_energy_odes_under_compression() {
    let (plasma, initial) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let original = initial.clone();
    let duration = 0.1;
    let before = plasma.indices().map(|i| initial.amounts()[i].value());
    let energies = [
        initial.heavy_internal_energy().value(),
        initial.electron_internal_energy().value(),
    ];
    let transport = FusionNetworkTransport {
        particles: Some(particle_source()),
        ..FusionNetworkTransport::default()
    };
    for ratio in [1.0, 0.7] {
        let geometry = FusionVolumeHistory::new(
            vec![
                FusionVolumeKnot {
                    time: Time::new(0.0),
                    volume: initial.volume(),
                },
                FusionVolumeKnot {
                    time: Time::new(duration),
                    volume: Volume::new(initial.volume().value() * ratio),
                },
            ],
            || false,
        )
        .unwrap();
        let result = advance_fusion_network_with_volume(
            &initial,
            &plasma,
            FusionNetworkChannels::NONE,
            heating(),
            &transport,
            &geometry,
            controls(),
            || false,
        )
        .unwrap();
        let expected_amounts = exact_amounts(before, duration);
        let expected = reference(before, energies, duration, ratio, 4096);
        let reference_half_resolution = reference(before, energies, duration, ratio, 2048);
        for (fine, coarse) in expected.into_iter().zip(reference_half_resolution) {
            close(fine, coarse, 2e-12, 1e-8);
        }
        let b = result.balance();
        let p = b.particles;
        for (slot, index) in plasma.indices().into_iter().enumerate() {
            close(
                result.state().amounts()[index].value(),
                expected_amounts[slot],
                8e-8,
                1e-18,
            );
            close(
                p.injected_ions_mol[slot],
                SUPPLY[slot] * duration,
                3e-14,
                1e-20,
            );
            close(
                p.exhausted_ions_mol[slot],
                before[slot] + SUPPLY[slot] * duration - expected_amounts[slot],
                1e-7,
                1e-18,
            );
        }
        let final_energy = [
            result.state().heavy_internal_energy().value(),
            result.state().electron_internal_energy().value(),
        ];
        let work = [b.boundary_work_on_heavy_j, b.boundary_work_on_electrons_j];
        for pool in 0..2 {
            close(final_energy[pool], expected[pool], 2e-7, 1e-6);
            close(
                p.exhausted_thermal_energy_j[pool],
                expected[pool + 2],
                2e-7,
                1e-6,
            );
            close(work[pool], expected[pool + 4], 2e-7, 1e-6);
            close(
                final_energy[pool] - energies[pool] + p.exhausted_thermal_energy_j[pool],
                work[pool] + p.injected_thermal_energy_j[pool],
                4e-12,
                1e-7,
            );
        }
        close(
            p.injected_electrons_mol,
            (0..5).map(|i| CHARGE[i] * p.injected_ions_mol[i]).sum(),
            3e-14,
            1e-20,
        );
        close(
            p.exhausted_electrons_mol,
            (0..5).map(|i| CHARGE[i] * p.exhausted_ions_mol[i]).sum(),
            3e-14,
            1e-20,
        );
        assert_eq!(b.reaction_extents_mol, [0.0; 4]);
        assert!(b.total_energy_residual_j.abs() < 1e-6);
        assert!(b.relative_baryon_residual.abs() < 1e-12);
        assert!(b.relative_charge_residual.abs() < 1e-12);
    }
    // Simultaneous supply/exhaust is noncommuting. Verify the actual fixed-step
    // method refines to the exact independent He4 inventory at second order.
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let mut c = controls();
        c.relative_tolerance = 0.1;
        c.absolute_temperature_tolerance_k = 1e6;
        c.energy_tolerance_j = 1e3;
        c.initial_step_s = duration / f64::from(steps);
        c.minimum_step_s = c.initial_step_s;
        c.maximum_step_s = c.initial_step_s;
        let result = advance_fusion_network(
            &initial,
            &plasma,
            FusionNetworkChannels::NONE,
            heating(),
            &transport,
            Time::new(duration),
            c,
            || false,
        )
        .unwrap();
        assert_eq!(result.balance().rejected_steps, 0);
        errors.push(
            (result.state().amounts()[plasma.indices()[4]].value()
                - exact_amounts(before, duration)[4])
                .abs(),
        );
    }
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    assert_eq!(initial, original);
}

#[test]
fn g3_particles_fueling_and_ash_exhaust_share_all_channels_and_thermal_operators() {
    let (plasma, initial) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let mut transport = transport(&initial);
    transport.particles = Some(particle_source());
    transport.energy_confinement = Some(
        FusionEnergyConfinement::new(
            Time::new(0.2),
            Time::new(0.13),
            Temperature::new(2.0 * KELVIN_PER_KEV),
            Temperature::new(KELVIN_PER_KEV),
        )
        .unwrap(),
    );
    let result = advance_fusion_network(
        &initial,
        &plasma,
        FusionNetworkChannels::ALL,
        mixed_heating(),
        &transport,
        Time::new(0.03),
        controls(),
        || false,
    )
    .unwrap();
    let b = result.balance();
    assert!(b.reaction_extents_mol.iter().all(|x| *x > 0.0));
    assert!(b.photon_energy_escaped_j > 0.0);
    assert!(b.electron_to_heavy_coulomb_energy_j != 0.0);
    assert!(b.heavy_transport_energy_out_j > 0.0 && b.electron_transport_energy_out_j > 0.0);
    assert!(b.particles.injected_ions_mol.iter().all(|x| *x > 0.0));
    assert!(b.particles.exhausted_ions_mol.iter().all(|x| *x > 0.0));
    close(
        result.state().internal_energy().value() - initial.internal_energy().value()
            + b.neutron_energy_escaped_j
            + b.photon_energy_escaped_j
            + b.heavy_transport_energy_out_j
            + b.electron_transport_energy_out_j
            + b.particles.exhausted_thermal_energy_j.iter().sum::<f64>(),
        b.nuclear_energy_released_j + b.particles.injected_thermal_energy_j.iter().sum::<f64>(),
        3e-12,
        1e-6,
    );
    for (slot, index) in plasma.indices().into_iter().enumerate() {
        let species = [
            NuclearSpeciesV1::Proton,
            NuclearSpeciesV1::Deuterium,
            NuclearSpeciesV1::Tritium,
            NuclearSpeciesV1::Helium3,
            NuclearSpeciesV1::Helium4,
        ][slot];
        let produced: f64 = FusionReactionV1::ALL
            .into_iter()
            .zip(b.reaction_extents_mol)
            .map(|(r, x)| f64::from(r.stoichiometric_coefficient(species)) * x)
            .sum();
        close(
            result.state().amounts()[index].value() - initial.amounts()[index].value(),
            produced + b.particles.injected_ions_mol[slot] - b.particles.exhausted_ions_mol[slot],
            3e-12,
            1e-18,
        );
    }
    assert_eq!(result.state().thermodynamics(), initial.thermodynamics());
}
