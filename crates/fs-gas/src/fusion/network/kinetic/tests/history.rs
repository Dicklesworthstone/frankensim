//! Source-history recording preserves the coupled kinetic and thermal solve.

use super::*;
use crate::fusion::FusionEnergyConfinement;
use crate::fusion::tests::history::check_history;
use crate::fusion::volume::FusionVolumeKnot;
use fs_qty::Volume;

#[test]
fn g3_kinetic_history_retains_neutron_moments_without_changing_alpha_evolution() {
    let (initial, alpha) = scenario(1e16, &grid(5000.0));
    let mut thermal = thermal_transport(&alpha, 0.15);
    thermal.energy_confinement = Some(
        FusionEnergyConfinement::new(
            Time::new(0.4),
            Time::new(0.3),
            Temperature::new(2.0 * KELVIN_PER_KEV),
            Temperature::new(KELVIN_PER_KEV),
        )
        .unwrap(),
    );
    let volume = FusionVolumeHistory::new(
        [(0.0, 1.0), (0.003, 0.95)]
            .into_iter()
            .map(|(time, volume)| FusionVolumeKnot {
                time: Time::new(time),
                volume: Volume::new(volume),
            })
            .collect(),
        || false,
    )
    .unwrap();
    let mut c = controls();
    c.relative_tolerance = 1e-5;
    let plain = advance_kinetic_fusion_with_volume(
        &initial,
        DT,
        heating(),
        &thermal,
        &alpha,
        &volume,
        c,
        || false,
    )
    .unwrap();
    c.record_neutron_history = true;
    let recorded = advance_kinetic_fusion_with_volume(
        &initial,
        DT,
        heating(),
        &thermal,
        &alpha,
        &volume,
        c,
        || false,
    )
    .unwrap();
    assert!(plain.neutron_history().is_none());
    assert_eq!(recorded.state(), plain.state());
    assert_eq!(recorded.balance(), plain.balance());
    assert_eq!(recorded.alpha_balance(), plain.alpha_balance());
    check_accounts(&recorded);
    check_history(
        recorded.neutron_history(),
        recorded.balance().accepted_steps,
        recorded.neutron_emission().unwrap(),
    );
    assert!(recorded.alpha_balance().born_alpha_amount_mol > 0.0);
    assert!(recorded.alpha_balance().deposited_electron_energy_j > 0.0);
    assert!(recorded.alpha_balance().boundary_work_on_alphas_j > 0.0);
    assert!(recorded.balance().heavy_transport_energy_out_j > 0.0);
}
