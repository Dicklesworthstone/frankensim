//! G3 thermal confinement acts on the Maxwellian pools, retaining fast alphas.

use super::*;
use crate::fusion::FusionEnergyConfinement;
use crate::fusion::volume::FusionVolumeKnot;
use fs_qty::Volume;

#[test]
fn g3_kinetic_confinement_couples_births_slowing_radiation_and_pressure_work() {
    let (initial, alpha) = scenario(1e16, &grid(5000.0));
    let original = initial.clone();
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
    let geometry = FusionVolumeHistory::new(
        [(0.0, 1.0), (0.015, 0.9), (0.03, 1.0)]
            .into_iter()
            .map(|(time, volume)| FusionVolumeKnot {
                time: Time::new(time),
                volume: Volume::new(volume),
            })
            .collect(),
        || false,
    )
    .unwrap();
    let moving = advance_kinetic_fusion_with_volume(
        &initial,
        DT,
        heating(),
        &thermal,
        &alpha,
        &geometry,
        controls(),
        || false,
    )
    .unwrap();
    let rigid = advance_kinetic_fusion(
        &initial,
        DT,
        heating(),
        &thermal,
        &alpha,
        geometry.duration(),
        controls(),
        || false,
    )
    .unwrap();
    for result in [&moving, &rigid] {
        let b = result.balance();
        let a = result.alpha_balance();
        assert!(b.heavy_transport_energy_out_j > 0.0 && b.electron_transport_energy_out_j > 0.0);
        assert!(b.photon_energy_escaped_j > 0.0);
        assert!(a.born_alpha_amount_mol > 0.0);
        assert!(a.deposited_electron_energy_j > 0.0);
        close(
            result.state().alpha_amount().value(),
            initial.alpha_amount().value() + a.born_alpha_amount_mol,
            3e-12,
            1e-20,
        );
        assert_eq!(result.thermal_transport(), &thermal);
        check_accounts(result);
    }
    assert_eq!(moving.state().volume(), initial.volume());
    assert_eq!(moving.state().energy_grid(), initial.energy_grid());
    assert!(moving.alpha_balance().boundary_work_on_alphas_j != 0.0);
    assert_eq!(rigid.total_boundary_work().value(), 0.0);
    assert_eq!(initial, original);
}
