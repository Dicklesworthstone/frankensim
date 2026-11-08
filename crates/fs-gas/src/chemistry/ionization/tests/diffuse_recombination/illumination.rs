//! G1/G3/G4: exact directional attenuation and original diffuse superposition.
use super::*;
use fs_thermochem::kinetics::atomic_ionization::photoionization::PhotonNumberFluxV1;

fn beam(energy_ev: f64, flux: f64, mu: f64) -> DiffuseRecombinationBeam {
    DiffuseRecombinationBeam {
        photons: PhotonNumberFluxV1 {
            photon_energy: Energy::new(energy_ev * EV),
            number_flux_per_m2_s: flux,
        },
        direction_cosine: mu,
    }
}

#[test]
fn g1_oblique_boundary_beams_resolve_spectral_shielding_and_excess_heat() {
    let (a, sa) = parcel(0.0, 10000.0, 1.0, 0.2);
    let (b, sb) = parcel(0.0, 20000.0, 2.0, 0.4);
    let cells = [
        DiffuseRecombinationCell { state: &a, source: &sa },
        DiffuseRecombinationCell { state: &b, source: &sb },
    ];
    let area = Area::new(SIGMA * NA);
    let beams = [beam(20.0, 2e20, 0.5), beam(30.0, 3e20, -0.8)];
    let result = evaluate_diffuse_recombination_with_illumination(
        &cells, area, &beams, controls(1e-7), || false,
    ).unwrap();
    assert_eq!(result.incident_beams(), &beams);
    assert_eq!(result.balance().emitted_photons_per_s, 0.0);
    assert_eq!(result.balance().spectral_panels, 0);
    let left = area.value() * 0.5 * 2e20;
    let right = area.value() * 0.8 * 3e20;
    let chi = sa.source().channel().ionization_energy_j_mol() / NA;
    for (c, (left_in, right_in, depth)) in [
        (left, right * (-2.0_f64 / 0.8).exp(), 1.0),
        (left * (-1.0_f64 / 0.5).exp(), right, 2.0),
    ].into_iter().enumerate() {
        let from_left = left_in * -(-depth / 0.5_f64).exp_m1();
        let from_right = right_in * -(-depth / 0.8_f64).exp_m1();
        let absorbed = from_left + from_right;
        let power = EV * (20.0 * from_left + 30.0 * from_right);
        let heat = (20.0 * EV - chi) * from_left + (30.0 * EV - chi) * from_right;
        let balance = result.cells()[c].balance();
        close(balance.incident_absorbed_photons_per_s, absorbed, 4e-14);
        close(balance.incident_absorbed_power_w, power, 4e-14);
        close(balance.net_ionization_rate_mol_per_s * NA, absorbed, 4e-14);
        close(balance.net_electron_power_w, heat, 4e-14);
        close(balance.net_gas_power_w, power, 4e-14);
        assert_eq!(balance.self_absorbed_photons_per_s, 0.0);
        assert_eq!(balance.nonlocal_absorbed_photons_per_s, 0.0);
    }
    close(result.interfaces()[0].right_photons_per_s, left, 1e-15);
    close(result.interfaces()[2].left_photons_per_s, right, 1e-15);
    close(result.interfaces()[2].right_photons_per_s, left * (-6.0_f64).exp(), 4e-14);
    close(result.interfaces()[0].left_photons_per_s, right * (-3.75_f64).exp(), 4e-14);
    let balance = result.balance();
    close(balance.incident_photons_per_s, left + right, 1e-15);
    assert!(balance.photon_residual_per_s.abs() < 1e-14 * balance.incident_photons_per_s);
    assert!(balance.energy_residual_w.abs() < 1e-14 * balance.incident_power_w);
}

#[test]
fn g3_incidence_superposes_without_recounting_spontaneous_recycling() {
    let (state, source) = parcel(0.2, 10000.0, 1.0, 1e5);
    let cells = [DiffuseRecombinationCell { state: &state, source: &source }];
    let area = Area::new(SIGMA * NA);
    let dark = evaluate_diffuse_recombination(&cells, area, controls(1e-7), || false).unwrap();
    let beams = [beam(25.0, 1e22, 0.6), beam(10.0, 2e22, -1.0)];
    let lit = evaluate_diffuse_recombination_with_illumination(
        &cells, area, &beams, controls(1e-7), || false,
    ).unwrap();
    let d = dark.cells()[0].balance();
    let l = lit.cells()[0].balance();
    let absorbed = area.value() * 0.6 * 1e22 * -(-0.8_f64 / 0.6).exp_m1();
    assert_eq!(l.emitted_photons_per_s, d.emitted_photons_per_s);
    assert_eq!(l.self_absorbed_photons_per_s, d.self_absorbed_photons_per_s);
    assert_eq!(l.captured_electron_power_w, d.captured_electron_power_w);
    close(l.incident_absorbed_photons_per_s, absorbed, 4e-14);
    close(l.net_gas_power_w, d.net_gas_power_w + 25.0 * EV * absorbed, 4e-14);
    close(l.net_ionization_rate_mol_per_s, d.net_ionization_rate_mol_per_s + absorbed / NA, 4e-14);
    // The subthreshold beam traverses the complete slab without heating it.
    close(lit.interfaces()[0].left_photons_per_s - dark.interfaces()[0].left_photons_per_s,
        area.value() * 2e22, 4e-14);
    assert_eq!(lit.cells()[0].state(), &state);
}

#[test]
fn g1_tiny_optical_depth_retains_absorption_below_transmission_roundoff() {
    let (state, source) = parcel(0.0, 10000.0, 1.0, 0.2);
    let area = Area::new(SIGMA * NA / 1e-20);
    let result = evaluate_diffuse_recombination_with_illumination(
        &[DiffuseRecombinationCell { state: &state, source: &source }],
        area, &[beam(20.0, 1.0, 1.0)], controls(1e-7), || false,
    ).unwrap();
    close(result.cells()[0].balance().incident_absorbed_photons_per_s, SIGMA * NA, 4e-14);
    assert!(result.cells()[0].balance().net_electron_power_w > 0.0);
    assert_eq!(result.interfaces()[0].right_photons_per_s, result.interfaces()[1].right_photons_per_s);
}

#[test]
fn g4_empty_columns_and_source_work_cancellation_admission() {
    let beams = [beam(20.0, 3e20, 0.5), beam(40.0, 2e20, -0.25)];
    let empty = evaluate_diffuse_recombination_with_illumination(
        &[], Area::new(2.0), &beams, controls(1e-7), || false,
    ).unwrap();
    assert_eq!(empty.balance().incident_photons_per_s, empty.balance().escaped_photons_per_s);
    assert_eq!(empty.balance().incident_power_w, empty.balance().escaped_power_w);
    assert_eq!(empty.balance().absorbed_power_w, 0.0);
    let (state, source) = parcel(0.0, 10000.0, 1.0, 0.2);
    let cells = [DiffuseRecombinationCell { state: &state, source: &source }];
    for bad in [beam(20.0, 1.0, 0.0), beam(20.0, 1.0, 1.1), beam(20.0, -1.0, 1.0), beam(f64::NAN, 1.0, 1.0)] {
        assert!(matches!(evaluate_diffuse_recombination_with_illumination(
            &cells, Area::new(1.0), &[bad], controls(1e-7), || false,
        ), Err(DiffuseRecombinationError::Invalid(_))));
    }
    let outside = DiffuseRecombinationBeam {
        photons: PhotonNumberFluxV1 {
            photon_energy: Energy::new(2.0 * source.source().maximum_energy().value()),
            number_flux_per_m2_s: 1.0,
        },
        direction_cosine: 1.0,
    };
    assert!(matches!(evaluate_diffuse_recombination_with_illumination(
        &cells, Area::new(1.0), &[outside], controls(1e-7), || false,
    ), Err(DiffuseRecombinationError::Absorption { index: 0, .. })));
    let mut limited = controls(1e-7);
    limited.maximum_work = 1;
    assert_eq!(evaluate_diffuse_recombination_with_illumination(
        &cells, Area::new(1.0), &beams, limited, || false,
    ).unwrap_err(), DiffuseRecombinationError::WorkBudget);
    assert_eq!(evaluate_diffuse_recombination_with_illumination(
        &cells, Area::new(1.0), &beams, controls(1e-7), || true,
    ).unwrap_err(), DiffuseRecombinationError::Cancelled);
}
