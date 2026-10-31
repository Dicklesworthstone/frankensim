use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1e-20),
        "{actual:e} != {expected:e}"
    );
}

fn membrane(lp: f64, species: &[(f64, f64)]) -> MembraneTransport {
    MembraneTransport::new(
        HydraulicPermeance::new(lp).unwrap(),
        species
            .iter()
            .map(|&(sigma, p)| {
                (
                    ReflectionCoefficient::new(sigma).unwrap(),
                    SolutePermeance::new(p).unwrap(),
                )
            })
            .collect(),
    )
    .unwrap()
}

#[test]
fn g1_millimolar_units_and_independent_particles() {
    close(
        ideal_dilute_osmotic_pressure_pa(298.15, &[100.0]).unwrap(),
        247_895.702_960_238_8,
    );
    // Explicit independent Na and Cl particles; no second dissociation factor.
    close(
        ideal_dilute_osmotic_pressure_pa(298.15, &[100.0, 100.0]).unwrap(),
        495_791.405_920_477_6,
    );
    assert_eq!(ideal_dilute_osmotic_pressure_pa(298.15, &[]).unwrap(), 0.0);
}

#[test]
fn g1_activity_pressure_and_dilute_roundtrip() {
    close(
        osmotic_pressure_from_water_activity_pa(298.15, 0.99, 18e-6).unwrap(),
        1_384_130.595_216_796_5,
    );
    for pi in [0.0, 1e3, 2e5, 1e8] {
        let activity = water_activity_from_osmotic_pressure(298.15, pi, 18e-6).unwrap();
        let reconstructed =
            osmotic_pressure_from_water_activity_pa(298.15, activity, 18e-6).unwrap();
        assert!((reconstructed - pi).abs() < 3e-8 + 1e-12 * pi);
    }
}

#[test]
fn g0_solution_dilution_conserves_amounts_and_refreshes_activity() {
    let solution = IdealDiluteSolution::new(298.15, 18e-6, 1e-3, vec![0.1, 0.05]).unwrap();
    let diluted = solution.with_solvent_volume_m3(2e-3).unwrap();
    assert_eq!(solution.solute_moles(), diluted.solute_moles());
    close(
        diluted.osmotic_pressure_pa().unwrap(),
        0.5 * solution.osmotic_pressure_pa().unwrap(),
    );
    assert!(diluted.water_activity().unwrap() > solution.water_activity().unwrap());
    let mut c = [0.0; 2];
    diluted.concentrations_into(&mut c).unwrap();
    assert_eq!(c, [50.0, 25.0]);
    let mut candidate = solution.clone();
    assert!(candidate.replace_inventories(1e-3, &[0.1, -1.0]).is_err());
    assert_eq!(candidate, solution);
    assert!(candidate.replace_inventories(0.0, &[0.1, 0.0]).is_err());
    assert_eq!(candidate, solution);
    candidate.replace_inventories(2e-3, &[0.02, 0.03]).unwrap();
    assert_eq!(candidate.solute_moles(), &[0.02, 0.03]);
}

#[test]
fn g1_forward_and_reverse_osmosis_have_physical_signs() {
    let law = membrane(1e-11, &[(1.0, 0.0)]);
    let mut solute = [9.0];
    let pi = ideal_dilute_osmotic_pressure_pa(300.0, &[100.0]).unwrap();
    let forward = law
        .flux_into(300.0, 0.0, 0.0, &[0.0], &[100.0], &mut solute)
        .unwrap();
    close(forward.volume_flux_m_s, 1e-11 * pi);
    assert_eq!(solute, [0.0]);
    let balanced = law
        .flux_into(300.0, 0.0, pi, &[0.0], &[100.0], &mut solute)
        .unwrap();
    assert_eq!(balanced.volume_flux_m_s, 0.0);
    let reverse = law
        .flux_into(300.0, 0.0, 2.0 * pi, &[0.0], &[100.0], &mut solute)
        .unwrap();
    close(reverse.volume_flux_m_s, -forward.volume_flux_m_s);
}

#[test]
fn g0_nonselective_interface_has_no_osmotic_force() {
    let law = membrane(1e-11, &[(0.0, 1e-7)]);
    let mut solute = [0.0];
    let flow = law
        .flux_into(300.0, 0.0, 0.0, &[10.0], &[0.0], &mut solute)
        .unwrap();
    assert_eq!(flow.volume_flux_m_s, 0.0);
    close(solute[0], 1e-6);
}

#[test]
fn g1_selective_mixture_uses_each_reflection_coefficient() {
    let law = membrane(2e-12, &[(1.0, 0.0), (0.0, 1e-6), (0.5, 1e-7)]);
    let mut out = [0.0; 3];
    let flow = law
        .flux_into(
            300.0,
            0.0,
            0.0,
            &[0.0, 400.0, 0.0],
            &[20.0, 0.0, 60.0],
            &mut out,
        )
        .unwrap();
    close(flow.volume_flux_m_s, 2e-12 * R * 300.0 * 50.0);
    assert_eq!(
        flow,
        law.water_flux(300.0, 0.0, 0.0, &[0.0, 400.0, 0.0], &[20.0, 0.0, 60.0])
            .unwrap()
    );
    assert_eq!(out[0], 0.0);
    assert!(out[1] > 0.0);
}

#[test]
fn g1_exponential_fit_matches_independent_closed_form() {
    let law = membrane(1e-8, &[(0.0, 2e-6)]);
    for u in [-1e-4, -1e-6, -1e-13, 0.0, 1e-13, 1e-6, 1e-4] {
        let mut out = [0.0];
        law.flux_into(300.0, u / 1e-8, 0.0, &[20.0], &[10.0], &mut out)
            .unwrap();
        let expected = if u == 0.0 {
            2e-6 * 10.0
        } else if u > 0.0 {
            u * (20.0 - (-u / 2e-6).exp() * 10.0) / -(-u / 2e-6).exp_m1()
        } else {
            u * (10.0 - (u / 2e-6).exp() * 20.0) / -(u / 2e-6).exp_m1()
        };
        close(out[0], expected);
        let (left, right) = law.solute_stream_coefficients(0, u).unwrap();
        assert!(left >= 0.0 && right >= 0.0);
        close(out[0], left * 20.0 - right * 10.0);
    }
}

#[test]
fn g0_zero_diffusion_uses_true_donor_and_perfect_rejection_stops_drag() {
    let law = membrane(1e-8, &[(0.0, 0.0), (1.0, 0.0)]);
    let mut out = [0.0; 2];
    let forward = law
        .flux_into(300.0, 100.0, 0.0, &[5.0, 1.0], &[0.0, 1.0], &mut out)
        .unwrap();
    close(out[0], 5.0 * forward.volume_flux_m_s);
    assert_eq!(out[1], 0.0);
    law.flux_into(300.0, -100.0, 0.0, &[5.0, 1.0], &[0.0, 1.0], &mut out)
        .unwrap();
    assert_eq!(out[0], 0.0);
    assert_eq!(out[1], 0.0);
}

#[test]
fn g3_orientation_reversal_reverses_every_flux() {
    let law = membrane(1e-11, &[(0.8, 1e-7), (0.1, 5e-6)]);
    let mut forward = [0.0; 2];
    let mut reverse = [0.0; 2];
    let a = law
        .flux_into(310.0, 1e5, 2e5, &[20.0, 80.0], &[60.0, 10.0], &mut forward)
        .unwrap();
    let b = law
        .flux_into(310.0, 2e5, 1e5, &[60.0, 10.0], &[20.0, 80.0], &mut reverse)
        .unwrap();
    assert_eq!(a.volume_flux_m_s, -b.volume_flux_m_s);
    assert_eq!(forward, reverse.map(|f| -f));
}

#[test]
fn g3_positive_concentration_dissipation_is_nonnegative() {
    let law = membrane(1e-11, &[(0.8, 1e-7), (0.1, 5e-6)]);
    for pressure in [-1e7, -1e5, 0.0, 1e5, 1e7] {
        for (a, b) in [(10.0, 50.0), (50.0, 10.0), (50.0, 50.0)] {
            let mut solute = [0.0; 2];
            let left = [a, 20.0];
            let right = [b, 60.0];
            let flux = law
                .flux_into(300.0, pressure, 0.0, &left, &right, &mut solute)
                .unwrap();
            let chemical = R
                * 300.0
                * solute
                    .iter()
                    .enumerate()
                    .map(|(i, js)| js * (left[i] / right[i]).ln())
                    .sum::<f64>();
            let solvent = flux.volume_flux_m_s
                * (pressure - R * 300.0 * ((left[0] - right[0]) + (left[1] - right[1])));
            assert!(solvent + chemical >= -1e-12);
        }
    }
}

#[test]
fn g0_invalid_or_overflowing_candidate_does_not_publish_flux() {
    let law = membrane(1e-11, &[(0.8, 1e-7)]);
    let mut out = [123.0];
    assert!(
        law.flux_into(300.0, 0.0, 0.0, &[f64::NAN], &[10.0], &mut out)
            .is_err()
    );
    assert_eq!(out, [123.0]);
    assert!(
        law.flux_into(300.0, 0.0, 0.0, &[f64::MAX], &[10.0], &mut out)
            .is_err()
    );
    assert_eq!(out, [123.0]);
    assert!(ReflectionCoefficient::new(1.1).is_err());
    assert!(SolutePermeance::new(-1.0).is_err());
    assert!(HydraulicPermeance::new(f64::INFINITY).is_err());
    assert!(osmotic_pressure_from_water_activity_pa(300.0, 0.0, 18e-6).is_err());
    assert!(water_activity_from_osmotic_pressure(300.0, 1e300, 18e-6).is_err());
    assert!(
        IdealDiluteSolution::new(300.0, 18e-6, 1.0, vec![0.0; MAX_OSMOTIC_SOLUTES + 1]).is_err()
    );
}

#[test]
fn g1_darcy_crosswalk_consumes_actual_fluid_viscosity_and_thickness() {
    let fluid =
        crate::fluid_transport::FluidTransportState::try_new(1000.0, 4200.0, 0.6, 1e-3).unwrap();
    let lp = HydraulicPermeance::from_darcy_slab(1e-18, 1e-4, &fluid).unwrap();
    close(lp.value(), 1e-11);
    let thicker = HydraulicPermeance::from_darcy_slab(1e-18, 2e-4, &fluid).unwrap();
    close(thicker.value(), 0.5 * lp.value());
}
