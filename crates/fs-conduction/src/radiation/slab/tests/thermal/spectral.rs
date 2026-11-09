use super::*;
use crate::radiation::spectral::planck::BlackbodyBand;

fn optical(depths: &[f64], index: f64, exterior: f64, high_k: f64) -> AbsorbingSlab {
    AbsorbingSlab::new_in_medium(
        depths
            .iter()
            .map(|&depth| {
                AbsorbingLayer::new(
                    1.0,
                    depth,
                    [100.0, high_k],
                    "authored finite-band absorption",
                )
                .unwrap()
                .with_refractive_index(index, "authored finite-band material index")
                .unwrap()
            })
            .collect(),
        exterior,
        "authored finite-band exterior index",
        SlabQuadrature::default(),
        || false,
    )
    .unwrap()
}

fn band(window: [f64; 2], optics: AbsorbingSlab, bath_k: [f64; 2]) -> SlabThermalBand {
    let window = BlackbodyBand::new(window[0], window[1]).unwrap();
    let index_squared = optics.exterior_refractive_index().powi(2);
    SlabThermalBand {
        window,
        optics,
        irradiation_w_m2: bath_k.map(|t| index_squared * window.evaluate(t).unwrap().power_w_m2),
    }
}

// Independent physical-wavelength Simpson integral of the dimensional Planck
// radiance. It does not use production cumulative moments or band differences.
fn planck_oracle(window: BlackbodyBand, temperature: f64) -> f64 {
    let h: f64 = 6.626_070_15e-34;
    let c: f64 = 299_792_458.0;
    let k: f64 = 1.380_649e-23;
    let sigma_si = 2.0 * core::f64::consts::PI.powi(5) * k.powi(4) / (15.0 * h.powi(3) * c.powi(2));
    let lo = window.lower_wavelength_m().ln();
    let width = window.upper_wavelength_m().ln() - lo;
    let panels = 4096;
    let mut sum = 0.0;
    for i in 0..=panels {
        let wavelength = (lo + width * i as f64 / panels as f64).exp();
        let x = h * c / (k * wavelength * temperature);
        let radiance = 2.0 * core::f64::consts::PI * h * c * c / (wavelength.powi(4) * x.exp_m1());
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * radiance;
    }
    sum * width / (3.0 * panels as f64) * SIGMA / sigma_si
}

#[test]
fn g1_spectral_absorption_uses_physical_planck_windows_and_original_latent_storage() {
    let opaque = optical(&[1000.0], 1.0, 1.0, 3000.0);
    let transparent = optical(&[0.0], 1.0, 1.0, 3000.0);
    let law = Caloric {
        capacity: 1000.0,
        datum: 1e9,
        latent: 20_000.0,
    };
    let old = state(&law, &opaque, &[law.energy(500.0) + 10_000.0], 0.01, 0.0, 0);
    let mut active = band([3e-6, 6e-6], opaque.clone(), [500.0; 2]);
    active.irradiation_w_m2 = active.irradiation_w_m2.map(|power| 2.0 * power);
    let mut clear = band([8e-6, 14e-6], transparent, [500.0; 2]);
    clear.irradiation_w_m2 = [1e6; 2];
    let expected_power = 0.02 * planck_oracle(active.window, 500.0);
    let source =
        SpectralSlabThermal::new(vec![active, clear], 0.01, &bindings(&law, 1), &old, || {
            Ok(())
        })
        .unwrap();
    let transport = source.endpoint_transport(&old, || Ok(())).unwrap();
    close(
        transport.heat.bulk.sources_w[&CellId(0)],
        expected_power,
        2e-10,
    );
    assert!(expected_power < 0.01 * exitance(500.0));
    let (next, _, audit) = old
        .advance_observed(&source, 0.05, config(), || Ok(()))
        .unwrap();
    close(audit.external_heat_j, 0.05 * expected_power, 1e-10);
    close(audit.stored_energy_change_j, 0.05 * expected_power, 1e-7);
    assert_eq!(next.temperature_k(CellId(0)).unwrap(), 500.0);
    assert_eq!(
        next.cells()[&CellId(0)].geometry,
        old.cells()[&CellId(0)].geometry
    );
    assert_eq!(
        source.bands()[0].optics.layers().as_ptr(),
        opaque.layers().as_ptr()
    );
    assert!(
        source.bands()[0].window.upper_wavelength_m()
            < source.bands()[1].window.lower_wavelength_m()
    );
}

fn two_bands(bath: [f64; 2]) -> Vec<SlabThermalBand> {
    vec![
        band(
            [1e-6, 4e-6],
            optical(&[1000.0, 1000.0], 1.5, 1.0, 1800.0),
            bath,
        ),
        band(
            [8e-6, 14e-6],
            optical(&[1000.0, 1000.0], 2.0, 2.0, 3000.0),
            bath,
        ),
    ]
}

#[test]
fn g3_spectral_pair_coalescing_conserves_native_heat_and_combines_true_band_stiffness() {
    let bands = two_bands([900.0, 300.0]);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = state(
        &law,
        &bands[0].optics,
        &[law.energy(900.0), law.energy(300.0)],
        0.01,
        0.5,
        1,
    );
    let expected_power: f64 = bands
        .iter()
        .map(|b| {
            0.01 * b.optics.refractive_index().unwrap().powi(2)
                * (planck_oracle(b.window, 900.0) - planck_oracle(b.window, 300.0))
        })
        .sum();
    let mut source =
        SpectralSlabThermal::new(bands.clone(), 0.01, &bindings(&law, 2), &old, || Ok(())).unwrap();
    let transport = source.endpoint_transport(&old, || Ok(())).unwrap();
    assert_eq!(transport.transfers.len(), 1); // two bands, ONE physical pair/budget slot
    assert!(transport.heat.bulk.sources_w.is_empty());
    close(transport.transfers[0].power_w, expected_power, 2e-9);
    let (next, _, audit) = old
        .advance_observed(&source, 1e-4, config(), || Ok(()))
        .unwrap();
    let energy = 1e-4 * (expected_power + 0.5 * 600.0);
    close(
        next.cells()[&CellId(0)].specific_h_j_kg,
        law.energy(900.0) - energy,
        2e-10,
    );
    close(
        next.cells()[&CellId(1)].specific_h_j_kg,
        law.energy(300.0) + energy,
        2e-10,
    );
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(next.links(), old.links());
    let mut decay = 0.0;
    let mut wrong_four_e_over_t = 0.0;
    for b in &bands {
        let high = b.optics.layers()[0].temperature_support_k()[1];
        let emission = b.window.evaluate(high).unwrap();
        let row = b.optics.exterior_refractive_index().powi(2)
            * b.optics.face_coupling()[0].iter().sum::<f64>()
            + b.optics.layer_coupling()[1];
        decay += 0.01 * row * emission.temperature_derivative_w_m2_k / law.capacity;
        wrong_four_e_over_t += 0.01 * row * 4.0 * emission.power_w_m2 / high / law.capacity;
    }
    assert!((decay / wrong_four_e_over_t - 1.0).abs() > 0.01);
    let proposed = 0.25;
    close(
        source
            .limit_timestep_s(&old, proposed, &mut || Ok(()))
            .unwrap(),
        proposed / (1.0 + proposed * decay / 0.8),
        1e-15,
    );
}

#[test]
fn g1_spectral_transport_drives_original_backward_euler_endpoint() {
    let bands = two_bands([550.0, 350.0]);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let power: f64 = bands
        .iter()
        .map(|b| {
            0.01 * b.optics.refractive_index().unwrap().powi(2)
                * b.window.difference(550.0, 350.0).unwrap()
        })
        .sum();
    let energy = power + 0.5 * 200.0;
    let old = state(
        &law,
        &bands[0].optics,
        &[law.energy(550.0) + energy, law.energy(350.0) - energy],
        0.01,
        0.5,
        1,
    );
    let source =
        SpectralSlabThermal::new(bands, 0.01, &bindings(&law, 2), &old, || Ok(())).unwrap();
    let step = BackwardEulerStep::new(&old, 1.0, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_transport(
            &[-energy, energy],
            |candidate| source.endpoint_transport(candidate, || Ok(())),
            || Ok(()),
        )
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-8);
    let (next, (), audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
    close(next.temperature_k(CellId(0)).unwrap(), 550.0, 1e-12);
    close(next.temperature_k(CellId(1)).unwrap(), 350.0, 1e-12);
    assert_eq!(audit.external_heat_j, 0.0);
    assert_eq!(next.links(), old.links());
}

#[test]
fn g3_spectral_window_subdivision_preserves_heat_and_exact_equilibrium() {
    let optics = optical(&[0.7, 1.2], 1.5, 1.3, 3000.0);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = state(&law, &optics, &[700_000.0, 900_000.0], 0.01, 0.0, 1);
    let broad = SpectralSlabThermal::new(
        vec![band([1e-6, 20e-6], optics.clone(), [300.0, 400.0])],
        0.01,
        &bindings(&law, 2),
        &old,
        || Ok(()),
    )
    .unwrap();
    let windows = [[1e-6, 3e-6], [3e-6, 8e-6], [8e-6, 20e-6]];
    let split = SpectralSlabThermal::new(
        windows
            .map(|w| band(w, optics.clone(), [300.0, 400.0]))
            .to_vec(),
        0.01,
        &bindings(&law, 2),
        &old,
        || Ok(()),
    )
    .unwrap();
    let a = broad.endpoint_transport(&old, || Ok(())).unwrap();
    let b = split.endpoint_transport(&old, || Ok(())).unwrap();
    for id in [CellId(0), CellId(1)] {
        close(
            a.heat.bulk.sources_w[&id],
            b.heat.bulk.sources_w[&id],
            2e-10,
        );
    }
    close(a.transfers[0].power_w, b.transfers[0].power_w, 2e-10);
    let same = state(&law, &optics, &[500_000.0; 2], 0.01, 0.0, 1);
    let equilibrium = SpectralSlabThermal::new(
        windows
            .map(|w| band(w, optics.clone(), [500.0; 2]))
            .to_vec(),
        0.01,
        &bindings(&law, 2),
        &same,
        || Ok(()),
    )
    .unwrap();
    let transport = equilibrium.endpoint_transport(&same, || Ok(())).unwrap();
    assert!(transport.heat.bulk.sources_w.is_empty());
    assert!(transport.transfers.is_empty());
    let (next, _, audit) = same
        .advance_observed(&equilibrium, 0.2, config(), || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), same.cells());
    assert_eq!(audit.external_heat_j, 0.0);
}

#[test]
fn g0_spectral_near_equilibrium_retains_positive_heat_in_adjacent_wavelength_windows() {
    let optics = optical(&[1000.0, 1000.0], 1.5, 1.0, 3000.0);
    let law = Caloric {
        capacity: 1.0,
        datum: 0.0,
        latent: 0.0,
    };
    let cold = f64::from_bits(500.0_f64.to_bits() - 1);
    let old = state(&law, &optics, &[500.0, cold], 0.01, 0.0, 1);
    let lower = 30e-6_f64;
    for window in [[1e-6, 20e-6], [lower, f64::from_bits(lower.to_bits() + 1)]] {
        let input = band(window, optics.clone(), [500.0, cold]);
        let expected = 0.01
            * 2.25
            * input
                .window
                .evaluate(500.0)
                .unwrap()
                .temperature_derivative_w_m2_k
            * (500.0 - cold);
        let source =
            SpectralSlabThermal::new(vec![input], 0.01, &bindings(&law, 2), &old, || Ok(()))
                .unwrap();
        let transport = source.endpoint_transport(&old, || Ok(())).unwrap();
        assert!(transport.heat.bulk.sources_w.is_empty());
        let heat = transport.transfers[0].power_w;
        assert!(heat > 0.0);
        assert!(
            (heat / expected - 1.0).abs() < 2e-13,
            "{heat:e} != {expected:e}"
        );
    }
}

#[test]
fn g2_spectral_overlap_geometry_support_union_budget_and_cancellation_refuse() {
    let optics = optical(&[1.0, 1.0], 1.0, 1.0, 3000.0);
    let law = Caloric {
        capacity: 1000.0,
        datum: 0.0,
        latent: 0.0,
    };
    let old = state(&law, &optics, &[300_000.0, 400_000.0], 0.01, 0.0, 1);
    let first = band([1e-6, 4e-6], optics.clone(), [300.0; 2]);
    let make = |bands| SpectralSlabThermal::new(bands, 0.01, &bindings(&law, 2), &old, || Ok(()));
    assert!(make(vec![]).is_err());
    assert!(make(vec![first.clone(); MAX_SLAB_SPECTRAL_BANDS + 1]).is_err());
    assert!(
        make(vec![
            first.clone(),
            band([3e-6, 6e-6], optics.clone(), [300.0; 2])
        ])
        .is_err()
    );
    let different = AbsorbingSlab::new(
        vec![AbsorbingLayer::new(2.0, 1.0, [100.0, 3000.0], "different thickness").unwrap(); 2],
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    assert!(
        make(vec![
            first.clone(),
            band([4e-6, 6e-6], different, [300.0; 2])
        ])
        .is_err()
    );
    let narrow_support = optical(&[1.0, 1.0], 1.0, 1.0, 350.0);
    assert!(matches!(
        make(vec![
            first.clone(),
            band([4e-6, 6e-6], narrow_support, [300.0; 2])
        ]),
        Err(DepositionError::MaterialDomain)
    ));
    let left = optical(&[1000.0, 1000.0, 0.0], 1.0, 1.0, 3000.0);
    let right = optical(&[0.0, 1000.0, 1000.0], 1.0, 1.0, 3000.0);
    let three = state(
        &law,
        &left,
        &[300_000.0, 400_000.0, 500_000.0],
        0.01,
        0.0,
        1,
    );
    assert!(matches!(
        SpectralSlabThermal::new(
            vec![
                band([1e-6, 4e-6], left, [300.0; 2]),
                band([4e-6, 8e-6], right, [300.0; 2])
            ],
            0.01,
            &bindings(&law, 3),
            &three,
            || Ok(())
        ),
        Err(DepositionError::Budget(_))
    ));
    let source = make(vec![first]).unwrap();
    let saved = old.cells().clone();
    let mut calls = 0;
    assert!(matches!(
        old.advance_observed(&source, 0.1, config(), || {
            calls += 1;
            if calls == 8 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(old.cells(), &saved);
    assert_eq!(old.time_s(), 0.0);
    assert!(source.endpoint_transport(&old, || Ok(())).is_ok());
}
