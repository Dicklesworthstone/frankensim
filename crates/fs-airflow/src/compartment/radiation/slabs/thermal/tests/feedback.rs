use super::*;

fn one_cell_surfaces(baths: [f64; 2]) -> [SlabThermalSurface; 4] {
    let mut result = bindings(baths);
    result[2] = SlabThermalSurface::Stored(CellId(1));
    result
}

fn response(
    thickness: f64,
    absorption: f64,
    index: f64,
) -> Result<AbsorbingLayer, ConductionError> {
    AbsorbingLayer::new(
        thickness,
        absorption,
        [250.0, 2200.0],
        "declared state-dependent absorption",
    )?
    .with_refractive_index(index, "declared state-dependent real refractive index")
}

// Independent single-layer ray sum and composite Simpson hemispherical integral.
// The two equal interfaces give A=(1-R)(1-t)/(1-R*t), separately for s and p.
fn absorptance(tau: f64, n: f64) -> f64 {
    const PANELS: usize = 4096;
    let mut total = 0.0;
    for j in 1..=PANELS {
        let mu = j as f64 / PANELS as f64;
        let internal = (1.0 - (1.0 - mu * mu) / (n * n)).sqrt();
        let loss = -(-tau / internal).exp_m1();
        let transmission = (-tau / internal).exp();
        let rs = ((mu - n * internal) / (mu + n * internal)).powi(2);
        let rp = ((n * mu - internal) / (n * mu + internal)).powi(2);
        let a = 0.5
            * ((1.0 - rs) * loss / (1.0 - rs * transmission)
                + (1.0 - rp) * loss / (1.0 - rp * transmission));
        let weight = if j == PANELS {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        total += weight * 2.0 * mu * a / (3.0 * PANELS as f64);
    }
    total
}

fn bisect(mut low: f64, mut high: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..64 {
        let mid = 0.5 * (low + high);
        if residual(mid) < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    0.5 * (low + high)
}

#[test]
fn implicit_refractive_feedback_matches_independent_ray_and_energy_root() {
    let initial = state(&[(0.2, 0.0, 1)], &[], 0.0);
    let enclosure = enclosure();
    let slabs = [slab(&[0.12])];
    let absorption = absorption();
    let optics = SlabThermalOptics::Gray {
        absorption: &absorption,
        slabs: &slabs,
        solve: RadiationSolve::default(),
    };
    let mut sampled = Vec::new();
    let result = enclosure
        .candidate_slab_thermal_with_feedback(
            &initial,
            Time::new(3.0),
            &one_cell_surfaces([700.0, 300.0]),
            &[vec![CellId(1)]],
            [Temperature::new(300.0); 2],
            optics,
            |query| {
                assert_eq!(query.cell.id, CellId(1));
                assert_eq!(query.band, None);
                assert_eq!(query.wavelengths, None);
                assert_eq!(query.reference, &slabs[0].optics.layers()[0]);
                let t = query.temperature.value();
                sampled.push(t);
                response(0.01, 12.0 + 0.3 * (t - 300.0), 1.3 + 0.001 * (t - 300.0))
            },
            black,
            SlabThermalSolve::default(),
            &mut enclosure.slab_thermal_workspace(&initial).unwrap(),
            || false,
        )
        .unwrap();
    let expected = bisect(300.0, 700.0, |t| {
        let a = absorptance(0.12 + 0.003 * (t - 300.0), 1.3 + 0.001 * (t - 300.0));
        200.0 * (t - 300.0) - 3.0 * a * (blackbody(700.0) + blackbody(300.0) - 2.0 * blackbody(t))
    });
    let t = result.thermal.temperature_k(CellId(1)).unwrap();
    near(t, expected, 2e-7);
    assert!(sampled.len() > 3 && sampled.iter().any(|x| *x > 310.0));
    let WallRadiation::Gray(report) = &result.report.radiation else {
        panic!("gray report");
    };
    near(
        report.slabs[0].optics.layers()[0].refractive_index(),
        1.3 + 0.001 * (t - 300.0),
        1e-14,
    );
    near(
        report.slabs[0].optics.layers()[0].absorption_m_inv(),
        12.0 + 0.3 * (t - 300.0),
        1e-12,
    );
    assert!(result.report.combined_energy_residual_j.abs() < 1e-5);
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, 0.0);
    assert_eq!(slabs[0].optics.layers()[0].refractive_index(), 1.0);
}

#[test]
fn latent_plateau_uses_native_enthalpy_for_optical_feedback() {
    let initial = state(&[(0.1, 20_000.0, 2)], &[], 0.0);
    let enclosure = enclosure();
    let slabs = [slab(&[0.4])];
    let absorption = absorption();
    let mut enthalpies = Vec::new();
    let result = enclosure
        .candidate_slab_thermal_with_feedback(
            &initial,
            Time::new(0.2),
            &one_cell_surfaces([500.0, 300.0]),
            &[vec![CellId(1)]],
            [Temperature::new(300.0); 2],
            SlabThermalOptics::Gray {
                absorption: &absorption,
                slabs: &slabs,
                solve: RadiationSolve::default(),
            },
            |query| {
                assert_eq!(query.temperature.value(), 310.0);
                assert_eq!(query.cell.material_id, MaterialId(2));
                let h = query.cell.specific_h_j_kg;
                enthalpies.push(h);
                response(0.01, 20.0 + 0.001 * h, 1.5 + 1e-6 * h)
            },
            black,
            SlabThermalSolve::default(),
            &mut enclosure.slab_thermal_workspace(&initial).unwrap(),
            || false,
        )
        .unwrap();
    let expected = bisect(20_000.0, 100_000.0, |h| {
        0.1 * (h - 20_000.0)
            - 0.2
                * absorptance(0.2 + 1e-5 * h, 1.5 + 1e-6 * h)
                * (blackbody(500.0) + blackbody(300.0) - 2.0 * blackbody(310.0))
    });
    near(
        result.thermal.cells()[&CellId(1)].specific_h_j_kg,
        expected,
        2e-6,
    );
    assert_eq!(result.thermal.temperature_k(CellId(1)).unwrap(), 310.0);
    assert!(enthalpies.iter().any(|h| *h > 20_100.0));
}

#[test]
fn spectral_feedback_keeps_real_bands_and_returns_actual_endpoint_sources() {
    let initial = state(&[(0.2, 0.0, 1)], &[], 0.0);
    let enclosure = enclosure();
    let ranges = [[1e-7, 4e-6], [4e-6, 1e-3]];
    let optics = SpectralOptics::new_with_slab_bands(
        ranges
            .map(|range| SpectralBand {
                wavelengths: range.map(Length::new),
                absorption: absorption(),
            })
            .to_vec(),
        [None; 2],
        (0..2)
            .map(|_| {
                (0..2)
                    .map(|_| {
                        BandSurfaceEmissivity::new(
                            1.0,
                            [250.0, 2200.0].map(Temperature::new),
                            "declared black spectral bath",
                        )
                        .unwrap()
                    })
                    .collect()
            })
            .collect(),
        vec![vec![slab(&[0.12])], vec![slab(&[0.4])]],
    )
    .unwrap();
    let mut sampled = [0; 2];
    let result = enclosure
        .candidate_slab_thermal_with_feedback(
            &initial,
            Time::new(0.5),
            &one_cell_surfaces([700.0, 300.0]),
            &[vec![CellId(1)]],
            [Temperature::new(300.0); 2],
            SlabThermalOptics::Spectral {
                optics: &optics,
                solve: SpectralSolve::default(),
            },
            |query| {
                let band = query.band.unwrap();
                assert_eq!(query.wavelengths.unwrap().map(Length::value), ranges[band]);
                sampled[band] += 1;
                response(
                    0.01,
                    10.0 + 30.0 * band as f64 + 0.1 * query.temperature.value(),
                    1.4 + 0.2 * band as f64,
                )
            },
            |_, _| panic!("spectral bands must not query gray cards"),
            SlabThermalSolve::default(),
            &mut enclosure.slab_thermal_workspace(&initial).unwrap(),
            || false,
        )
        .unwrap();
    let t = result.thermal.temperature_k(CellId(1)).unwrap();
    let WallRadiation::Spectral(report) = &result.report.radiation else {
        panic!("spectral report");
    };
    let mut expected_power = 0.0;
    for (band, range) in ranges.iter().enumerate() {
        let layer = &report.optics.slab_bands().unwrap()[band][0].optics.layers()[0];
        let kappa = 10.0 + 30.0 * band as f64 + 0.1 * t;
        near(layer.absorption_m_inv(), kappa, 1e-12);
        let black = BlackbodyBand::new(range[0], range[1]).unwrap();
        let b = [700.0, 300.0, t].map(|x| black.evaluate(x).unwrap().power_w_m2);
        expected_power +=
            absorptance(0.01 * kappa, 1.4 + 0.2 * band as f64) * (b[0] + b[1] - 2.0 * b[2]);
    }
    near(200.0 * (t - 300.0), 0.5 * expected_power, 2e-6);
    assert_eq!(sampled[0], sampled[1]);
    assert!(sampled[0] > 1);
    assert_eq!(report.factorizations, 2);
    assert_eq!(report.wavelength_evaluations, 0);
    assert_eq!(report.optics.bands(), optics.bands());
    assert_eq!(
        optics.slab_bands().unwrap()[0][0].optics.layers()[0].absorption_m_inv(),
        12.0
    );
}

#[test]
fn feedback_refuses_geometry_domain_and_cancellation_without_advancing_original_state() {
    let initial = state(&[(0.2, 0.0, 1)], &[], 0.0);
    let original_cells = initial.cells().clone();
    let enclosure = enclosure();
    let slabs = [slab(&[0.4])];
    let absorption = absorption();
    let optics = SlabThermalOptics::Gray {
        absorption: &absorption,
        slabs: &slabs,
        solve: RadiationSolve::default(),
    };
    let mut work = enclosure.slab_thermal_workspace(&initial).unwrap();
    for mode in 0..4 {
        let result = enclosure.candidate_slab_thermal_with_feedback(
            &initial,
            Time::new(0.2),
            &one_cell_surfaces([500.0, 300.0]),
            &[vec![CellId(1)]],
            [Temperature::new(300.0); 2],
            optics,
            |query| match mode {
                0 => response(0.02, 40.0, 1.5),
                1 => AbsorbingLayer::new(0.01, 40.0, [400.0, 500.0], "out-of-domain fixture"),
                3 => Err(ConductionError::Cancelled {
                    stage: "layer optical source",
                    at: 0,
                }),
                _ => Ok(query.reference.clone()),
            },
            black,
            SlabThermalSolve::default(),
            &mut work,
            || mode == 2,
        );
        assert!(matches!(
            result,
            Err(SlabThermalError::Invalid(_) | SlabThermalError::Cancelled)
        ));
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.cells(), &original_cells);
    }
    let identity = enclosure
        .candidate_slab_thermal_with_feedback(
            &initial,
            Time::new(0.2),
            &one_cell_surfaces([500.0, 300.0]),
            &[vec![CellId(1)]],
            [Temperature::new(300.0); 2],
            optics,
            |query| Ok(query.reference.clone()),
            black,
            SlabThermalSolve::default(),
            &mut work,
            || false,
        )
        .unwrap();
    let fixed = enclosure
        .candidate_slab_thermal(
            &initial,
            Time::new(0.2),
            &one_cell_surfaces([500.0, 300.0]),
            &[vec![CellId(1)]],
            [Temperature::new(300.0); 2],
            optics,
            black,
            SlabThermalSolve::default(),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(identity.thermal.cells(), fixed.thermal.cells());
    assert_eq!(identity.report, fixed.report);
}
