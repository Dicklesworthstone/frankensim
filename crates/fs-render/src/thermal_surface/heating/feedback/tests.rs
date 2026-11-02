use super::*;
use crate::conductor::{
    CauchyCoatingLayer, CoatingIorSample, CoherentAbsorbingCoating, CoherentDielectricCoating,
    CoherentPeriodicCoating,
};
use crate::dielectric::CauchyIor;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, Event, StepConfig,
    ThermalLaw, TimedEvent,
};

fn source(id: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([id; 32]), ConductorDataStatus::Representative).unwrap()
}

fn metal(eta: f64, k: f64) -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        source(11),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-14, "{actual} != {expected}");
}

#[test]
fn constitutive_midpoint_uses_complex_index_before_original_fresnel() {
    let left = metal(0.2, 3.0);
    let right = metal(4.2, 1.0);
    let law = SurfaceConductorLaw::try_temperature(source(12), &[(300.0, left), (700.0, right)], 2)
        .unwrap();
    let midpoint = law.at_temperature(500.0).unwrap();
    for sample in midpoint.samples() {
        close(sample.eta(), 2.2);
        close(sample.k(), 2.0);
    }
    let actual = midpoint.power_balance(550.0, 1.0, 1.0).unwrap();
    let independent = 4.0 * 2.2 / (3.2_f64.powi(2) + 4.0);
    close(actual.absorbed[0], independent);
    let intensity_mean = 0.5
        * (left.power_balance(550.0, 1.0, 1.0).unwrap().absorbed[0]
            + right.power_balance(550.0, 1.0, 1.0).unwrap().absorbed[0]);
    assert!((actual.absorbed[0] - intensity_mean).abs() > 0.1);
    assert_eq!(midpoint.source(), source(12));
    assert_eq!(law.at_temperature(300.0).unwrap().samples(), left.samples());
    assert_eq!(
        law.at_temperature(700.0).unwrap().samples(),
        right.samples()
    );
    assert!(law.at_specific_enthalpy(500.0, 500.0).is_err());
    for outside in [299.0, 701.0, f64::NAN, f64::INFINITY] {
        assert!(law.at_temperature(outside).is_err());
    }

    let constant =
        SurfaceConductorLaw::try_temperature(source(12), &[(300.0, left), (700.0, left)], 2)
            .unwrap();
    assert_eq!(
        constant.at_temperature(300.0).unwrap(),
        constant.at_temperature(541.7).unwrap()
    );
    assert!(Arc::ptr_eq(&law.rows, &law.clone().rows));
}

#[test]
fn subnormal_complex_index_interpolation_retains_small_positive_contributions() {
    let smallest = f64::from_bits(1);
    let next = f64::from_bits(2);
    let law = SurfaceConductorLaw::try_temperature(
        source(12),
        &[
            (300.0, metal(smallest, next)),
            (700.0, metal(next, smallest)),
        ],
        2,
    )
    .unwrap();
    for (temperature, eta_bits, k_bits) in [(400.0, 1, 2), (500.0, 2, 2), (600.0, 2, 1)] {
        for sample in law.at_temperature(temperature).unwrap().samples() {
            // Midpoint is 1.5 minimum-subnormal units, rounding to the even
            // significand 2. Independently rounded convex products lose it.
            assert_eq!(sample.eta().to_bits(), eta_bits);
            assert_eq!(sample.k().to_bits(), k_bits);
        }
    }
}

#[test]
fn unchanged_transparent_absorbing_and_periodic_films_keep_original_response() {
    let transparent = CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(1.45).unwrap(),
        90e-9,
        source(13),
    )
    .unwrap();
    let absorbing = CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 2.1, 0.2).unwrap()
        }),
        45e-9,
        source(14),
    )
    .unwrap();
    let periodic = CoherentPeriodicCoating::try_new(
        &[
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(2.1).unwrap(), 0.05, 45e-9)
                .unwrap(),
            CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.45).unwrap(), 0.0, 90e-9)
                .unwrap(),
        ],
        2,
        source(15),
    )
    .unwrap();
    for kind in 0..3 {
        let attach = |optics: ConductorOptics| match kind {
            0 => optics.with_coating(transparent),
            1 => optics.with_absorbing_coating(absorbing),
            _ => optics.with_periodic_coating(periodic),
        };
        let left = attach(metal(0.2, 3.0));
        let law = SurfaceConductorLaw::try_temperature(
            source(12),
            &[(300.0, left), (700.0, attach(metal(4.2, 1.0)))],
            2,
        )
        .unwrap();
        let resolved = law.at_temperature(500.0).unwrap();
        let expected = attach(metal(2.2, 2.0));
        assert_eq!(resolved.coating(), expected.coating());
        assert_eq!(resolved.absorbing_coating(), expected.absorbing_coating());
        assert_eq!(resolved.periodic_coating(), expected.periodic_coating());
        for wavelength in [480.0, 550.0, 730.0] {
            for cosine in [0.43, 1.0] {
                let actual = resolved.power_balance(wavelength, 1.0, cosine).unwrap();
                let reference = expected.power_balance(wavelength, 1.0, cosine).unwrap();
                for polarization in 0..2 {
                    close(
                        actual.reflected[polarization],
                        reference.reflected[polarization],
                    );
                    close(
                        actual.absorbed[polarization],
                        reference.absorbed[polarization],
                    );
                }
            }
        }
        assert!(
            SurfaceConductorLaw::try_temperature(
                source(12),
                &[(300.0, left), (700.0, metal(4.2, 1.0))],
                2,
            )
            .is_err()
        );
    }
}

#[test]
fn admission_retains_sources_and_domains_without_upgrading_or_extrapolating() {
    let rows = [(300.0, metal(0.2, 3.0)), (700.0, metal(4.2, 1.0))];
    let law = SurfaceConductorLaw::try_temperature(source(12), &rows, 2).unwrap();
    assert_eq!(law.row_count(), 2);
    assert_eq!(law.temperature_bounds_k(), (300.0, 700.0));
    assert_eq!(law.caloric_identity(), None);
    assert_eq!(law.specific_enthalpy_bounds_j_kg(), None);
    assert_eq!(
        law.identity(),
        SurfaceConductorLaw::try_temperature(source(12), &rows, 8)
            .unwrap()
            .identity()
    );
    let measured = ConductorSource::try_new(
        ContentHash([16; 32]),
        ConductorDataStatus::CallerAssertedMeasured,
    )
    .unwrap();
    assert!(SurfaceConductorLaw::try_temperature(measured, &rows, 2).is_err());
    let measured_rows = rows.map(|(t, optics)| {
        (
            t,
            ConductorOptics::try_new(optics.samples(), measured).unwrap(),
        )
    });
    assert!(SurfaceConductorLaw::try_temperature(measured, &measured_rows, 2).is_ok());
    let mut shifted_grid = rows;
    let mut samples = shifted_grid[1].1.samples();
    samples[1] = ConductorIorSample::try_new(440.0, 4.2, 1.0).unwrap();
    shifted_grid[1].1 = ConductorOptics::try_new(samples, source(11)).unwrap();
    assert!(SurfaceConductorLaw::try_temperature(source(12), &shifted_grid, 2).is_err());
    assert!(matches!(
        SurfaceConductorLaw::try_temperature(source(12), &rows, 1),
        Err(SurfaceHeatingError::Budget)
    ));
    for coordinates in [
        [700.0, 300.0],
        [300.0, 300.0],
        [0.0, 700.0],
        [300.0, f64::NAN],
    ] {
        let invalid = [(coordinates[0], rows[0].1), (coordinates[1], rows[1].1)];
        assert!(SurfaceConductorLaw::try_temperature(source(12), &invalid, 2).is_err());
    }
    for (caloric, bounds, coordinates) in [
        (ContentHash([0; 32]), (500.0, 500.0), [-100.0, 0.0]),
        (ContentHash([21; 32]), (0.0, 500.0), [-100.0, 0.0]),
        (ContentHash([21; 32]), (500.0, 499.0), [-100.0, 0.0]),
        (ContentHash([21; 32]), (500.0, 500.0), [-f64::MAX, f64::MAX]),
    ] {
        assert!(
            SurfaceConductorLaw::try_enthalpy(
                source(12),
                caloric,
                bounds,
                &[(coordinates[0], rows[0].1), (coordinates[1], rows[1].1)],
                2,
            )
            .is_err()
        );
    }
    let changed_source = [
        (300.0, rows[0].1),
        (
            700.0,
            ConductorOptics::try_new(rows[1].1.samples(), source(17)).unwrap(),
        ),
    ];
    for different in [
        SurfaceConductorLaw::try_temperature(source(13), &rows, 2).unwrap(),
        SurfaceConductorLaw::try_temperature(source(12), &changed_source, 2).unwrap(),
        SurfaceConductorLaw::try_temperature(
            source(12),
            &[(300.0, rows[0].1), (700.0, metal(4.3, 1.0))],
            2,
        )
        .unwrap(),
        SurfaceConductorLaw::try_enthalpy(
            source(12),
            ContentHash([21; 32]),
            (300.0, 700.0),
            &rows,
            2,
        )
        .unwrap(),
    ] {
        assert_ne!(law.identity(), different.identity());
    }
}

struct Plateau {
    identity: Option<[u8; 32]>,
    basis: ThermalStorageBasis,
}
impl ThermalLaw for Plateau {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.identity
    }
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(-200.0..=200.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(500.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
}

fn native(law: &Plateau, h: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(1),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(7),
                mass_kg: 1.0,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 500.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(7), law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 1,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn native_enthalpy_resolves_latent_states_only_with_original_caloric_datum() {
    let rows = [(-100.0, metal(0.2, 3.0)), (0.0, metal(4.2, 1.0))];
    let law = SurfaceConductorLaw::try_enthalpy(
        source(12),
        ContentHash([21; 32]),
        (500.0, 500.0),
        &rows,
        2,
    )
    .unwrap();
    assert_eq!(law.caloric_identity(), Some(ContentHash([21; 32])));
    assert_eq!(law.specific_enthalpy_bounds_j_kg(), Some((-100.0, 0.0)));
    assert_eq!(law.temperature_bounds_k(), (500.0, 500.0));
    assert!(law.at_temperature(500.0).is_err());
    for (h, t) in [
        (-101.0, 500.0),
        (1.0, 500.0),
        (-50.0, 501.0),
        (f64::NAN, 500.0),
    ] {
        assert!(law.at_specific_enthalpy(h, t).is_err());
    }
    let caloric = Plateau {
        identity: Some([21; 32]),
        basis: ThermalStorageBasis::Enthalpy,
    };
    let a = native(&caloric, -100.0);
    let b = native(&caloric, -50.0);
    assert_eq!(
        a.temperature_k(CellId(1)).unwrap(),
        b.temperature_k(CellId(1)).unwrap()
    );
    let first = law.resolve(&a, CellId(1), MaterialId(7)).unwrap();
    let midpoint = law.resolve(&b, CellId(1), MaterialId(7)).unwrap();
    assert_ne!(first, midpoint);
    assert_eq!(midpoint, law.at_specific_enthalpy(-50.0, 500.0).unwrap());
    close(midpoint.complex_ior(550.0).unwrap().eta(), 2.2);
    assert!(law.resolve(&b, CellId(1), MaterialId(8)).is_err());
    assert!(law.resolve(&b, CellId(2), MaterialId(7)).is_err());
    assert!(
        law.resolve(&native(&caloric, 1.0), CellId(1), MaterialId(7))
            .is_err()
    );
    for wrong in [
        Plateau {
            identity: Some([22; 32]),
            basis: ThermalStorageBasis::Enthalpy,
        },
        Plateau {
            identity: None,
            basis: ThermalStorageBasis::Enthalpy,
        },
        Plateau {
            identity: Some([21; 32]),
            basis: ThermalStorageBasis::InternalEnergy,
        },
    ] {
        assert!(
            law.resolve(&native(&wrong, -50.0), CellId(1), MaterialId(7))
                .is_err()
        );
    }
    let unsupported = SurfaceConductorLaw::try_enthalpy(
        source(12),
        ContentHash([21; 32]),
        (501.0, 600.0),
        &rows,
        2,
    )
    .unwrap();
    assert!(matches!(
        unsupported.resolve(&b, CellId(1), MaterialId(7)),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(1),
            ..
        })
    ));
}
