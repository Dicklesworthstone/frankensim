use super::*;
use crate::conductor::{ConductorIorSample, ConductorOptics};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, Event, StepConfig,
    ThermalLaw, TimedEvent,
};

fn source(id: u8) -> ConductorSource {
    ConductorSource::try_new(ContentHash([id; 32]), ConductorDataStatus::Representative).unwrap()
}

fn film(eta: f64, k: f64) -> CoherentAbsorbingCoating {
    CoherentAbsorbingCoating::try_new(
        core::array::from_fn(|i| {
            CoatingIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        }),
        90e-9,
        source(31),
    )
    .unwrap()
}

fn substrate() -> ConductorOptics {
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }),
        source(32),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 2e-14, "{actual} != {expected}");
}

// G0/G3: constitutive interpolation acts on n/k before coherent field work.
#[test]
fn state_resolves_complex_index_before_coherent_absorption() {
    let rows = [
        (300.0, film(1.5, 0.0)),
        (500.0, film(2.5, 0.4)),
        (900.0, film(3.5, 0.8)),
    ];
    let law = SurfaceCoatingLaw::try_temperature(source(33), &rows, 3).unwrap();
    let coating = law.at_temperature(700.0).unwrap();
    for sample in coating.samples() {
        close(sample.eta(), 3.0);
        close(sample.k(), 0.6);
    }
    assert_eq!(coating.thickness_m(), rows[0].1.thickness_m());
    assert_eq!(coating.source(), source(33));
    let base = substrate();
    let absorption = |coating| {
        base.with_absorbing_coating(coating)
            .layer_power_balance(550.0, 1.0, 0.7)
            .unwrap()
            .layer_absorbed()[0]
    };
    let actual = absorption(coating);
    let expected = absorption(film(3.0, 0.6));
    for pol in 0..2 {
        close(actual[pol], expected[pol]);
        let mixed_intensity = 0.5 * (absorption(rows[1].1)[pol] + absorption(rows[2].1)[pol]);
        assert!((actual[pol] - mixed_intensity).abs() > 1e-4);
    }
    for (coordinate, endpoint) in rows {
        assert_eq!(
            law.at_temperature(coordinate).unwrap().samples(),
            endpoint.samples()
        );
    }
    let constant = SurfaceCoatingLaw::try_temperature(
        source(33),
        &[(300.0, rows[0].1), (900.0, rows[0].1)],
        2,
    )
    .unwrap();
    assert_eq!(
        constant.at_temperature(300.0).unwrap(),
        constant.at_temperature(631.25).unwrap()
    );
    let cloned = law.clone();
    let Model::Rows { rows, .. } = &law.model else {
        panic!("authored coating rows")
    };
    let Model::Rows { rows: other, .. } = &cloned.model else {
        panic!("cloned authored coating rows")
    };
    assert!(Arc::ptr_eq(rows, other));
}

// G0: exact transparent endpoints and representable weak losses stay distinct.
#[test]
fn transparent_and_subnormal_extinction_survive_state_interpolation() {
    let law = SurfaceCoatingLaw::try_temperature(
        source(33),
        &[
            (300.0, film(f64::from_bits(1), 0.0)),
            (700.0, film(f64::from_bits(2), f64::from_bits(3))),
        ],
        2,
    )
    .unwrap();
    for (temperature, eta_bits, k_bits) in [
        (300.0, 1, 0),
        (400.0, 1, 1),
        (500.0, 2, 2),
        (600.0, 2, 2),
        (700.0, 2, 3),
    ] {
        for sample in law.at_temperature(temperature).unwrap().samples() {
            assert_eq!(sample.eta().to_bits(), eta_bits);
            assert_eq!(sample.k().to_bits(), k_bits);
        }
    }
    let weak = SurfaceCoatingLaw::try_temperature(
        source(33),
        &[(300.0, film(1.5, 0.0)), (700.0, film(1.5, 1e-20))],
        2,
    )
    .unwrap();
    let resolved = weak.at_temperature(500.0).unwrap();
    assert_eq!(resolved.samples()[0].k(), 5e-21);
    let work = substrate()
        .with_absorbing_coating(resolved)
        .layer_power_balance(550.0, 1.0, 1.0)
        .unwrap();
    assert!(
        work.layer_absorbed()[0]
            .iter()
            .all(|&a| a > 0.0 && a < 1e-19)
    );
    let transparent = substrate()
        .with_absorbing_coating(weak.at_temperature(300.0).unwrap())
        .layer_power_balance(550.0, 1.0, 1.0)
        .unwrap();
    assert_eq!(transparent.layer_absorbed()[0], [0.0; 2]);
}

#[test]
fn admission_preserves_fixed_geometry_sources_and_domains() {
    let rows = [(300.0, film(1.5, 0.0)), (700.0, film(3.5, 0.8))];
    let law = SurfaceCoatingLaw::try_temperature(source(33), &rows, 2).unwrap();
    assert_eq!(law.row_count(), 2);
    assert_eq!(law.temperature_bounds_k(), (300.0, 700.0));
    assert_eq!(law.caloric_identity(), None);
    assert_eq!(law.specific_enthalpy_bounds_j_kg(), None);
    assert_eq!(
        law.identity(),
        SurfaceCoatingLaw::try_temperature(source(33), &rows, 10)
            .unwrap()
            .identity()
    );
    assert!(matches!(
        SurfaceCoatingLaw::try_temperature(source(33), &rows, 1),
        Err(SurfaceHeatingError::Budget)
    ));
    assert!(SurfaceCoatingLaw::try_temperature(source(33), &[], 2).is_err());
    assert!(SurfaceCoatingLaw::try_temperature(source(33), &rows[..1], 2).is_err());
    for coordinates in [
        [700.0, 300.0],
        [300.0, 300.0],
        [0.0, 700.0],
        [300.0, f64::NAN],
    ] {
        let invalid = [(coordinates[0], rows[0].1), (coordinates[1], rows[1].1)];
        assert!(SurfaceCoatingLaw::try_temperature(source(33), &invalid, 2).is_err());
    }
    for outside in [299.0, 701.0, f64::NAN, f64::INFINITY] {
        assert!(law.at_temperature(outside).is_err());
    }
    assert!(law.at_specific_enthalpy(500.0, 500.0).is_err());
    let mut shifted_samples = rows[1].1.samples();
    shifted_samples[1] = CoatingIorSample::try_new(440.0, 3.5, 0.8).unwrap();
    for changed in [
        CoherentAbsorbingCoating::try_new(shifted_samples, 90e-9, source(31)).unwrap(),
        CoherentAbsorbingCoating::try_new(rows[1].1.samples(), 91e-9, source(31)).unwrap(),
    ] {
        assert!(
            SurfaceCoatingLaw::try_temperature(
                source(33),
                &[(300.0, rows[0].1), (700.0, changed)],
                2,
            )
            .is_err()
        );
    }
    let measured = ConductorSource::try_new(
        ContentHash([34; 32]),
        ConductorDataStatus::CallerAssertedMeasured,
    )
    .unwrap();
    assert!(SurfaceCoatingLaw::try_temperature(measured, &rows, 2).is_err());
    let measured_rows = rows.map(|(x, coating)| {
        (
            x,
            CoherentAbsorbingCoating::try_new(coating.samples(), 90e-9, measured).unwrap(),
        )
    });
    assert!(SurfaceCoatingLaw::try_temperature(measured, &measured_rows, 2).is_ok());
    let changed_source = [
        rows[0],
        (
            700.0,
            CoherentAbsorbingCoating::try_new(rows[1].1.samples(), 90e-9, source(35)).unwrap(),
        ),
    ];
    for different in [
        SurfaceCoatingLaw::try_temperature(source(34), &rows, 2).unwrap(),
        SurfaceCoatingLaw::try_temperature(source(33), &changed_source, 2).unwrap(),
        SurfaceCoatingLaw::try_temperature(
            source(33),
            &[(300.0, rows[0].1), (700.0, film(3.5, 0.7))],
            2,
        )
        .unwrap(),
        SurfaceCoatingLaw::try_enthalpy(
            source(33),
            ContentHash([36; 32]),
            (300.0, 700.0),
            &rows,
            2,
        )
        .unwrap(),
    ] {
        assert_ne!(law.identity(), different.identity());
    }
    for (identity, bounds, coordinates) in [
        (ContentHash([0; 32]), (500.0, 500.0), [-100.0, 0.0]),
        (ContentHash([36; 32]), (0.0, 500.0), [-100.0, 0.0]),
        (ContentHash([36; 32]), (500.0, 499.0), [-100.0, 0.0]),
        (ContentHash([36; 32]), (500.0, 500.0), [-f64::MAX, f64::MAX]),
    ] {
        assert!(
            SurfaceCoatingLaw::try_enthalpy(
                source(33),
                identity,
                bounds,
                &[(coordinates[0], rows[0].1), (coordinates[1], rows[1].1)],
                2,
            )
            .is_err()
        );
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

// G0/G3: equal-temperature latent states retain their own film absorption.
#[test]
fn native_enthalpy_resolves_latent_states_only_with_exact_caloric_binding() {
    let rows = [(-100.0, film(1.5, 0.0)), (0.0, film(3.5, 0.8))];
    let law = SurfaceCoatingLaw::try_enthalpy(
        source(33),
        ContentHash([36; 32]),
        (500.0, 500.0),
        &rows,
        2,
    )
    .unwrap();
    assert_eq!(law.caloric_identity(), Some(ContentHash([36; 32])));
    assert_eq!(law.specific_enthalpy_bounds_j_kg(), Some((-100.0, 0.0)));
    assert_eq!(law.temperature_bounds_k(), (500.0, 500.0));
    assert!(law.at_temperature(500.0).is_err());
    for (h, temperature) in [
        (-101.0, 500.0),
        (1.0, 500.0),
        (-50.0, 501.0),
        (f64::NAN, 500.0),
    ] {
        assert!(law.at_specific_enthalpy(h, temperature).is_err());
    }
    let caloric = Plateau {
        identity: Some([36; 32]),
        basis: ThermalStorageBasis::Enthalpy,
    };
    let first = native(&caloric, -100.0);
    let middle = native(&caloric, -50.0);
    assert_eq!(
        first.temperature_k(CellId(1)).unwrap(),
        middle.temperature_k(CellId(1)).unwrap()
    );
    let first_film = law.resolve(&first, CellId(1), MaterialId(7)).unwrap();
    let midpoint = law.resolve(&middle, CellId(1), MaterialId(7)).unwrap();
    assert_eq!(first_film.samples()[0].k(), 0.0);
    assert_eq!(midpoint, law.at_specific_enthalpy(-50.0, 500.0).unwrap());
    close(midpoint.samples()[0].eta(), 2.5);
    close(midpoint.samples()[0].k(), 0.4);
    let work = substrate()
        .with_absorbing_coating(midpoint)
        .layer_power_balance(550.0, 1.0, 1.0)
        .unwrap();
    assert!(work.layer_absorbed()[0].iter().all(|&a| a > 0.0));
    assert!(law.resolve(&middle, CellId(1), MaterialId(8)).is_err());
    assert!(law.resolve(&middle, CellId(2), MaterialId(7)).is_err());
    assert!(
        law.resolve(&native(&caloric, 1.0), CellId(1), MaterialId(7))
            .is_err()
    );
    for wrong in [
        Plateau {
            identity: Some([37; 32]),
            basis: ThermalStorageBasis::Enthalpy,
        },
        Plateau {
            identity: None,
            basis: ThermalStorageBasis::Enthalpy,
        },
        Plateau {
            identity: Some([36; 32]),
            basis: ThermalStorageBasis::InternalEnergy,
        },
    ] {
        assert!(
            law.resolve(&native(&wrong, -50.0), CellId(1), MaterialId(7))
                .is_err()
        );
    }
    let unsupported = SurfaceCoatingLaw::try_enthalpy(
        source(33),
        ContentHash([36; 32]),
        (501.0, 600.0),
        &rows,
        2,
    )
    .unwrap();
    assert!(matches!(
        unsupported.resolve(&middle, CellId(1), MaterialId(7)),
        Err(SurfaceHeatingError::TemperatureDomain {
            thermal_cell: CellId(1),
            ..
        })
    ));
}

#[test]
fn temperature_law_reads_native_temperature_with_either_storage_basis() {
    let law = SurfaceCoatingLaw::try_temperature(
        source(33),
        &[(400.0, film(1.5, 0.0)), (600.0, film(3.5, 0.8))],
        2,
    )
    .unwrap();
    for basis in [
        ThermalStorageBasis::Enthalpy,
        ThermalStorageBasis::InternalEnergy,
    ] {
        let caloric = Plateau {
            identity: None,
            basis,
        };
        let native = native(&caloric, -75.0);
        let resolved = law.resolve(&native, CellId(1), MaterialId(7)).unwrap();
        assert_eq!(resolved, law.at_temperature(500.0).unwrap());
        assert_eq!(
            law.coordinate_at(&native, CellId(1), MaterialId(7))
                .unwrap(),
            500.0
        );
    }
}
