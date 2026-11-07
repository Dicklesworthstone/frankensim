use super::*;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, StepConfig, ThermalLaw,
    ThermalStorageBasis, TimedEvent,
};

struct CaloricSource {
    identity: Option<[u8; 32]>,
    basis: ThermalStorageBasis,
    datum: f64,
}
impl ThermalLaw for CaloricSource {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.identity
    }
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let h = h - self.datum;
        if !h.is_finite() || !(0.0..=200.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(500.0 + (h - 100.0).max(0.0))
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static CALORIC: CaloricSource = CaloricSource {
    identity: Some([2; 32]),
    basis: ThermalStorageBasis::Enthalpy,
    datum: 0.0,
};
static SHIFTED: CaloricSource = CaloricSource {
    identity: Some([3; 32]),
    basis: ThermalStorageBasis::Enthalpy,
    datum: -1000.0,
};
static INTERNAL: CaloricSource = CaloricSource {
    identity: Some([2; 32]),
    basis: ThermalStorageBasis::InternalEnergy,
    datum: 0.0,
};
static UNIDENTIFIED: CaloricSource = CaloricSource {
    identity: None,
    basis: ThermalStorageBasis::Enthalpy,
    datum: 0.0,
};

fn simulation(law: &'static CaloricSource, h: f64, source_w: f64) -> DepositionSimulation<'static> {
    let simulation = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(7),
                segment_id: 17,
                body_id: 1,
                material_id: MaterialId(5),
                mass_kg: 1.0,
                specific_h_j_kg: law.datum + h,
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
                    source_w,
                },
            }),
        }],
        &[(MaterialId(5), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 32,
        },
    )
    .unwrap();
    simulation.advance_to(0.0, config(), || Ok(())).unwrap().0
}

fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1e-10,
    }
}

fn segment() -> BoundSegment {
    BoundSegment {
        path: BeamCellSegment {
            cell_index: 3,
            density_kg_m3: 1.0,
            length_m: 1.0,
        },
        thermal_cell: CellId(7),
        table: 0,
        material_id: Some(MaterialId(5)),
    }
}

fn table(datum: f64, identity: u8, max_h: f64, values: [f64; 2]) -> EnthalpyAbsorptionTable {
    EnthalpyAbsorptionTable::try_new(
        ContentHash([1; 32]),
        ContentHash([identity; 32]),
        (500.0, 500.0),
        &[
            (datum, MassAbsorptionSpectrum::constant(values[0]).unwrap()),
            (
                datum + max_h,
                MassAbsorptionSpectrum::constant(values[1]).unwrap(),
            ),
        ],
        2,
    )
    .unwrap()
}

fn compile(table: &EnthalpyAbsorptionTable) -> CompiledAbsorption<1> {
    SourceTable::Enthalpy(table)
        .compile(
            &[550.0],
            &mut DomainHasher::new("test.direct-beam.enthalpy"),
            &mut || Ok(()),
        )
        .unwrap()
}

#[test]
fn g1_compiled_beam_table_tracks_native_latent_storage() {
    let compiled = compile(&table(0.0, 2, 100.0, [0.2, 1.2]));
    let original = simulation(&CALORIC, 10.0, 10.0);
    let (before, h) = compiled.at(&original, &segment(), 500.0).unwrap();
    assert!((before[0] - 0.3).abs() < 1e-15);
    assert_eq!(h, Some(10.0));
    let (heated, audit) = original.advance_to(2.0, config(), || Ok(())).unwrap();
    assert_eq!(heated.temperature_k(CellId(7)).unwrap(), 500.0);
    assert_eq!(audit.external_heat_j, 20.0);
    let (after, h) = compiled.at(&heated, &segment(), 500.0).unwrap();
    assert!((after[0] - 0.5).abs() < 1e-15);
    assert_eq!(h, Some(30.0));
    assert_eq!(original.cells()[&CellId(7)].specific_h_j_kg, 10.0);
}

#[test]
fn g3_negative_caloric_datum_preserves_coefficients_without_losing_the_state_key() {
    let original = simulation(&CALORIC, 10.0, 10.0);
    let shifted = simulation(&SHIFTED, 10.0, 10.0);
    let a = compile(&table(0.0, 2, 100.0, [0.2, 1.2]));
    let b = compile(&table(-1000.0, 3, 100.0, [0.2, 1.2]));
    for time in [0.0, 0.5, 2.0] {
        let a_state = original.advance_to(time, config(), || Ok(())).unwrap().0;
        let b_state = shifted.advance_to(time, config(), || Ok(())).unwrap().0;
        let (ka, ha) = a.at(&a_state, &segment(), 500.0).unwrap();
        let (kb, hb) = b.at(&b_state, &segment(), 500.0).unwrap();
        assert_eq!(ka, kb);
        assert_eq!(ha.unwrap(), hb.unwrap() + 1000.0);
    }
}

#[test]
fn g0_native_material_identity_storage_and_temperature_domains_remain_distinct() {
    let compiled = compile(&table(0.0, 2, 100.0, [0.2, 1.2]));
    let original = simulation(&CALORIC, 10.0, 0.0);
    for law in [&SHIFTED, &INTERNAL, &UNIDENTIFIED] {
        let state = simulation(law, 10.0, 0.0);
        assert!(matches!(
            compiled.check(&state, &segment(), 500.0),
            Err(OpticalFeedbackError::Invalid(_))
        ));
    }
    for material_id in [None, Some(MaterialId(6))] {
        let mut wrong = segment();
        wrong.material_id = material_id;
        assert!(matches!(
            compiled.at(&original, &wrong, 500.0),
            Err(OpticalFeedbackError::Invalid(_))
        ));
    }
    let hot = simulation(&CALORIC, 110.0, 0.0);
    assert!(matches!(
        compiled.check(&hot, &segment(), hot.temperature_k(CellId(7)).unwrap()),
        Err(OpticalFeedbackError::TemperatureDomain { .. })
    ));
}

#[test]
fn g4_transparent_beam_still_refuses_out_of_domain_native_enthalpy() {
    let compiled = compile(&table(0.0, 2, 5.0, [0.0; 2]));
    let original = simulation(&CALORIC, 1.0, 10.0);
    assert_eq!(compiled.at(&original, &segment(), 500.0).unwrap().0, [0.0]);
    let (heated, _) = original.advance_to(0.5, config(), || Ok(())).unwrap();
    assert_eq!(heated.temperature_k(CellId(7)).unwrap(), 500.0);
    assert!(matches!(
        compiled.check(&heated, &segment(), 500.0),
        Err(OpticalFeedbackError::EnthalpyDomain {
            thermal_cell: CellId(7),
            specific_h_j_kg: 6.0,
            ..
        })
    ));
    assert_eq!(original.cells()[&CellId(7)].specific_h_j_kg, 1.0);
}

#[test]
fn g0_temperature_parity_spectral_support_and_compilation_cancellation() {
    let original = simulation(&CALORIC, 10.0, 0.0);
    let h = table(0.0, 2, 100.0, [0.4; 2]);
    let t = TemperatureAbsorptionTable::try_new(
        ContentHash([4; 32]),
        &[
            (450.0, MassAbsorptionSpectrum::constant(0.4).unwrap()),
            (550.0, MassAbsorptionSpectrum::constant(0.4).unwrap()),
        ],
        2,
    )
    .unwrap();
    let thermal = SourceTable::Temperature(&t)
        .compile(
            &[550.0],
            &mut DomainHasher::new("test.direct-beam.temperature"),
            &mut || Ok(()),
        )
        .unwrap();
    let (h_coefficients, h_key) = compile(&h).at(&original, &segment(), 500.0).unwrap();
    let (t_coefficients, t_key) = thermal.at(&original, &segment(), 500.0).unwrap();
    assert_eq!(h_coefficients, t_coefficients);
    assert_eq!(h_key, Some(10.0));
    assert_eq!(t_key, None);
    assert!(matches!(
        SourceTable::Enthalpy(&h).compile(
            &[1000.0],
            &mut DomainHasher::new("test.direct-beam.off-band"),
            &mut || Ok(()),
        ),
        Err(OpticalFeedbackError::Transport(
            VolumeTransportError::Wavelength
        ))
    ));
    let mut polls = 0;
    assert!(matches!(
        SourceTable::Enthalpy(&h).compile(
            &[550.0],
            &mut DomainHasher::new("test.direct-beam.cancel"),
            &mut || {
                polls += 1;
                if polls == 2 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(polls, 2);
}
