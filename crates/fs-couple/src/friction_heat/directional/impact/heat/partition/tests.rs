//! Exact partition helper with original finite storage; analytical caloric fixtures.
use super::*;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, ThermalStorageBasis, TimedEvent,
};

#[derive(Clone, Copy)]
enum Caloric {
    Cp,
    Cv,
    Latent,
}
impl ThermalLaw for Caloric {
    fn storage_basis(&self) -> ThermalStorageBasis {
        if matches!(self, Self::Cv) {
            ThermalStorageBasis::InternalEnergy
        } else {
            ThermalStorageBasis::Enthalpy
        }
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([*self as u8 + 1; 32])
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        if !energy.is_finite() || !(0.0..=100_000.0).contains(&energy) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(match self {
            // Cp=T and Cv=T/2, with independent original H/U coordinates.
            Self::Cp => (90_000.0 + 2.0 * energy).sqrt(),
            Self::Cv => (90_000.0 + 4.0 * energy).sqrt(),
            Self::Latent => 300.0 + (energy - 100.0).max(0.0),
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static LAWS: [Caloric; 3] = [Caloric::Cp, Caloric::Cv, Caloric::Latent];
const IDS: [CellId; 3] = [CellId(11), CellId(22), CellId(33)];

fn thermal() -> DepositionSimulation<'static> {
    let events = IDS
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: [0.02, 0.04, 2.0][i],
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    let laws: Vec<_> = LAWS
        .iter()
        .enumerate()
        .map(|(i, law)| (MaterialId(i as u64), law as &dyn ThermalLaw))
        .collect();
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 3,
            max_links: 0,
            max_events: 3,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(
        7.0,
        StepConfig {
            max_dt_s: 7.0,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a:e} != {b:e}");
}

#[test]
fn original_nonlinear_cp_cv_and_latent_inventories_receive_one_declared_loss() {
    let start = thermal();
    let fractions = BTreeMap::from([(IDS[0], 0.25), (IDS[1], 0.5), (IDS[2], 0.125)]);
    let partition = ImpactHeatPartition::new(&start, &fractions, 1e-9, 1e-9, || Ok(())).unwrap();
    let (next, receipt) = partition.deposit(&start, 1200.0, || Ok(())).unwrap();
    close(next.temperature_k(IDS[0]).unwrap(), 120_000.0_f64.sqrt());
    close(next.temperature_k(IDS[1]).unwrap(), 150_000.0_f64.sqrt());
    assert_eq!(next.temperature_k(IDS[2]).unwrap(), 300.0);
    assert_eq!(next.cells()[&IDS[2]].specific_h_j_kg, 75.0);
    close(receipt.thermal.stored_energy_change_j, 1050.0);
    close(receipt.thermal.external_heat_j, 1050.0);
    assert_eq!(receipt.unallocated_energy_j, 150.0);
    assert_eq!(receipt.partition_energy_residual_j, 0.0);
    assert_eq!(receipt.thermal.residual_j, 0.0);
    assert_eq!((receipt.time_s, next.time_s()), (7.0, 7.0));
    assert_eq!((receipt.thermal.events, receipt.thermal.substeps), (0, 0));
    for id in IDS {
        let old = &start.cells()[&id];
        let mut expected = old.clone();
        expected.specific_h_j_kg = next.cells()[&id].specific_h_j_kg;
        assert_eq!(next.cells()[&id], expected);
        assert_eq!(
            next.material_storage_basis(old.material_id).unwrap(),
            start.material_storage_basis(old.material_id).unwrap()
        );
        assert_eq!(
            next.material_caloric_identity(old.material_id).unwrap(),
            start.material_caloric_identity(old.material_id).unwrap()
        );
    }
}

#[test]
fn fractions_cannot_overspend_even_below_one_ulp_and_zero_heat_preserves_the_receiver() {
    let start = thermal();
    for fractions in [
        BTreeMap::from([(IDS[0], -0.1)]),
        BTreeMap::from([(IDS[0], f64::NAN)]),
        BTreeMap::from([(IDS[0], f64::INFINITY)]),
        BTreeMap::from([(IDS[0], 1.01)]),
        BTreeMap::from([(IDS[0], 1.0), (IDS[1], f64::MIN_POSITIVE)]),
        BTreeMap::from([(CellId(999), 0.0)]),
    ] {
        assert!(ImpactHeatPartition::new(&start, &fractions, 0.0, 0.0, || Ok(())).is_err());
    }
    for fractions in [BTreeMap::new(), BTreeMap::from([(IDS[0], 0.0)])] {
        let (next, receipt) = ImpactHeatPartition::new(&start, &fractions, 0.0, 0.0, || Ok(()))
            .unwrap()
            .deposit(&start, 23.0, || Ok(()))
            .unwrap();
        assert_eq!(next.cells(), start.cells());
        assert_eq!(receipt.unallocated_energy_j, 23.0);
        assert_eq!(receipt.thermal.external_heat_j, 0.0);
        assert_eq!(next.time_s(), start.time_s());
    }
    let fractions = BTreeMap::from([(IDS[0], 0.5)]);
    assert!(matches!(
        ImpactHeatPartition::new(&start, &fractions, 0.0, 0.0, || Ok(()))
            .unwrap()
            .deposit(&start, f64::from_bits(1), || Ok(())),
        Err(DepositionError::Arithmetic)
    ));
}

#[test]
fn late_cancellation_and_caloric_refusal_preserve_original_storage_and_retry() {
    let start = thermal();
    let fractions = BTreeMap::from([(IDS[0], 0.25), (IDS[1], 0.5)]);
    let prepare = || ImpactHeatPartition::new(&start, &fractions, 1e-9, 1e-9, || Ok(())).unwrap();
    let mut calls = 0;
    let (reference, receipt) = prepare()
        .deposit(&start, 1200.0, || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut retries = 0;
    let refused = prepare().deposit(&start, 1200.0, || {
        retries += 1;
        if retries == calls {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(refused, Err(DepositionError::Cancelled)));
    assert!(matches!(
        prepare().deposit(&start, 1e12, || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    assert_eq!(start.time_s(), 7.0);
    assert!(
        start
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 0.0)
    );
    let (retry, retried) = prepare().deposit(&start, 1200.0, || Ok(())).unwrap();
    assert_eq!(reference.cells(), retry.cells());
    assert_eq!(receipt.thermal, retried.thermal);
    assert_eq!(receipt.unallocated_energy_j, retried.unallocated_energy_j);
}
