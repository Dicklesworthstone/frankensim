//! Manufactured composite exercised through the existing explicit timeline and
//! backward-Euler residual. The simulation's supplied mass remains invariant.

use fs_blake3::ContentHash;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, composite::PcmComposite,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, TimedEvent,
    backward_euler::{BackwardEulerStep, StepCoefficients},
};

fn constituents() -> (EquilibriumEnthalpyPhaseCurve, EquilibriumEnthalpyPhaseCurve) {
    let knot = |h, t, f, rho| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: rho,
    };
    (
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([1; 32]),
            SolidLiquidPhase::Solid,
            vec![
                knot(0.0, 300.0, 0.0, 8000.0),
                knot(200_000.0, 500.0, 0.0, 8000.0),
            ],
        )
        .unwrap(),
        EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([2; 32]),
            vec![
                knot(0.0, 300.0, 0.0, 1000.0),
                knot(200_000.0, 400.0, 0.0, 1000.0),
                knot(400_000.0, 400.0, 1.0, 800.0),
                knot(650_000.0, 500.0, 1.0, 800.0),
            ],
        )
        .unwrap(),
    )
}

fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.25,
        energy_tolerance_j: 1e-7,
    }
}

fn initial<'a>(law: &'a PcmComposite<'_>) -> DepositionSimulation<'a> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(1),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: 2.0,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.0, 0.0, 1.0],
                    volume_m3: law.volume_m3(0.0, 2.0).unwrap(),
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 100_000.0,
                },
            }),
        }],
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 128,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn timeline_crosses_latent_plateau_without_creating_mass_or_double_counting_heat() {
    let (matrix, pcm) = constituents();
    let law = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    let original = initial(&law);
    let reference_volume = original.cells()[&CellId(1)].geometry.volume_m3;
    let mut accepted = original.clone();
    let mut previous_time = 0.0;
    for (time, temperature, total_liquid_fraction) in [
        (2.5, 400.0, 0.0),
        (3.0, 400.0, 0.125),
        (3.5, 400.0, 0.25),
        (6.25, 500.0, 0.25),
    ] {
        let (next, audit) = accepted.advance_to(time, config(), || Ok(())).unwrap();
        let cell = &next.cells()[&CellId(1)];
        let state = law
            .values_at_specific_enthalpy(cell.specific_h_j_kg)
            .unwrap();
        assert!((next.temperature_k(CellId(1)).unwrap() - temperature).abs() < 1e-10);
        assert!((state.liquid_mass_fraction - total_liquid_fraction).abs() < 1e-12);
        assert!((cell.specific_h_j_kg - time * 50_000.0).abs() < 1e-8);
        assert_eq!(next.mass_kg(), 2.0);
        assert_eq!(cell.geometry.volume_m3, reference_volume);
        assert!((audit.external_heat_j - 100_000.0 * (time - previous_time)).abs() < 1e-7);
        assert!(audit.residual_j.abs() < 1e-7);
        if total_liquid_fraction > 0.0 {
            assert!(law.volume_m3(cell.specific_h_j_kg, cell.mass_kg).unwrap() > reference_volume);
        }
        previous_time = time;
        accepted = next;
    }
    assert_eq!(original.cells()[&CellId(1)].specific_h_j_kg, 0.0);
    assert_eq!(original.mass_kg(), accepted.mass_kg());
    assert!((accepted.energy_j() - 625_000.0).abs() < 1e-7);
    assert!(accepted.advance_to(7.0, config(), || Ok(())).is_err());
    assert_eq!(accepted.temperature_k(CellId(1)).unwrap(), 500.0);
}

#[test]
fn native_implicit_residual_accepts_latent_storage_and_original_entropy() {
    let (matrix, pcm) = constituents();
    let law = PcmComposite::try_new(&matrix, &pcm, 0.25, 0.0).unwrap();
    let old = initial(&law);
    let step = BackwardEulerStep::new(&old, 3.0, || Ok(())).unwrap();
    let trial = step
        .evaluate_with_entropy(
            &[300_000.0],
            |_| {
                Ok::<_, DepositionError>(StepCoefficients {
                    sources_w: Default::default(),
                    conductances_w_k: Default::default(),
                    payload: (),
                })
            },
            || Ok(()),
        )
        .unwrap();
    assert!(trial.maximum_cell_residual_j() < 1e-8);
    let entropy = trial.entropy_audit().unwrap();
    let expected_entropy = 2.0 * (1250.0 * (400.0_f64 / 300.0).ln() + 25_000.0 / 400.0);
    assert!((entropy.thermal_storage_change_j_k - expected_entropy).abs() < 1e-9);
    let (next, (), audit) = trial.accept(1e-8, 1e-8, || Ok(())).unwrap();
    assert_eq!(next.cells()[&CellId(1)].specific_h_j_kg, 150_000.0);
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 400.0);
    assert_eq!(next.mass_kg(), 2.0);
    assert_eq!(
        law.values_at_specific_enthalpy(150_000.0)
            .unwrap()
            .liquid_mass_fraction,
        0.125
    );
    assert!(audit.residual_j.abs() < 1e-8);
    assert_eq!(audit.external_heat_j, 300_000.0);
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 0.0);
    let refused = step.evaluate(
        &[700_000.0],
        |_| {
            Ok::<_, DepositionError>(StepCoefficients {
                sources_w: Default::default(),
                conductances_w_k: Default::default(),
                payload: (),
            })
        },
        || Ok(()),
    );
    assert!(refused.is_err());
    assert_eq!(old.energy_j(), 0.0);
}
