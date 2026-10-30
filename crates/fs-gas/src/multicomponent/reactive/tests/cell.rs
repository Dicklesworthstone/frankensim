use super::*;
use crate::multicomponent::{NormalCell, NormalPrimitive};

fn prepared(heat: f64, reference: f64, method: IsochoricKineticsMethod) -> CellChemistry {
    let mechanism = mechanism(heat, 5.0, reference);
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let controls = IsochoricKineticsSolve {
        method,
        ..chemistry()
    };
    CellChemistry::new(gas, mechanism, controls, || false).unwrap()
}

fn planar(gas: &MixtureGas, velocity: [f64; 2]) -> [f64; 5] {
    let mut q = [0.0; 5];
    gas.normal_conserved_into(
        NormalPrimitive {
            partial_densities: &[1.0, 0.0],
            normal_velocity: velocity[0],
            transverse_velocity: velocity[1],
            pressure: 1.0,
        },
        &mut q,
    )
    .unwrap();
    q
}

fn primitive(gas: &MixtureGas, q: &[f64; 5]) -> super::super::super::Primitive {
    gas.normal_primitive(
        NormalCell {
            partial_densities: &q[..2],
            normal_momentum: q[2],
            transverse_momentum: q[3],
            energy: q[4],
        },
        || false,
    )
    .unwrap()
}

#[test]
fn g1_planar_chemical_heat_matches_original_reactor_and_keeps_both_momenta() {
    for method in [
        IsochoricKineticsMethod::BackwardEuler,
        IsochoricKineticsMethod::Sdirk2,
    ] {
        let mut source = prepared(40.0, 0.0, method);
        let mut q = planar(source.gas(), [0.7, -1.1]);
        let before = q;
        let volume = Volume::new(0.2);
        let duration = Time::new(0.02);
        let exact = advance_adiabatic(
            source.mechanism(),
            &[(
                SpeciesId::new("A").unwrap(),
                Amount::new(volume.value() / R),
            )],
            volume,
            Temperature::new(1.0),
            duration,
            source.controls(),
            || false,
        )
        .unwrap();
        let workspace = source.workspace_bytes();
        let mut report = ChemistryReport::default();
        source
            .advance_planar(&mut q, volume, duration, 17, &mut report, || false)
            .unwrap();
        assert_eq!(&q[2..4], &before[2..4]);
        assert_eq!(source.workspace_bytes(), workspace);
        assert!(workspace > 0);
        let p = primitive(source.gas(), &q);
        assert!(p.temperature > 1.0);
        assert!((p.temperature - exact.temperature().value()).abs() < 2e-9);
        assert!((q[1] - exact.amounts()[1].value() * R / volume.value()).abs() < 2e-10);
        let physical_change = q[4] - before[4]
            + (q[0] - before[0]) * source.gas().energy_references()[0]
            + (q[1] - before[1]) * source.gas().energy_references()[1];
        assert!(physical_change.abs() < 1e-13);
        assert_eq!(report.cell_solves, 1);
        assert!(report.entropy_change_j_k > 0.0);
        assert!(report.maximum_reconstruction_residual_j < 1e-9);
    }
}

#[test]
fn g3_shared_chemical_update_respects_axial_limit_rotation_and_reference_gauge() {
    let mut baseline: Option<(f64, f64)> = None;
    for reference in [0.0, 1e6] {
        for velocity in [[0.7, 0.0], [0.7, -1.1], [-1.1, 0.7]] {
            let mut source = prepared(40.0, reference, IsochoricKineticsMethod::Sdirk2);
            let mut q = planar(source.gas(), velocity);
            let original = q;
            let mut report = ChemistryReport::default();
            source
                .advance_planar(
                    &mut q,
                    Volume::new(0.2),
                    Time::new(0.02),
                    0,
                    &mut report,
                    || false,
                )
                .unwrap();
            let p = primitive(source.gas(), &q);
            match baseline {
                Some((product, temperature)) => {
                    assert!((q[1] - product).abs() < 1e-8);
                    assert!((p.temperature - temperature).abs() < 1e-8);
                }
                None => baseline = Some((q[1], p.temperature)),
            }
            if velocity[1] == 0.0 {
                let mut axial = [original[0], original[1], original[2], original[4]];
                let mut axial_report = ChemistryReport::default();
                source
                    .advance_axial(
                        &mut axial,
                        Volume::new(0.2),
                        Time::new(0.02),
                        0,
                        &mut axial_report,
                        || false,
                    )
                    .unwrap();
                assert_eq!(axial, [q[0], q[1], q[2], q[4]]);
                assert_eq!(axial_report, report);
            }
        }
    }
}

#[test]
fn g4_planar_source_publishes_cell_and_report_together_after_final_cancellation_poll() {
    let mut source = prepared(40.0, 0.0, IsochoricKineticsMethod::Sdirk2);
    let original = planar(source.gas(), [0.7, -1.1]);
    let mut expected = original;
    let mut expected_report = ChemistryReport::default();
    let mut count = 0;
    source
        .advance_planar(
            &mut expected,
            Volume::new(0.2),
            Time::new(0.02),
            8,
            &mut expected_report,
            || {
                count += 1;
                false
            },
        )
        .unwrap();
    let mut q = original;
    let mut report = ChemistryReport::default();
    let mut current = 0;
    let result = source.advance_planar(
        &mut q,
        Volume::new(0.2),
        Time::new(0.02),
        8,
        &mut report,
        || {
            current += 1;
            current == count
        },
    );
    assert_eq!(result, Err(ReactiveError::Cancelled));
    assert_eq!(current, count);
    assert_eq!(q, original);
    assert_eq!(report, ChemistryReport::default());
    source
        .advance_planar(
            &mut q,
            Volume::new(0.2),
            Time::new(0.02),
            8,
            &mut report,
            || false,
        )
        .unwrap();
    assert_eq!(q, expected);
    assert_eq!(report, expected_report);

    let mut overflow = ChemistryReport {
        cell_solves: usize::MAX,
        ..ChemistryReport::default()
    };
    let old_report = overflow;
    let old_q = q;
    assert!(
        source
            .advance_planar(
                &mut q,
                Volume::new(0.2),
                Time::new(0.02),
                8,
                &mut overflow,
                || false
            )
            .is_err()
    );
    assert_eq!(q, old_q);
    assert_eq!(overflow, old_report);
}

#[test]
fn g0_cell_chemistry_refuses_wrong_sources_dimensions_and_caloric_paths() {
    let mut source = prepared(40.0, 0.0, IsochoricKineticsMethod::BackwardEuler);
    let mut q = planar(source.gas(), [0.7, -1.1]);
    let original = q;
    let mut report = ChemistryReport::default();
    assert!(
        CellChemistry::new(
            source.gas().clone(),
            mechanism(41.0, 5.0, 0.0),
            chemistry(),
            || false
        )
        .is_err()
    );
    assert!(
        CellChemistry::new(
            source.gas().clone(),
            source.mechanism().clone(),
            IsochoricKineticsSolve {
                maximum_temperature_k: 51.0,
                ..chemistry()
            },
            || false
        )
        .is_err()
    );
    assert!(
        source
            .advance_planar(
                &mut q[..4],
                Volume::new(0.2),
                Time::new(0.02),
                0,
                &mut report,
                || false
            )
            .is_err()
    );
    assert!(
        source
            .advance_planar(
                &mut q,
                Volume::new(0.0),
                Time::new(0.02),
                0,
                &mut report,
                || false
            )
            .is_err()
    );
    assert_eq!(q, original);
    assert_eq!(report, ChemistryReport::default());

    let mut narrow = CellChemistry::new(
        source.gas().clone(),
        source.mechanism().clone(),
        IsochoricKineticsSolve {
            maximum_temperature_k: 1.01,
            maximum_steps: 50,
            ..chemistry()
        },
        || false,
    )
    .unwrap();
    assert!(
        narrow
            .advance_planar(
                &mut q,
                Volume::new(0.2),
                Time::new(0.02),
                0,
                &mut report,
                || false
            )
            .is_err()
    );
    assert_eq!(q, original);
    assert_eq!(report, ChemistryReport::default());
}
