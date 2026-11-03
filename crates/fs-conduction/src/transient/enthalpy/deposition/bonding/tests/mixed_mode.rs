use super::*;
use fs_material::healing::strength::{MixedModeStrengthConfig, MixedModeStrengthKnot};

#[test]
fn actual_deposition_exposure_drives_mode_dependent_failure_and_preserves_area() {
    let mobility = mobility();
    let map = StrengthMap::new_mixed_mode(
        "declared-interlayer-modes".into(),
        "test-mobility".into(),
        vec![
            MixedModeStrengthKnot {
                dose: 0.0,
                normal_strength_pa: 0.0,
                tangential_strength_pa: 0.0,
                mode_i_fracture_energy_j_m2: 0.0,
                mode_ii_fracture_energy_j_m2: 0.0,
            },
            MixedModeStrengthKnot {
                dose: 2.0,
                normal_strength_pa: 2e5,
                tangential_strength_pa: 3e5,
                mode_i_fracture_energy_j_m2: 100.0,
                mode_ii_fracture_energy_j_m2: 240.0,
            },
        ],
        MixedModeStrengthConfig {
            normal_stiffness_pa_m: 2e9,
            tangential_stiffness_pa_m: 5e8,
            bk_exponent: 2.0,
            compression_stiffness_pa_m: 4e9,
            maximum_separation_m: 0.01,
        },
    )
    .unwrap();
    let whole = simulation(schedule(500.0, 0.0), &mobility, &[7], 1000);
    let split = simulation(
        vec![
            birth(1, 500.0, 0.0),
            birth(2, 500.0, 0.0),
            event(0.0, Event::Contact(contact(7, 0.0001, 0.1))),
            event(0.0, Event::Contact(contact(8, 0.0001, 0.1))),
        ],
        &mobility,
        &[7, 8],
        1000,
    );
    let (whole, audit, _) = whole.advance_to(1.0, config(0.25), || Ok(())).unwrap();
    let (split, split_audit, _) = split.advance_to(1.0, config(0.25), || Ok(())).unwrap();
    assert_eq!(whole.thermal().cells(), split.thermal().cells());
    // Splitting a contact adds a declaration event, but changes no physical
    // mass, heat or stored energy.
    assert_eq!(
        [
            audit.incoming_mass_kg,
            audit.incoming_energy_j,
            audit.external_heat_j,
            audit.stored_energy_change_j,
            audit.residual_j
        ],
        [
            split_audit.incoming_mass_kg,
            split_audit.incoming_energy_j,
            split_audit.external_heat_j,
            split_audit.stored_energy_change_j,
            split_audit.residual_j
        ]
    );
    let original = whole.freeze_bond(7, &map).unwrap();
    assert_eq!(original.interface_id, 7);
    assert_eq!((original.a, original.b), (CellId(1), CellId(2)));
    assert_eq!(original.area_m2, 0.0002);
    close(original.bond.dose, 1.0, 1e-14);
    let law = original.bond.law.unwrap();
    let p = law.mixed_mode_parameters().unwrap();
    close(p.normal_strength_pa, 1e5, 1e-8);
    close(p.tangential_strength_pa, 1.5e5, 1e-8);
    for (axis, strength, fracture) in [(0, 1e5, 50.0), (1, 1.5e5, 120.0)] {
        let mut jump = [0.0; 3];
        jump[axis] = 1.01 * 2.0 * fracture / strength;
        let step = law.bonded_point().trial(jump, 1e-9).unwrap();
        assert!(step.point.failed());
        let work = original.area_m2 * step.energy.work;
        close(work, original.area_m2 * fracture, 1e-12);
        let split_work: f64 = [7, 8]
            .into_iter()
            .map(|id| {
                let patch = split.freeze_bond(id, &map).unwrap();
                let law = patch.bond.law.unwrap();
                patch.area_m2 * law.bonded_point().trial(jump, 1e-9).unwrap().energy.work
            })
            .sum();
        close(work, split_work, 1e-12);
    }
}
