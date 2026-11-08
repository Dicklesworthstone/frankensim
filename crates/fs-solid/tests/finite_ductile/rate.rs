//! G1/G3: original equilibrium, physical intervals and momentum retain kinetics.
use super::*;
use fs_material::ductile::evolution::DuctileJ2Conditions;
use fs_material::plastic::finite::rate::FiniteOverstressParameters;
use fs_solid::hyper3::plastic::ductile::dynamics::{
    FiniteDuctileMotionControls, FiniteDuctileTetMotion,
};
use fs_solid::hyper3::plastic::ductile::path::FiniteDuctilePathControls;

const VOLUME: f64 = 1.0 / 6.0;

fn kinetic_model(width: f64) -> FiniteDuctileJ2<TabulatedDuctileLocus> {
    model(width, 1.0).with_overstress(FiniteOverstressParameters {
        viscosity_pa_s: 100.0,
        stress_scale_pa: 3.0,
        exponent: 1.0,
    }).unwrap()
}

#[test]
fn force_controlled_continuation_uses_retained_equilibrium_and_each_physical_subinterval() {
    with_cx(|_, cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = kinetic_model(bands[0].band_width_m());
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry, &materials, &bands, TetAssemblyBudget::standard(), cx,
        ).unwrap();
        let target = stretch(0.007);
        let full_supports = affine(&NODES, target, [0.0; 3], 0);
        let supports: Vec<_> = full_supports.into_iter().filter(|&(dof, _)| dof != 3).collect();
        let expected = material.trial(&target, &material.initial_state(), DuctileJ2Conditions {
            duration_s: conditions().duration_s, temperature_k: 300.0,
        }).unwrap();
        let force = VOLUME * expected.first_piola_pa[0];
        let loaded = solid.advance(
            &supports, &[(3, force)], [0.0; 3], conditions(),
            IncrementSettings::default(), controls(), cx,
        ).unwrap();
        near(loaded.solid.displacement_m()[1][0], target[0] - 1.0);
        near(loaded.solid.element_states()[0].native_state().alpha, expected.state.native_state().alpha);
        let old = &loaded.solid.element_states()[0];
        assert!(material.native().validate_committed(old.deformation(), old.native_state(), 1e-8).is_err());
        material.validate_committed(old, 1e-8).unwrap();
        let duration = 0.05;
        let proposed_hold = material.stress(old.deformation(), old, DuctileJ2Conditions {
            duration_s: duration, temperature_k: 300.0,
        }).unwrap();
        // A relaxation proposal is NOT the accepted force state. Checking
        // initial equilibrium through a trial would incorrectly refuse here.
        assert!((VOLUME * proposed_hold[0] - force).abs() > 1e-4);
        let settings = FiniteDuctilePathControls {
            initial_fraction: 0.25,
            minimum_fraction: 1.0 / 4096.0,
            maximum_fraction: 0.25,
            maximum_accepted_steps: 32,
            maximum_trials: 64,
            maximum_damage_increment: 0.1,
            initial_equilibrium_tolerance_n: 1e-7,
            work: controls(),
        };
        let path = loaded.solid.advance_path(
            &supports, &[(3, [force, force])], [[0.0; 3]; 2],
            DuctileStepConditions { duration_s: duration, ..conditions() },
            IncrementSettings::default(), settings, cx,
        ).unwrap();
        assert_eq!(path.records.len(), 4);
        near(path.records.iter().map(|record| record.duration_s).sum(), duration);
        assert!(path.endpoint.solid.displacement_m()[1][0] > loaded.solid.displacement_m()[1][0]);
        assert!(path.endpoint.solid.element_states()[0].damage() > old.damage());
        assert!(path.energy.material.plastic_dissipation_j > 0.0);
        assert!(path.energy.material.damage_release_j > 0.0);
        near(path.energy.balance_error_j, 0.0);

        // Replaying the published physical partition advances each history
        // exactly once; Newton trials and initial inspection consume no time.
        let mut replay = loaded.solid;
        for record in &path.records {
            replay = replay.advance(
                &supports, &[(3, force)], [0.0; 3],
                DuctileStepConditions { duration_s: record.duration_s, ..conditions() },
                IncrementSettings::default(), controls(), cx,
            ).unwrap().solid;
        }
        assert_eq!(path.endpoint.solid.element_states(), replay.element_states());
        assert_eq!(path.endpoint.solid.displacement_m(), replay.displacement_m());
    });
}

#[test]
fn moving_and_held_specimens_keep_rate_heat_separate_from_reference_inertia() {
    with_cx(|_, cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = kinetic_model(bands[0].band_width_m());
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry, &materials, &bands, TetAssemblyBudget::standard(), cx,
        ).unwrap();
        let motion = FiniteDuctileTetMotion::try_new(solid, &[12.0], &[[0.0; 3]; 4], cx).unwrap();
        let motion_controls = FiniteDuctileMotionControls {
            material: controls(), maximum_numerical_kinetic_loss_j: 1.0,
        };
        let supports = affine(&NODES, stretch(0.02), [0.0; 3], 0);
        let first = motion.advance(
            &supports, &[], [0.0; 3], conditions(),
            IncrementSettings::default(), motion_controls, cx,
        ).unwrap();
        let f = first.motion.solid().element_deformation()[0];
        let point_conditions = DuctileJ2Conditions {
            duration_s: conditions().duration_s, temperature_k: 300.0,
        };
        let expected = material.trial(&f, &material.initial_state(), point_conditions).unwrap();
        assert_eq!(&expected.state, &first.motion.solid().element_states()[0]);
        near(first.material_energy.plastic_dissipation_j, VOLUME * expected.energy.plastic_dissipation_j_m3);
        near(first.material_energy.damage_release_j, VOLUME * expected.energy.damage_release_j_m3);
        assert!(first.inertia.numerical_loss_j > 0.0);
        near(first.balance_error_j, 0.0);
        let retained = first.motion.solid().element_states()[0].clone();
        let held = first.motion.advance(
            &supports, &[], [0.0; 3], conditions(),
            IncrementSettings::default(), motion_controls, cx,
        ).unwrap();
        let relaxed = material.trial(&f, &retained, point_conditions).unwrap();
        assert_eq!(&relaxed.state, &held.motion.solid().element_states()[0]);
        assert!(held.motion.solid().element_states()[0].damage() > retained.damage());
        assert!(held.material_energy.plastic_dissipation_j > 0.0);
        assert!(held.material_energy.damage_release_j > 0.0);
        near(held.material_energy.internal_endpoint_work_j, 0.0);
        near(held.material_energy.plastic_dissipation_j, VOLUME * relaxed.energy.plastic_dissipation_j_m3);
        near(held.inertia.kinetic_change_j + held.inertia.numerical_loss_j, 0.0);
        near(held.balance_error_j, 0.0);
        assert_eq!(first.motion.solid().element_states()[0], retained);
        assert_eq!(held.motion.velocity_mps(), &[[0.0; 3]; 4]);
        near(held.motion.reference_mass_kg(), 2.0);
    });
}
