//! Native inertia and irreversible tensile failure share one source-bound owner.
//! These synthetic coupons check the stated fixed-plane small-strain model.
use super::*;

use fs_material::brittle::source::dynamics::resolve_brittle_dynamic_material;
use fs_material::state_point::DENSITY_PROPERTY;
use fs_qty::Density;
use fs_solid::brittle::solid::dynamics::{
    BrittleTetDynamicControls, BrittleTetDynamicError, BrittleTetMotion,
};

fn dynamic_card(density_kg_m3: f64) -> MaterialCard {
    // At 300 K this existing source fixture has E = YOUNG, ft = STRENGTH
    // and nu = 0, giving an independent one-coordinate tensile oracle.
    let base = source_card(STRENGTH, "dynamic reference inventory", true);
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            value: PropertyValue::Scalar {
                value: density_kg_m3,
                dims: Density::DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic reference density in kg/m^3".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn dynamic_material(
    source: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
) -> ResolvedBrittleCrackBand {
    let query = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1.0, 0.0, 0.0]).unwrap();
    resolve_brittle_dynamic_material(
        source,
        &query,
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.03, 1e-7),
    )
    .unwrap()
}

fn dynamic_controls() -> BrittleTetDynamicControls {
    BrittleTetDynamicControls {
        work: controls(),
        kinetic_work_tolerance_j: 1e-10,
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}

fn uniform_x_force(total_n: f64) -> Vec<(usize, f64)> {
    (0..4).map(|node| (3 * node, total_n / 4.0)).collect()
}

#[test]
fn source_mass_preserves_free_momentum_and_scales_force_acceleration() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let initial_velocity = [[0.2, -0.1, 0.05]; 4];
        let dt = 1e-6;
        let total_force = 2.0;
        for density in [1000.0, 2000.0] {
            let source = dynamic_card(density);
            let materials = [dynamic_material(&source, &geometry)];
            let solid =
                BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                    .unwrap();
            let mass = density * solid.band_geometry()[0].reference_volume_m3();
            let initial_kinetic =
                0.5 * mass * (0.2_f64.powi(2) + 0.1_f64.powi(2) + 0.05_f64.powi(2));
            let old =
                BrittleTetMotion::try_from_solid(&solid, &initial_velocity, 1.25, 4, cx).unwrap();
            assert_eq!(old.reference_density_kg_m3(), &[density]);
            close(old.reference_mass_kg(), mass, 1e-13);
            close(old.kinetic_energy_j(), initial_kinetic, 1e-13);
            assert!(std::ptr::eq(
                old.solid().materials()[0].law(),
                materials[0].law()
            ));

            let coast = old
                .advance(
                    dt,
                    &[],
                    &[],
                    IncrementSettings::default(),
                    dynamic_controls(),
                    cx,
                )
                .unwrap();
            for (velocity, displacement) in coast
                .motion
                .velocity_mps()
                .iter()
                .zip(&coast.motion.solid().solution().displacement_m)
            {
                for axis in 0..3 {
                    close(velocity[axis], initial_velocity[0][axis], 2e-8);
                    close(displacement[axis], dt * initial_velocity[0][axis], 2e-8);
                }
            }
            close(coast.motion.kinetic_energy_j(), initial_kinetic, 2e-8);
            assert!(coast.energy.inertia.numerical_loss_j < 1e-18);
            assert!(coast.energy.fracture_dissipation_j.abs() < 1e-12);
            assert!(coast.energy.total_balance_error_j.abs() < 1e-10);

            let forced = coast
                .motion
                .advance(
                    dt,
                    &[],
                    &uniform_x_force(total_force),
                    IncrementSettings::default(),
                    dynamic_controls(),
                    cx,
                )
                .unwrap();
            let change = dt * total_force / mass;
            let expected_velocity = [initial_velocity[0][0] + change, -0.1, 0.05];
            for (velocity, displacement) in forced
                .motion
                .velocity_mps()
                .iter()
                .zip(&forced.motion.solid().solution().displacement_m)
            {
                for axis in 0..3 {
                    close(velocity[axis], expected_velocity[axis], 2e-8);
                    close(
                        displacement[axis],
                        dt * (initial_velocity[0][axis] + expected_velocity[axis]),
                        2e-8,
                    );
                }
            }
            close(
                forced.energy.dead_load_endpoint_work_j,
                total_force * dt * expected_velocity[0],
                2e-8,
            );
            close(
                forced.energy.inertia.numerical_loss_j,
                0.5 * mass * change.powi(2),
                2e-7,
            );
            assert!(
                forced.motion.solid().solution().element_strain[0]
                    .iter()
                    .all(|strain| strain.abs() < 1e-12)
            );
            assert_eq!(
                forced.motion.solid().solution().element_states[0].damage(),
                0.0
            );
            assert!(forced.energy.fracture_dissipation_j.abs() < 1e-12);
            assert!(forced.energy.total_balance_error_j.abs() < 1e-10);
            assert_eq!(forced.motion.accepted_steps(), 2);
            assert_eq!(forced.motion.remaining_steps(), 2);
            close(forced.motion.time_s(), 1.25 + 2.0 * dt, 1e-14);
            assert_eq!(old.velocity_mps(), &initial_velocity);
            assert_eq!(old.accepted_steps(), 0);
            assert_eq!(old.time_s(), 1.25);
            assert_eq!(old.solid().solution().displacement_m, vec![[0.0; 3]; 4]);
        }
    });
}

#[test]
fn native_inertia_crosses_the_static_tensile_peak_and_keeps_fracture_history() {
    with_cx(false, |cx| {
        let length = 0.03;
        let nodes = points(length);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let density = 2000.0;
        let source = dynamic_card(density);
        let materials = [dynamic_material(&source, &geometry)];
        let solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let band = &solid.band_geometry()[0];
        let volume = band.reference_volume_m3();
        let area = volume / length;
        let onset = STRENGTH / YOUNG;
        let failure = 2.0 * FRACTURE_WORK / (STRENGTH * band.band_width_m());
        let force = 1.1 * STRENGTH * area;
        let supports: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let settings = IncrementSettings {
            max_newton_iterations: 8,
            ..IncrementSettings::default()
        };
        assert!(
            solid
                .advance(&supports, &[(3, force)], settings, controls(), cx)
                .is_err()
        );

        let dt = 1e-5;
        // The free coordinate sees M33 = rho*V/10, not lumped nodal mass.
        // Its softening tangent is negative, but inertia makes the total
        // tangent positive and gives a unique admissible dynamic endpoint.
        let modal_mass = density * volume / 10.0;
        let intercept = area * STRENGTH * failure / (failure - onset);
        let softening_stiffness = -area * STRENGTH / (length * (failure - onset));
        let displacement = (force - intercept) / (modal_mass / (dt * dt) + softening_stiffness);
        let strain = displacement / length;
        assert!(strain > onset && strain < failure);
        let stress = STRENGTH * (failure - strain) / (failure - onset);
        let material_work = volume
            * (0.5 * STRENGTH * onset
                + STRENGTH / (failure - onset)
                    * (failure * (strain - onset) - 0.5 * (strain * strain - onset * onset)));
        let stored = volume * 0.5 * stress * strain;
        let physical_loss = material_work - stored;
        let kinetic = 0.5 * modal_mass * (displacement / dt).powi(2);
        let old = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 4, cx).unwrap();
        let loaded = old
            .advance(
                dt,
                &supports,
                &[(3, force)],
                settings,
                dynamic_controls(),
                cx,
            )
            .unwrap();
        close(
            loaded.motion.solid().solution().displacement_m[1][0],
            displacement,
            2e-8,
        );
        close(loaded.motion.velocity_mps()[1][0], displacement / dt, 2e-8);
        close(
            loaded.motion.solid().solution().element_strain[0][0],
            strain,
            2e-8,
        );
        close(
            loaded.motion.solid().solution().element_stress[0][0],
            stress,
            2e-8,
        );
        assert!(loaded.motion.solid().solution().element_states[0].damage() > 0.0);
        close(loaded.energy.stored_after_j, stored, 2e-8);
        close(loaded.energy.material_work_j, material_work, 2e-8);
        close(loaded.energy.fracture_dissipation_j, physical_loss, 2e-8);
        close(loaded.energy.inertia.kinetic_after_j, kinetic, 2e-8);
        close(loaded.energy.inertia.numerical_loss_j, kinetic, 2e-8);
        close(
            loaded.energy.internal_endpoint_work_j,
            volume * stress * strain,
            2e-8,
        );
        close(
            loaded.energy.dead_load_endpoint_work_j,
            force * displacement,
            2e-8,
        );
        close(
            loaded.energy.constitutive_endpoint_defect_j,
            volume * stress * strain - material_work,
            2e-8,
        );
        assert_eq!(loaded.energy.support_endpoint_work_j, 0.0);
        assert!(loaded.energy.total_balance_error_j.abs() < 1e-9);
        assert!(loaded.energy.fracture_dissipation_j > 0.0);
        assert!(loaded.energy.inertia.numerical_loss_j > 0.0);
        close(
            loaded.element_energy[0].dissipated * volume,
            physical_loss,
            2e-8,
        );
        assert_eq!(old.solid().solution().element_states[0].damage(), 0.0);

        // Fix the attained displacement to stop motion without new opening.
        // Lost kinetic energy is numerical damping; damage cannot heal or
        // acquire extra fracture work from the velocity change alone.
        let hold: Vec<_> = loaded
            .motion
            .solid()
            .solution()
            .displacement_m
            .as_flattened()
            .iter()
            .copied()
            .enumerate()
            .collect();
        let stopped = loaded
            .motion
            .advance(dt, &hold, &[], settings, dynamic_controls(), cx)
            .unwrap();
        assert_eq!(stopped.motion.velocity_mps(), &[[0.0; 3]; 4]);
        close(stopped.energy.inertia.numerical_loss_j, kinetic, 2e-8);
        assert!(stopped.energy.fracture_dissipation_j.abs() < 1e-12);
        assert_eq!(
            stopped.motion.solid().solution().element_states,
            loaded.motion.solid().solution().element_states
        );
        close(
            stopped.motion.solid().fracture_dissipation_j(),
            physical_loss,
            2e-8,
        );
        close(stopped.motion.solid().stored_energy_j(), stored, 2e-8);
        close(stopped.motion.reference_mass_kg(), density * volume, 1e-13);
        assert!(stopped.energy.total_balance_error_j.abs() < 1e-9);
        let resumed = BrittleTetMotion::try_from_solid(
            stopped.motion.solid(),
            stopped.motion.velocity_mps(),
            stopped.motion.time_s(),
            2,
            cx,
        )
        .unwrap();
        assert_eq!(
            resumed.solid().solution(),
            stopped.motion.solid().solution()
        );
        close(
            resumed.solid().fracture_dissipation_j(),
            physical_loss,
            2e-8,
        );
    });
}

#[test]
fn moving_support_work_includes_consistent_inertia_and_applied_nodal_loads() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = dynamic_card(2000.0);
        let materials = [dynamic_material(&source, &geometry)];
        let solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        let old = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 1, cx).unwrap();
        let dt = 1e-4;
        let travel = 1e-5;
        let force = 3.0;
        let supports: Vec<_> = (0..12)
            .map(|dof| (dof, if dof % 3 == 0 { travel } else { 0.0 }))
            .collect();
        let step = old
            .advance(
                dt,
                &supports,
                &uniform_x_force(force),
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            )
            .unwrap();
        let kinetic = 0.5 * old.reference_mass_kg() * (travel / dt).powi(2);
        close(step.energy.inertia.kinetic_after_j, kinetic, 1e-12);
        close(step.energy.inertia.numerical_loss_j, kinetic, 1e-12);
        close(step.energy.inertia.inertial_work_j, 2.0 * kinetic, 1e-12);
        close(step.energy.dead_load_endpoint_work_j, force * travel, 1e-12);
        close(
            step.energy.support_endpoint_work_j,
            2.0 * kinetic - force * travel,
            1e-12,
        );
        assert!(step.energy.support_endpoint_work_j < 0.0);
        assert!(step.energy.internal_endpoint_work_j.abs() < 1e-12);
        assert!(step.energy.fracture_dissipation_j.abs() < 1e-12);
        assert!(step.energy.total_balance_error_j.abs() < 1e-12);
        assert!(
            step.motion.solid().solution().element_strain[0]
                .iter()
                .all(|strain| strain.abs() < 1e-14)
        );
        for velocity in step.motion.velocity_mps() {
            close(velocity[0], travel / dt, 1e-12);
            assert_eq!(velocity[1..], [0.0; 2]);
        }
    });
}

#[test]
fn dynamic_admission_late_energy_refusal_and_cancellation_preserve_the_accepted_owner() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let no_density_card = card(STRENGTH, "mechanical-only source");
        let no_density_materials = [material(&no_density_card, &geometry, 0)];
        let no_density_solid = BrittleTetSolid::try_new(
            &geometry,
            &no_density_materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            BrittleTetMotion::try_from_solid(&no_density_solid, &[[0.0; 3]; 4], 0.0, 1, cx),
            Err(BrittleTetDynamicError::Source { element: 0, .. })
        ));

        let source = dynamic_card(2000.0);
        let materials = [dynamic_material(&source, &geometry)];
        let solid =
            BrittleTetSolid::try_new(&geometry, &materials, TetAssemblyBudget::standard(), cx)
                .unwrap();
        for time in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], time, 1, cx).is_err());
        }
        assert!(BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 3], 0.0, 1, cx).is_err());
        assert!(BrittleTetMotion::try_from_solid(&solid, &[[f64::NAN; 3]; 4], 0.0, 1, cx).is_err());
        assert!(BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 0, cx).is_err());
        let old = BrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 4], 0.0, 1, cx).unwrap();
        let snapshot = old.solid().solution().clone();
        let forces = uniform_x_force(2.0);
        let dt = 1e-6;
        for invalid_dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                old.advance(
                    invalid_dt,
                    &[],
                    &forces,
                    IncrementSettings::default(),
                    dynamic_controls(),
                    cx
                )
                .is_err()
            );
        }
        assert!(
            old.advance(
                dt,
                &[(12, 0.0)],
                &[],
                IncrementSettings::default(),
                dynamic_controls(),
                cx
            )
            .is_err()
        );
        assert!(
            old.advance(
                dt,
                &[],
                &[(12, 1.0)],
                IncrementSettings::default(),
                dynamic_controls(),
                cx
            )
            .is_err()
        );
        let strict = BrittleTetDynamicControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..dynamic_controls()
        };
        assert!(matches!(
            old.advance(dt, &[], &forces, IncrementSettings::default(), strict, cx),
            Err(BrittleTetDynamicError::Solid(BrittleTetError::EnergyGate {
                account: "backward-Euler numerical kinetic loss",
                ..
            }))
        ));
        with_cx(true, |cancelled| {
            assert!(matches!(
                old.advance(
                    dt,
                    &[],
                    &forces,
                    IncrementSettings::default(),
                    dynamic_controls(),
                    cancelled
                ),
                Err(BrittleTetDynamicError::Solid(BrittleTetError::Cancelled))
            ));
        });
        assert_eq!(old.solid().solution(), &snapshot);
        assert_eq!(old.velocity_mps(), &[[0.0; 3]; 4]);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.accepted_steps(), 0);
        assert_eq!(old.kinetic_energy_j(), 0.0);
        assert_eq!(old.solid().fracture_dissipation_j(), 0.0);

        let accepted = old
            .advance(
                dt,
                &[],
                &forces,
                IncrementSettings::default(),
                dynamic_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(accepted.motion.accepted_steps(), 1);
        assert_eq!(accepted.motion.remaining_steps(), 0);
        assert_eq!(accepted.motion.time_s(), dt);
        let accepted_solution = accepted.motion.solid().solution().clone();
        let accepted_velocity = accepted.motion.velocity_mps().to_vec();
        assert!(matches!(
            accepted.motion.advance(
                dt,
                &[],
                &forces,
                IncrementSettings::default(),
                dynamic_controls(),
                cx
            ),
            Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: 1,
                maximum_steps: 1
            })
        ));
        assert_eq!(accepted.motion.solid().solution(), &accepted_solution);
        assert_eq!(accepted.motion.velocity_mps(), &accepted_velocity);
        assert_eq!(accepted.motion.time_s(), dt);
        assert_eq!(old.solid().solution(), &snapshot);
    });
}

#[path = "dynamics/temperature.rs"]
mod temperature;

#[path = "dynamics/time.rs"]
mod time;
