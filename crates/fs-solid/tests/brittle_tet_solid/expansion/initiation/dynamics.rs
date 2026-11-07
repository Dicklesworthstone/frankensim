//! Native midpoint thermal motion: independent scalar force/work and P1 mass.
//! Manufactured coefficients verify this split model, not a material population.

use super::*;
use fs_material::brittle::source::expansion::dynamics::resolve_thermal_brittle_dynamic_material;
use fs_material::state_point::DENSITY_PROPERTY;
use fs_qty::Density;
use fs_solid::brittle::solid::dynamics::BrittleTetDynamicError;
use fs_solid::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use fs_solid::brittle::solid::expansion::dynamics::midpoint::ThermalBrittleTetMidpointEnergy;
use fs_solid::brittle::solid::expansion::initiation::dynamics::{
    ThermalTensileElasticTetMidpointStep, ThermalTensileElasticTetMotion,
    ThermalTensileTetDynamicError, ThermalTensileTetMidpointProbe,
};

fn moving_card() -> MaterialCard {
    let base = weakening_card();
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 2000.0), (400.0, 1000.0)],
                dims: Density::DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                400.0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic current density; original reference inventory stays fixed"
                    .into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn moving_reference(
    source: &MaterialCard,
    geometry: &TetLoadGeometry3<'_>,
) -> ResolvedThermalBrittle {
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1.0, 0.0, 0.0]).unwrap();
    resolve_thermal_brittle_dynamic_material(
        source,
        &temperature(300.0),
        &temperature(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.03, 1e-7),
    )
    .unwrap()
}

fn midpoint_controls() -> BrittleTetMidpointControls {
    BrittleTetMidpointControls {
        material_energy_tolerance_j: 1e-9,
        equilibrium_work_tolerance_j: 1e-9,
        kinetic_work_tolerance_j: 1e-9,
        total_energy_tolerance_j: 1e-9,
        difference_step_strain: Some(1e-7),
    }
}

fn moving_elastic<'a>(
    probe: ThermalTensileTetMidpointProbe<'a>,
) -> ThermalTensileElasticTetMidpointStep<'a> {
    match probe {
        ThermalTensileTetMidpointProbe::Elastic(step) => step,
        other => panic!("expected elastic thermal motion, got {other:?}"),
    }
}

fn balances(energy: &ThermalBrittleTetMidpointEnergy) {
    for defect in [
        energy.strain_split_error_j,
        energy.material_balance_error_j,
        energy.mean_work_error_j,
        energy.equilibrium_work_residual_j,
        energy.inertia.balance_error_j,
        energy.total_balance_error_j,
    ] {
        assert!(defect.abs() < 1e-9, "{energy:?}");
    }
    assert_eq!(energy.inertia.numerical_loss_j, 0.0);
}

/// M_ij = rho V (1 + delta_ij)/20 for one linear tetrahedron.
fn p1_kinetic(velocity: &[[f64; 3]], density: f64, volume: f64) -> f64 {
    let mut sum = 0.0;
    for axis in 0..3 {
        let total: f64 = velocity.iter().map(|v| v[axis]).sum();
        let squares: f64 = velocity.iter().map(|v| v[axis] * v[axis]).sum();
        sum += total * total + squares;
    }
    density * volume * sum / 40.0
}

fn free_y_supports(nodes: &[[f64; 3]], free: f64) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 7)
        .map(|dof| (dof, free * nodes[dof / 3][dof % 3]))
        .collect()
}

#[test]
fn thermal_midpoint_elastic_motion_and_cooling_onset_keep_original_mass_and_exact_work() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let card = moving_card();
        let reference = [moving_reference(&card, &geometry)];
        let warm = [reference[0].at_temperature(&card, 304.0).unwrap()];
        let virgin = ThermalTensileElasticTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let h = 1e-5;
        let length = 0.02;
        let volume = geometry.element_volumes_m3()[0];
        let mass = 2000.0 * volume / 10.0;
        let stiffness = YOUNG * volume / (length * length);
        let inertia = 2.0 * mass / (h * h);
        let peak_strain = STRENGTH / YOUNG;
        let e0 = 0.2 * peak_strain;
        let q0 = e0 * length;
        let supports = free_y_supports(&nodes, 0.0);
        let preload = elastic(
            virgin
                .probe(
                    &reference,
                    &supports,
                    &[(7, stiffness * q0)],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let mut velocity = [[0.0; 3]; 4];
        velocity[2][1] = 0.02;
        let initial = ThermalTensileElasticTetMotion::try_from_solid(
            &preload.solid,
            &velocity,
            2.0,
            2,
            initiation_settings(),
            cx,
        )
        .unwrap();
        let free = 1e-5 * 4.0 + 0.5e-7 * 4.0 * 4.0;
        close(warm[0].free_linear_strain(), free, 2e-12);
        let e1 = 0.3 * peak_strain;
        let q1 = (free + e1) * length;
        let next_young = 0.98 * YOUNG;
        let mean_stress = 0.5 * next_young * (e0 + e1);
        let force = mean_stress * volume / length + inertia * (q1 - q0 - h * velocity[2][1]);
        let warm_supports = free_y_supports(&nodes, free);
        let accepted = moving_elastic(
            initial
                .probe_midpoint(
                    &warm,
                    h,
                    &warm_supports,
                    &[(7, force)],
                    IncrementSettings::default(),
                    midpoint_controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let motion = &accepted.motion;
        close(motion.solution().displacement_m[2][1], q1, 2e-8);
        close(
            motion.velocity_mps()[2][1],
            2.0 * (q1 - q0) / h - velocity[2][1],
            2e-8,
        );
        close(accepted.mean_element_stress_pa[0][1], mean_stress, 2e-8);
        close(
            motion.solution().element_stress[0][1],
            next_young * e1,
            2e-8,
        );
        close(
            accepted.energy.parameter_work_j,
            0.5 * (next_young - YOUNG) * e0 * e0 * volume,
            2e-8,
        );
        close(
            accepted.energy.mechanical_work_j,
            mean_stress * (q1 - q0) * volume / length,
            2e-8,
        );
        close(
            accepted.energy.eigenstrain_work_j,
            -mean_stress * free * volume,
            2e-8,
        );
        close(
            motion.kinetic_energy_j(),
            p1_kinetic(motion.velocity_mps(), 2000.0, volume),
            2e-8,
        );
        assert_eq!(motion.reference_density_kg_m3(), &[2000.0]);
        close(
            motion.materials()[0]
                .band()
                .reference_density_kg_m3()
                .unwrap(),
            1960.0,
            2e-12,
        );
        assert_eq!(motion.accepted_steps(), 1);
        assert_eq!(motion.remaining_steps(), 1);
        assert_eq!(accepted.energy.fracture_dissipation_j, 0.0);
        balances(&accepted.energy);

        let before = motion.solution().clone();
        // An out-of-domain native solve must leave the warm motion reusable.
        assert!(
            motion
                .probe_midpoint(
                    &reference,
                    h,
                    &supports,
                    &[(7, 1e9)],
                    IncrementSettings::default(),
                    midpoint_controls(),
                    None,
                    cx,
                )
                .is_err()
        );
        assert_eq!(motion.solution(), &before);
        assert_eq!(motion.accepted_steps(), 1);

        let q_peak = peak_strain * length;
        let mechanical_start = e1 * length;
        let target = 2.0 * q_peak;
        let band = TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap();
        let q_failure = 2.0 * FRACTURE_WORK * length / (STRENGTH * band.band_width_m());
        let peak_force = STRENGTH * volume / length;
        // Independent exact integral of the next-temperature bilinear traction.
        let strain_work = 0.5 * stiffness * (q_peak * q_peak - mechanical_start * mechanical_start)
            + peak_force / (q_failure - q_peak)
                * (q_failure * (target - q_peak) - 0.5 * (target * target - q_peak * q_peak));
        let mean_material_force = strain_work / (target - mechanical_start);
        let total_delta = target - q1;
        let force = mean_material_force + inertia * (total_delta - h * motion.velocity_mps()[2][1]);
        let forces = [(7, 0.6 * force), (7, 0.4 * force)];
        let prediction = match motion
            .probe_midpoint(
                &reference,
                h,
                &supports,
                &forces,
                IncrementSettings::default(),
                midpoint_controls(),
                None,
                cx,
            )
            .unwrap()
        {
            ThermalTensileTetMidpointProbe::Onset(onset) => onset,
            other => panic!("expected cooling/load tensile onset, got {other:?}"),
        };
        assert_eq!(prediction.initial_solution(), &before);
        assert_eq!(prediction.mean_nodal_forces_n(), &forces);
        assert_eq!(prediction.prescribed_end_m(), supports.as_slice());
        assert_eq!(
            prediction.physical_interval_s(),
            [motion.time_s(), motion.time_s() + h]
        );
        close(
            prediction.events()[0].as_ref().unwrap().native().normal()[1].abs(),
            1.0,
            2e-12,
        );
        let predictor = (force + inertia * (q1 + h * motion.velocity_mps()[2][1])
            - 0.5 * stiffness * mechanical_start)
            / (inertia + 0.5 * stiffness);
        close(
            prediction.elastic_predictor().displacement_m[2][1],
            predictor,
            2e-8,
        );
        close(
            prediction.elastic_predictor().element_stress[0][1],
            YOUNG * predictor / length,
            2e-8,
        );
        assert!(predictor > q_peak);
        assert!(0.5 * (mechanical_start + predictor) < q_peak);
        let bands = prediction.resolve_bands(cx).unwrap();
        let handoff = prediction.handoff(&bands.initial, 1e-9, cx).unwrap();
        assert_eq!(handoff.motion.velocity_mps(), motion.velocity_mps());
        assert_eq!(handoff.motion.kinetic_energy_j(), motion.kinetic_energy_j());
        assert_eq!(
            handoff.motion.reference_mass_kg(),
            motion.reference_mass_kg()
        );
        assert_eq!(handoff.motion.reference_density_kg_m3(), &[2000.0]);
        close(
            handoff.motion.solid().materials()[0]
                .band()
                .reference_density_kg_m3()
                .unwrap(),
            1960.0,
            2e-12,
        );
        assert_eq!(
            handoff.motion.solid().solution().displacement_m,
            before.displacement_m
        );
        assert_eq!(
            handoff.motion.solid().solution().element_strain,
            before.element_strain
        );
        assert_eq!(handoff.motion.time_s(), motion.time_s());
        assert_eq!(handoff.motion.accepted_steps(), 1);
        assert_eq!(handoff.motion.remaining_steps(), 1);
        assert_eq!(handoff.motion.solid().fracture_dissipation_j(), 0.0);
        with_cx(true, |cancelled| {
            assert!(
                prediction
                    .replay(&bands, midpoint_controls(), cancelled)
                    .is_err()
            );
        });
        let replay = prediction.replay(&bands, midpoint_controls(), cx).unwrap();
        close(
            replay.motion.solid().solution().displacement_m[2][1],
            target,
            2e-7,
        );
        close(
            replay.motion.velocity_mps()[2][1],
            2.0 * total_delta / h - motion.velocity_mps()[2][1],
            2e-7,
        );
        close(replay.energy.constitutive_strain_work_j, strain_work, 2e-7);
        close(
            replay.energy.mechanical_work_j,
            mean_material_force * total_delta,
            2e-7,
        );
        close(
            replay.energy.eigenstrain_work_j,
            mean_material_force * free * length,
            2e-7,
        );
        close(
            replay.energy.parameter_work_j,
            0.5 * (YOUNG - next_young) * e1 * e1 * volume,
            2e-7,
        );
        close(
            replay.motion.kinetic_energy_j(),
            p1_kinetic(replay.motion.velocity_mps(), 2000.0, volume),
            2e-8,
        );
        assert!(replay.energy.fracture_dissipation_j > 0.0);
        assert_eq!(replay.energy.parameter_fracture_dissipation_j, 0.0);
        assert_eq!(replay.motion.accepted_steps(), 2);
        assert_eq!(replay.motion.remaining_steps(), 0);
        assert_eq!(replay.motion.time_s(), motion.time_s() + h);
        balances(&replay.energy);
        let repeated = prediction.replay(&bands, midpoint_controls(), cx).unwrap();
        assert_eq!(
            repeated.motion.solid().solution(),
            replay.motion.solid().solution()
        );
        assert_eq!(repeated.energy, replay.energy);
        assert_eq!(motion.solution(), &before);
        assert!(matches!(
            replay.motion.advance_midpoint(
                &bands.next,
                h,
                &supports,
                &forces,
                IncrementSettings::default(),
                midpoint_controls(),
                cx,
            ),
            Err(BrittleTetDynamicError::StepBudget { .. })
        ));
    });
}

#[test]
fn thermal_midpoint_coefficient_activation_preserves_old_translation_and_charges_native_work() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let card = moving_card();
        let reference = [moving_reference(&card, &geometry)];
        let hot = [reference[0].at_temperature(&card, 400.0).unwrap()];
        let virgin = ThermalTensileElasticTetSolid::try_new(
            &geometry,
            &reference,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let strain = 0.75 * STRENGTH / YOUNG;
        let targets = affine_targets(&nodes, [0.0, strain, 0.0]);
        let solid = elastic(
            virgin
                .probe(
                    &reference,
                    &targets,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    controls(),
                    None,
                    cx,
                )
                .unwrap(),
        )
        .solid;
        let velocity = [0.1, -0.2, 0.05];
        let initial = ThermalTensileElasticTetMotion::try_from_solid(
            &solid,
            &[velocity; 4],
            3.0,
            2,
            initiation_settings(),
            cx,
        )
        .unwrap();
        let h = 1e-4;
        let translated: Vec<_> = targets
            .iter()
            .map(|&(dof, u)| (dof, u + h * velocity[dof % 3]))
            .collect();
        let accepted = moving_elastic(
            initial
                .probe_midpoint(
                    &reference,
                    h,
                    &translated,
                    &[],
                    IncrementSettings::default(),
                    midpoint_controls(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let motion = &accepted.motion;
        let before = motion.solution().clone();
        let free = 1e-5 * 100.0 + 0.5e-7 * 100.0 * 100.0;
        close(hot[0].free_linear_strain(), free, 2e-12);
        let expanded: Vec<_> = translated
            .iter()
            .map(|&(dof, u)| {
                (
                    dof,
                    u + h * velocity[dof % 3] + free * nodes[dof / 3][dof % 3],
                )
            })
            .collect();
        let forces = [(0, 1.0)];
        // Early activation must not hide malformed forces behind a later replay.
        assert!(matches!(
            motion.probe_midpoint(
                &hot,
                h,
                &expanded,
                &[(12, 1.0)],
                IncrementSettings::default(),
                midpoint_controls(),
                None,
                cx,
            ),
            Err(ThermalTensileTetDynamicError::Initiation(
                ThermalTensileTetInitiationError::Invalid(_)
            ))
        ));
        let activation = match motion
            .probe_midpoint(
                &hot,
                h,
                &expanded,
                &forces,
                IncrementSettings::default(),
                midpoint_controls(),
                None,
                cx,
            )
            .unwrap()
        {
            ThermalTensileTetMidpointProbe::Activation(activation) => activation,
            other => panic!("expected held-strain coefficient activation, got {other:?}"),
        };
        assert_eq!(activation.initial_solution(), &before);
        assert_eq!(activation.mean_nodal_forces_n(), &forces);
        assert_eq!(activation.prescribed_end_m(), expanded.as_slice());
        assert_eq!(
            activation.physical_interval_s(),
            [motion.time_s(), motion.time_s() + h]
        );
        close(
            activation.activations()[0]
                .as_ref()
                .unwrap()
                .native()
                .normal()[1]
                .abs(),
            1.0,
            2e-12,
        );
        let bands = activation.resolve_bands(cx).unwrap();
        let handoff = activation.handoff(&bands.initial, 1e-9, cx).unwrap();
        assert_eq!(handoff.motion.velocity_mps(), motion.velocity_mps());
        assert_eq!(handoff.motion.kinetic_energy_j(), motion.kinetic_energy_j());
        assert_eq!(
            handoff.motion.reference_mass_kg(),
            motion.reference_mass_kg()
        );
        assert_eq!(handoff.motion.time_s(), motion.time_s());
        assert_eq!(handoff.motion.accepted_steps(), 1);
        assert_eq!(handoff.motion.remaining_steps(), 1);
        assert_eq!(handoff.motion.solid().fracture_dissipation_j(), 0.0);
        let replay = activation.replay(&bands, midpoint_controls(), cx).unwrap();
        let volume = geometry.element_volumes_m3()[0];
        let hot_young = 0.5 * YOUNG;
        let hot_strength = 0.25 * STRENGTH;
        let failure = 2.0 * FRACTURE_WORK
            / (hot_strength * replay.motion.solid().band_geometry()[0].band_width_m());
        let stress = hot_strength * (failure - strain) / (failure - hot_strength / hot_young);
        let induced_fracture = 0.5 * (hot_young * strain * strain - stress * strain) * volume;
        close(
            replay.energy.parameter_work_j,
            0.5 * (hot_young - YOUNG) * strain * strain * volume,
            2e-7,
        );
        close(
            replay.energy.parameter_fracture_dissipation_j,
            induced_fracture,
            2e-7,
        );
        close(replay.energy.fracture_dissipation_j, induced_fracture, 2e-7);
        close(
            replay.energy.stored_after_j,
            0.5 * stress * strain * volume,
            2e-7,
        );
        close(
            replay.energy.mechanical_work_j,
            stress * free * volume,
            2e-7,
        );
        close(
            replay.energy.eigenstrain_work_j,
            -stress * free * volume,
            2e-7,
        );
        assert!(replay.energy.constitutive_strain_work_j.abs() < 1e-10);
        let expected_velocity: Vec<[f64; 3]> = nodes
            .iter()
            .map(|xyz| core::array::from_fn(|axis| velocity[axis] + 2.0 * free * xyz[axis] / h))
            .collect();
        for (actual, expected) in replay.motion.velocity_mps().iter().zip(&expected_velocity) {
            for axis in 0..3 {
                close(actual[axis], expected[axis], 2e-8);
            }
        }
        close(
            replay.motion.kinetic_energy_j(),
            p1_kinetic(&expected_velocity, 2000.0, volume),
            2e-8,
        );
        assert_eq!(replay.motion.reference_density_kg_m3(), &[2000.0]);
        assert_eq!(
            replay.motion.solid().materials()[0]
                .band()
                .reference_density_kg_m3()
                .unwrap(),
            1000.0
        );
        assert_eq!(replay.motion.accepted_steps(), 2);
        assert_eq!(replay.motion.remaining_steps(), 0);
        assert_eq!(motion.solution(), &before);
        assert_eq!(motion.accepted_steps(), 1);
        balances(&replay.energy);
    });
}
