//! Actual tensile layer geometry, native mechanical loss and finite calorimetry.
use super::*;
use fs_solid::cohesive::dynamics::CohesiveMotion;

const THICKNESS: f64 = 0.1;
const DT: f64 = 0.01;

fn deformation() -> CohesiveLayerDeformation {
    CohesiveLayerDeformation::new(0.1, 0.1).unwrap()
}

fn normal_system<'a>(
    elastic: &'a dyn LinearOp,
    source: &'a ResolvedBilinearCohesive,
) -> CohesiveSystem<'a> {
    let patch = CohesivePatch::from_interface(
        7,
        1.0,
        source,
        [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
        vec![JumpTerm {
            dof: 0,
            local: [1.0, 0.0, 0.0],
        }],
    )
    .unwrap();
    CohesiveSystem::new(
        elastic,
        vec![patch],
        Budget {
            max_dofs: 1,
            max_patches: 1,
            max_jump_terms: 1,
            max_krylov_basis_entries: 128,
        },
    )
    .unwrap()
}

fn loaded<'s, 'a>(
    joint: &'s CohesiveSystem<'a>,
    mass: &'s Scalar,
    gap: f64,
    force: f64,
    velocity: f64,
) -> CohesiveMotion<'s, 'a> {
    joint
        .bonded_state()
        .solve_over_time_controlled(
            &[force],
            &[[gap, 0.0, 0.0]],
            1.0,
            dynamics().solve,
            &mut || false,
        )
        .unwrap()
        .state
        .into_dynamics(mass, vec![velocity], 0.0, dynamics())
        .unwrap()
}

// Independent physical-distance network: current layer length, unchanged face
// resistance and heat capacities. This does not use the production route helper.
fn stretched_reference(stretch: f64, heat: f64) -> ([f64; 2], [f64; 3]) {
    let conductances: [f64; 2] =
        std::array::from_fn(|i| 1.0 / (0.5 / K[i] + 0.5 * THICKNESS * stretch / K[2] + R[i]));
    let weights: [f64; 2] =
        std::array::from_fn(|i| DT * conductances[i] / (CAPACITY[i] + DT * conductances[i]));
    let rise = heat / (CAPACITY[2] + CAPACITY[0] * weights[0] + CAPACITY[1] * weights[1]);
    (
        conductances,
        [
            300.0 + weights[0] * rise,
            300.0 + weights[1] * rise,
            300.0 + rise,
        ],
    )
}

#[test]
fn normal_tensile_motion_uses_current_thickness_and_original_viscous_heat() {
    let (source, card) = material_and_card(true, Some(ADHESIVE));
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = normal_system(&elastic, &source);
    let gap = 0.0008;
    let motion = loaded(&joint, &mass, gap, (1000.0 + 2.0) * gap, 0.01);
    close(motion.state().displacements_m()[0], 0.0, 1e-13);
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let delta = (2.0 * mass.0 * 0.01 / DT - 1000.0 * gap)
        / (2.0 * mass.0 / (DT * DT) + 0.5 * (1000.0 + 1000.0) + 2.0 / DT);
    let viscous = 2.0 * delta * delta / DT;
    let stretch = 1.0 + (gap + delta) / THICKNESS;
    let next = propose_deformable_cohesive_layer_heat_with_deformation(
        &motion,
        &[0.0],
        &[[gap, 0.0, 0.0]],
        DT,
        binding(&route, &card),
        deformation(),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| Ok::<_, &'static str>(transport()),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    close(
        next.mechanical.motion.state().displacements_m()[0],
        delta,
        1e-12,
    );
    close(
        next.mechanical.motion.velocity_m_s()[0],
        2.0 * delta / DT - 0.01,
        1e-10,
    );
    close(next.heat.layer_normal_stretch, stretch, 1e-12);
    close(next.heat.dissipated_j, viscous, 1e-12);
    let (conductances, expected) = stretched_reference(stretch, viscous);
    for i in 0..2 {
        close(next.heat.conductances_w_k[i], conductances[i], 1e-13);
    }
    for i in 0..3 {
        let id = CellId(i as u64);
        close(next.thermal.temperature_k(id).unwrap(), expected[i], 2e-10);
        assert_eq!(
            next.thermal.cells()[&id].mass_kg,
            thermal.cells()[&id].mass_kg
        );
        assert_eq!(
            next.thermal.cells()[&id].geometry,
            thermal.cells()[&id].geometry
        );
    }
    assert_eq!(next.mechanical.motion.state().points()[0].damage(), 0.0);
    close(next.thermal.energy_j() - thermal.energy_j(), viscous, 1e-11);
    close(next.combined_energy_residual_j, 0.0, 2e-10);
    assert!(
        next.heat
            .entropy
            .as_ref()
            .unwrap()
            .network_conduction_generation_j_k
            >= 0.0
    );
    // The old fixed-reference contact entry has not silently acquired this model.
    let calls = Cell::new(0);
    let refused = propose_deformable_cohesive_layer_heat(
        &motion,
        &[0.0],
        &[[gap, 0.0, 0.0]],
        DT,
        binding(&route, &card),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| {
            calls.set(calls.get() + 1);
            Ok::<_, &'static str>(transport())
        },
        &CancelGate::new_clock_free(),
    );
    assert!(matches!(refused, Err(CohesiveHeatError::Invalid(_))));
    assert_eq!(calls.get(), 0);
    assert_eq!(motion.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn zero_normal_extension_recovers_original_shear_heating_exactly() {
    let (source, card) = material_and_card(true, Some(ADHESIVE));
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], true);
    let motion = joint
        .bonded_state()
        .into_dynamics(&mass, vec![0.02], 0.0, dynamics())
        .unwrap();
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let original = propose_deformable_cohesive_layer_heat(
        &motion,
        &[0.0],
        &[[0.0; 3]],
        DT,
        binding(&route, &card),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| Ok::<_, &'static str>(transport()),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let stretched = propose_deformable_cohesive_layer_heat_with_deformation(
        &motion,
        &[0.0],
        &[[0.0; 3]],
        DT,
        binding(&route, &card),
        CohesiveLayerDeformation::new(0.0, 0.1).unwrap(),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| Ok::<_, &'static str>(transport()),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(stretched.heat.layer_normal_stretch, 1.0);
    assert!(stretched.heat.dissipated_j > 0.0);
    assert_eq!(stretched.thermal.cells(), original.thermal.cells());
    assert_eq!(
        stretched.heat.conductances_w_k,
        original.heat.conductances_w_k
    );
    assert_eq!(stretched.heat.dissipated_j, original.heat.dissipated_j);
    assert_eq!(
        stretched.mechanical.motion.state().points(),
        original.mechanical.motion.state().points()
    );
    assert_eq!(
        stretched.mechanical.motion.velocity_m_s(),
        original.mechanical.motion.velocity_m_s()
    );
}

#[test]
fn progressive_tensile_damage_keeps_fracture_and_thermal_area_distinct() {
    let (source, card) = material_and_card(true, Some(ADHESIVE));
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], false);
    let motion = loaded(&joint, &mass, 0.0005, 0.0, 0.0);
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let next = propose_deformable_cohesive_layer_heat_with_deformation(
        &motion,
        &[0.0],
        &[[0.003, 0.0, 0.0]],
        DT,
        binding(&route, &card),
        deformation(),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| Ok::<_, &'static str>(transport()),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    close(
        next.mechanical.motion.state().points()[0].damage(),
        8.0 / 9.0,
        1e-12,
    );
    let survival_integral = 0.0005 + (4.0 / 3.0) * 0.001 * 3.0_f64.ln() - 0.002 / 3.0;
    let viscous = 2.0 * (0.003 - 0.0005) / DT * survival_integral;
    let bond_work =
        0.5 * 1000.0 * (0.001_f64.powi(2) - 0.0005_f64.powi(2)) + 0.5 * (1.0 + 1.0 / 3.0) * 0.002;
    let storage = 0.5 * 1000.0 * ((1.0 / 9.0) * 0.003_f64.powi(2) - 0.0005_f64.powi(2));
    close(next.heat.dissipated_j, viscous, 1e-11);
    close(next.fracture_work_j, bond_work - storage, 1e-11);
    close(next.thermal.energy_j() - thermal.energy_j(), viscous, 1e-10);
    let (conductances, expected) = stretched_reference(1.03, viscous);
    for i in 0..2 {
        close(next.heat.conductances_w_k[i], conductances[i], 1e-13);
    }
    for i in 0..3 {
        close(
            next.thermal.temperature_k(CellId(i as u64)).unwrap(),
            expected[i],
            2e-10,
        );
    }
    close(next.combined_energy_residual_j, 0.0, 2e-10);
    assert!(next.fracture_work_j > viscous);
}

#[test]
fn rupture_and_strain_refusals_leave_original_motion_and_heat_unmodified() {
    let (source, card) = material_and_card(true, Some(ADHESIVE));
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &source, &[7], false);
    let motion = loaded(&joint, &mass, 0.0005, 0.0, 0.0);
    let thermal = layer_network();
    let cells = thermal.cells().clone();
    let point = motion.state().points()[0];
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let calls = Cell::new(0);
    for (target, normal_limit, shear_limit) in [
        ([0.0045, 0.0, 0.0], 0.1, 0.1), // Full bond rupture at admissible strain.
        ([0.0008, 0.0, 0.0], 0.006, 0.1), // Endpoint normal strain.
        ([0.0005, 0.0002, 0.0], 0.1, 0.001), // Endpoint shear strain.
        ([0.0003, 0.0, 0.0], 0.004, 0.1), // Initial, not merely endpoint, strain.
    ] {
        let refused = propose_deformable_cohesive_layer_heat_with_deformation(
            &motion,
            &[0.0],
            &[target],
            DT,
            binding(&route, &card),
            CohesiveLayerDeformation::new(normal_limit, shear_limit).unwrap(),
            &thermal,
            layer_controls(),
            &[0.0; 3],
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok::<_, &'static str>(transport())
            },
            &CancelGate::new_clock_free(),
        );
        assert!(matches!(refused, Err(CohesiveHeatError::Invalid(_))));
        assert_eq!(motion.state().points()[0], point);
        assert_eq!(motion.state().displacements_m(), &[0.0]);
        assert_eq!(motion.velocity_m_s(), &[0.0]);
        assert_eq!(motion.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells(), &cells);
    }
    assert_eq!(calls.get(), 0);
}

#[test]
fn source_temperature_weakening_preserves_tensile_history_and_finite_layer_heat() {
    let (cold, card) =
        super::temperature::temperature_material("manufactured tensile thermal layer");
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let elastic = Scalar(1000.0);
    let mass = Scalar(10.0);
    let joint = system(&elastic, &cold, &[7], false);
    let hot_joint = joint
        .with_temperature_sources_controlled(&[&hot], "T", &mut || false)
        .unwrap();
    let motion = loaded(&joint, &mass, 0.0006, 0.0, 0.0);
    let point = motion.state().points()[0];
    let thermal = layer_network();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let next = propose_deformable_cohesive_layer_heat_with_temperature_and_deformation(
        &motion,
        &[0.0],
        &[[0.0005, 0.0, 0.0]],
        DT,
        super::temperature::source_temperature(&hot_joint),
        binding(&route, &card),
        deformation(),
        &thermal,
        layer_controls(),
        &[0.0; 3],
        |_, _, _| Ok::<_, &'static str>(transport()),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let heated = &next.heated;
    let viscous = (17.0 / 27.0) * 0.0001_f64.powi(2) / DT;
    close(
        heated.mechanical.motion.state().points()[0].damage(),
        10.0 / 27.0,
        1e-12,
    );
    close(next.temperature.parameter_work_j, -0.00009, 1e-13);
    close(next.fracture_work_j, 1.0 / 30000.0, 1e-13);
    close(heated.fracture_work_j, 0.0, 1e-13);
    close(heated.heat.dissipated_j, viscous, 1e-12);
    close(
        heated.thermal.energy_j() - thermal.energy_j(),
        viscous,
        1e-11,
    );
    let (conductances, expected) = stretched_reference(1.005, viscous);
    for i in 0..2 {
        close(heated.heat.conductances_w_k[i], conductances[i], 1e-13);
    }
    for i in 0..3 {
        close(
            heated.thermal.temperature_k(CellId(i as u64)).unwrap(),
            expected[i],
            2e-10,
        );
    }
    close(next.combined_energy_residual_j, 0.0, 2e-10);
    close(next.summed_energy_residual_j, 0.0, 2e-10);
    assert_eq!(motion.state().points()[0], point);
    assert_eq!(motion.state().points()[0].damage(), 0.0);
    assert_eq!(motion.time_s(), 0.0);
}
