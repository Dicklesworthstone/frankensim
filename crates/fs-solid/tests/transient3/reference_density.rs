//! Analytical graded-tetrahedron inertia through native transient operations.
use super::*;
use fs_blake3::ContentHash;
use fs_scenario::FrameTree;
use fs_solid::loads3::fields::AccelerationField3;
use fs_solid::loads3::history::{TetLoadHistory3, VolumeHistory3};
use fs_solid::loads3::transient::adjoint::support::NodalMotionCotangent3;
use fs_solid::loads3::transient::events::MaterialRebindConfig3;

// Unit reference tet: V=1/6, mean density=6 and total mass=1 kg.
// Its graded mass is different from a uniform 6 kg/m³ specimen.
const RHO: [f64; 4] = [2.0, 4.0, 8.0, 10.0];
const FRAME: ContentHash = ContentHash([119; 32]);
const ROW_MASS: [f64; 4] = [13.0 / 60.0, 7.0 / 30.0, 4.0 / 15.0, 17.0 / 60.0];
const COLUMN_1: [f64; 4] = [1.0 / 24.0, 4.0 / 45.0, 1.0 / 20.0, 19.0 / 360.0];

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 3e-10 * (1.0 + expected.abs()),
        "{actual:.17e} != {expected:.17e}");
}

fn cfg() -> TransientConfig3 {
    TransientConfig3 {
        solve: TetStaticSolveConfig {
            relative_residual_tolerance: 1e-12,
            ..Default::default()
        },
        ..config()
    }
}

#[test]
fn free_body_history_and_checkpoint_retain_the_original_density_inventory() {
    with_cx(false, |cx| {
        // Deliberately inconsistent card density proves inertia and gravity
        // both come from the explicitly authored spatial reference field.
        let material = TetElasticMaterial::try_new(600.0, 60.0, 0.25, [119; 32]).unwrap();
        let p = problem(&material, &[]);
        let g = geometry(cx);
        let frames = FrameTree::new();
        let acceleration = [0.2, -0.4, 0.6];
        let velocity = [0.1, -0.05, 0.2];
        let history = TetLoadHistory3 {
            surfaces: &[],
            volume: VolumeHistory3::FixedAcceleration(AccelerationField3::Uniform(acceleration)),
            nodal_forces_n: &[],
            frames: &frames,
        };
        let mut authored = RHO;
        let mut run = history.begin_transient_with_nodal_reference_density(
            &g, &p, &authored, None, 0.0, &ZERO, &[velocity; 4], cfg(), cx).unwrap();
        authored.fill(600.0);
        assert_eq!(run.nodal_reference_density_kg_m3(), Some(RHO.as_slice()));
        close(run.reference_mass_kg(), 1.0);
        for _ in 0..4 {
            let step = run.advance_history(&history, 0.01, cx).unwrap();
            close(step.balance_error_j, 0.0);
        }
        let t = run.state().time_s;
        let motion = run.nodal_motion(cx).unwrap();
        for node in 0..4 {
            for c in 0..3 {
                close(motion[0][node][c], velocity[c] * t + 0.5 * acceleration[c] * t * t);
                close(motion[1][node][c], velocity[c] + acceleration[c] * t);
                close(motion[2][node][c], acceleration[c]);
            }
        }
        close(run.state().kinetic_energy_j,
            0.5 * (0..3).map(|c| (velocity[c] + acceleration[c] * t).powi(2)).sum::<f64>());
        let response = run.recover_response(FRAME, 0.01, cx).unwrap();
        close(response.free_force_residual_n, 0.0);
        for (node, mass) in ROW_MASS.into_iter().enumerate() {
            for c in 0..3 {
                close(response.inertia_force_n[node][c], mass * acceleration[c]);
            }
        }
        let mut fork = run.clone();
        let checkpoint = run.into_checkpoint(cx).unwrap();
        let mut rebound = checkpoint.bind(&p, cx).unwrap();
        assert_eq!(rebound.nodal_reference_density_kg_m3(), Some(RHO.as_slice()));
        assert_eq!(fork.advance_history(&history, 0.01, cx).unwrap(),
            rebound.advance_history(&history, 0.01, cx).unwrap());
        assert_eq!(fork.state(), rebound.state());

        let schedule = fs_solid::loads3::history::dynamics::events::EventLoadHistory3::new(
            history, 32, cx).unwrap();
        let mut scheduled = schedule.begin_transient_with_nodal_reference_density(
            &g, &p, &RHO, None, 0.0, &ZERO, &ZERO, cfg(), cx).unwrap();
        scheduled.advance_to(0.02, 0.02, cx).unwrap();
        for node in scheduled.run().nodal_motion(cx).unwrap()[2].iter() {
            for c in 0..3 { close(node[c], acceleration[c]); }
        }
    });
}

#[test]
fn oscillator_impulse_and_prestressed_pencil_share_exact_graded_mass() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(600.0, 60.0, 0.0, [119; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let p = problem(&material, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let mut displacement = ZERO;
        displacement[1][0] = 0.001;
        let mut run = load.begin_transient_with_nodal_reference_density(
            &p, &RHO, None, 0.0, &displacement, &ZERO, cfg(), cx).unwrap();
        // The original cell has K_ff=10 N/m, but graded M_ff=4/45 kg.
        // Replacing the density by its mean would incorrectly give M_ff=.1.
        close(run.state().acceleration[0], -0.01 / COLUMN_1[1]);
        close(run.state().total_energy_j, 5e-6);
        let response = run.recover_response(FRAME, 0.01, cx).unwrap();
        for node in 0..4 {
            close(response.inertia_force_n[node][0],
                COLUMN_1[node] * run.state().acceleration[0]);
        }
        let tangent = run.assemble_prestressed_constitutive_tangent(
            p.materials, FRAME, 0.01, cx).unwrap();
        close(tangent.pencil().mass.row(0).1[0], COLUMN_1[1]);
        close(tangent.pencil().total_mass_kg, 1.0);

        let impulse = 0.04;
        let report = run.apply_nodal_impulse(&[(3, impulse)], Default::default(), cx).unwrap();
        let delta_v = impulse / COLUMN_1[1];
        close(run.state().velocity[0], delta_v);
        for node in 0..4 {
            close(report.nodal_momentum_change_kg_m_s[node][0], COLUMN_1[node] * delta_v);
            if node != 1 {
                close(report.support_impulse_kg_m_s[node][0], COLUMN_1[node] * delta_v);
            }
        }
        close(report.kinetic_change_j, 0.5 * impulse * delta_v);
        close(report.applied_impulse_work_j, report.kinetic_change_j);
        run.advance(&load, 0.002, cx).unwrap();
    });
}

#[test]
fn fully_prescribed_body_retains_graded_reactions_damping_and_support_work() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(600.0, 60.0, 0.25, [119; 32]).unwrap();
        let fixed: Vec<_> = (0..12).collect();
        let p = problem(&material, &fixed);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let acceleration = [0.2, -0.4, 0.6];
        let velocity = [0.1, -0.05, 0.2];
        let eta = 0.4;
        let mut run = load.begin_transient_with_nodal_reference_density_and_supports(
            &p, &RHO, None, 0.0, &ZERO, &[velocity; 4], &[acceleration; 4],
            TransientConfig3 { mass_damping_per_s: eta, ..cfg() }, cx).unwrap();
        assert!(run.free_dofs().is_empty());
        close(run.reference_mass_kg(), 1.0);
        let initial_energy = run.state().kinetic_energy_j;
        let supports: Vec<_> = (0..12).map(|d| (d, acceleration[d % 3])).collect();
        let step = run.advance_with_support_accelerations(&load, &supports, 0.02, cx).unwrap();
        close(step.support_work_j,
            run.state().kinetic_energy_j - initial_energy + step.damping_loss_j);
        let response = run.recover_response(FRAME, 0.01, cx).unwrap();
        for (node, mass) in ROW_MASS.into_iter().enumerate() {
            for c in 0..3 {
                let speed = velocity[c] + 0.02 * acceleration[c];
                close(response.inertia_force_n[node][c], mass * acceleration[c]);
                close(response.damping_force_n[node][c], eta * mass * speed);
                close(response.support_reactions_n[node][c], mass * (acceleration[c] + eta * speed));
            }
        }
        // With every acceleration prescribed, scaling the spatial inventory
        // scales the complete reaction directly, including mass damping.
        let gradient = run.reaction_pullback(&[(2, 1.0)], Default::default(), cx).unwrap();
        close(gradient.material.element_scales[0][1], response.support_reactions_n[0][2]);
        close(gradient.material.rayleigh[0], ROW_MASS[0] * (velocity[2] + 0.02 * acceleration[2]));
    });
}

#[test]
fn material_and_thermal_refresh_preserve_explicit_reference_inventory() {
    with_cx(false, |cx| {
        let original = TetElasticMaterial::try_new(600.0, 60.0, 0.0, [119; 32]).unwrap();
        let replacement = TetElasticMaterial::try_new(300.0, 120.0, 0.0, [120; 32]).unwrap();
        let heated = TetElasticMaterial::try_new(150.0, 240.0, 0.0, [121; 32]).unwrap();
        let p = problem(&original, &[]);
        let next = problem(&replacement, &[]);
        let hot = problem(&heated, &[]);
        let g = geometry(cx);
        let load = g.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let velocity = [[0.1, -0.05, 0.2]; 4];
        let run = load.begin_transient_with_nodal_reference_density(
            &p, &RHO, None, 0.0, &ZERO, &velocity, cfg(), cx).unwrap();
        let (rebound, event) = run.rebind_material_candidate(
            &next, None, &load, FRAME, MaterialRebindConfig3 {
                maximum_absolute_work_j: 1e-12,
                absolute_momentum_tolerance_kg_m_s: 1e-12,
                maximum_strain: 0.01,
            }, cx).unwrap();
        assert_eq!(event.reference_mass_kg, [1.0; 2]);
        assert_eq!(run.state().velocity, rebound.state().velocity);
        close(event.kinetic_work_j, 0.0);
        let thermal = TetThermalStrainState::try_new(301.0, 300.0, [0.0; 6],
            heated.material_state_identity, FRAME, 0.01).unwrap();
        let (refreshed, step) = rebound.advance_thermal_candidate(
            &hot, TetThermalStrainField::Uniform(&thermal), &load, 0.01, FRAME, 0.01, cx).unwrap();
        assert_eq!(refreshed.nodal_reference_density_kg_m3(), Some(RHO.as_slice()));
        close(refreshed.reference_mass_kg(), 1.0);
        close(refreshed.state().kinetic_energy_j, run.state().kinetic_energy_j);
        close(step.physical_balance_error_j, 0.0);
        for node in refreshed.nodal_motion(cx).unwrap()[1].iter() {
            for c in 0..3 { close(node[c], velocity[0][c]); }
        }
    });
}

#[test]
fn graded_support_and_density_adjoints_match_independent_acceleration_derivatives() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(600.0, 60.0, 0.0, [119; 32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|d| *d != 3).collect();
        let p = problem(&material, &fixed);
        let g = geometry(cx);
        let force = 0.03;
        let load = g.assemble(&[], VolumeLoad3::None, &[(3, force)], cx).unwrap();
        let run_at = |density_scale: f64, support_acceleration: f64| {
            let density = RHO.map(|value| value * density_scale);
            let mut prescribed = ZERO;
            prescribed[0][0] = support_acceleration;
            load.begin_transient_with_nodal_reference_density_and_supports(
                &p, &density, None, 0.0, &ZERO, &ZERO, &prescribed, cfg(), cx).unwrap()
        };
        let run = run_at(1.0, 0.2);
        close(run.state().acceleration[0], (force - COLUMN_1[0] * 0.2) / COLUMN_1[1]);
        let adjoint = run.prescribed_adjoint(Default::default(), Default::default(), cx).unwrap();
        let mut seed = NodalMotionCotangent3 {
            q: ZERO.to_vec(), v: ZERO.to_vec(), a: ZERO.to_vec(),
        };
        seed.a[1][0] = 1.0;
        let pullback = adjoint.pullback_initialization(&seed, cx).unwrap();
        close(pullback.initial.a[0][0], -COLUMN_1[0] / COLUMN_1[1]);
        let epsilon = 1e-6;
        close(pullback.initial.a[0][0],
            (run_at(1.0, 0.2 + epsilon).state().acceleration[0]
                - run_at(1.0, 0.2 - epsilon).state().acceleration[0]) / (2.0 * epsilon));
        let density = run.material_scale_pullback(&pullback.force, Default::default(), cx).unwrap();
        close(density.element_scales[0][1], -force / COLUMN_1[1]);
        close(density.element_scales[0][1],
            (run_at(1.0 + epsilon, 0.2).state().acceleration[0]
                - run_at(1.0 - epsilon, 0.2).state().acceleration[0]) / (2.0 * epsilon));
    });
}
