//! Material substitution -> draft/attitude equilibrium -> native floating motion.
//!
//! Run with the authored demonstration inputs:
//! `cargo run -p fs-mbd-material --features hydrostatic-equilibrium --example floating_equilibrium`
//! Optional arguments are LEFT_DENSITY RIGHT_DENSITY LIQUID_DENSITY [kg/m3].
//! For example, `-- 600 400 1000` reverses the denser side and the equilibrium tilt.
//! Diagnostics go to stderr; the native motion time series is CSV on stdout.

mod source;

use std::error::Error;

use fs_ga::{Se3, So3, So3Tangent, Twist, Vec3};
use fs_material::liquid::ResolvedLiquidState;
use fs_mbd::{Gravity, RigidBodyState};
use fs_mbd_material::hydrostatic::{
    DisplacementMesh, FloatingEquilibrium, FloatingEquilibriumControls, HydrostaticLoad,
    HydrostaticMedium, HydrostaticStabilityControls, HydrostaticStepControls,
    analyze_hydrostatic_stability, solve_floating_equilibrium,
};
use fs_mbd_material::principal::PrincipalBody;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn gravity() -> Result<Gravity> {
    Ok(Gravity::new(fs_mbd::Vec3::new(0.0, 0.0, -10.0))?)
}

fn stability_controls(roll_direction_world: fs_mbd::Vec3) -> HydrostaticStabilityControls {
    HydrostaticStabilityControls {
        roll_direction_world,
        length_scale_m: 1.0,
        heave_step_m: 1e-4,
        tilt_step_radians: 1e-4,
        maximum_force_defect_n: 1e-5,
        maximum_torque_defect_nm: 1e-5,
        maximum_relative_stiffness_defect: 1e-4,
        neutral_stiffness_tolerance_j: 1e-5,
    }
}

fn initialize(
    body: &PrincipalBody<'_>,
    medium: &HydrostaticMedium<'_, ResolvedLiquidState>,
    mesh: &DisplacementMesh,
) -> Result<FloatingEquilibrium> {
    // Both draft and attitude deliberately start away from equilibrium. The
    // solve holds the actual material COM's lateral position, not the authored
    // geometric origin; the principal frame conversion remains the sole owner.
    let reference = body.state_from_kinematics(
        Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(0.08, -0.10, 0.0)))?,
            Vec3::new(0.0, 0.0, 0.15),
        )?,
        Twist::zero(),
    )?;
    Ok(solve_floating_equilibrium(
        medium,
        mesh,
        body.frame(),
        reference,
        FloatingEquilibriumControls {
            roll_direction_world: fs_mbd::Vec3::new(1.0, 0.0, 0.0),
            maximum_force_defect_n: 1e-6,
            maximum_torque_defect_nm: 1e-6,
            length_scale_m: 1.0,
            difference_step: 1e-6,
            maximum_heave_offset_m: 1.0,
            maximum_tilt_radians: 1.2,
            maximum_newton_iterations: 32,
            maximum_evaluations: 2048,
        },
        || false,
    )?)
}

/// Signed inclinations of the authored body's +x/+y axes toward gravity.
/// This readback stays meaningful if the free yaw gauge changes; it is not an
/// Euler-angle decomposition of the integrator's principal-frame quaternion.
fn downhill_angles(pose: Se3) -> Result<[f64; 2]> {
    let x = pose.rotation().rotate(Vec3::new(1.0, 0.0, 0.0))?;
    let y = pose.rotation().rotate(Vec3::new(0.0, 1.0, 0.0))?;
    Ok([(-x.z).atan2(x.x.hypot(x.y)), (-y.z).atan2(y.x.hypot(y.y))])
}

fn observe(
    body: &PrincipalBody<'_>,
    state: RigidBodyState,
    load: HydrostaticLoad,
    time_s: f64,
    cumulative_energy_defect_j: f64,
) -> Result<()> {
    let pose = body.world_from_body(state)?;
    let angles = downhill_angles(pose)?;
    println!(
        "{time_s:.3},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{cumulative_energy_defect_j:.9}",
        pose.translation().z,
        state.center_of_mass_velocity_world(body.properties())?.z,
        angles[0].to_degrees(),
        angles[1].to_degrees(),
        load.displaced().volume_m3(),
        load.displaced_mass_kg(),
    );
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let densities = match args.as_slice() {
        [] => [400.0, 600.0, 1000.0],
        [left, right, liquid] => [left.parse()?, right.parse()?, liquid.parse()?],
        _ => return Err("usage: floating_equilibrium [LEFT_DENSITY_KG_M3 RIGHT_DENSITY_KG_M3 LIQUID_DENSITY_KG_M3]".into()),
    };
    let (source, liquid) = source::resolve(densities[0], densities[1], densities[2])?;
    let body = source.principal_body(1e-12, || false)?;
    let mesh = DisplacementMesh::cuboid(source::SIDES_M, || false)?;
    let medium = HydrostaticMedium::new(&liquid, gravity()?, Vec3::new(0.0, 0.0, 0.0))?;
    let equilibrium = initialize(&body, &medium, &mesh)?;
    let pose = body.world_from_body(equilibrium.state)?;
    let angles = downhill_angles(pose)?;

    eprintln!("Synthetic material data; 2 x 2 x 1 m box with two disjoint solid halves.");
    eprintln!(
        "left_density_kg_m3={}; right_density_kg_m3={}; liquid_density_kg_m3={}",
        densities[0], densities[1], densities[2]
    );
    eprintln!(
        "body_mass_kg={:.9}; displaced_fluid_mass_kg={:.9}; submerged_volume_m3={:.9}",
        body.properties().mass(),
        equilibrium.load.displaced_mass_kg(),
        equilibrium.load.displaced().volume_m3()
    );
    eprintln!(
        "material_com_body_m={:?}; equilibrium_com_world_m={:?}; center_of_buoyancy_world_m={:?}",
        source.inertia().center_of_mass(),
        equilibrium.state.pose().position_world(),
        equilibrium.load.center_of_buoyancy_world()
    );
    eprintln!(
        "envelope_origin_height_m={:.9}; x_axis_downhill_deg={:.9}; y_axis_downhill_deg={:.9}",
        pose.translation().z,
        angles[0].to_degrees(),
        angles[1].to_degrees()
    );
    eprintln!(
        "net_force_world_n={:?}; torque_world_nm={:?}; load_evaluations={}; iterations={}",
        equilibrium.net_force_world_n,
        equilibrium.torque_world_nm,
        equilibrium.evaluations,
        equilibrium.iterations
    );
    eprintln!(
        "Default 400/600/1000 case: independent torque balance gives 5*tan(theta)+4*tan(theta)^3=1.2, theta=12.965625921 deg."
    );

    let stability = analyze_hydrostatic_stability(
        &medium,
        &mesh,
        body.frame(),
        equilibrium.state,
        stability_controls(fs_mbd::Vec3::new(1.0, 0.0, 0.0)),
        || false,
    )?;
    eprintln!(
        "Local heave/roll/pitch stability={:?}; horizontal drift and yaw are excluded from this classification.",
        stability.classification
    );
    eprintln!(
        "restoring_stiffness_j={:?}; coordinates=(heave/1m, roll_rad, pitch_rad)",
        stability.restoring_stiffness_j
    );
    eprintln!(
        "restoring_eigenvalues_j={:?}; estimated_resolution_j={:.9}; relative_refinement_defect={:.3e}; relative_symmetry_defect={:.3e}; relative_one_sided_defect={:.3e}; load_evaluations={}",
        stability.eigenvalues_j,
        stability.estimated_resolution_j,
        stability.relative_refinement_defect,
        stability.relative_symmetry_defect,
        stability.relative_one_sided_defect,
        stability.evaluations
    );

    // A declared initial-condition change tests the real native restoring
    // motion after the equilibrium solve. This is not physical relaxation:
    // the prescribed hydrostatic reservoir includes no waves, added mass,
    // viscosity, drag or damping. Oscillations need not decay.
    let perturbation = So3::exp(So3Tangent::new(Vec3::new(0.006, -0.01, 0.0)))?;
    let mut state = body.state_from_kinematics(
        Se3::from_parts(
            perturbation.compose(pose.rotation())?,
            pose.translation() + Vec3::new(0.0, 0.0, 0.01),
        )?,
        Twist::zero(),
    )?;
    eprintln!(
        "Native motion starts from +0.01 m height and a small attitude perturbation; no physical damping is modeled."
    );
    println!(
        "time_s,envelope_origin_height_m,com_vertical_velocity_m_s,x_axis_downhill_deg,y_axis_downhill_deg,submerged_volume_m3,displaced_fluid_mass_kg,cumulative_numerical_energy_defect_j"
    );
    observe(
        &body,
        state,
        medium.load(&mesh, &body, state, || false)?,
        0.0,
        0.0,
    )?;
    let mut cumulative = 0.0;
    for index in 1..=1000 {
        let step = medium.step(
            &mesh,
            &body,
            state,
            HydrostaticStepControls {
                duration_s: 0.001,
                maximum_energy_defect_j: 1e-3,
            },
            || false,
        )?;
        state = step.dynamics.state_after;
        cumulative += step.energy_defect_j;
        if index % 50 == 0 {
            observe(
                &body,
                state,
                step.load_after,
                f64::from(index) * 0.001,
                cumulative,
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_mbd_material::hydrostatic::HydrostaticStability;

    #[test]
    fn swapping_material_halves_reverses_solved_tilt_without_changing_displaced_mass() -> Result<()>
    {
        for (left, right, sign) in [(400.0, 600.0, 1.0), (600.0, 400.0, -1.0)] {
            let (source, liquid) = source::resolve(left, right, 1000.0)?;
            let body = source.principal_body(1e-12, || false)?;
            let mesh = DisplacementMesh::cuboid(source::SIDES_M, || false)?;
            let medium = HydrostaticMedium::new(&liquid, gravity()?, Vec3::new(0.0, 0.0, 0.0))?;
            let solved = initialize(&body, &medium, &mesh)?;
            let pose = body.world_from_body(solved.state)?;
            let angles = downhill_angles(pose)?;
            assert!((solved.load.displaced_mass_kg() - 2000.0).abs() < 1e-7);
            assert!((angles[0] - sign * 0.22629286190286405).abs() < 1e-8);
            assert!(angles[1].abs() < 1e-8);

            // Align pitch with the actual body y axis, independently of the
            // free yaw gauge chosen by the two-tilt initialization chart.
            let x = pose.rotation().rotate(Vec3::new(1.0, 0.0, 0.0))?;
            let stability = analyze_hydrostatic_stability(
                &medium,
                &mesh,
                body.frame(),
                solved.state,
                stability_controls(fs_mbd::Vec3::new(x.x, x.y, 0.0)),
                || false,
            )?;
            assert_eq!(stability.classification, HydrostaticStability::Restoring);
            let k = stability.restoring_stiffness_j;
            let free_heave_pitch = k[2][2] - k[2][0] * k[0][2] / k[0][0];
            // With heave minimized, direct box integration gives
            // E(theta)/(rho*g) = 2/3 sec(theta) - 1/6 cos(theta)
            //                   - sign*0.2 sin(theta).
            // Its second derivative at the independent cubic root is this
            // stiffness. This checks the mixed COM-heave/tilt terms too.
            assert!(
                (free_heave_pitch - 9639.26346922437).abs() < 1e-3,
                "free-heave pitch stiffness {free_heave_pitch}"
            );
        }
        Ok(())
    }
}
