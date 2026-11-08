//! Original graded inventory through the native loaded-shell vibration pencil.
//!
//! Density determines the retained reference inertia. Nonlinear membrane/beam
//! equilibrium, tethers, tangent stiffness and stability remain owned by the
//! original preload path. Prescribed dead loads are explicit: changing density
//! does not silently change those loads or invent a gravity field. Small
//! vibrations retain the original moderate-rotation/reference-frame scope.

use super::{
    PlateSection, PreloadOptions, PreloadProblem, PreloadedShell, ShellBeam, ShellMesh,
    ShellSection, ShellSupport, ShellTether, equilibrate_sections_with_assembly,
    equilibrate_sections_with_load_assembly,
};
use crate::shell::stiffened::{
    BeamMass, assemble_stiffened_shell_materials_with_reference_density,
    assemble_stiffened_shell_materials_with_reference_density_and_body_acceleration,
    assemble_stiffened_shell_with_reference_density,
    assemble_stiffened_shell_with_reference_density_and_body_acceleration,
};

/// Solve supported nonlinear equilibrium with the original graded reference mass.
///
/// Supply one P1 density triple per shell facet and one endpoint density pair
/// per bonded beam, in their ORIGINAL input order. The shared assembly uses
/// consistent panel inertia and the selected beam inertia, including eccentric
/// translation/rotation coupling. The complete assembled mass survives the
/// nonlinear equilibrium solve and accompanies its accepted tangent for modes.
///
/// `tethers` may be empty. Their original conservative potential, force and
/// geometric stiffness participate in equilibrium, without added mass. Source
/// coverage, load/work/deformation bounds, the injected nonlinear solver and
/// native force-balance/stability acceptance are the same as ordinary preload.
/// No density resampling, averaging, extra mass assembly or replacement solver
/// occurs during continuation. The inputs remain immutable on every refusal.
///
/// # Errors
/// Original supported-preload refusals and invalid density coverage/values or
/// unrepresentable graded inertia. Density is admitted before invoking `solve`.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_stiffened_shell_with_reference_density(
    mesh: &ShellMesh,
    sections: &[PlateSection],
    reference_corner_density_kg_m3: &[[f64; 3]],
    boundary: &[usize],
    support: ShellSupport,
    beams: &[ShellBeam],
    reference_beam_density_kg_m3: &[[f64; 2]],
    mass: BeamMass,
    loads: &[f64],
    tethers: &[ShellTether],
    options: PreloadOptions,
    solve: &mut impl FnMut(&PreloadProblem<'_>, Vec<f64>, usize) -> Result<(Vec<f64>, usize), String>,
) -> Result<PreloadedShell, String> {
    equilibrate_sections_with_assembly(
        mesh,
        sections,
        boundary,
        support,
        beams,
        loads,
        tethers,
        options,
        solve,
        || {
            assemble_stiffened_shell_with_reference_density(
                mesh,
                sections,
                reference_corner_density_kg_m3,
                boundary,
                support,
                beams,
                reference_beam_density_kg_m3,
                mass,
            )
        },
    )
}

/// Independent membrane/bending sections with the same graded preload path.
///
/// Density remains an explicit original inventory; it does not infer elastic
/// coupling or change the source sections. All other physical semantics and
/// refusals are those of [`equilibrate_stiffened_shell_with_reference_density`].
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_shell_materials_with_reference_density(
    mesh: &ShellMesh,
    sections: &[ShellSection],
    reference_corner_density_kg_m3: &[[f64; 3]],
    boundary: &[usize],
    support: ShellSupport,
    beams: &[ShellBeam],
    reference_beam_density_kg_m3: &[[f64; 2]],
    mass: BeamMass,
    loads: &[f64],
    tethers: &[ShellTether],
    options: PreloadOptions,
    solve: &mut impl FnMut(&PreloadProblem<'_>, Vec<f64>, usize) -> Result<(Vec<f64>, usize), String>,
) -> Result<PreloadedShell, String> {
    equilibrate_sections_with_assembly(
        mesh,
        sections,
        boundary,
        support,
        beams,
        loads,
        tethers,
        options,
        solve,
        || {
            assemble_stiffened_shell_materials_with_reference_density(
                mesh,
                sections,
                reference_corner_density_kg_m3,
                boundary,
                support,
                beams,
                reference_beam_density_kg_m3,
                mass,
            )
        },
    )
}

/// Solve graded reinforced-shell equilibrium under uniform body acceleration.
///
/// `acceleration_m_s2` is an explicit physical acceleration in the reference
/// mesh's Cartesian axes, for example authored gravity. The shared assembly
/// projects rigid translation through each original panel/beam mass BEFORE
/// support restriction. This retains consistent-mass contributions from fixed
/// nodes and the generalized moments of eccentric reinforcement. Multiplying
/// only the reduced mass by free translation would lose those contributions.
///
/// The complete body force is added once to `loads`, then the original preload
/// owner applies support restriction, continuation, tethers and acceptance.
/// Supply only additional dead forces in `loads`; do not include this same
/// weight again. One assembly supplies both the body load and retained inertia.
///
/// This is a dead load under the native moderate-rotation/reference-frame
/// kinematics. It does not infer a gravity field, update a finite-rotation
/// follower load, or include centripetal/angular-acceleration terms.
///
/// # Errors
/// Nonfinite acceleration, invalid original inventory or supports, arithmetic
/// overflow, or any native equilibrium/work/stability refusal. Invalid input
/// cannot invoke `solve` or publish a partial loaded model.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_stiffened_shell_with_reference_density_and_body_acceleration(
    mesh: &ShellMesh,
    sections: &[PlateSection],
    reference_corner_density_kg_m3: &[[f64; 3]],
    boundary: &[usize],
    support: ShellSupport,
    beams: &[ShellBeam],
    reference_beam_density_kg_m3: &[[f64; 2]],
    mass: BeamMass,
    acceleration_m_s2: [f64; 3],
    loads: &[f64],
    tethers: &[ShellTether],
    options: PreloadOptions,
    solve: &mut impl FnMut(&PreloadProblem<'_>, Vec<f64>, usize) -> Result<(Vec<f64>, usize), String>,
) -> Result<PreloadedShell, String> {
    admit_acceleration(acceleration_m_s2)?;
    equilibrate_sections_with_load_assembly(
        mesh,
        sections,
        boundary,
        support,
        beams,
        loads,
        tethers,
        options,
        solve,
        || {
            let (model, body) =
                assemble_stiffened_shell_with_reference_density_and_body_acceleration(
                    mesh,
                    sections,
                    reference_corner_density_kg_m3,
                    boundary,
                    support,
                    beams,
                    reference_beam_density_kg_m3,
                    mass,
                )?;
            let force = body.full_load(acceleration_m_s2)?;
            Ok((model, force))
        },
    )
}

/// Independent membrane/bending sections under the same complete body load.
///
/// Physical scope, additive dead loads and admission are identical to
/// [`equilibrate_stiffened_shell_with_reference_density_and_body_acceleration`].
/// The original independent constitutive sections remain authoritative.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_shell_materials_with_reference_density_and_body_acceleration(
    mesh: &ShellMesh,
    sections: &[ShellSection],
    reference_corner_density_kg_m3: &[[f64; 3]],
    boundary: &[usize],
    support: ShellSupport,
    beams: &[ShellBeam],
    reference_beam_density_kg_m3: &[[f64; 2]],
    mass: BeamMass,
    acceleration_m_s2: [f64; 3],
    loads: &[f64],
    tethers: &[ShellTether],
    options: PreloadOptions,
    solve: &mut impl FnMut(&PreloadProblem<'_>, Vec<f64>, usize) -> Result<(Vec<f64>, usize), String>,
) -> Result<PreloadedShell, String> {
    admit_acceleration(acceleration_m_s2)?;
    equilibrate_sections_with_load_assembly(
        mesh,
        sections,
        boundary,
        support,
        beams,
        loads,
        tethers,
        options,
        solve,
        || {
            let (model, body) =
                assemble_stiffened_shell_materials_with_reference_density_and_body_acceleration(
                    mesh,
                    sections,
                    reference_corner_density_kg_m3,
                    boundary,
                    support,
                    beams,
                    reference_beam_density_kg_m3,
                    mass,
                )?;
            let force = body.full_load(acceleration_m_s2)?;
            Ok((model, force))
        },
    )
}

fn admit_acceleration(acceleration_m_s2: [f64; 3]) -> Result<(), String> {
    if acceleration_m_s2.iter().any(|value| !value.is_finite()) {
        return Err("uniform shell body acceleration must be finite".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests;
