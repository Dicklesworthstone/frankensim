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
};
use crate::shell::stiffened::{
    BeamMass, assemble_stiffened_shell_materials_with_reference_density,
    assemble_stiffened_shell_with_reference_density,
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

#[cfg(test)]
mod tests;
