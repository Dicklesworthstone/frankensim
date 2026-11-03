//! Density-weighted pressure correction with an authored volume constraint.

use super::{CancelGate, ChtError, FlowFace, FlowField, ProjectionReport, VoxelDomain};
use super::{finite, finite_positive, poll, poll_every};

mod forcing;
pub use forcing::{
    DensityForceControls, DensityForceInput, DensityForcePrediction, DensityForceReport,
    ForcedDensityProjection, ForcedDensityProjectionReport,
};
mod reconstruction;
pub use reconstruction::{
    DensityMomentumReconstruction, DensityMomentumReconstructionControls,
    DensityMomentumReconstructionReport,
};

/// Physical step and independent readback budgets for a density projection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DensityProjectionControls {
    /// Duration represented by the pressure impulse, seconds; positive.
    pub duration_s: f64,
    /// Relative true-residual tolerance of the native row-scaled Krylov solve.
    pub relative_tolerance: f64,
    /// Maximum absolute per-cell `net_outflow - target_outflow`, m^3/s.
    pub maximum_volume_defect_m3_s: f64,
    /// Maximum absolute defect of the discrete pressure-work identity, joules.
    pub maximum_energy_defect_j: f64,
}

/// Independent flux and mechanical-energy readback of the returned field.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityProjectionReport {
    /// Native solve evidence; its divergence fields remain *raw* net outflow.
    pub projection: ProjectionReport,
    /// Largest absolute per-cell volume-constraint defect before correction.
    pub max_volume_defect_before_m3_s: f64,
    /// Largest absolute per-cell volume-constraint defect after correction.
    pub max_volume_defect_after_m3_s: f64,
    /// Kinetic energy on active staggered face dual volumes before correction.
    pub kinetic_energy_before_j: f64,
    /// Kinetic energy on the same frozen dual masses after correction.
    pub kinetic_energy_after_j: f64,
    /// Signed kinetic-energy increment, evaluated without subtracting totals.
    pub kinetic_energy_change_j: f64,
    /// `dt sum(p_c target_outflow_c)`, joules, using correction pressure.
    pub pressure_source_work_j: f64,
    /// `dt sum(p_c Q_out)` on held `Fixed` faces, joules; subtracted from
    /// source work. The combined work is invariant to a sealed-component gauge.
    pub fixed_boundary_pressure_work_j: f64,
    /// Nonnegative `sum(M_face delta_u^2 / 2)`, joules. This numerical
    /// projection loss is reported separately and is never converted to heat.
    pub projection_loss_j: f64,
    /// `delta_K + projection_loss - source_work + fixed_boundary_work`, joules.
    pub energy_defect_j: f64,
}

/// Corrected native face fluxes, pressure impulse and accepted work evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityProjection {
    /// Actual corrected volumetric fluxes, m^3/s, in the original native layout.
    pub field: FlowField,
    /// Cell pressure correction in Pa; zero in solids and at `Free` boundary
    /// planes. A sealed component pins one fluid cell. This is not an absolute
    /// thermodynamic pressure or an equation-of-state update.
    pub pressure_correction_pa: Vec<f64>,
    /// Independent readback admitted before the result was returned.
    pub report: DensityProjectionReport,
}

pub(super) struct FrozenDensity<'a> {
    density_kg_m3: &'a [f64],
    pub(super) target_outflow_m3_s: &'a [f64],
    duration_s: f64,
    pub(super) uniform_density: Option<f64>,
}

impl FrozenDensity<'_> {
    fn face_mass(&self, domain: &VoxelDomain, c: usize, f: usize, nb: Option<usize>) -> f64 {
        face_dual_mass(domain, self.density_kg_m3, c, f, nb)
    }

    pub(super) fn conductance(
        &self,
        domain: &VoxelDomain,
        c: usize,
        f: usize,
        nb: Option<usize>,
    ) -> f64 {
        let area = domain.face_area(c, f / 2);
        self.duration_s * area * (area / self.face_mass(domain, c, f, nb))
    }
}

fn face_dual_mass(
    domain: &VoxelDomain,
    density_kg_m3: &[f64],
    c: usize,
    f: usize,
    nb: Option<usize>,
) -> f64 {
    let axis = f / 2;
    let half_c = 0.5 * domain.widths(c)[axis];
    let half_nb = nb.map_or(0.0, |n| 0.5 * domain.widths(n)[axis]);
    let mass_per_area = density_kg_m3[c] * half_c + nb.map_or(0.0, |n| density_kg_m3[n] * half_nb);
    domain.face_area(c, axis) * mass_per_area
}

impl FlowField {
    /// Project onto a prescribed per-cell volume flux with frozen physical
    /// density. `target_outflow_m3_s[c]` is the desired **net outward** volume
    /// flux, not a fractional divergence or a mass rate. Positive values
    /// require local expansion/injection and negative values contraction/removal.
    ///
    /// The native projection uses face dual mass
    /// `M = A (rho_left h_left/2 + rho_right h_right/2)` and pressure impulse
    /// `delta_Q = dt A^2 (p_left - p_right) / M`. Open boundary faces use the
    /// adjacent half-cell mass. Thus a sharp density jump uses the actual two
    /// half-cell inertias, including on a graded grid. Homogeneous density uses
    /// the original geometric potential solve and exactly retains its face
    /// arithmetic before converting potential to pressure.
    ///
    /// `Wall` and solid faces must already carry zero flux. `Fixed` boundary
    /// fluxes remain fixed; `Free` boundary correction pressure is zero. Every
    /// connected component without a free face must independently balance its
    /// total source with its prescribed boundary outflow. No mean is removed
    /// to conceal an incompatible source. Source entries in solids must be
    /// zero; density entries there are ignored. The supplied domain must be
    /// the same geometry/occupancy used to construct this field.
    ///
    /// This is a pressure/continuity component, not a density or phase time
    /// integrator. The caller owns the frozen density, any associated mass and
    /// momentum source, thermal storage, and source/boundary work. In a
    /// constant-density material a volume source implies a mass source
    /// `rho * target_outflow`; closed-material phase expansion instead needs a
    /// compatible evolving-density/enthalpy owner. Existing fixed-density
    /// thermal and phase entry points retain their original admission gates.
    /// Numerical projection loss never supplies latent or sensible heat.
    ///
    /// # Errors
    /// Invalid dimensions, non-finite data, nonpositive fluid density or time,
    /// blocked-face flux, solid-cell source, component incompatibility, native
    /// solver refusal, exceeded volume/energy readback budget, or cancellation.
    /// The original field is unchanged on every refusal.
    pub fn project_density_volume_source(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        density_kg_m3: &[f64],
        target_outflow_m3_s: &[f64],
        controls: DensityProjectionControls,
        gate: &CancelGate,
    ) -> Result<DensityProjection, ChtError> {
        let density = self.admit_density_projection(
            domain,
            faces,
            density_kg_m3,
            target_outflow_m3_s,
            controls,
            gate,
        )?;
        let mut field = self.clone();
        let (projection, cells, potential) = field.project_with_density(
            domain,
            faces,
            controls.relative_tolerance,
            Some(&density),
            gate,
        )?;
        let mut pressure = vec![0.0; domain.cell_count()];
        for (row, (&c, phi)) in cells.iter().zip(potential).enumerate() {
            poll_every(gate, row)?;
            pressure[c] = density
                .uniform_density
                .map_or(phi, |rho| phi * (rho / controls.duration_s));
            finite("flow.pressure_correction_pa", pressure[c])?;
        }
        let report = audit(
            self, &field, domain, faces, &density, &pressure, projection, controls, gate,
        )?;
        poll(gate)?;
        Ok(DensityProjection {
            field,
            pressure_correction_pa: pressure,
            report,
        })
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn admit_density_projection<'a>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        density_kg_m3: &'a [f64],
        target_outflow_m3_s: &'a [f64],
        controls: DensityProjectionControls,
        gate: &CancelGate,
    ) -> Result<FrozenDensity<'a>, ChtError> {
        poll(gate)?;
        finite_positive("flow.projection_duration_s", controls.duration_s)?;
        finite_positive("flow.projection_tolerance", controls.relative_tolerance)?;
        nonnegative(
            "flow.maximum_volume_defect_m3_s",
            controls.maximum_volume_defect_m3_s,
        )?;
        nonnegative(
            "flow.maximum_energy_defect_j",
            controls.maximum_energy_defect_j,
        )?;
        let uniform_density = self.admit_density_field(
            domain,
            faces,
            density_kg_m3,
            Some(target_outflow_m3_s),
            gate,
        )?;
        let density = FrozenDensity {
            density_kg_m3,
            target_outflow_m3_s,
            duration_s: controls.duration_s,
            uniform_density,
        };
        if let Some(reference) = uniform_density {
            finite_positive("flow.pressure_scale", reference / controls.duration_s)?;
        }
        for c in 0..domain.cell_count() {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                continue;
            }
            for f in 0..6 {
                let nb = domain.neighbor(c, f);
                if active_face(domain, faces, f, nb) {
                    finite_positive(
                        "flow.face_dual_mass_kg",
                        density.face_mass(domain, c, f, nb),
                    )?;
                    finite_positive(
                        "flow.pressure_conductance",
                        density.conductance(domain, c, f, nb),
                    )?;
                }
            }
        }
        Ok(density)
    }

    /// Shared layout, occupancy and frozen-density admission. A force predictor
    /// has no prescribed continuity source yet and does not allocate a dummy one.
    fn admit_density_field(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        density_kg_m3: &[f64],
        target_outflow_m3_s: Option<&[f64]>,
        gate: &CancelGate,
    ) -> Result<Option<f64>, ChtError> {
        poll(gate)?;
        let [nx, ny, nz] = domain.dims();
        if self.dims != domain.dims()
            || self.dx != domain.min_width()
            || self.fx.len() != (nx + 1) * ny * nz
            || self.fy.len() != nx * (ny + 1) * nz
            || self.fz.len() != nx * ny * (nz + 1)
            || density_kg_m3.len() != domain.cell_count()
            || target_outflow_m3_s.is_some_and(|source| source.len() != domain.cell_count())
        {
            return Err(refuse(
                "flow.density_projection_shape",
                "field, density and source must match the original domain",
            ));
        }
        let mut reference = None;
        let mut homogeneous = true;
        for c in 0..domain.cell_count() {
            poll_every(gate, c)?;
            if let Some(source) = target_outflow_m3_s {
                finite("flow.target_outflow_m3_s", source[c])?;
            }
            if domain.is_fluid(c) {
                finite_positive("flow.density_kg_m3", density_kg_m3[c])?;
                homogeneous &= *reference.get_or_insert(density_kg_m3[c]) == density_kg_m3[c];
            } else if target_outflow_m3_s.is_some_and(|source| source[c] != 0.0) {
                return Err(refuse(
                    "flow.solid_volume_source",
                    "a solid cell cannot carry a fluid volume source",
                ));
            }
            finite("flow.initial_net_outflow", self.net_outflow(domain, c))?;
            for f in 0..6 {
                let q = self.outward(domain, c, f);
                finite("flow.initial_face_flux", q)?;
                let blocked = !domain.is_fluid(c)
                    || domain
                        .neighbor(c, f)
                        .map_or(faces[f] == FlowFace::Wall, |n| !domain.is_fluid(n));
                if blocked && q != 0.0 {
                    return Err(refuse(
                        "flow.blocked_face_flux",
                        "wall and solid faces must carry exactly zero flux",
                    ));
                }
            }
        }
        let reference = reference.ok_or_else(|| {
            refuse(
                "flow.fluid_cells",
                "projection needs at least one fluid cell",
            )
        })?;
        Ok(homogeneous.then_some(reference))
    }
}

fn active_face(domain: &VoxelDomain, faces: [FlowFace; 6], f: usize, nb: Option<usize>) -> bool {
    match nb {
        Some(n) => f % 2 == 1 && domain.is_fluid(n),
        None => faces[f] != FlowFace::Wall,
    }
}

#[derive(Default)]
struct Sum {
    value: f64,
    compensation: f64,
}

impl Sum {
    fn add(&mut self, value: f64) {
        let corrected = value - self.compensation;
        let next = self.value + corrected;
        self.compensation = (next - self.value) - corrected;
        self.value = next;
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn audit(
    before: &FlowField,
    after: &FlowField,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    density: &FrozenDensity<'_>,
    pressure: &[f64],
    projection: ProjectionReport,
    controls: DensityProjectionControls,
    gate: &CancelGate,
) -> Result<DensityProjectionReport, ChtError> {
    let mut max_before = 0.0f64;
    let mut max_after = 0.0f64;
    let mut initial = Sum::default();
    let mut final_energy = Sum::default();
    let mut change = Sum::default();
    let mut loss = Sum::default();
    let mut source_work = Sum::default();
    let mut boundary_work = Sum::default();
    for c in 0..domain.cell_count() {
        poll_every(gate, c)?;
        if !domain.is_fluid(c) {
            continue;
        }
        let source = density.target_outflow_m3_s[c];
        let defect_before = before.net_outflow(domain, c) - source;
        let defect_after = after.net_outflow(domain, c) - source;
        finite("flow.volume_defect_before", defect_before)?;
        finite("flow.volume_defect_after", defect_after)?;
        max_before = max_before.max(defect_before.abs());
        max_after = max_after.max(defect_after.abs());
        source_work.add(controls.duration_s * pressure[c] * source);
        for f in 0..6 {
            let nb = domain.neighbor(c, f);
            if !active_face(domain, faces, f, nb) {
                continue;
            }
            let area = domain.face_area(c, f / 2);
            let mass = density.face_mass(domain, c, f, nb);
            let q0 = before.outward(domain, c, f);
            let q1 = after.outward(domain, c, f);
            let u0 = q0 / area;
            let u1 = q1 / area;
            let du = (q1 - q0) / area;
            initial.add(0.5 * mass * u0 * u0);
            final_energy.add(0.5 * mass * u1 * u1);
            change.add(mass * (u0 + 0.5 * du) * du);
            loss.add(0.5 * mass * du * du);
            if nb.is_none() && faces[f] == FlowFace::Fixed {
                boundary_work.add(controls.duration_s * pressure[c] * q1);
            }
        }
    }
    for value in [
        initial.value,
        final_energy.value,
        change.value,
        loss.value,
        source_work.value,
        boundary_work.value,
    ] {
        finite("flow.projection_energy", value)?;
    }
    let energy_defect = change.value + loss.value - source_work.value + boundary_work.value;
    finite("flow.projection_energy_defect", energy_defect)?;
    if max_after > controls.maximum_volume_defect_m3_s {
        return Err(refuse(
            "flow.volume_defect_budget",
            &format!(
                "actual per-cell volume defect {max_after:e} m^3/s exceeds {}",
                controls.maximum_volume_defect_m3_s
            ),
        ));
    }
    if energy_defect.abs() > controls.maximum_energy_defect_j {
        return Err(refuse(
            "flow.energy_defect_budget",
            &format!(
                "pressure-work defect {energy_defect:e} J exceeds {}",
                controls.maximum_energy_defect_j
            ),
        ));
    }
    Ok(DensityProjectionReport {
        projection,
        max_volume_defect_before_m3_s: max_before,
        max_volume_defect_after_m3_s: max_after,
        kinetic_energy_before_j: initial.value,
        kinetic_energy_after_j: final_energy.value,
        kinetic_energy_change_j: change.value,
        pressure_source_work_j: source_work.value,
        fixed_boundary_pressure_work_j: boundary_work.value,
        projection_loss_j: loss.value,
        energy_defect_j: energy_defect,
    })
}

fn nonnegative(field: &'static str, value: f64) -> Result<(), ChtError> {
    finite(field, value)?;
    if value < 0.0 {
        Err(refuse(field, "must be nonnegative"))
    } else {
        Ok(())
    }
}

fn refuse(field: &'static str, reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field,
        reason: reason.into(),
    }
}

#[cfg(test)]
mod tests;
