//! Density-weighted pressure correction with an authored volume constraint.

use super::{
    CancelGate, ChtError, FlowFace, FlowField, ProjectionPotential, ProjectionReport, VoxelDomain,
};
use super::{finite, finite_positive, poll, poll_every};

mod feedback;
pub use feedback::{
    MaterialDensityMomentumControls, MaterialDensityMomentumInput, MaterialDensityMomentumReport,
    MaterialDensityMomentumStep,
};
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

/// Borrowed reservoir correction pressures in the native staggered face layout.
///
/// Each slice has exactly the length of its corresponding [`FlowField`] flux
/// buffer. `Some(p)` prescribes a finite face-constant pressure in Pa at a
/// fluid-exposed `Free` boundary face; `None` retains zero there. All interior,
/// solid, `Wall` and `Fixed` slots must be `None`, including when the authored
/// pressure would be zero. Thus a profile cannot also prescribe a held velocity.
///
/// Storage is caller-owned and can be reused across steps. Admission checks
/// the three existing face counts before indexing or cloning the flow field;
/// no profile-sized allocation or callback occurs in the projection. The
/// caller owns face averaging and spatial refinement of a continuous profile.
#[derive(Debug, Clone, Copy)]
pub struct StaggeredBoundaryPressure<'a> {
    /// X-normal face pressures: `(z * ny + y) * (nx + 1) + i`, `0 <= i <= nx`.
    pub x_pa: &'a [Option<f64>],
    /// Y-normal face pressures: `(z * (ny + 1) + j) * nx + x`, `0 <= j <= ny`.
    pub y_pa: &'a [Option<f64>],
    /// Z-normal face pressures: `(k * ny + y) * nx + x`, `0 <= k <= nz`.
    pub z_pa: &'a [Option<f64>],
}

/// Open-face reservoir correction pressure for composed native flow owners.
/// A uniform plane array needs no staggered profile allocation. A borrowed
/// profile retains the original spatial head and native face layout. Both use
/// the same pressure projection, half-cell inertia and mechanical work account;
/// neither supplies a thermodynamic reservoir state.
#[derive(Debug, Clone, Copy)]
pub enum ReservoirPressure<'a> {
    /// Uniform pressures in Pa on `[x-, x+, y-, y+, z-, z+]`; `None` means zero.
    /// A supplied value must be finite and its plane must be `Free`.
    Planes([Option<f64>; 6]),
    /// Separate pressure at each original staggered boundary face.
    Profile(StaggeredBoundaryPressure<'a>),
}

impl From<[Option<f64>; 6]> for ReservoirPressure<'_> {
    fn from(planes: [Option<f64>; 6]) -> Self {
        Self::Planes(planes)
    }
}

impl<'a> From<StaggeredBoundaryPressure<'a>> for ReservoirPressure<'a> {
    fn from(profile: StaggeredBoundaryPressure<'a>) -> Self {
        Self::Profile(profile)
    }
}

impl<'a> ReservoirPressure<'a> {
    fn admit(
        self,
        faces: [FlowFace; 6],
        gate: &CancelGate,
    ) -> Result<BoundaryPressure<'a>, ChtError> {
        poll(gate)?;
        match self {
            Self::Planes(planes) => {
                for (f, pressure) in planes.iter().enumerate() {
                    if let Some(pressure) = pressure {
                        finite("flow.boundary_pressure_pa", *pressure)?;
                        if faces[f] != FlowFace::Free {
                            return Err(refuse(
                                "flow.boundary_pressure_pa",
                                "pressure is only admissible on a Free face",
                            ));
                        }
                    }
                }
                Ok(BoundaryPressure::Planes(planes.map(|p| p.unwrap_or(0.0))))
            }
            // Original projection admission checks the native layout and roles
            // before indexing a profile or allocating its corrected field.
            Self::Profile(profile) => Ok(BoundaryPressure::Profile(profile)),
        }
    }
}

impl StaggeredBoundaryPressure<'_> {
    fn value(self, axis: usize, slot: usize) -> Option<f64> {
        [self.x_pa, self.y_pa, self.z_pa][axis][slot]
    }

    fn admit(
        self,
        field: &FlowField,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        poll(gate)?;
        if [self.x_pa, self.y_pa, self.z_pa]
            .iter()
            .zip(field.face_arrays())
            .any(|(pressure, flux)| pressure.len() != flux.len())
        {
            return Err(refuse(
                "flow.boundary_pressure_shape",
                "pressure buffers must match all three native staggered face layouts",
            ));
        }
        for c in 0..domain.cell_count() {
            poll_every(gate, c)?;
            for f in 0..6 {
                let nb = domain.neighbor(c, f);
                // Check every native slot once, including faces of solids.
                if nb.is_some() && f % 2 == 0 {
                    continue;
                }
                let (slot, _) = field.slot(domain, c, f);
                if let Some(pressure) = self.value(f / 2, slot) {
                    finite("flow.boundary_pressure_pa", pressure)?;
                    if nb.is_some() || !domain.is_fluid(c) || faces[f] != FlowFace::Free {
                        return Err(refuse(
                            "flow.boundary_pressure_role",
                            "pressure is only admissible on a fluid-exposed Free exterior face",
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
enum BoundaryPressure<'a> {
    Planes([f64; 6]),
    Profile(StaggeredBoundaryPressure<'a>),
}

impl BoundaryPressure<'_> {
    fn value(self, f: usize, slot: usize) -> f64 {
        match self {
            Self::Planes(values) => values[f],
            Self::Profile(profile) => profile.value(f / 2, slot).unwrap_or(0.0),
        }
    }
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
    /// `dt sum(p_c target_outflow_c)`, joules. Work uses pressure relative to
    /// the first open reservoir in each fluid component (zero for sealed
    /// components), so a common pressure offset cannot amplify solver error.
    pub pressure_source_work_j: f64,
    /// `dt sum(p_c Q_out)` on held `Fixed` faces, joules; subtracted from
    /// source work. The combined work is invariant to a sealed-component gauge.
    pub fixed_boundary_pressure_work_j: f64,
    /// `dt sum(p_boundary Q_out)` on `Free` faces, joules, in the same
    /// component gauge as source work. Negative values supply mechanical
    /// energy from the reservoirs; positive values deliver it to them.
    pub open_boundary_pressure_work_j: f64,
    /// Nonnegative `sum(M_face delta_u^2 / 2)`, joules. This numerical
    /// projection loss is reported separately and is never converted to heat.
    pub projection_loss_j: f64,
    /// `delta_K + projection_loss - source_work + fixed_boundary_work
    /// + open_boundary_work`, joules, in the common component gauge.
    pub energy_defect_j: f64,
}

/// Corrected native face fluxes, pressure impulse and accepted work evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct DensityProjection {
    /// Actual corrected volumetric fluxes, m^3/s, in the original native layout.
    pub field: FlowField,
    /// Cell pressure correction in Pa, in the authored boundary-pressure
    /// gauge; zero in solids. The original entry prescribes zero on `Free`
    /// boundary planes. A sealed component pins one fluid cell. This is not
    /// an absolute thermodynamic pressure or an equation-of-state update.
    pub pressure_correction_pa: Vec<f64>,
    /// Independent readback admitted before the result was returned.
    pub report: DensityProjectionReport,
}

// Fresh composed stages retain the native relative solution until their
// impulse/work readback is complete. Reuse the solver's original row/component
// maps and values; subtracting restored absolute pressures can lose a small
// driving head under a large reservoir reference. No second pressure map or
// solve is needed, and no private state is added to the public result.
struct PressureReadback {
    potential: ProjectionPotential,
    uniform_pressure_scale: Option<f64>,
}

impl PressureReadback {
    fn cell_pressure(&self, cell: usize) -> f64 {
        let phi = self.potential.values[self.potential.row_of[cell]];
        self.uniform_pressure_scale.map_or(phi, |scale| phi * scale)
    }

    fn difference(&self, cell: usize, neighbor: Option<usize>, boundary_pa: f64) -> f64 {
        let other = neighbor.map_or_else(
            || {
                let row = self.potential.row_of[cell];
                boundary_pa - self.potential.reference_pressure_pa[self.potential.component[row]]
            },
            |neighbor| self.cell_pressure(neighbor),
        );
        self.cell_pressure(cell) - other
    }
}

pub(super) struct FrozenDensity<'a> {
    density_kg_m3: &'a [f64],
    pub(super) target_outflow_m3_s: &'a [f64],
    duration_s: f64,
    pub(super) uniform_density: Option<f64>,
    boundary_pressure: BoundaryPressure<'a>,
}

impl FrozenDensity<'_> {
    pub(super) fn boundary_pressure_pa(&self, f: usize, slot: usize) -> f64 {
        self.boundary_pressure.value(f, slot)
    }

    pub(super) fn boundary_potential(&self, f: usize, slot: usize, reference_pa: f64) -> f64 {
        let relative_pressure = self.boundary_pressure_pa(f, slot) - reference_pa;
        self.uniform_density.map_or(relative_pressure, |rho| {
            relative_pressure / (rho / self.duration_s)
        })
    }

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

pub(super) fn face_dual_mass(
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
        self.project_density_volume_source_with_boundary_pressure(
            domain,
            faces,
            [None; 6],
            density_kg_m3,
            target_outflow_m3_s,
            controls,
            gate,
        )
    }

    /// Project using prescribed reservoir correction pressures, in Pa, on
    /// `Free` domain planes ordered `[x-, x+, y-, y+, z-, z+]`. `None` retains
    /// the original zero pressure. A supplied pressure on a `Wall` or `Fixed`
    /// face is refused: those faces already prescribe normal velocity.
    ///
    /// The face pressure acts at the boundary plane through the actual
    /// adjacent half-cell inertia. Pressure differences can drive either
    /// inflow or outflow, without an artificial outflow-only clamp. The solve
    /// and work audit remove a common reservoir offset independently in each
    /// connected fluid component; pressure readback restores the authored
    /// gauge. Consequently atmospheric or other large reference pressures do
    /// not hide a small pressure-driven acceleration in the solver tolerance.
    /// Pressures are uniform on each authored domain plane. Use
    /// [`Self::project_density_volume_source_with_pressure_profile`] for a
    /// spatially varying head; reservoir thermodynamic states belong to the caller.
    ///
    /// All density, source, compatibility and rollback requirements of
    /// [`Self::project_density_volume_source`] continue to apply. This entry
    /// applies a pressure impulse for `duration_s`; it does not solve a steady
    /// viscous duct or update a thermodynamic EOS.
    ///
    /// # Errors
    /// A non-finite or velocity-overconstraining boundary pressure, or any
    /// refusal documented by [`Self::project_density_volume_source`].
    #[allow(clippy::too_many_arguments)]
    pub fn project_density_volume_source_with_boundary_pressure(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        boundary_pressure_pa: [Option<f64>; 6],
        density_kg_m3: &[f64],
        target_outflow_m3_s: &[f64],
        controls: DensityProjectionControls,
        gate: &CancelGate,
    ) -> Result<DensityProjection, ChtError> {
        let boundary_pressure =
            ReservoirPressure::Planes(boundary_pressure_pa).admit(faces, gate)?;
        self.project_density_volume_source_with_reservoirs(
            domain,
            faces,
            boundary_pressure,
            density_kg_m3,
            target_outflow_m3_s,
            controls,
            gate,
        )
    }

    /// Project using a spatial reservoir-pressure profile on the open faces.
    ///
    /// Each native boundary slot acts through its own physical half-cell mass.
    /// A hydrostatic pressure profile can therefore balance gravity along a
    /// vertical opening; it is not replaced by one pressure for the whole plane.
    /// The first exposed reservoir slot supplies the reference independently
    /// in each connected fluid component. This preserves small differences
    /// even when components share a plane but have different large references.
    /// Reservoir work is `dt sum((p_face - p_reference) Q_out)` over the actual
    /// returned boundary fluxes, with the same reference used for source and
    /// held-boundary work. The original pressure/volume/energy owner is reused.
    ///
    /// For hydrostatic rest, first use [`Self::predict_density_body_force`]
    /// with the same density and duration, then project its returned field
    /// with the corresponding pressure head. Pressure is a mechanical impulse;
    /// this does not evolve reservoir density, thermodynamic state, or energy.
    /// All source, geometry and rollback requirements of
    /// [`Self::project_density_volume_source`] apply.
    ///
    /// # Errors
    /// Incorrect native face counts, a non-finite pressure, a pressure on an
    /// interior/solid/held face, or any native projection refusal. The input
    /// field and caller-owned profile are unchanged on every refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn project_density_volume_source_with_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        density_kg_m3: &[f64],
        target_outflow_m3_s: &[f64],
        controls: DensityProjectionControls,
        gate: &CancelGate,
    ) -> Result<DensityProjection, ChtError> {
        self.project_density_volume_source_with_reservoirs(
            domain,
            faces,
            BoundaryPressure::Profile(profile),
            density_kg_m3,
            target_outflow_m3_s,
            controls,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn project_density_volume_source_with_reservoirs<'a>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        boundary_pressure: BoundaryPressure<'a>,
        density_kg_m3: &'a [f64],
        target_outflow_m3_s: &'a [f64],
        controls: DensityProjectionControls,
        gate: &CancelGate,
    ) -> Result<DensityProjection, ChtError> {
        self.project_density_volume_source_with_pressure_readback(
            domain,
            faces,
            boundary_pressure,
            density_kg_m3,
            target_outflow_m3_s,
            controls,
            gate,
        )
        .map(|(projection, _)| projection)
    }

    #[allow(clippy::too_many_arguments)]
    fn project_density_volume_source_with_pressure_readback<'a>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        boundary_pressure: BoundaryPressure<'a>,
        density_kg_m3: &'a [f64],
        target_outflow_m3_s: &'a [f64],
        controls: DensityProjectionControls,
        gate: &CancelGate,
    ) -> Result<(DensityProjection, PressureReadback), ChtError> {
        let mut density = self.admit_density_projection(
            domain,
            faces,
            density_kg_m3,
            target_outflow_m3_s,
            controls,
            gate,
        )?;
        if let BoundaryPressure::Profile(profile) = boundary_pressure {
            profile.admit(self, domain, faces, gate)?;
        }
        density.boundary_pressure = boundary_pressure;
        let mut field = self.clone();
        let (projection, potential) = field.project_with_density(
            domain,
            faces,
            controls.relative_tolerance,
            Some(&density),
            gate,
        )?;
        let mut pressure = vec![0.0; domain.cell_count()];
        for (row, &c) in potential.cells.iter().enumerate() {
            poll_every(gate, row)?;
            let phi = potential.values[row];
            pressure[c] = density
                .uniform_density
                .map_or(phi, |rho| phi * (rho / controls.duration_s));
            finite("flow.pressure_correction_pa", pressure[c])?;
        }
        let report = audit(
            self, &field, domain, faces, &density, &pressure, &potential, projection, controls,
            gate,
        )?;
        for (row, &c) in potential.cells.iter().enumerate() {
            poll_every(gate, row)?;
            pressure[c] += potential.reference_pressure_pa[potential.component[row]];
            finite("flow.pressure_correction_pa", pressure[c])?;
        }
        poll(gate)?;
        Ok((
            DensityProjection {
                field,
                pressure_correction_pa: pressure,
                report,
            },
            PressureReadback {
                potential,
                uniform_pressure_scale: density
                    .uniform_density
                    .map(|rho| rho / controls.duration_s),
            },
        ))
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
            boundary_pressure: BoundaryPressure::Planes([0.0; 6]),
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
    potential: &ProjectionPotential,
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
    let mut open_work = Sum::default();
    for (row, &c) in potential.cells.iter().enumerate() {
        poll_every(gate, row)?;
        let reference_pa = potential.reference_pressure_pa[potential.component[row]];
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
            } else if nb.is_none() && faces[f] == FlowFace::Free {
                let (slot, _) = after.slot(domain, c, f);
                open_work.add(
                    controls.duration_s
                        * (density.boundary_pressure_pa(f, slot) - reference_pa)
                        * q1,
                );
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
        open_work.value,
    ] {
        finite("flow.projection_energy", value)?;
    }
    let energy_defect =
        change.value + loss.value - source_work.value + boundary_work.value + open_work.value;
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
        open_boundary_pressure_work_j: open_work.value,
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
