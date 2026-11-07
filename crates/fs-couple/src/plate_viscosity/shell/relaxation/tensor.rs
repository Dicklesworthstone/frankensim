//! Full source tensor memory on the original shear-deformable shell.
//!
//! Only thickness-normal traction is eliminated. All five physical retained
//! strains, including in-plane/transverse-shear coupling, use the complete
//! history-dependent constraint; branches are never condensed independently.

use super::*;
use crate::render::schedule::force::coupled::{ModalConnection, ModalPlaneStressPoint};
use crate::tensor_material::{PlaneStressPoint, PlaneStressRelaxation, TensorMaterialError};
use fs_material::visco_tensor::coordinates::{
    RelaxationCoordinateBudget, RelaxationCoordinates, RelaxationMemoryMap,
};
use fs_material::visco_tensor::source::ResolvedTensorMaterial;
use fs_material::visco_tensor::{TensorProjectionBudget, TensorProny};
use fs_math::det;
use fs_plate::shell::mitc::strain::{PreparedElementStrain, integration_rule};
use fs_plate::shell::mitc::{MitcShellModel, ShearShellSection};
use std::ops::Range;
use std::sync::Arc;

fn tensor_error(error: TensorMaterialError) -> AcousticRealizeError {
    match error {
        TensorMaterialError::Cancelled => AcousticRealizeError::Cancelled,
        other => owner(other),
    }
}

fn plate_error(error: fs_plate::PlateError) -> AcousticRealizeError {
    match error {
        fs_plate::PlateError::Cancelled => AcousticRealizeError::Cancelled,
        other => owner(other),
    }
}

/// One homogeneous source layer about the original facet reference surface.
/// The explicit proper rotation maps material axes to the facet frame.
///
/// The caller's positive shear correction is a work-conjugate congruence:
/// gamma_yz and gamma_zx enter the original law multiplied by sqrt(kappa).
/// This applies the correction once to shear work and consistently to mixed
/// strain terms. It defines an equivalent-single-layer shell, not a resolved
/// three-dimensional through-thickness shear-stress distribution.
#[derive(Clone, Debug)]
pub struct TensorShellMaterial {
    source: Arc<ResolvedTensorMaterial>,
    law: TensorProny,
    faces_m: [f64; 2],
    shear_correction: f64,
}

impl TensorShellMaterial {
    pub fn new(
        source: Arc<ResolvedTensorMaterial>,
        faces_m: [f64; 2],
        material_to_facet: [[f64; 3]; 3],
        shear_correction: f64,
    ) -> Result<Self, AcousticRealizeError> {
        let h = faces_m[1] - faces_m[0];
        if faces_m.iter().any(|v| !v.is_finite())
            || !h.is_finite()
            || h <= 0.0
            || !shear_correction.is_finite()
            || shear_correction <= 0.0
        {
            return Err(invalid(
                "tensor shell requires ordered finite faces and explicit positive shear correction",
            ));
        }
        let law = source.rotated_law(material_to_facet).map_err(owner)?;
        Ok(Self {
            source,
            law,
            faces_m,
            shear_correction,
        })
    }

    pub fn source(&self) -> &ResolvedTensorMaterial {
        &self.source
    }
    pub const fn material_law(&self) -> &TensorProny {
        &self.law
    }
    pub const fn faces_m(&self) -> [f64; 2] {
        self.faces_m
    }
    pub const fn shear_correction(&self) -> f64 {
        self.shear_correction
    }

    fn moments(&self) -> Result<[f64; 3], AcousticRealizeError> {
        let h = self.faces_m[1] - self.faces_m[0];
        let center = self.faces_m[0] + 0.5 * h;
        let moments = [h, h * center, h * (center * center + h * h / 12.0)];
        if moments.iter().any(|v| !v.is_finite()) || moments[2] <= 0.0 {
            return Err(invalid(
                "tensor shell thickness moments are unrepresentable",
            ));
        }
        Ok(moments)
    }
}

/// Explicit cold projection limits. Every original material history survives;
/// the selected time runtime separately admits its complete augmented state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TensorShellRelaxationLimits {
    pub max_modes: usize,
    pub max_material_points: usize,
    /// Evolved shared coordinates, after exact equal-time history sharing.
    pub max_memory_states: usize,
    pub max_projection_terms: usize,
}

/// An original area/thickness quadrature point and all of its source history.
#[derive(Clone, Debug)]
pub struct TensorShellMaterialPoint {
    facet: usize,
    barycentric: [f64; 3],
    thickness_coordinate_m: f64,
    memory: Range<usize>,
    memory_map: Option<RelaxationMemoryMap>,
    image: PlaneStressRelaxation,
}

impl TensorShellMaterialPoint {
    pub const fn facet(&self) -> usize {
        self.facet
    }
    pub const fn barycentric(&self) -> [f64; 3] {
        self.barycentric
    }
    pub const fn thickness_coordinate_m(&self) -> f64 {
        self.thickness_coordinate_m
    }
    pub fn memory_range(&self) -> Range<usize> {
        self.memory.clone()
    }
    pub const fn image(&self) -> &PlaneStressRelaxation {
        &self.image
    }
}

/// Complete source-bound Kelvin–Maxwell loss on the original supported MITC basis.
/// Equilibrium elasticity and inertia are validated against that model and are
/// counted once. No numerical drilling or fitted modal decay is material heat.
#[derive(Clone, Debug)]
pub struct TensorShellRelaxation {
    facets: Vec<TensorShellMaterial>,
    points: Vec<TensorShellMaterialPoint>,
    modal: ModalShellRelaxation,
    shared: RelaxationCoordinates,
    viscosity: ModalPlateViscosity,
    equilibrium_omega_squared: Vec<f64>,
    band_rad_s: (f64, f64),
}

impl TensorShellRelaxation {
    /// Integrate the original seven-point area rule and two thickness Gauss
    /// points. Membrane/bending cross work, both transverse shears and all
    /// branch-history coupling survive on the original condensed trial space.
    /// The unit-volume constitutive change of coordinates is prepared once per
    /// facet, then exactly projected. The existing equal-time coordinate owner
    /// shares driven histories without pole merging or a rank cutoff. Every
    /// original point history remains reconstructible for observation.
    ///
    /// A homogeneous source must reproduce each complete original equilibrium
    /// section and all three reference inertia moments at every facet corner.
    /// Sources with equilibrium in-plane/shear cross terms unsupported by the
    /// current native section refuse; their Maxwell increments may contain
    /// arbitrary cross terms. The mass basis is checked, never re-normalized.
    pub fn from_mitc_shell(
        model: &MitcShellModel,
        facets: &[TensorShellMaterial],
        modes: &[Vec<f64>],
        limits: TensorShellRelaxationLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let (mesh, directors, sections) = model.original_geometry();
        let n = modes.len();
        let count = facets
            .len()
            .checked_mul(14)
            .ok_or_else(|| invalid("tensor shell material point extent overflows"))?;
        if !(1..=64).contains(&limits.max_modes)
            || n == 0
            || n > limits.max_modes
            || !(1..=4096).contains(&limits.max_material_points)
            || limits.max_memory_states > 8192
            || limits.max_projection_terms == 0
            || facets.len() != mesh.tris.len()
            || facets.is_empty()
            || count > limits.max_material_points
            || modes.iter().any(|v| v.len() != model.pencil().free)
        {
            return Err(invalid(
                "tensor shell needs complete original facets, modes and bounded point/history/work capacities",
            ));
        }
        let states = facets
            .iter()
            .try_fold(0usize, |sum, facet| {
                // A thickness-normal dashpot retains one additional physical
                // strain history, even for a source without Maxwell branches.
                let normal = usize::from(
                    facet
                        .law
                        .kv_pa_s()
                        .is_some_and(|d| d[2].iter().any(|v| *v != 0.0)),
                );
                facet
                    .law
                    .branches()
                    .len()
                    .checked_mul(6)
                    .and_then(|r| r.checked_add(normal))
                    .and_then(|r| r.checked_mul(14))
                    .and_then(|r| sum.checked_add(r))
            })
            .ok_or_else(|| invalid("tensor shell material history extent overflows"))?;
        let setup = facets.len() as u128
            * (PreparedElementStrain::PREPARATION_WORK + n * PreparedElementStrain::REAL_FIELD_WORK)
                as u128
            + model.pencil().m.nnz() as u128 * n as u128
            + 16 * model.pencil().k.nnz() as u128 * n as u128
            + 17 * model.pencil().free as u128 * (n as u128).pow(2)
            + mesh.nodes.len() as u128 * 6 * n as u128
            + count as u128 * n as u128 * 256
            + count as u128 * (n as u128).pow(2)
            + 64 * (n as u128).pow(3);
        if setup >= limits.max_projection_terms as u128 {
            return Err(invalid(
                "tensor shell original geometry and modal pencil exceed projection work",
            ));
        }
        let work = limits.max_projection_terms - setup as usize;
        let factor_work = work / 3;
        if states.checked_mul(n).is_none_or(|v| v > factor_work) {
            return Err(invalid(
                "tensor shell original histories exceed their explicit projection work",
            ));
        }
        let point_work = (work - factor_work) / 2 / (count + facets.len());
        if point_work == 0 {
            return Err(invalid(
                "tensor shell constitutive projection requires positive work",
            ));
        }
        validate_mass_basis(model.pencil(), modes, gate)?;
        let equilibrium_omega_squared = project_original_stiffness(model, modes, gate)?;
        let physical = modes
            .iter()
            .map(|v| {
                poll(gate)?;
                model.expand(v).map_err(plate_error)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut rows = Vec::with_capacity(states);
        let mut times = Vec::with_capacity(states);
        let mut points = Vec::with_capacity(count);
        let mut viscosity = vec![fs_math::dd::Dd::ZERO; n * n];
        let mut band = (0.0_f64, f64::MAX);
        let rule = integration_rule();
        for (e, facet) in facets.iter().enumerate() {
            poll(gate)?;
            let unit: [[f64; 5]; 5] =
                core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
            let constitutive = PlaneStressPoint::new_shell(
                &facet.law,
                &unit,
                1.0,
                TensorProjectionBudget {
                    max_coordinates: 6,
                    max_scalar_terms: point_work,
                },
                gate,
            )
            .map_err(tensor_error)?
            .viscoelastic_image(point_work, gate)
            .map_err(tensor_error)?;
            validate_section(
                model,
                e,
                &sections[e],
                facet,
                constitutive.equilibrium_stiffness(),
            )?;
            let source_band = facet.law.band_rad_s();
            band = (band.0.max(source_band.0), band.1.min(source_band.1));
            let observer =
                PreparedElementStrain::new(mesh, e, directors, &sections[e], &mut || {
                    gate.is_requested()
                })
                .map_err(plate_error)?;
            let tri = mesh.tris[e];
            let fields = physical
                .iter()
                .map(|u| {
                    observer
                        .real_field(
                            &core::array::from_fn(|i| u[6 * tri[i / 6] + i % 6]),
                            &mut || gate.is_requested(),
                        )
                        .map_err(plate_error)
                })
                .collect::<Result<Vec<_>, _>>()?;
            let area = mesh.facet(e).map_err(plate_error)?.area_m2;
            let h = facet.faces_m[1] - facet.faces_m[0];
            let center = facet.faces_m[0] + 0.5 * h;
            let offset = 0.5 * h / det::sqrt(3.0);
            let shear = det::sqrt(facet.shear_correction);
            for &(lambda, weight) in &rule {
                let strains = fields
                    .iter()
                    .map(|field| field.at_barycentric(lambda).map_err(plate_error))
                    .collect::<Result<Vec<_>, _>>()?;
                for z in [center - offset, center + offset] {
                    poll(gate)?;
                    let columns: Vec<[f64; 5]> = strains
                        .iter()
                        .map(|s| {
                            [
                                s.membrane[0] + z * s.gradient_per_m[0],
                                s.membrane[1] + z * s.gradient_per_m[1],
                                s.membrane[2] + z * s.gradient_per_m[2],
                                shear * s.transverse_shear[1],
                                shear * s.transverse_shear[0],
                            ]
                        })
                        .collect();
                    let image = constitutive
                        .project_shell(
                            &columns,
                            area * weight * h * 0.5,
                            TensorProjectionBudget {
                                max_coordinates: n + 1,
                                max_scalar_terms: point_work,
                            },
                            gate,
                        )
                        .map_err(tensor_error)?;
                    let start = rows.len();
                    for (total, &value) in viscosity.iter_mut().zip(image.feedthrough_viscosity()) {
                        *total = *total + fs_math::dd::Dd::from_f64(value);
                    }
                    rows.extend_from_slice(image.coupling_rows());
                    times.extend(
                        image
                            .relaxation_rates_per_s()
                            .iter()
                            .map(|rate| rate.recip()),
                    );
                    points.push(TensorShellMaterialPoint {
                        facet: e,
                        barycentric: lambda,
                        thickness_coordinate_m: z,
                        memory: start..rows.len(),
                        memory_map: None,
                        image,
                    });
                }
            }
        }
        if band.0 != 0.0 || band.1 < band.0 || rows.len() != states || points.len() != count {
            return Err(invalid(
                "tensor shell needs complete zero-to-transient source support and original histories",
            ));
        }
        if equilibrium_omega_squared
            .iter()
            .any(|&value| det::sqrt(value) > band.1)
        {
            return Err(invalid(
                "original shell modal frequencies exceed source support",
            ));
        }
        let shared = RelaxationCoordinates::new(
            &rows,
            &times,
            n,
            RelaxationCoordinateBudget {
                maximum_states: limits.max_memory_states,
                maximum_scalar_terms: factor_work,
            },
            || gate.is_requested(),
        )
        .map_err(|error| match error {
            fs_material::visco_tensor::TensorPronyError::Cancelled => {
                AcousticRealizeError::Cancelled
            }
            other => owner(other),
        })?;
        // Reuse the original-unit compensated passive-work owner. The full
        // spatial matrix retains cross-modal loss, including both shears;
        // no diagonal decay or viscosity sampled at a modal frequency is used.
        let viscosity = ModalPlateViscosity::from_psd_modal_matrix(
            viscosity.into_iter().map(|value| value.to_f64()).collect(),
            n,
            1.0,
            gate,
        )?;
        for point in &mut points {
            poll(gate)?;
            point.memory_map = Some(shared.memory_map(point.memory.clone()).map_err(owner)?);
        }
        rows = shared.coupling_rows().to_vec();
        times = shared.relaxation_times_s().to_vec();
        poll(gate)?;
        Ok(Self {
            facets: facets.to_vec(),
            points,
            band_rad_s: band,
            shared,
            viscosity,
            equilibrium_omega_squared,
            modal: ModalShellRelaxation {
                projected: Vec::new(),
                rows,
                times,
                modes: n,
                limits: PlateRelaxationLimits {
                    max_modes: limits.max_modes,
                    max_memory_states: limits.max_memory_states,
                    max_projection_terms: limits.max_projection_terms,
                },
            },
        })
    }

    pub fn facets(&self) -> &[TensorShellMaterial] {
        &self.facets
    }
    pub fn material_points(&self) -> &[TensorShellMaterialPoint] {
        &self.points
    }
    pub const fn mode_count(&self) -> usize {
        self.modal.modes
    }
    pub fn coupling_rows_per_s(&self) -> &[Vec<f64>] {
        &self.modal.rows
    }
    pub fn relaxation_times_s(&self) -> &[f64] {
        &self.modal.times
    }
    pub const fn band_rad_s(&self) -> (f64, f64) {
        self.band_rad_s
    }
    pub const fn shared_coordinates(&self) -> &RelaxationCoordinates {
        &self.shared
    }
    /// Complete instantaneous reduced C [1/s], with its original passive
    /// channels. Normal-stress relaxation remains in the material histories.
    pub const fn kelvin_viscosity(&self) -> &ModalPlateViscosity {
        &self.viscosity
    }

    /// Include these dashpots in the original coupled system's construction,
    /// then pass their exact contiguous range to
    /// [`Self::attach_to_with_kelvin_connections`]. The reference resistance is
    /// only a channel normalization; the full source matrix is unchanged.
    pub fn kelvin_connections(&self, component: usize) -> Vec<ModalConnection> {
        self.viscosity.connections(component)
    }
    /// Original diagonal entries of Phi^T K Phi [1/s²], including the
    /// original numerical elastic regularization, with no material loss added.
    pub fn equilibrium_omega_squared(&self) -> &[f64] {
        &self.equilibrium_omega_squared
    }

    fn validate_component(
        &self,
        component: &fs_modal::time::ModalAcousticTimeModel,
    ) -> Result<(), AcousticRealizeError> {
        if component.modes().len() != self.mode_count()
            || component.modes().iter().enumerate().any(|(i, mode)| {
                let actual = mode.angular_frequency_rad_s * mode.angular_frequency_rad_s;
                let expected = self.equilibrium_omega_squared[i];
                !actual.is_finite()
                    || (actual - expected).abs() > 1e-8 * expected
                    || component.damping_rate_per_s(i) != 0.0
            })
        {
            return Err(invalid(
                "tensor shell target must retain its original equilibrium frequencies and zero intrinsic damping",
            ));
        }
        Ok(())
    }

    /// Reconstruct every original local relaxation coordinate from the shared
    /// runtime state. `memory_range` addresses this original layout, not a slice
    /// of the evolved shared state.
    pub fn point_memory(
        &self,
        point: usize,
        memory: &[f64],
    ) -> Result<Vec<f64>, AcousticRealizeError> {
        let point = self
            .points
            .get(point)
            .ok_or_else(|| invalid("unknown tensor shell material point"))?;
        point
            .memory_map
            .as_ref()
            .ok_or_else(|| invalid("tensor shell original history map is absent"))?
            .reconstruct(memory)
            .map_err(owner)
    }

    /// Existing pHS displacement-only branches. Callers of this low-level
    /// bridge must also include `kelvin_viscosity()` in their physical
    /// resistance and retain source-band and accepted-strain admission.
    pub fn branches(
        &self,
        first_mode: usize,
        total_modes: usize,
        gate: &CancelGate,
    ) -> Result<Vec<RelaxationBranch>, AcousticRealizeError> {
        self.modal.branches(first_mode, total_modes, gate)
    }

    /// Reconstruct the complete original tensor at one integration point.
    /// The returned shear is the documented energy-equivalent corrected strain.
    /// This enforces the source domain at the point; it is not a whole-facet or
    /// between-step maximum-strain certificate.
    pub fn point_strain(
        &self,
        point: usize,
        displacement: &[f64],
        memory: &[f64],
    ) -> Result<fs_material::Voigt, AcousticRealizeError> {
        if memory.len() != self.modal.rows.len() || memory.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "tensor shell observation needs its complete finite history",
            ));
        }
        let history = self.point_memory(point, memory)?;
        let point = self
            .points
            .get(point)
            .ok_or_else(|| invalid("unknown tensor shell material point"))?;
        point
            .image
            .strain(displacement, &history)
            .map_err(tensor_error)
    }

    /// Attach to the existing cold exact modal owner, including every original
    /// quadrature source-domain check. Previous material memories are preserved.
    /// The original runtime admits the augmented stiffness/clock/work capacity;
    /// its complete linear frequency bound must fit every source's band.
    ///
    /// The component's equilibrium frequencies and zero intrinsic modal damping
    /// are checked. Spatial mode signs and orientation within degenerate modal
    /// subspaces remain caller obligations: the generic runtime has no shapes.
    /// This is one bounded source-point binding, as required by
    /// the existing coupled owner. Point checks do not certify unsampled spatial
    /// or temporal strain peaks. No heat is inferred from numerical energy loss.
    pub fn attach_to(
        &self,
        system: CoupledModalSystem,
        component: usize,
        initial: InitialPlateMemory,
        gate: &CancelGate,
    ) -> Result<CoupledModalSystem, AcousticRealizeError> {
        if self.viscosity.channel_count() != 0 {
            return Err(invalid(
                "tensor shell Kelvin channels require an explicit original connection range",
            ));
        }
        self.attach_to_with_kelvin_connections(system, component, 0..0, initial, gate)
    }

    /// Bind complete source memory after constructing the system with
    /// `kelvin_connections(component)`. Verify every original full modal
    /// column and coefficient before attaching any history, so missing,
    /// rescaled or unrelated dashpots cannot silently discard source loss.
    /// Kelvin work remains in the runtime's connection heat account; internal
    /// normal and Maxwell relaxation use its separate material heat account.
    ///
    /// With retained normal viscosity, a nonzero displaced body's zero source
    /// history is not the transformed image's zero history. Such an unrelaxed
    /// request refuses; relaxed initialization or a zero displacement is exact.
    pub fn attach_to_with_kelvin_connections(
        &self,
        system: CoupledModalSystem,
        component: usize,
        kelvin_connections: Range<usize>,
        initial: InitialPlateMemory,
        gate: &CancelGate,
    ) -> Result<CoupledModalSystem, AcousticRealizeError> {
        poll(gate)?;
        self.validate_component(
            system
                .components()
                .get(component)
                .ok_or_else(|| invalid("tensor shell component is absent"))?,
        )?;
        let start = system
            .components()
            .get(..component)
            .ok_or_else(|| invalid("tensor shell component is absent"))?
            .iter()
            .map(|m| m.modes().len())
            .sum::<usize>();
        if kelvin_connections.start > kelvin_connections.end
            || kelvin_connections.len() != self.viscosity.channel_count()
            || system
                .constant_loss_connections()
                .iter()
                .any(|anchor| kelvin_connections.contains(&anchor.connection))
        {
            return Err(invalid(
                "tensor shell requires every original Kelvin channel without constant-loss substitution",
            ));
        }
        for (channel, connection) in kelvin_connections.enumerate() {
            poll(gate)?;
            let actual = system
                .connection_column(connection)
                .ok_or_else(|| invalid("tensor shell Kelvin connection is absent"))?;
            let law = system
                .connection_coefficients(connection)
                .ok_or_else(|| invalid("tensor shell Kelvin law is absent"))?;
            self.validate_kelvin_channel(
                channel,
                start,
                actual,
                law.stiffness_n_m,
                law.damping_n_s_m,
                law.rest_extension_m,
            )?;
        }
        if initial == InitialPlateMemory::Unrelaxed
            && self
                .points
                .iter()
                .any(|point| point.image.point().normal_strain_coordinates() != 0)
            && system.components()[component]
                .states()
                .iter()
                .any(|state| state.displacement_m_sqrt_kg != 0.0)
        {
            return Err(invalid(
                "displaced shell with normal Kelvin history requires relaxed source initialization",
            ));
        }
        let memory_start = system.material_memory().len();
        let system = self.modal.attach_to(system, component, initial, gate)?;
        if system.linear_frequency_upper_bound_rad_s() > self.band_rad_s.1 {
            return Err(invalid(
                "complete coupled spectrum exceeds tensor shell source support",
            ));
        }
        let bindings = self
            .points
            .iter()
            .map(|point| ModalPlaneStressPoint {
                modes: start..start + self.mode_count(),
                memories: memory_start..memory_start + self.modal.rows.len(),
                memory_map: point.memory_map.clone(),
                // A conservative correction also bounds unscaled transverse
                // engineering strains when kappa is below one.
                strain_scale: (1.0 / det::sqrt(self.facets[point.facet].shear_correction)).max(1.0),
                image: point.image.clone(),
            })
            .collect();
        system
            .with_plane_stress_material_points(bindings, gate)
            .map_err(owner)
    }

    fn validate_kelvin_channel(
        &self,
        channel: usize,
        first_mode: usize,
        actual: &[f64],
        stiffness: f64,
        damping: f64,
        rest_extension: f64,
    ) -> Result<(), AcousticRealizeError> {
        let expected = self
            .viscosity
            .channel_shapes()
            .get(channel)
            .ok_or_else(|| invalid("unknown tensor shell Kelvin channel"))?;
        let end = first_mode
            .checked_add(self.mode_count())
            .ok_or_else(|| invalid("tensor shell Kelvin mode extent overflows"))?;
        if end > actual.len()
            || stiffness != 0.0
            || rest_extension != 0.0
            || damping != self.viscosity.reference_resistance_n_s_m()
            || actual.iter().enumerate().any(|(i, &value)| {
                value
                    != if (first_mode..end).contains(&i) {
                        expected[i - first_mode]
                    } else {
                        0.0
                    }
            })
        {
            return Err(invalid(
                "tensor shell Kelvin connection must retain its complete original source work",
            ));
        }
        Ok(())
    }
}

// The scalar-oscillator runtime requires a diagonal ORIGINAL equilibrium
// pencil, not merely an arbitrary mass-orthonormal basis. Keep its values for
// attachment admission; compensated sparse products avoid erasing weak work.
fn project_original_stiffness(
    model: &MitcShellModel,
    modes: &[Vec<f64>],
    gate: &CancelGate,
) -> Result<Vec<f64>, AcousticRealizeError> {
    use fs_math::dd::Dd;
    let n = modes.len();
    let mut projected = vec![0.0; n * n];
    let mut weighted = vec![Dd::ZERO; model.pencil().free];
    for (j, mode) in modes.iter().enumerate() {
        poll(gate)?;
        for (i, value) in weighted.iter_mut().enumerate() {
            if i % 256 == 0 {
                poll(gate)?;
            }
            let (columns, values) = model.pencil().k.row(i);
            *value = columns
                .iter()
                .zip(values)
                .fold(Dd::ZERO, |sum, (&column, &coefficient)| {
                    sum + Dd::from_f64(coefficient) * Dd::from_f64(mode[column])
                });
        }
        for (i, shape) in modes.iter().enumerate() {
            poll(gate)?;
            projected[n * i + j] = shape
                .iter()
                .zip(&weighted)
                .fold(Dd::ZERO, |sum, (&x, &force)| sum + Dd::from_f64(x) * force)
                .to_f64();
        }
    }
    let diagonal: Vec<_> = (0..n).map(|i| projected[n * i + i]).collect();
    if projected.iter().any(|v| !v.is_finite())
        || diagonal.iter().any(|v| *v < 0.0)
        || (0..n).any(|i| {
            (0..n).any(|j| {
                i != j
                    && projected[n * i + j].abs()
                        > 1e-8 * det::sqrt(diagonal[i]) * det::sqrt(diagonal[j])
            })
        })
    {
        return Err(invalid(
            "tensor shell basis must diagonalize the original nonnegative equilibrium stiffness",
        ));
    }
    Ok(diagonal)
}

fn validate_section(
    model: &MitcShellModel,
    e: usize,
    section: &ShearShellSection,
    source: &TensorShellMaterial,
    c: &[f64],
) -> Result<(), AcousticRealizeError> {
    let moments = source.moments()?;
    let (a, b, d, shear) = section.coupled_matrices();
    let close = |actual: f64, expected: f64, scale: f64| {
        actual.is_finite()
            && expected.is_finite()
            && scale.is_finite()
            && scale > 0.0
            && (actual - expected).abs() <= 1e-9 * scale
    };
    for i in 0..3 {
        for j in 0..3 {
            let scale = det::sqrt(c[5 * i + i]) * det::sqrt(c[5 * j + j]);
            if !close(a[3 * i + j], moments[0] * c[5 * i + j], moments[0] * scale)
                || !close(
                    b[3 * i + j],
                    moments[1] * c[5 * i + j],
                    det::sqrt(moments[0]) * det::sqrt(moments[2]) * scale,
                )
                || !close(d[3 * i + j], moments[2] * c[5 * i + j], moments[2] * scale)
            {
                return Err(invalid(
                    "tensor shell source does not reproduce the original equilibrium A/B/D section",
                ));
            }
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            let (u, v) = (4 - i, 4 - j);
            let factor = moments[0] * source.shear_correction;
            let scale = factor * det::sqrt(c[5 * u + u]) * det::sqrt(c[5 * v + v]);
            if !close(shear[i][j], factor * c[5 * u + v], scale) {
                return Err(invalid(
                    "tensor shell source does not reproduce original transverse-shear work",
                ));
            }
        }
    }
    for i in 0..3 {
        for j in 3..5 {
            if !close(
                c[5 * i + j],
                0.0,
                det::sqrt(c[5 * i + i]) * det::sqrt(c[5 * j + j]),
            ) {
                return Err(invalid(
                    "original shell equilibrium cannot represent source in-plane/shear coupling",
                ));
            }
        }
    }
    let density = source.source.density_kg_m3();
    let expected = moments.map(|v| density * v);
    let scales = [
        expected[0],
        det::sqrt(expected[0]) * det::sqrt(expected[2]),
        expected[2],
    ];
    for inertia in model
        .reference_inertia_moments(e)
        .ok_or_else(|| invalid("original shell inertia is absent"))?
    {
        let actual = [
            inertia.areal_mass_kg_m2(),
            inertia.mass_first_moment_kg_m(),
            inertia.rotary_inertia_kg(),
        ];
        if (0..3).any(|i| !close(actual[i], expected[i], scales[i])) {
            return Err(invalid(
                "tensor shell source must reproduce every original reference inertia moment",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
