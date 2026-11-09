//! Common-diffusivity mixing on accepted, fixed carrier inventories.
//!
//! The scalar law `j_k = -rho D grad(Y_k)` transports each original constituent
//! mass with the SAME diffusivity. The binary entry retains complementary Fick
//! fluxes, including its trace-safe arithmetic. A common diffusivity is also
//! appropriate for an explicitly authored turbulent gradient-diffusion model;
//! it is not a general multicomponent molecular diffusion law. D is prescribed,
//! never inferred from viscosity or a constituent index by this kernel.
//! Half-cell resistances give one shared flux across a spatial coefficient jump.
//! Each face has an exact positive exponential flow; forward and reverse half
//! sweeps are a symmetric second-order splitting for frozen linear coefficients.
//! A large stable step is not an accuracy certificate: refine physical substeps.
//!
//! The mass-only entry is an isothermal fixed-carrier constitutive component.
//! The enthalpy extension transports original species partial enthalpies; its
//! gas adapter reconstructs temperature and material properties from the same
//! inventories. Neither entry evolves pressure or mixture volume.
//! Unequal molecular species diffusivities, nonideal thermodynamic factors and
//! Soret/barodiffusion require their own original constitutive owners.
//!
//! Fick mass-frame convention: COMSOL Chemical Reaction Engineering, eq. 3-86,
//! https://doc.comsol.com/6.3/doc/com.comsol.help.chem/chem_ug_chemsptrans.08.241.html.
//! Symmetric splitting: Strang (1968), https://doi.org/10.1137/0705041.

use super::{
    CancelGate, ChtError, MAX_MATERIAL_SPECIES, MaterialAdvectionStep, Sum, nonnegative, poll,
    poll_every, ratio,
};
use crate::conjugate::VoxelDomain;
use fs_math::det;

mod boundary;
pub use boundary::{
    AdaptiveBoundarySpeciesDiffusionStep, BoundarySpeciesDiffusionStep, SpeciesBoundaryDiffusionReport,
    SpeciesBoundaryTransfer, SpeciesDiffusionBoundary,
};
mod enthalpy;
pub use enthalpy::{SpeciesEnthalpyDiffusionReport, SpeciesEnthalpyDiffusionStep};
pub use enthalpy::boundary::adaptive::{
    AdaptiveBoundarySpeciesEnthalpyDiffusionStep, BoundarySpeciesEnthalpyAccuracy,
};
pub use enthalpy::boundary::{
    BoundarySpeciesEnthalpyDiffusionStep, SpeciesBoundaryEnthalpyDiffusionReport,
    SpeciesBoundaryEnthalpyTransfer,
};
pub use enthalpy::adaptive::{
    AdaptiveSpeciesEnthalpyDiffusionControls, AdaptiveSpeciesEnthalpyDiffusionReport,
    AdaptiveSpeciesEnthalpyDiffusionStep, SpeciesEnthalpyDiffusionModel,
};
mod gas;
pub use gas::{
    BoundaryThermochemicalSpeciesDiffusionStep, GasSpeciesEnthalpyBoundarySource,
    GasSpeciesDiffusivity, GasSpeciesEnthalpyField, GasSpeciesEnthalpySource,
    ThermochemicalSpeciesDiffusionReport, ThermochemicalSpeciesDiffusionStep,
};

/// Explicit diffusion work, storage and independent constituent-mass budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeciesDiffusionControls<const N: usize> {
    /// Equal physical substeps, each a forward/reverse pair of half sweeps.
    /// Positive; refine this count to resolve splitting error.
    pub substeps: usize,
    /// Bound on two times the active interior AND authored reservoir face count
    /// times substeps. Each geometric
    /// pair update evolves N constituents, so component work is bounded by
    /// N times this value; N must be in 1..=MAX_MATERIAL_SPECIES.
    pub maximum_pair_updates: usize,
    /// Conservative payload bound for cell scratch, output and prepared faces.
    /// Existing carrier storage and allocator overhead are additional.
    pub maximum_storage_bytes: usize,
    /// Per-species local inventory-representation defects plus the absolute
    /// mass-balance defect [kg], including authored reservoir transfers when
    /// present. Trace species retain their own allowance.
    pub maximum_mass_defect_kg: [f64; N],
    /// Sum of absolute initial and final per-cell partition defects [kg],
    /// measured against the SAME immutable accepted carrier mass. Boundary
    /// entries also charge every map's absolute summed carrier transfer.
    pub maximum_partition_defect_kg: f64,
}

/// Physical readback of one complete diffusion interval.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesDiffusionReport<const N: usize> {
    pub duration_s: f64,
    pub substeps: usize,
    /// Faces with positive conductance, including authored reservoir faces in
    /// the boundary entries and faces with exactly zero advective mass flux.
    pub diffusive_faces: usize,
    pub pair_updates: usize,
    /// Checked N times pair_updates: constituent pair maps evaluated.
    pub component_updates: usize,
    /// Largest half-map relaxation: 1-exp(-G*dt_half*(1/M_l+1/M_r))
    /// internally, or 1-exp(-G*dt_half/M) at a fixed reservoir. This measures
    /// stiffness, not the temporal error of a sweep.
    pub maximum_pair_relaxation: f64,
    pub mass_before_kg: [f64; N],
    pub mass_after_kg: [f64; N],
    pub mass_change_kg: [f64; N],
    pub absolute_mass_defect_kg: [f64; N],
    pub initial_partition_defect_kg: f64,
    pub endpoint_partition_defect_kg: f64,
}

/// All authoritative constituent inventories, with no carrier or input mutation.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesDiffusionStep<const N: usize> {
    pub species_mass_kg: Vec<[f64; N]>,
    pub report: SpeciesDiffusionReport<N>,
}

/// Automatic refinement of the existing positive, conservative face sweeps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveSpeciesDiffusionControls<const N: usize> {
    /// Original physical budgets. `substeps` is the initial coarse resolution;
    /// the pair-update limit covers ALL coarse and fine solves, including
    /// rejected resolutions. Storage includes the retained comparison state.
    pub solver: SpeciesDiffusionControls<N>,
    /// Positive absolute mass-fraction tolerance for EACH constituent. A trace
    /// species does not inherit the dominant species' absolute allowance.
    pub absolute_mass_fraction_tolerance: [f64; N],
    /// Nonnegative relative allowance on max(initial, fine) mass fraction.
    pub relative_tolerance: f64,
    /// Maximum fine solves. Each doubles the previous coarse substep count.
    pub maximum_refinements: usize,
}

/// Complete diffusion interval accepted by a temporal refinement estimate.
#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveSpeciesDiffusionStep<const N: usize> {
    /// Original fine-grid-in-time inventories and independent physical account.
    /// No Richardson extrapolation, clipping or renormalization is applied.
    pub step: SpeciesDiffusionStep<N>,
    /// Per-species maximum over cells of |Y_fine - Y_coarse| / 3. This is a
    /// second-order temporal estimate, not a certified global error bound.
    pub estimated_mass_fraction_error: [f64; N],
    /// Largest estimate divided by its absolute-plus-relative allowance.
    pub maximum_error_ratio: f64,
    /// Fine solves performed, including resolutions that were not accepted.
    /// Zero when there are no conducting faces.
    pub refinements: usize,
    /// Aggregate geometric pair work, including the initial coarse solve.
    pub total_pair_updates: usize,
    /// Aggregate constituent pair work, including discarded coarse solutions.
    pub total_component_updates: usize,
}

/// Original binary Fick controls, with unchanged geometric work semantics.
pub type BinarySpeciesDiffusionControls = SpeciesDiffusionControls<2>;
/// Original binary Fick mass and partition account.
pub type BinarySpeciesDiffusionReport = SpeciesDiffusionReport<2>;
/// Original binary Fick constituent result.
pub type BinarySpeciesDiffusionStep = SpeciesDiffusionStep<2>;

#[derive(Clone, Copy)]
struct Face {
    cells: [usize; 2],
    /// H*(1-exp(-G*dt_half/H)), H = M_l*M_r/(M_l+M_r) [kg].
    exchange_capacity_kg: f64,
}

#[derive(Clone, Copy)]
struct Inventory<const N: usize> {
    mass: [Sum; N],
    increment: [Sum; N],
}

trait DiffusionObserver<const N: usize> {
    fn storage_bytes(&self, _: usize, _: usize) -> Option<usize> { Some(0) }
    fn prepare_storage(&mut self, _: usize, _: usize) {}
    fn prepare(&mut self, _: [usize; 2], _: usize) -> Result<(), ChtError> { Ok(()) }
    fn transferred(&mut self, _: usize, _: [f64; N]) -> Result<(), ChtError> { Ok(()) }
    fn boundary_faces(&mut self, _: &CancelGate) -> Result<usize, ChtError> { Ok(0) }
    fn prepare_boundaries(&mut self, _: f64, _: &CancelGate) -> Result<(), ChtError> { Ok(()) }
    fn boundary_half_sweep(
        &mut self, _: &mut [Inventory<N>], _: bool, _: &CancelGate,
    ) -> Result<(), ChtError> { Ok(()) }
    fn boundary_mass_transfer(&self) -> [f64; N] { [0.0; N] }
    fn boundary_carrier_defect(&self) -> f64 { 0.0 }
    fn maximum_boundary_relaxation(&self) -> f64 { 0.0 }
}

impl<const N: usize> DiffusionObserver<N> for () {}

fn refuse(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.species_diffusion",
        reason: reason.into(),
    }
}

fn positive(value: f64) -> Result<f64, ChtError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(refuse(
            "positive diffusion geometry or coefficient is unrepresentable",
        ))
    }
}

impl MaterialAdvectionStep {
    /// The shared native half-cell resistance, after carrier/domain admission.
    /// Both supplied diffusivities must be strictly positive.
    pub(in crate::conjugate) fn diffusion_face_conductance(
        &self,
        domain: &VoxelDomain,
        cells: [usize; 2],
        axis: usize,
        diffusivity_m2_s: [f64; 2],
    ) -> Result<f64, ChtError> {
        let [c, next] = cells;
        let rho_left = positive(self.mass_kg[c] / domain.volume(c))?;
        let rho_right = positive(self.mass_kg[next] / domain.volume(next))?;
        let left =
            positive((0.5 * domain.widths(c)[axis]) / positive(rho_left * diffusivity_m2_s[0])?)?;
        let right = positive(
            (0.5 * domain.widths(next)[axis]) / positive(rho_right * diffusivity_m2_s[1])?,
        )?;
        let largest = left.max(right);
        positive((domain.face_area(c, axis) / largest) / (1.0 + left.min(right) / largest))
    }

    /// Diffuse two constituent masses on this step's accepted carrier field.
    ///
    /// Supply post-advection masses, e.g. from advect_species::<2>, or the
    /// previous diffusion endpoint. Their sum must partition this SAME accepted
    /// carrier mass. The original public carrier mass is checked for edits first.
    /// The domain must be the original accepted geometry; its dimensions and
    /// fluid mask are checked, while identical physical widths remain the caller's
    /// obligation as in the original pressure/material path.
    ///
    /// Each fluid cell supplies one nonnegative D [m²/s], held fixed throughout
    /// this interval. Density is the accepted mass divided by native cell volume.
    /// For neighboring cells, G=A/(h_l/(2*rho_l*D_l)+h_r/(2*rho_r*D_r)) [kg/s].
    /// This series resistance supports graded Cartesian cells and spatial
    /// coefficient jumps. Zero D blocks the adjacent face; solids and all
    /// exterior faces are impermeable. No boundary concentration is invented.
    ///
    /// Paired transfers update BOTH constituent inventories; the second is never
    /// reconstructed as carrier minus first, and neither is clipped/normalized.
    /// Each pair approaches its mass-weighted equilibrium without overshoot.
    /// Reusable prepared rates keep exponentials outside the sweep loop; work is
    /// O(cells + faces*substeps), storage O(cells + faces), independent of stiffness.
    /// Physical substep and mesh refinement remain necessary for accuracy.
    ///
    /// No thermodynamic heat account is supplied by anonymous species indices.
    /// This does not promote an authored D into a sourced material law or turn
    /// binary Fick transport into arbitrary multicomponent diffusion.
    /// Failure/cancellation publishes no partial state and changes no input.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_binary_species(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; 2]],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        controls: BinarySpeciesDiffusionControls,
        gate: &CancelGate,
    ) -> Result<BinarySpeciesDiffusionStep, ChtError> {
        self.diffuse_species(
            domain,
            species_mass_kg,
            diffusivity_m2_s,
            duration_s,
            controls,
            gate,
        )
    }

    /// Mix 1..=MAX_MATERIAL_SPECIES original inventories with one common D.
    ///
    /// Every species obeys `j_k = -rho D grad(Y_k)` on this accepted carrier.
    /// A partition of unity therefore has zero net diffusive carrier flux.
    /// Each constituent retains its OWN inventory, positive pair relaxation,
    /// and mass-defect budget. No constituent is selected as a remainder, and
    /// no negative mass is clipped or composition normalized. The binary case
    /// retains its complementary trace-safe face arithmetic exactly.
    ///
    /// Equal diffusivity represents an explicit common turbulent mixing law,
    /// or a molecular model whose equal-diffusivity assumptions are justified
    /// separately. It cannot represent arbitrary molecular diffusion matrices.
    /// The caller supplies D in m²/s, for example an independently declared
    /// `nu_t / Sc_t`; this kernel does not choose a turbulent Schmidt number.
    ///
    /// Geometry, frozen-density, closed-boundary, second-order splitting and
    /// failure-atomic semantics match [`Self::diffuse_binary_species`]. Shared
    /// face rates are prepared once: work is O(N*(cells + faces*substeps)), and
    /// storage is O(N*cells + faces). Both geometric and component counts use
    /// checked arithmetic. The independent carrier-partition account includes
    /// all N constituent inventories, including trace species.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_species<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        controls: SpeciesDiffusionControls<N>,
        gate: &CancelGate,
    ) -> Result<SpeciesDiffusionStep<N>, ChtError> {
        self.diffuse_species_checked(
            domain,
            species_mass_kg,
            diffusivity_m2_s,
            duration_s,
            controls,
            gate,
            |_| poll(gate),
        )
    }

    /// Resolve frozen-coefficient diffusion with automatic temporal refinement.
    ///
    /// Compare the complete interval at n and 2n physical substeps, retaining
    /// the finer conservative solution. A rejected fine solution becomes the
    /// next coarse comparison; it is not recomputed. Acceptance requires every
    /// cell and constituent to satisfy
    /// `|Y_fine-Y_coarse|/3 <= atol[k] + rtol*max(Y_initial,Y_fine)`.
    /// All solves start from the SAME supplied original inventories, so this
    /// comparison measures error over the requested interval, not an
    /// accumulated sum of unrelated local tolerances.
    ///
    /// The factor three follows the second-order symmetric face composition
    /// in its asymptotic regime. This is an estimate, not a rigorous bound;
    /// spatial, constitutive and advection-splitting errors are not estimated.
    /// Physical mass/partition budgets still apply independently to every
    /// candidate. Work limits include rejected refinements, and peak payload
    /// includes the simultaneously retained coarse result. Failure or
    /// cancellation publishes nothing and changes no input.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_species_adaptive<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        controls: AdaptiveSpeciesDiffusionControls<N>,
        gate: &CancelGate,
    ) -> Result<AdaptiveSpeciesDiffusionStep<N>, ChtError> {
        self.diffuse_species_adaptive_checked(
            domain, species_mass_kg, diffusivity_m2_s, duration_s, controls, gate,
            |_| poll(gate),
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn diffuse_species_adaptive_checked<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        controls: AdaptiveSpeciesDiffusionControls<N>,
        gate: &CancelGate,
        checkpoint: impl FnMut(usize) -> Result<(), ChtError>,
    ) -> Result<AdaptiveSpeciesDiffusionStep<N>, ChtError> {
        self.refine_species(
            domain, species_mass_kg, controls, gate, checkpoint, 0,
            |solver| Ok((self.diffuse_species(
                domain, species_mass_kg, diffusivity_m2_s, duration_s, solver, gate,
            )?, ())),
            |_, _| Ok(0.0),
        ).map(|(step, ())| step)
    }

    // Shared complete-interval refinement. Extra retained payload and transfer
    // comparisons let boundary owners use the same acceptance/work semantics.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn refine_species<const N: usize, T>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        controls: AdaptiveSpeciesDiffusionControls<N>,
        gate: &CancelGate,
        mut checkpoint: impl FnMut(usize) -> Result<(), ChtError>,
        extra_retained_bytes: usize,
        mut solve: impl FnMut(SpeciesDiffusionControls<N>)
            -> Result<(SpeciesDiffusionStep<N>, T), ChtError>,
        mut compare_transfers: impl FnMut(&T, &T) -> Result<f64, ChtError>,
    ) -> Result<(AdaptiveSpeciesDiffusionStep<N>, T), ChtError> {
        poll(gate)?;
        if N == 0 || N > MAX_MATERIAL_SPECIES || controls.maximum_refinements == 0 {
            return Err(refuse("adaptive diffusion component or refinement count is invalid"));
        }
        nonnegative(controls.relative_tolerance)?;
        for tolerance in controls.absolute_mass_fraction_tolerance {
            positive(tolerance)?;
        }
        // The original solver already budgets its own scratch, output, faces
        // and fixed arrays. Add the retained coarse output and this controller's
        // fixed result/estimate state before admitting either solve.
        let retained_bytes = domain.cell_count()
            .checked_mul(core::mem::size_of::<[f64; N]>())
            .and_then(|n| n.checked_add(core::mem::size_of::<SpeciesDiffusionStep<N>>()))
            .and_then(|n| n.checked_add(core::mem::size_of::<AdaptiveSpeciesDiffusionStep<N>>()))
            .and_then(|n| n.checked_add(4 * core::mem::size_of::<[f64; N]>() + 512))
            .and_then(|n| n.checked_add(extra_retained_bytes))
            .ok_or_else(|| refuse("adaptive diffusion storage extent overflows"))?;
        let solver_bytes = controls.solver.maximum_storage_bytes.checked_sub(retained_bytes)
            .ok_or_else(|| refuse("adaptive diffusion comparison storage budget exhausted"))?;
        let mut solver = SpeciesDiffusionControls {
            maximum_storage_bytes: solver_bytes,
            ..controls.solver
        };
        let (mut coarse, mut coarse_transfers) = solve(solver)?;
        let mut total_pair_updates = coarse.report.pair_updates;
        if coarse.report.diffusive_faces == 0 {
            poll(gate)?;
            return Ok((AdaptiveSpeciesDiffusionStep {
                step: coarse,
                estimated_mass_fraction_error: [0.0; N],
                maximum_error_ratio: 0.0,
                refinements: 0,
                total_pair_updates,
                total_component_updates: 0,
            }, coarse_transfers));
        }
        for refinement in 1..=controls.maximum_refinements {
            checkpoint(refinement)?;
            solver.substeps = solver.substeps.checked_mul(2)
                .ok_or_else(|| refuse("adaptive diffusion substep extent overflows"))?;
            solver.maximum_pair_updates = controls.solver.maximum_pair_updates
                .checked_sub(total_pair_updates)
                .ok_or_else(|| refuse("adaptive diffusion cumulative work budget exhausted"))?;
            let fine_pairs = coarse.report.diffusive_faces.checked_mul(2)
                .and_then(|n| n.checked_mul(solver.substeps))
                .ok_or_else(|| refuse("adaptive diffusion work extent overflows"))?;
            let next_total_pairs = total_pair_updates.checked_add(fine_pairs)
                .ok_or_else(|| refuse("adaptive diffusion cumulative work extent overflows"))?;
            let total_component_updates = next_total_pairs.checked_mul(N)
                .ok_or_else(|| refuse("adaptive diffusion component work extent overflows"))?;
            let (fine, fine_transfers) = solve(solver)?;
            total_pair_updates = next_total_pairs;
            let mut estimated_mass_fraction_error = [0.0_f64; N];
            let mut maximum_error_ratio = compare_transfers(&coarse_transfers, &fine_transfers)?;
            for c in 0..domain.cell_count() {
                poll_every(gate, c)?;
                if !domain.is_fluid(c) {
                    continue;
                }
                for (k, estimate) in estimated_mass_fraction_error.iter_mut().enumerate() {
                    // Subtract original masses first: dividing separately can
                    // erase a small but representable concentration change.
                    let error = (fine.species_mass_kg[c][k] - coarse.species_mass_kg[c][k])
                        .abs() / self.mass_kg[c] / 3.0;
                    let scale = ratio(species_mass_kg[c][k], self.mass_kg[c])?
                        .max(ratio(fine.species_mass_kg[c][k], self.mass_kg[c])?);
                    let tolerance = positive(controls.absolute_mass_fraction_tolerance[k]
                        + controls.relative_tolerance * scale)?;
                    nonnegative(error)?;
                    *estimate = estimate.max(error);
                    // An overflowing ratio is a valid rejection, never a
                    // reason to publish an insufficiently resolved state.
                    maximum_error_ratio = maximum_error_ratio.max(error / tolerance);
                }
            }
            if maximum_error_ratio <= 1.0 {
                poll(gate)?;
                return Ok((AdaptiveSpeciesDiffusionStep {
                    step: fine,
                    estimated_mass_fraction_error,
                    maximum_error_ratio,
                    refinements: refinement,
                    total_pair_updates,
                    total_component_updates,
                }, fine_transfers));
            }
            coarse = fine;
            coarse_transfers = fine_transfers;
        }
        Err(refuse("adaptive diffusion temporal tolerance is unresolved at its refinement limit"))
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn diffuse_species_checked<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        controls: SpeciesDiffusionControls<N>,
        gate: &CancelGate,
        mut checkpoint: impl FnMut(usize) -> Result<(), ChtError>,
    ) -> Result<SpeciesDiffusionStep<N>, ChtError> {
        self.diffuse_species_observed(
            domain, species_mass_kg, diffusivity_m2_s, duration_s, controls, gate,
            &mut checkpoint, &mut (),
        )
    }

    /// Both entries share the original mass solve. A thermal consumer prepares
    /// one immutable face image and observes its actual signed constituent
    /// transfers; it never reconstructs fluxes from rounded endpoint masses.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn diffuse_species_observed<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        controls: SpeciesDiffusionControls<N>,
        gate: &CancelGate,
        mut checkpoint: impl FnMut(usize) -> Result<(), ChtError>,
        observer: &mut impl DiffusionObserver<N>,
    ) -> Result<SpeciesDiffusionStep<N>, ChtError> {
        poll(gate)?;
        let transport = self.accepted_mass_transport(gate)?;
        let count = domain.cell_count();
        if N == 0
            || N > MAX_MATERIAL_SPECIES
            || self.projection.field.dims() != domain.dims()
            || count != self.mass_kg.len()
            || species_mass_kg.len() != count
            || diffusivity_m2_s.len() != count
            || controls.substeps == 0
            || controls.substeps as u64 > (1u64 << 53)
        {
            return Err(refuse(
                "diffusion component count, dimensions or substeps are invalid",
            ));
        }
        positive(duration_s)?;
        for value in controls
            .maximum_mass_defect_kg
            .into_iter()
            .chain([controls.maximum_partition_defect_kg])
        {
            nonnegative(value)?;
        }
        let half_dt = positive(0.5 * (duration_s / controls.substeps as f64))?;
        let mut initial_partition = Sum::default();
        let mut face_count = 0usize;
        for c in 0..count {
            poll_every(gate, c)?;
            nonnegative(diffusivity_m2_s[c])?;
            for &mass in &species_mass_kg[c] {
                nonnegative(mass)?;
            }
            if domain.is_fluid(c) != (transport.row_of[c] != usize::MAX) {
                return Err(refuse(
                    "diffusion fluid mask differs from the accepted carrier",
                ));
            }
            if !domain.is_fluid(c) {
                if species_mass_kg[c] != [0.0; N] || diffusivity_m2_s[c] != 0.0 {
                    return Err(refuse("solid cells cannot carry species or diffusivity"));
                }
                continue;
            }
            positive(domain.volume(c))?;
            let mut partition = Sum::default();
            partition.add(-self.mass_kg[c]);
            for &mass in &species_mass_kg[c] {
                ratio(mass, self.mass_kg[c])?;
                partition.add(mass);
            }
            initial_partition.add(partition.value.abs());
            for side in [1, 3, 5] {
                if let Some(next) = domain.neighbor(c, side)
                    && domain.is_fluid(next)
                    && diffusivity_m2_s[c] > 0.0
                    && diffusivity_m2_s[next] > 0.0
                {
                    face_count = face_count
                        .checked_add(1)
                        .ok_or_else(|| refuse("diffusion face extent overflows"))?;
                }
            }
        }
        nonnegative(initial_partition.value)?;
        if initial_partition.value > controls.maximum_partition_defect_kg {
            return Err(refuse("species do not partition the accepted carrier mass"));
        }
        let boundary_count = observer.boundary_faces(gate)?;
        let all_face_count = face_count.checked_add(boundary_count)
            .ok_or_else(|| refuse("diffusion boundary face extent overflows"))?;
        let pair_updates = all_face_count
            .checked_mul(2)
            .and_then(|n| n.checked_mul(controls.substeps))
            .ok_or_else(|| refuse("diffusion work extent overflows"))?;
        let component_updates = pair_updates
            .checked_mul(N)
            .ok_or_else(|| refuse("diffusion component work extent overflows"))?;
        if pair_updates > controls.maximum_pair_updates {
            return Err(refuse("diffusion pair-update budget exhausted"));
        }
        let fixed_bytes = core::mem::size_of::<SpeciesDiffusionStep<N>>()
            .checked_add(4 * core::mem::size_of::<[Sum; N]>())
            .and_then(|n| n.checked_add(core::mem::size_of::<[f64; N]>()))
            .ok_or_else(|| refuse("diffusion fixed storage extent overflows"))?;
        let bytes = count
            .checked_mul(core::mem::size_of::<Inventory<N>>() + core::mem::size_of::<[f64; N]>())
            .and_then(|n| n.checked_add(face_count.checked_mul(core::mem::size_of::<Face>())?))
            .and_then(|n| n.checked_add(fixed_bytes.max(2048)))
            .and_then(|n| n.checked_add(observer.storage_bytes(count, face_count)?))
            .ok_or_else(|| refuse("diffusion storage extent overflows"))?;
        if bytes > controls.maximum_storage_bytes {
            return Err(refuse("diffusion storage budget exhausted"));
        }
        observer.prepare_storage(count, face_count);
        observer.prepare_boundaries(half_dt, gate)?;
        let mut faces = Vec::with_capacity(face_count);
        let mut maximum_pair_relaxation = observer.maximum_boundary_relaxation();
        for c in 0..count {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) || diffusivity_m2_s[c] == 0.0 {
                continue;
            }
            for side in [1, 3, 5] {
                let Some(next) = domain.neighbor(c, side).filter(|&n| domain.is_fluid(n)) else {
                    continue;
                };
                if diffusivity_m2_s[next] == 0.0 {
                    continue;
                }
                let axis = side / 2;
                let conductance = self.diffusion_face_conductance(
                    domain,
                    [c, next],
                    axis,
                    [diffusivity_m2_s[c], diffusivity_m2_s[next]],
                )?;
                let small = self.mass_kg[c].min(self.mass_kg[next]);
                let large = self.mass_kg[c].max(self.mass_kg[next]);
                let capacity = positive(small / (1.0 + small / large))?;
                let mut exponent = (conductance / capacity) * half_dt;
                if !exponent.is_normal() {
                    // An infinite rate-time product means full pair equilibration;
                    // balancing the factors also retains finite scaled products.
                    exponent =
                        det::exp(det::ln(conductance) - det::ln(capacity) + det::ln(half_dt));
                }
                let relaxation = -det::expm1(-exponent);
                if !relaxation.is_finite() || !(0.0..=1.0).contains(&relaxation) {
                    return Err(refuse("diffusion relaxation is unrepresentable"));
                }
                maximum_pair_relaxation = maximum_pair_relaxation.max(relaxation);
                observer.prepare([c, next], axis)?;
                faces.push(Face {
                    cells: [c, next],
                    exchange_capacity_kg: positive(capacity * relaxation)?,
                });
            }
        }
        let mut inventories = Vec::with_capacity(count);
        for (c, &mass) in species_mass_kg.iter().enumerate() {
            poll_every(gate, c)?;
            inventories.push(Inventory {
                mass: mass.map(|value| Sum {
                    value,
                    compensation: 0.0,
                }),
                increment: [Sum::default(); N],
            });
        }
        // With no conducting faces, no repeated empty sweeps consume a huge
        // caller-specified substep count.
        if all_face_count > 0 {
            for substep in 0..controls.substeps {
                poll(gate)?;
                observer.boundary_half_sweep(&mut inventories, false, gate)?;
                for (index, face) in faces.iter().enumerate() {
                    poll_every(gate, index)?;
                    observer.transferred(index, exchange(face, &self.mass_kg, &mut inventories)?)?;
                }
                for (index, face) in faces.iter().enumerate().rev() {
                    poll_every(gate, index)?;
                    observer.transferred(index, exchange(face, &self.mass_kg, &mut inventories)?)?;
                }
                observer.boundary_half_sweep(&mut inventories, true, gate)?;
                checkpoint(substep + 1)?;
            }
        }
        let mut before = [Sum::default(); N];
        let mut after = [Sum::default(); N];
        let mut change = [Sum::default(); N];
        let mut defect = [Sum::default(); N];
        let mut endpoint_partition = Sum::default();
        let mut result = Vec::with_capacity(count);
        for c in 0..count {
            poll_every(gate, c)?;
            let endpoint = inventories[c].mass.map(|m| m.value);
            let mut partition = Sum::default();
            partition.add(-self.mass_kg[c]);
            for k in 0..N {
                nonnegative(endpoint[k])?;
                before[k].add(species_mass_kg[c][k]);
                after[k].add(endpoint[k]);
                let delta = endpoint[k] - species_mass_kg[c][k];
                change[k].add(delta);
                defect[k].add((delta - inventories[c].increment[k].value).abs());
                partition.add(endpoint[k]);
            }
            endpoint_partition.add(partition.value.abs());
            result.push(endpoint);
        }
        let boundary_mass = observer.boundary_mass_transfer();
        let absolute_mass_defect_kg =
            core::array::from_fn(|k| defect[k].value + (change[k].value - boundary_mass[k]).abs());
        for k in 0..N {
            for value in [before[k].value, after[k].value, absolute_mass_defect_kg[k]] {
                nonnegative(value)?;
            }
            if absolute_mass_defect_kg[k] > controls.maximum_mass_defect_kg[k] {
                return Err(refuse("diffusion exceeds its constituent-mass budget"));
            }
        }
        let partition_defect = initial_partition.value + endpoint_partition.value
            + observer.boundary_carrier_defect();
        nonnegative(partition_defect)?;
        if partition_defect > controls.maximum_partition_defect_kg {
            return Err(refuse(
                "diffusion does not preserve the accepted mass partition",
            ));
        }
        poll(gate)?;
        Ok(SpeciesDiffusionStep {
            species_mass_kg: result,
            report: SpeciesDiffusionReport {
                duration_s,
                substeps: controls.substeps,
                diffusive_faces: all_face_count,
                pair_updates,
                component_updates,
                maximum_pair_relaxation,
                mass_before_kg: before.map(|s| s.value),
                mass_after_kg: after.map(|s| s.value),
                mass_change_kg: change.map(|s| s.value),
                absolute_mass_defect_kg,
                initial_partition_defect_kg: initial_partition.value,
                endpoint_partition_defect_kg: endpoint_partition.value,
            },
        })
    }
}

fn exchange<const N: usize>(
    face: &Face,
    carrier: &[f64],
    inventory: &mut [Inventory<N>],
) -> Result<[f64; N], ChtError> {
    let [left, right] = face.cells;
    let mut transfers = [0.0; N];
    if N == 2 {
        let a_left = ratio(inventory[left].mass[0].value, carrier[left])?;
        let a_right = ratio(inventory[right].mass[0].value, carrier[right])?;
        let b_left = ratio(inventory[left].mass[1].value, carrier[left])?;
        let b_right = ratio(inventory[right].mass[1].value, carrier[right])?;
        // Retain the original binary path exactly. Complementary gradients
        // are identical, so the rarer constituent resolves a trace gradient
        // even when the majority's two fractions both round to one.
        let gap = if a_left + a_right <= b_left + b_right {
            a_left - a_right
        } else {
            b_right - b_left
        };
        let transferred = face.exchange_capacity_kg * gap;
        if !transferred.is_finite() || (gap != 0.0 && transferred == 0.0) {
            return Err(refuse("nonzero diffusion transfer is unrepresentable"));
        }
        transfers[0] = transferred;
        transfers[1] = -transferred;
    } else {
        // All species see one identical pair map. Keeping every gradient and
        // inventory independent preserves trace transport and permutation
        // symmetry without choosing a majority/remainder constituent.
        for (k, transferred) in transfers.iter_mut().enumerate() {
            let gap = ratio(inventory[left].mass[k].value, carrier[left])?
                - ratio(inventory[right].mass[k].value, carrier[right])?;
            *transferred = face.exchange_capacity_kg * gap;
            if !transferred.is_finite() || (gap != 0.0 && *transferred == 0.0) {
                return Err(refuse("nonzero diffusion transfer is unrepresentable"));
            }
        }
    }
    for (cell, sign) in [(left, -1.0), (right, 1.0)] {
        for (k, &transfer) in transfers.iter().enumerate() {
            let delta = sign * transfer;
            inventory[cell].mass[k].add(delta);
            inventory[cell].increment[k].add(delta);
            nonnegative(inventory[cell].mass[k].value)?;
            if !inventory[cell].mass[k].compensation.is_finite()
                || !inventory[cell].increment[k].value.is_finite()
                || !inventory[cell].increment[k].compensation.is_finite()
            {
                return Err(refuse("diffusion inventory arithmetic overflowed"));
            }
        }
    }
    Ok(transfers)
}

#[cfg(test)]
mod tests;
