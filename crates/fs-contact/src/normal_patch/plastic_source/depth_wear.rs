//! Depth-dependent indentation hardness exposed by actual material removal.
//!
//! Each body retains its own original H(depth) test and Archard denominator.
//! Actual cumulative BODY recession supplies exposed material depth; relative
//! plastic indentation is never split between the bodies. The local pressure
//! cap is min(H_A, H_B) on the original native contact grid.
//!
//! `advance_wear` holds beginning pressure, k and local H during removal. The
//! integrated variants instead solve integral H(z) dz = k p ds through the
//! original depth profile, still at beginning pressure and fixed k. Recession
//! and cap exposure then use the original native held-approach, maintained-force
//! or finite loading-spring boundary. Both work receipts are returned; signed
//! parameter work is not heat. Only the complete successor is publishable.
//! Refine physical slip increments when pressure changes appreciably: exact
//! depth integration alone does not integrate pressure/cap feedback. No layered
//! elastic kernel, coating failure, debris transport, hardness-to-yield
//! conversion, or annealing is inferred. Elasticity remains homogeneous; explicit
//! temperature continuations refresh original sources between wear intervals.
//! A profile prepared from complete native J2/cyclic history also retains its
//! admitted preconditioning alpha throughout the declared material column.
//! Contact removal changes only exposed depth; no local indentation-to-alpha
//! conversion or transfer of bulk residual stress/plastic work is performed.

use std::sync::Arc;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::PropertyKey;
use fs_material::hardness::depth::{
    DepthHardnessError, EXPOSED_MATERIAL_DEPTH_AXIS, ResolvedDepthHardness,
};
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_material::state_point::{
    POISSON_RATIO_PROPERTY, ResolvedInterfaceStatePoint, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, QuantitySpec};
use fs_tribo::{ArchardLaw, LocalHardnessArchardSample, TriboError, WearAdvance, WearState};

use super::{
    BoundNormalInterfaceState, FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticError,
    FiniteGapPlasticLinearSolver, FiniteGapPlasticModel, FiniteGapPlasticReceipt,
    FiniteGapPlasticState, PlasticContactControls, PlasticContactEnvelope,
    SourcePlasticContactError, admit_native_candidate, validate_envelope, validate_fixed_normal,
    within,
};
use crate::normal_patch::{
    FiniteGapContactError, FiniteGapPlasticContinuation, FiniteGapPlasticContinuationReceipt,
    FiniteGapPlasticDriveError, FiniteGapPlasticLoadError,
};

mod control;
pub mod integration;
/// Prescribed independent body/interface temperatures with retained wear.
pub mod temperature;
pub use control::{
    DepthHardnessNormalDriveAdvance, DepthHardnessNormalLoadAdvance, DepthHardnessWearDriveAdvance,
    DepthHardnessWearLoadAdvance,
};

#[derive(Clone, Debug, PartialEq)]
struct Sources {
    normal: BoundNormalInterfaceState,
    coefficients: ResolvedInterfaceStatePoint,
    coefficient_keys: [PropertyKey; 2],
    k: [f64; 2],
    hardness: [ResolvedDepthHardness; 2],
    envelope: PlasticContactEnvelope,
    maximum_recession_m: [f64; 2],
    identity: ContentHash,
}

/// Complete accepted original-source contact, per-body surface and wear state.
/// All fields are private; ordinary steps return a complete immutable candidate.
/// Prepared source lines and native compliance kernels are shared by successors.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceDepthHardnessWear {
    sources: Arc<Sources>,
    model: FiniteGapPlasticModel,
    contact: FiniteGapPlasticState,
    volumes: [WearState; 2],
    recession_m: [Vec<f64>; 2],
}

#[derive(Clone, Debug, PartialEq)]
/// Native normal advance and the complete successor at unchanged exposed depths.
pub struct DepthHardnessNormalAdvance {
    /// Replace the accepted owner only after accepting this complete candidate.
    pub candidate: SourceDepthHardnessWear,
    /// Original native contact and plastic-work response.
    pub receipt: FiniteGapPlasticReceipt,
}

/// Actual removal and exposed-depth cap change, with both original work accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthHardnessWearAdvance {
    /// Complete successor contact, hardness field, profiles and inventories.
    pub candidate: SourceDepthHardnessWear,
    /// One original Archard inventory receipt for each ordered body.
    pub wear: [WearAdvance; 2],
    /// Representable actual material removal at each original body/grid cell [m].
    pub recession_increment_m: [Vec<f64>; 2],
    /// Recession at the interval's original caps and requested normal boundary.
    pub recession: FiniteGapPlasticContinuationReceipt,
    /// Subsequent native cap refresh at the newly exposed material depths.
    pub hardness: FiniteGapPlasticContinuationReceipt,
}

#[derive(Clone, Debug, PartialEq)]
/// Refusal leaves the complete accepted owner unchanged.
pub enum DepthHardnessWearError {
    /// Original depth support, typed source or hardness value was refused.
    Hardness(DepthHardnessError),
    /// Native contact, original applicability, geometry or cancellation refusal.
    Contact(SourcePlasticContactError),
    /// Original Archard inventory admission failed.
    Tribology(TriboError),
    /// Prescribed force is unattainable, ambiguous, mismatched or unresolved.
    Load(FiniteGapPlasticLoadError),
    /// Finite loading-spring equilibrium or predecessor balance was refused.
    Drive(FiniteGapPlasticDriveError),
    /// A source ownership, physical or arithmetic requirement was not met.
    Invalid(&'static str),
}
impl From<DepthHardnessError> for DepthHardnessWearError {
    fn from(error: DepthHardnessError) -> Self {
        Self::Hardness(error)
    }
}
impl From<SourcePlasticContactError> for DepthHardnessWearError {
    fn from(error: SourcePlasticContactError) -> Self {
        Self::Contact(error)
    }
}
impl From<FiniteGapPlasticError> for DepthHardnessWearError {
    fn from(error: FiniteGapPlasticError) -> Self {
        Self::Contact(error.into())
    }
}
impl From<TriboError> for DepthHardnessWearError {
    fn from(error: TriboError) -> Self {
        Self::Tribology(error)
    }
}
impl core::fmt::Display for DepthHardnessWearError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "exposed-depth hardness/contact/wear: {self:?}")
    }
}
impl std::error::Error for DepthHardnessWearError {}

struct PreparedWear {
    recession_m: [Vec<f64>; 2],
    increment_m: [Vec<f64>; 2],
    relative_increment_m: Vec<f64>,
    next_caps_pa: Vec<f64>,
    volumes: [WearState; 2],
    wear: [WearAdvance; 2],
}

#[derive(Clone, Copy)]
enum RemovalRule {
    BeginningHardness,
    IntegratedProfile,
}

impl SourceDepthHardnessWear {
    /// Bind original dry interface, homogeneous elasticity, separate measured
    /// H(depth) profiles and fixed k_A/k_B. The pristine surface starts at depth
    /// zero in both bodies; normal loading does not change either exposed depth.
    /// Profiles made by `ResolvedDepthHardness::from_j2_accepted` or
    /// `from_cyclic_accepted` retain their native preconditioning and exact
    /// selected elastic claims. That accepted alpha must represent each complete
    /// declared material column; this reduced contact law does not evolve it.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        normal: BoundNormalInterfaceState,
        coefficients: ResolvedInterfaceStatePoint,
        coefficient_keys: [PropertyKey; 2],
        hardness: [ResolvedDepthHardness; 2],
        grid: FiniteGapGrid,
        undeformed_gap_m: Vec<f64>,
        controls: PlasticContactControls,
        envelope: PlasticContactEnvelope,
        maximum_recession_m: [f64; 2],
        linear_solver: FiniteGapPlasticLinearSolver,
        cx: &Cx<'_>,
    ) -> Result<Self, DepthHardnessWearError> {
        checkpoint(cx)?;
        let (k, selected) = admit_depth_sources(
            &normal,
            &coefficients,
            &coefficient_keys,
            &hardness,
            envelope,
            maximum_recession_m,
        )?;
        let count = undeformed_gap_m.len();
        let cap = hardness[0]
            .hardness_pa_at(0.0)?
            .min(hardness[1].hardness_pa_at(0.0)?);
        let model = FiniteGapPlasticModel::new_with_linear_solver(
            FiniteGapPlasticConfig {
                grid,
                undeformed_gap_m,
                reduced_modulus_pa: normal.elastic().reduced_modulus_pa(),
                hardness_pa: vec![cap; count],
                maximum_active_set_iterations: controls.maximum_active_set_iterations,
                complementarity_tolerance_m: controls.complementarity_tolerance_m,
                boundary_clearance_cells: controls.boundary_clearance_cells,
            },
            linear_solver,
            cx,
        )?;
        let identity = depth_source_identity(
            &normal,
            &coefficients,
            &hardness,
            &selected,
            &model,
            envelope,
            maximum_recession_m,
        );
        let contact = model.initial_state();
        let source = Self {
            sources: Arc::new(Sources {
                normal,
                coefficients,
                coefficient_keys,
                k,
                hardness,
                envelope,
                maximum_recession_m,
                identity,
            }),
            model,
            contact,
            volumes: [WearState::default(); 2],
            recession_m: [vec![0.0; count], vec![0.0; count]],
        };
        checkpoint(cx)?;
        Ok(source)
    }

    /// Identity of original sources, reference geometry and physical bounds.
    #[must_use]
    pub fn source_identity(&self) -> ContentHash {
        self.sources.identity
    }
    /// Original fixed dry interface and ordered homogeneous elastic sources.
    #[must_use]
    pub fn normal_source(&self) -> &BoundNormalInterfaceState {
        &self.sources.normal
    }
    /// Complete original interface bundle selecting both fixed wear coefficients.
    #[must_use]
    pub fn coefficient_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.sources.coefficients
    }
    /// Exact coefficient keys in the original body order.
    #[must_use]
    pub fn coefficient_keys(&self) -> &[PropertyKey; 2] {
        &self.sources.coefficient_keys
    }
    /// Independently selected original indentation tests in body order.
    #[must_use]
    pub fn hardness_profiles(&self) -> &[ResolvedDepthHardness; 2] {
        &self.sources.hardness
    }
    /// Current native model, including spatial caps at accepted exposed depths.
    #[must_use]
    pub const fn model(&self) -> &FiniteGapPlasticModel {
        &self.model
    }
    /// Complete original native mechanical history at accepted exposed depths.
    #[must_use]
    pub const fn contact_state(&self) -> &FiniteGapPlasticState {
        &self.contact
    }
    /// Original per-body exposed material depth, advanced only by actual removal.
    #[must_use]
    pub fn recession_m(&self) -> [&[f64]; 2] {
        [&self.recession_m[0], &self.recession_m[1]]
    }
    /// Original cumulative removed material inventory for each body [m³].
    #[must_use]
    pub fn volumes_m3(&self) -> [f64; 2] {
        self.volumes.map(|state| state.volume_m3())
    }

    /// Native normal loading at the currently exposed hardness field; no wear.
    pub fn advance_normal(
        &self,
        approach_m: f64,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessNormalAdvance, DepthHardnessWearError> {
        checkpoint(cx)?;
        within(
            "approach [m]",
            approach_m,
            self.sources.envelope.maximum_approach_m,
        )?;
        let native = self.model.advance(&self.contact, approach_m, cx)?;
        self.admit(&self.model, &native.candidate, &self.recession_m, cx)?;
        let mut candidate = self.clone();
        candidate.contact = native.candidate;
        Ok(DepthHardnessNormalAdvance {
            candidate,
            receipt: native.receipt,
        })
    }

    /// Remove actual unsigned local sliding material, expose new hardness and
    /// re-equilibrate at the same approach. All candidates start from this
    /// accepted history; a source/geometry/arithmetic/cancellation refusal leaves
    /// both bodies and every original mechanical history untouched.
    pub fn advance_wear(
        &self,
        actual_unsigned_slip_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessWearAdvance, DepthHardnessWearError> {
        self.advance_wear_at_approach(actual_unsigned_slip_m, RemovalRule::BeginningHardness, cx)
    }

    /// Integrate each original H(depth) profile through the removed material,
    /// then equilibrate recession and exposed hardness at the held approach.
    ///
    /// Each cell solves `integral H(z) dz = k * p_begin * |ds|`, including all
    /// crossed original affine segments. The interval mean of that same profile
    /// supplies the original Archard inventory kernel. Beginning pressure and k
    /// remain held, so slip refinement is still required for pressure feedback.
    /// Source support, actual depth/volume resolution and all mechanical bounds
    /// must admit the complete result before any successor can be published.
    pub fn advance_wear_integrated(
        &self,
        actual_unsigned_slip_m: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessWearAdvance, DepthHardnessWearError> {
        self.advance_wear_at_approach(actual_unsigned_slip_m, RemovalRule::IntegratedProfile, cx)
    }

    fn advance_wear_at_approach(
        &self,
        actual_unsigned_slip_m: &[f64],
        rule: RemovalRule,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessWearAdvance, DepthHardnessWearError> {
        let mut prepared = self.prepare_wear(actual_unsigned_slip_m, rule, cx)?;
        let removed = self
            .model
            .recede(&self.contact, &prepared.relative_increment_m, cx)?;
        let model = self.model_after_removal(&mut prepared, &removed, cx)?;
        let refreshed = model.refresh_from(&self.model, &removed.candidate, cx)?;
        self.finish_wear(prepared, model, removed, refreshed, cx)
    }

    #[allow(clippy::too_many_lines)]
    fn prepare_wear(
        &self,
        actual_unsigned_slip_m: &[f64],
        rule: RemovalRule,
        cx: &Cx<'_>,
    ) -> Result<PreparedWear, DepthHardnessWearError> {
        checkpoint(cx)?;
        let pressure = self.contact.pressure_pa();
        if actual_unsigned_slip_m.len() != pressure.len()
            || actual_unsigned_slip_m
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "finite nonnegative actual local sliding on the original grid required",
            ));
        }
        let grid = self.model.config().grid;
        let area = grid.cell_width_m * grid.cell_depth_m;
        let mut recession = self.recession_m.clone();
        let mut increment = [vec![0.0; pressure.len()], vec![0.0; pressure.len()]];
        let mut relative = vec![0.0; pressure.len()];
        let mut hardness = vec![[0.0; 2]; pressure.len()];
        let mut next_caps = vec![0.0; pressure.len()];
        let mut measures = vec![0.0; pressure.len()];
        let mut geometric = [0.0; 2];
        let mut corrections = [0.0; 2];
        let mut roundoff = [0.0; 2];
        for (i, (&p, &slip)) in pressure.iter().zip(actual_unsigned_slip_m).enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            measures[i] = (p * area) * slip;
            if !measures[i].is_finite() {
                return Err(invalid("finite pressure-weighted local sliding required"));
            }
            let mut next_h = [0.0; 2];
            for side in 0..2 {
                let old = self.recession_m[side][i];
                let profile = &self.sources.hardness[side];
                let (total, depth) = match rule {
                    RemovalRule::BeginningHardness => {
                        hardness[i][side] = profile.hardness_pa_at(old)?;
                        let depth = (self.sources.k[side] * (p / hardness[i][side])) * slip;
                        (old + depth, depth)
                    }
                    RemovalRule::IntegratedProfile => {
                        let active = self.sources.k[side] > 0.0 && p > 0.0 && slip > 0.0;
                        let integral = if active {
                            (self.sources.k[side] * p) * slip
                        } else {
                            0.0
                        };
                        if !integral.is_finite() || (active && integral == 0.0) {
                            return Err(invalid(
                                "depth-hardness integral lies outside arithmetic resolution",
                            ));
                        }
                        let advance = profile.advance_by_integral(old, integral)?;
                        hardness[i][side] = advance.path_mean_hardness_pa;
                        (advance.depth_m, advance.increment_m)
                    }
                };
                let retained = total - old;
                if !depth.is_finite() || !total.is_finite() || (depth > 0.0 && retained == 0.0) {
                    return Err(invalid(
                        "Archard recession must be finite and representable in its original body profile",
                    ));
                }
                within(
                    "body cumulative recession [m]",
                    total,
                    self.sources.maximum_recession_m[side],
                )?;
                next_h[side] = self.sources.hardness[side].hardness_pa_at(total)?;
                let volume = retained * area;
                if !volume.is_finite() || (retained > 0.0 && volume == 0.0) {
                    return Err(invalid(
                        "wear increment lies outside volume arithmetic resolution",
                    ));
                }
                add(&mut geometric[side], &mut corrections[side], volume);
                roundoff[side] += 4.0 * f64::EPSILON * ((old * area).abs() + (total * area).abs());
                increment[side][i] = retained;
                recession[side][i] = total;
            }
            next_caps[i] = next_h[0].min(next_h[1]);
            relative[i] = increment[0][i] + increment[1][i];
            if !relative[i].is_finite() {
                return Err(invalid("finite relative removal required"));
            }
        }
        let mut volumes = self.volumes;
        let interface = self.sources.normal.elastic().interface().interface();
        let mut advance_volume = |side| -> Result<WearAdvance, DepthHardnessWearError> {
            let mut interrupted = false;
            let samples = (0..pressure.len()).scan((), |_, i| {
                if i % 64 == 0 && cx.checkpoint().is_err() {
                    interrupted = true;
                    return None;
                }
                Some(LocalHardnessArchardSample {
                    wear_coefficient: self.sources.k[side],
                    pressure_weighted_sliding_j: measures[i],
                    hardness_pa: hardness[i][side],
                })
            });
            let result = ArchardLaw::advance_local_hardness_pressure_weighted_sliding(
                interface,
                &mut volumes[side],
                samples,
            );
            if interrupted {
                return Err(cancelled());
            }
            Ok(result?)
        };
        let wear = [advance_volume(0)?, advance_volume(1)?];
        for side in 0..2 {
            let native = wear[side].increment_m3;
            let retained = volumes[side].volume_m3() - self.volumes[side].volume_m3();
            let arithmetic =
                roundoff[side] + 16.0 * f64::EPSILON * native.abs().max(geometric[side].abs());
            let inventory = arithmetic
                + 4.0
                    * f64::EPSILON
                    * (self.volumes[side].volume_m3().abs() + volumes[side].volume_m3().abs());
            if !arithmetic.is_finite()
                || !inventory.is_finite()
                || !geometric[side].is_finite()
                || (native > 0.0 && retained == 0.0)
                || (native - geometric[side]).abs() > arithmetic
                || (retained - geometric[side]).abs() > inventory
            {
                return Err(invalid(
                    "original Archard inventory and body recession disagree beyond roundoff",
                ));
            }
        }
        Ok(PreparedWear {
            recession_m: recession,
            increment_m: increment,
            relative_increment_m: relative,
            next_caps_pa: next_caps,
            volumes,
            wear,
        })
    }

    fn model_after_removal(
        &self,
        prepared: &mut PreparedWear,
        removed: &FiniteGapPlasticContinuation,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticModel, DepthHardnessWearError> {
        self.admit(&self.model, &removed.candidate, &prepared.recession_m, cx)?;
        Ok(self.model.with_material_properties(
            self.model.config().reduced_modulus_pa,
            std::mem::take(&mut prepared.next_caps_pa),
            cx,
        )?)
    }

    fn finish_wear(
        &self,
        prepared: PreparedWear,
        model: FiniteGapPlasticModel,
        removed: FiniteGapPlasticContinuation,
        refreshed: FiniteGapPlasticContinuation,
        cx: &Cx<'_>,
    ) -> Result<DepthHardnessWearAdvance, DepthHardnessWearError> {
        self.admit(&model, &refreshed.candidate, &prepared.recession_m, cx)?;
        checkpoint(cx)?;
        Ok(DepthHardnessWearAdvance {
            candidate: Self {
                sources: Arc::clone(&self.sources),
                model,
                contact: refreshed.candidate,
                volumes: prepared.volumes,
                recession_m: prepared.recession_m,
            },
            wear: prepared.wear,
            recession_increment_m: prepared.increment_m,
            recession: removed.receipt,
            hardness: refreshed.receipt,
        })
    }

    fn admit(
        &self,
        model: &FiniteGapPlasticModel,
        state: &FiniteGapPlasticState,
        recession: &[Vec<f64>; 2],
        cx: &Cx<'_>,
    ) -> Result<(), DepthHardnessWearError> {
        let radius = admit_native_candidate(model, self.sources.envelope, state, cx)?;
        for side in 0..2 {
            let maximum = recession[side].iter().copied().fold(0.0_f64, f64::max);
            let remaining = self.sources.envelope.half_space_depth_m[side] - maximum;
            if remaining <= 0.0 {
                return Err(invalid(
                    "wear exhausted the original homogeneous elastic depth",
                ));
            }
            within(
                "patch radius / remaining homogeneous elastic depth",
                radius / remaining,
                self.sources.envelope.maximum_patch_to_depth[side],
            )?;
        }
        checkpoint(cx)
    }
}

// One source admission for initial binding and prescribed-temperature refresh.
fn admit_depth_sources(
    normal: &BoundNormalInterfaceState,
    coefficients: &ResolvedInterfaceStatePoint,
    coefficient_keys: &[PropertyKey; 2],
    hardness: &[ResolvedDepthHardness; 2],
    envelope: PlasticContactEnvelope,
    maximum_recession_m: [f64; 2],
) -> Result<([f64; 2], [usize; 2]), DepthHardnessWearError> {
    validate_fixed_normal(normal)?;
    validate_envelope(envelope)?;
    let original = normal.interface_state();
    if coefficients.card_identity() != original.card_identity()
        || coefficients.surface_materials() != original.surface_materials()
        || coefficients.query_point() != original.query_point()
        || coefficients.axis_quantities() != original.axis_quantities()
    {
        return Err(invalid(
            "wear coefficients must retain the original ordered interface and fixed query",
        ));
    }
    let surfaces = [
        normal.elastic().surface_a_state(),
        normal.elastic().surface_b_state(),
    ];
    let mut k = [0.0; 2];
    let mut selected = [0; 2];
    for side in 0..2 {
        let surface = surfaces[side];
        let bulk = surface.resolved();
        let history = hardness[side].accepted_history();
        let owned = |name: &str| {
            name == EXPOSED_MATERIAL_DEPTH_AXIS || (history.is_some() && name == ALPHA)
        };
        let h = hardness[side].reference_source().resolved();
        if h.card_identity() != bulk.card_identity()
            || h.material() != bulk.material()
            || !h
                .query_point()
                .iter()
                .filter(|(name, _)| !owned(name))
                .eq(bulk.query_point().iter())
            || !h
                .axis_quantities()
                .iter()
                .filter(|(name, _)| !owned(name))
                .eq(bulk.axis_quantities().iter())
            || bulk
                .properties()
                .iter()
                .any(|property| !property.is_static_at_fixed_external_conditions())
        {
            return Err(invalid(
                "depth hardness must retain each original ordered material and fixed environment",
            ));
        }
        if let Some(history) = history {
            let mechanical = history.mechanical_source();
            if mechanical.card_identity() != bulk.card_identity()
                || mechanical.material() != bulk.material()
                || mechanical.query_point() != bulk.query_point()
                || mechanical.axis_quantities() != bulk.axis_quantities()
                || h.state_coordinate(ALPHA) != Some(history.alpha())
                || h.axis_quantities().get(ALPHA) != Some(&QuantitySpec::dimensional(Dims::NONE))
                || surface.elastic_state().derivation().is_some()
            {
                return Err(invalid(
                    "preconditioned depth hardness requires its admitted native alpha and original canonical elastic source",
                ));
            }
            for (name, value) in [
                (YOUNG_MODULUS_PROPERTY, surface.young_modulus_pa()),
                (POISSON_RATIO_PROPERTY, surface.poisson_ratio()),
            ] {
                let selected = mechanical
                    .property(name)
                    .ok_or_else(|| invalid("missing original native elastic source role"))?;
                if Some(selected) != bulk.property(name)
                    || !selected.is_static_at_fixed_external_conditions()
                    || selected.value_si().to_bits() != value.to_bits()
                {
                    return Err(invalid(
                        "preconditioned profile and half-space must retain identical original selected elastic claims",
                    ));
                }
            }
        }
        let coefficient = coefficients
            .property_by_key(&coefficient_keys[side])
            .ok_or_else(|| invalid("missing exact original wear coefficient"))?;
        if coefficient_keys[side].quantity() != QuantitySpec::dimensional(Dims::NONE)
            || !coefficient.value_si().is_finite()
            || coefficient.value_si() < 0.0
            || !coefficient.is_static_at_fixed_external_conditions()
            || !maximum_recession_m[side].is_finite()
            || maximum_recession_m[side] < 0.0
            || maximum_recession_m[side] > hardness[side].maximum_depth_m()
        {
            return Err(invalid(
                "fixed nonnegative k and bounded original per-body depth support required",
            ));
        }
        k[side] = coefficient.value_si();
        selected[side] = coefficients
            .properties()
            .iter()
            .position(|value| value == coefficient)
            .ok_or_else(|| invalid("wear coefficient must be retained by its original bundle"))?;
    }
    Ok((k, selected))
}

fn depth_source_identity(
    normal: &BoundNormalInterfaceState,
    coefficients: &ResolvedInterfaceStatePoint,
    hardness: &[ResolvedDepthHardness; 2],
    selected: &[usize; 2],
    model: &FiniteGapPlasticModel,
    envelope: PlasticContactEnvelope,
    maximum_recession_m: [f64; 2],
) -> ContentHash {
    let mut hash = DomainHasher::new("org.frankensim.fs-contact.source-depth-hardness-wear.v1");
    hash.update(normal.identity().as_bytes());
    hash.update(coefficients.identity().as_bytes());
    hash.update(model.identity().as_bytes());
    for side in 0..2 {
        hash.update(&(selected[side] as u64).to_le_bytes());
        hash.update(hardness[side].identity().as_bytes());
    }
    for value in [
        envelope.maximum_approach_m,
        envelope.maximum_plastic_depth_m,
        envelope.half_space_depth_m[0],
        envelope.half_space_depth_m[1],
        envelope.maximum_patch_to_depth[0],
        envelope.maximum_patch_to_depth[1],
        maximum_recession_m[0],
        maximum_recession_m[1],
    ] {
        hash.update(&value.to_bits().to_le_bytes());
    }
    hash.finalize()
}

fn invalid(reason: &'static str) -> DepthHardnessWearError {
    DepthHardnessWearError::Invalid(reason)
}
fn cancelled() -> DepthHardnessWearError {
    FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled).into()
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), DepthHardnessWearError> {
    cx.checkpoint().map_err(|_| cancelled())
}
fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let corrected = value - *correction;
    let next = *sum + corrected;
    *correction = (next - *sum) - corrected;
    *sum = next;
}
