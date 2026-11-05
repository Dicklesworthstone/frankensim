//! Original material hardness driving irreversible finite-gap indentation.
//!
//! This is a quasistatic surface pressure-saturation approximation. The ordered
//! source pair supplies the elastic modulus, and the softer surface's original
//! test-specific pressure-valued indentation hardness supplies the pressure cap.
//! Permanent depth describes the relative surface gap; it is not arbitrarily divided
//! between the bodies. No yield-stress or empirical hardness-scale conversion
//! is inferred.
//!
//! Material state and external conditions are fixed for each owner. Explicit
//! thermal and wear continuations retain the original irreversible history.

/// Original temperature-dependent material properties with retained history.
pub mod temperature;
/// Local Archard recession and its effect on the original contact pressure.
pub mod wear;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::hardness::IndentationHardnessStatePoint;

use crate::interface_binding::BoundNormalInterfaceState;

use super::{
    FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticError, FiniteGapPlasticModel,
    FiniteGapPlasticReceipt, FiniteGapPlasticState,
};

/// Numerical controls for the bounded half-space contact solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlasticContactControls {
    /// Maximum active-set iterations in one candidate solve.
    pub maximum_active_set_iterations: usize,
    /// Maximum displacement complementarity residual, in metres.
    pub complementarity_tolerance_m: f64,
    /// Required unloaded rows and columns at the grid boundary.
    pub boundary_clearance_cells: usize,
}

/// Explicit applicability bounds for this surface approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlasticContactEnvelope {
    /// Maximum signed approach relative to the original gap datum, in metres.
    pub maximum_approach_m: f64,
    /// Maximum accumulated permanent relative indentation, in metres.
    pub maximum_plastic_depth_m: f64,
    /// Homogeneous half-space depth for each original surface, in metres.
    pub half_space_depth_m: [f64; 2],
    /// Maximum contact footprint radius divided by each homogeneous depth.
    pub maximum_patch_to_depth: [f64; 2],
}

/// Source, physical-domain, or native numerical refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SourcePlasticContactError {
    /// A source or parameter does not satisfy the declared model.
    Invalid(&'static str),
    /// A candidate exceeds a declared physical applicability bound.
    OutsideApplicability {
        /// Physical quantity being checked.
        quantity: &'static str,
        /// Candidate value.
        value: f64,
        /// Declared upper bound.
        limit: f64,
    },
    /// The original finite-gap solver refused the candidate.
    Contact(FiniteGapPlasticError),
}

impl From<FiniteGapPlasticError> for SourcePlasticContactError {
    fn from(value: FiniteGapPlasticError) -> Self {
        Self::Contact(value)
    }
}

impl core::fmt::Display for SourcePlasticContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source hardness-limited contact: {self:?}")
    }
}

impl std::error::Error for SourcePlasticContactError {}

/// Fixed-state original material sources and their finite-gap indentation law.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticContact {
    normal: BoundNormalInterfaceState,
    hardness: [IndentationHardnessStatePoint; 2],
    model: FiniteGapPlasticModel,
    envelope: PlasticContactEnvelope,
    identity: ContentHash,
}

/// Immutable history sealed to the complete ordered original source pair.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticContactState {
    owner: ContentHash,
    contact: FiniteGapPlasticState,
}

impl SourcePlasticContactState {
    /// Read the original native contact history without changing its owner.
    pub const fn contact_state(&self) -> &FiniteGapPlasticState {
        &self.contact
    }
}

/// A complete candidate; the caller adopts it only after its own coupled audit.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePlasticContactAdvance {
    /// Identity of the complete ordered original sources and physical envelope.
    pub source_identity: ContentHash,
    /// Proposed successor history. The predecessor remains unchanged.
    pub candidate: SourcePlasticContactState,
    /// Original pressure, elastic energy, permanent depth, and plastic work.
    pub receipt: FiniteGapPlasticReceipt,
}

impl SourcePlasticContact {
    /// Bind the exact original surface cards, conditions, and indentation tests.
    pub fn new(
        normal: BoundNormalInterfaceState,
        hardness: [IndentationHardnessStatePoint; 2],
        grid: FiniteGapGrid,
        undeformed_gap_m: Vec<f64>,
        controls: PlasticContactControls,
        envelope: PlasticContactEnvelope,
        cx: &Cx<'_>,
    ) -> Result<Self, SourcePlasticContactError> {
        if normal.adhesion_energy_j_per_m2() != 0.0
            || normal
                .elastic()
                .interface()
                .card()
                .context()
                .third_body
                .is_some()
        {
            return Err(SourcePlasticContactError::Invalid(
                "pressure saturation requires a nonadhesive dry interface without a third body",
            ));
        }
        if normal
            .interface_state()
            .properties()
            .iter()
            .any(|property| !property.is_static_at_fixed_external_conditions())
        {
            return Err(SourcePlasticContactError::Invalid(
                "pressure saturation requires static interface data at fixed external conditions",
            ));
        }
        let original = [
            normal.elastic().surface_a_state().resolved(),
            normal.elastic().surface_b_state().resolved(),
        ];
        for (source, surface) in hardness.iter().zip(original) {
            let resolved = source.resolved();
            if resolved.card_identity() != surface.card_identity()
                || resolved.material() != surface.material()
                || resolved.query_point() != surface.query_point()
                || resolved.axis_quantities() != surface.axis_quantities()
                || !source.property().is_static_at_fixed_external_conditions()
                || surface
                    .properties()
                    .iter()
                    .any(|property| !property.is_static_at_fixed_external_conditions())
            {
                return Err(SourcePlasticContactError::Invalid(
                    "hardness must retain the original ordered surface, state, and indentation test",
                ));
            }
        }
        if !envelope.maximum_approach_m.is_finite()
            || envelope.maximum_approach_m <= 0.0
            || !envelope.maximum_plastic_depth_m.is_finite()
            || envelope.maximum_plastic_depth_m < 0.0
            || envelope
                .half_space_depth_m
                .iter()
                .chain(envelope.maximum_patch_to_depth.iter())
                .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(SourcePlasticContactError::Invalid(
                "finite positive approach, depths, and ratios; nonnegative permanent depth bound",
            ));
        }
        let cap = hardness[0].hardness_pa().min(hardness[1].hardness_pa());
        let hardness_pa = vec![cap; undeformed_gap_m.len()];
        let model = FiniteGapPlasticModel::new(
            FiniteGapPlasticConfig {
                grid,
                undeformed_gap_m,
                reduced_modulus_pa: normal.elastic().reduced_modulus_pa(),
                hardness_pa,
                maximum_active_set_iterations: controls.maximum_active_set_iterations,
                complementarity_tolerance_m: controls.complementarity_tolerance_m,
                boundary_clearance_cells: controls.boundary_clearance_cells,
            },
            cx,
        )?;
        let mut hasher = DomainHasher::new("org.frankensim.fs-contact.source-plastic-gap.v1");
        hasher.update(normal.identity().as_bytes());
        for source in &hardness {
            hasher.update(source.identity().as_bytes());
        }
        hasher.update(model.identity().as_bytes());
        for value in [
            envelope.maximum_approach_m,
            envelope.maximum_plastic_depth_m,
            envelope.half_space_depth_m[0],
            envelope.half_space_depth_m[1],
            envelope.maximum_patch_to_depth[0],
            envelope.maximum_patch_to_depth[1],
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
        let identity = hasher.finalize();
        Ok(Self {
            normal,
            hardness,
            model,
            envelope,
            identity,
        })
    }

    /// Complete ordered normal material/interface source.
    pub const fn normal_source(&self) -> &BoundNormalInterfaceState {
        &self.normal
    }

    /// Complete indentation-hardness sources in original surface order.
    pub const fn hardness_sources(&self) -> &[IndentationHardnessStatePoint; 2] {
        &self.hardness
    }

    /// The original native, precomputed compliance owner.
    pub const fn model(&self) -> &FiniteGapPlasticModel {
        &self.model
    }

    /// Complete source identity, including tests, conditions, and envelope.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Undeformed, unloaded history bound to these original sources.
    pub fn initial_state(&self) -> SourcePlasticContactState {
        SourcePlasticContactState {
            owner: self.identity,
            contact: self.model.initial_state(),
        }
    }

    /// Propose a quasistatic step without mutating accepted physical history.
    pub fn advance(
        &self,
        before: &SourcePlasticContactState,
        approach_m: f64,
        cx: &Cx<'_>,
    ) -> Result<SourcePlasticContactAdvance, SourcePlasticContactError> {
        if before.owner != self.identity {
            return Err(SourcePlasticContactError::Invalid(
                "history belongs to different original material sources",
            ));
        }
        within("approach [m]", approach_m, self.envelope.maximum_approach_m)?;
        let advance = self.model.advance(&before.contact, approach_m, cx)?;
        self.admit_candidate(&advance.candidate, cx)?;
        Ok(SourcePlasticContactAdvance {
            source_identity: self.identity,
            candidate: SourcePlasticContactState {
                owner: self.identity,
                contact: advance.candidate,
            },
            receipt: advance.receipt,
        })
    }

    // Mechanical, thermal and wear candidates face the same physical domain.
    // A continuation cannot evade a late depth or footprint refusal.
    fn admit_candidate(
        &self,
        candidate: &FiniteGapPlasticState,
        cx: &Cx<'_>,
    ) -> Result<(), SourcePlasticContactError> {
        if candidate.owner_identity() != self.model.identity() {
            return Err(FiniteGapPlasticError::ForeignState.into());
        }
        within(
            "approach [m]",
            candidate.approach_m(),
            self.envelope.maximum_approach_m,
        )?;
        let maximum_depth = candidate
            .plastic_depth_m()
            .iter()
            .copied()
            .fold(0.0_f64, f64::max);
        within(
            "permanent relative indentation [m]",
            maximum_depth,
            self.envelope.maximum_plastic_depth_m,
        )?;

        // Include the farthest corner of each loaded cell. An area-equivalent
        // disk would underestimate separated or elongated contact footprints.
        let grid = self.model.config().grid;
        let mut patch_radius = 0.0_f64;
        for (index, pressure) in candidate.pressure_pa().iter().enumerate() {
            if *pressure > 0.0 {
                let x = (index % grid.cells_x) as f64 + 0.5 - 0.5 * grid.cells_x as f64;
                let y = (index / grid.cells_x) as f64 + 0.5 - 0.5 * grid.cells_y as f64;
                let corner_x = (x.abs() + 0.5) * grid.cell_width_m;
                let corner_y = (y.abs() + 0.5) * grid.cell_depth_m;
                patch_radius = patch_radius.max(corner_x.hypot(corner_y));
            }
        }
        for side in 0..2 {
            within(
                "patch radius / homogeneous depth",
                patch_radius / self.envelope.half_space_depth_m[side],
                self.envelope.maximum_patch_to_depth[side],
            )?;
        }
        cx.checkpoint()
            .map_err(|_| FiniteGapPlasticError::Contact(super::FiniteGapContactError::Cancelled))?;
        Ok(())
    }
}

fn within(quantity: &'static str, value: f64, limit: f64) -> Result<(), SourcePlasticContactError> {
    if !value.is_finite() {
        return Err(SourcePlasticContactError::Invalid(quantity));
    }
    if value > limit {
        return Err(SourcePlasticContactError::OutsideApplicability {
            quantity,
            value,
            limit,
        });
    }
    Ok(())
}
