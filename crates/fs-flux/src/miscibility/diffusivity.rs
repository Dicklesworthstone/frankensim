//! Original-source kinetics for the existing conservative miscibility field.
//!
//! A single explicitly assigned binary Maxwell--Stefan source replaces the
//! authored constant face diffusivities, while geometry, thermodynamics,
//! inventories and elapsed time survive. At each canonical pair stage the
//! original source is sampled at the actual temperature, both current endpoint
//! compositions and their arithmetic mean. That positive face coefficient is
//! held during the scalar convex-split solve. This is first-order lagged
//! nonlinear kinetics, not a fully implicit D(T,x) scheme. The existing adaptive
//! timeline refreshes it on every coarse and fine stage.
//!
//! The same source survives prescribed-temperature changes, closed adiabatic
//! evolution and advection/diffusion composition. Every published endpoint is
//! checked against its original source support; rejected trials publish nothing.
//! Kinetic changes do not alter free energy, supply heat or shift equilibrium.
//! Temperature gradients, thermodiffusion and arbitrary heterogeneous chemical
//! potentials require their own constitutive models and are not inferred here.
//!
//! Material queries and transient receipts allocate separately from the original
//! numerical workspace. Each sweep makes O(cells + faces) bounded queries;
//! scalar root iterations reuse the stage coefficient without further queries.

use super::{
    BinaryMiscibilityError, BinaryMiscibilityNetwork, BinaryMiscibilityState, Result,
    face_mobility, poll,
};
use fs_material::liquid::solution::BinaryLiquidError;
use fs_material::liquid::solution::diffusion::MaterialBinaryDiffusivity;
use std::sync::Arc;

/// Atomic kinetic binding and the unchanged original physical field.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityDiffusivityCandidate {
    network: BinaryMiscibilityNetwork,
    state: BinaryMiscibilityState,
}

impl BinaryMiscibilityDiffusivityCandidate {
    pub const fn network(&self) -> &BinaryMiscibilityNetwork {
        &self.network
    }
    pub const fn state(&self) -> &BinaryMiscibilityState {
        &self.state
    }
    pub fn into_parts(self) -> (BinaryMiscibilityNetwork, BinaryMiscibilityState) {
        (self.network, self.state)
    }
}

impl BinaryMiscibilityNetwork {
    /// Bind a source to this SAME accepted field without resetting its clock,
    /// composition, paired inventories or compensated low words. The source's
    /// ordered species must match the thermodynamic law exactly. It replaces
    /// every authored constant D; a closed geometric face still has zero area.
    /// Explicit zero source D closes only that stage's molecular kinetics.
    ///
    /// Face diffusivities exposed by `faces()` become snapshots at this binding
    /// (and later temperature rebindings). Actual transport always re-evaluates
    /// the source, so a snapshot is never a mutable last-state coefficient cache.
    /// The returned network owns the new state. Prepared advection must be
    /// explicitly rebound with `rebind_advection` before further transport.
    pub fn bind_diffusivity_source(
        &self,
        initial: &BinaryMiscibilityState,
        source: MaterialBinaryDiffusivity,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityDiffusivityCandidate> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        if source.species() != self.model().thermodynamics().species() {
            return Err(BinaryMiscibilityError::Invalid(
                "diffusivity source and miscibility model require the same ordered species",
            ));
        }
        let mut owner = (*self.owner).clone();
        let mut faces = owner.faces.to_vec();
        for face in &mut faces {
            let [i, j] = face.physical.cells;
            let resolved = source
                .resolve_face(
                    self.model().temperature(),
                    [initial.fractions[i], initial.fractions[j]],
                    &mut cancelled,
                )
                .map_err(source_error)?;
            face.physical.diffusivity_m2_s = resolved.diffusivity_m2_s;
            face.mobility_mol2_j_s = face_mobility(
                self.model(),
                face.physical.area_m2 / face.physical.distance_m,
                resolved.diffusivity_m2_s,
            )?;
        }
        owner.faces = faces.into();
        owner.diffusivity_source = Some(Arc::new(source));
        let network = Self {
            owner: Arc::new(owner),
        };
        let mut state = initial.clone();
        state.owner = Arc::clone(&network.owner);
        network.qualify_diffusivity_sources(&state, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(BinaryMiscibilityDiffusivityCandidate { network, state })
    }

    /// Original immutable molecular source, with inspectable pinned claims and
    /// held physical context. None denotes the authored constant-D path.
    pub fn diffusivity_source(&self) -> Option<&MaterialBinaryDiffusivity> {
        self.owner.diffusivity_source.as_deref()
    }
}

pub(super) fn source_error(error: BinaryLiquidError) -> BinaryMiscibilityError {
    if error == BinaryLiquidError::Cancelled {
        BinaryMiscibilityError::Cancelled
    } else {
        BinaryMiscibilityError::DiffusivitySource(error)
    }
}
