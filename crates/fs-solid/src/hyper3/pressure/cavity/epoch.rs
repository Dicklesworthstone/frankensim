//! A physical gas reference bound to an accepted deformed cavity geometry.
//!
//! Heat exchange may change entropy at fixed current volume. Starting the next
//! mechanical isentrope from that parcel requires no fictitious expansion back
//! to the undeformed solid. This owner admits the actual gas/geometry pair;
//! original material owners still supply every caloric and acoustic operation.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::energy::entropy::{GasEntropyError, GasEntropySolve};
use fs_material::thermochemical_gas::energy::parcel::mechanics::GasMechanicalLoad;
use fs_material::thermochemical_gas::energy::parcel::{ClosedGasParcel, GasVolumeStep};

use super::{AdiabaticCavityGas, ClosedTetCavity, HyperTetError, checkpoint, invalid};

/// Shared immutable geometric anchor; Newton gas clones never copy the mesh.
#[derive(Debug)]
struct Anchor {
    cavity: Arc<ClosedTetCavity>,
    positions_m: Vec<[f64; 3]>,
    volume_m3: f64,
}

/// One actual gas parcel and its admitted accepted-geometry reference.
///
/// Construct this after a physical gas update at the current cavity volume.
/// Neither gas entropy nor energy is mapped to the undeformed volume. The
/// supplied prepared cavity is retained by identity; a different chamber or a
/// separately cloned preparation cannot reuse this binding.
///
/// The sealed gas interface delegates to the original ideal or Peng--Robinson
/// parcel, including phase admission, source calorics and acoustic stiffness.
/// Only compact anchor positions and cavity identity are added, shared through
/// an Arc across every Newton clone. Ordinary unwrapped gas references keep
/// their original strict undeformed-volume admission.
///
/// Accepted loads contain output-only epochs. Continue with the original
/// reference for the same isentrope. To begin a new physical interval after
/// heating, extract the accepted parcel with [`Self::into_gas`], perform the
/// original gas update, and explicitly bind it at accepted solid displacement
/// with [`Self::try_new`]. Even an output at unchanged volume is not silently
/// promoted to a new source reference. An already admitted reference may use
/// [`Self::with_gas_at_anchor`] for disposable heat probes at that same geometry.
/// This type performs no heat exchange, gas momentum integration or coupled
/// energy acceptance itself.
#[derive(Clone, Debug)]
pub struct CavityGasEpoch<G = ClosedGasParcel> {
    gas: G,
    anchor: Arc<Anchor>,
    reference: bool,
}

impl<G: AdiabaticCavityGas> CavityGasEpoch<G> {
    /// Bind a physical parcel to this exact accepted complete-solid geometry.
    ///
    /// The existing geometry owner admits every solid element and the positive
    /// cavity volume. The gas volume must match the represented accepted volume
    /// exactly. This method neither changes the parcel nor performs an entropy
    /// solve; all source identity, inventory and caloric data remain original.
    ///
    /// # Errors
    /// Invalid/deformed mesh, a different gas volume, or cancellation. Refusal
    /// publishes no geometry-bound reference and mutates no borrowed state.
    pub fn try_new(
        cavity: Arc<ClosedTetCavity>,
        accepted_displacement_m: &[f64],
        gas: G,
        cx: &Cx<'_>,
    ) -> Result<Self, HyperTetError> {
        checkpoint(cx)?;
        let (positions_m, volume_m3) =
            cavity.positions_and_volume(accepted_displacement_m, cx)?;
        if gas.volume_m3().to_bits() != volume_m3.to_bits() {
            return Err(invalid(
                "cavity gas epoch requires the exact accepted physical gas volume",
            ));
        }
        checkpoint(cx)?;
        Ok(Self {
            gas,
            anchor: Arc::new(Anchor {
                cavity,
                positions_m,
                volume_m3,
            }),
            reference: true,
        })
    }


    /// Explicitly bind another physical parcel at this same admitted anchor.
    ///
    /// Use this for source-owned heating probes at fixed accepted geometry.
    /// The new parcel must have exactly the anchor volume and this epoch must
    /// still be a reference, not an output from a mechanical evaluation. The
    /// immutable geometry and exact cavity identity are shared without another
    /// complete-mesh admission or compact-position copy.
    ///
    /// This operation does not infer heat, source continuity or inventory
    /// transfer. The caller supplies the original physical gas owner's admitted
    /// replacement and retains its heat/source accounting. A new geometry still
    /// requires [`Self::try_new`].
    ///
    /// # Errors
    /// An output-only epoch, a different represented volume, or cancellation.
    /// Refusal leaves this admitted reference unchanged.
    pub fn with_gas_at_anchor(
        &self,
        gas: G,
        cx: &Cx<'_>,
    ) -> Result<Self, HyperTetError> {
        checkpoint(cx)?;
        if !self.reference || gas.volume_m3().to_bits() != self.anchor.volume_m3.to_bits() {
            return Err(invalid(
                "same-anchor gas update requires a reference epoch and exact accepted volume",
            ));
        }
        checkpoint(cx)?;
        Ok(Self {
            gas,
            anchor: Arc::clone(&self.anchor),
            reference: true,
        })
    }

    /// Original immutable physical gas parcel, including the complete source.
    pub const fn gas(&self) -> &G {
        &self.gas
    }

    /// Recover the physical gas for source-owned heating or a fresh geometry bind.
    pub fn into_gas(self) -> G {
        self.gas
    }

    /// The exact immutable prepared cavity retained by this epoch.
    pub fn cavity(&self) -> &Arc<ClosedTetCavity> {
        &self.anchor.cavity
    }

    /// Actual accepted reference volume; it need not be the undeformed volume.
    pub fn anchor_volume_m3(&self) -> f64 {
        self.anchor.volume_m3
    }
}

impl<G: AdiabaticCavityGas> AdiabaticCavityGas for CavityGasEpoch<G> {
    type State = G::State;

    fn volume_m3(&self) -> f64 {
        self.gas.volume_m3()
    }

    fn reference_positions_m<'a>(
        &'a self,
        cavity: &'a ClosedTetCavity,
    ) -> Result<&'a [[f64; 3]], HyperTetError> {
        if !self.reference
            || !std::ptr::eq(self.anchor.cavity.as_ref(), cavity)
            || self.gas.volume_m3().to_bits() != self.anchor.volume_m3.to_bits()
        {
            return Err(invalid(
                "cavity gas epoch requires its exact admitted reference and cavity identity",
            ));
        }
        Ok(&self.anchor.positions_m)
    }

    fn at_volume_isentropic(
        &self,
        volume_m3: f64,
        controls: GasEntropySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, GasVolumeStep), GasEntropyError> {
        let (gas, step) = self.gas.at_volume_isentropic(volume_m3, controls, cancelled)?;
        Ok((
            Self {
                gas,
                anchor: Arc::clone(&self.anchor),
                reference: false,
            },
            step,
        ))
    }

    fn mechanical_load_with_limit<'a>(
        &'a self,
        volume_m3: f64,
        gradient: &'a [f64],
        exterior_pressure_pa: f64,
        maximum_coordinates: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GasMechanicalLoad<'a, Self::State>, GasEntropyError> {
        self.gas.mechanical_load_with_limit(
            volume_m3,
            gradient,
            exterior_pressure_pa,
            maximum_coordinates,
            cancelled,
        )
    }

    fn stored_energy_change_to(
        &self,
        after: &Self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, GasEnergyError> {
        self.gas.stored_energy_change_to(&after.gas, cancelled)
    }
}
