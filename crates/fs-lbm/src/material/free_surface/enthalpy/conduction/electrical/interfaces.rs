//! Pinned original interface sources on the moving-liquid electrical timeline.
//!
//! The caller supplies physical state coordinates and explicit connectivity.
//! Resistance always comes from refreshing the original complete source bundle;
//! one existing current solve owns its heat, motional voltage and Lorentz force.

use super::*;
use crate::conjugate::enthalpy::transport::interfaces::{
    ElectricalInterfaceResponse, ElectricalInterfaceState, SwitchingInterfaceEpoch,
    SwitchingInterfaceSample,
};
use fs_matdb::{InterfaceSystemCard, QueryPoint};
use fs_material::electrical_interface::ResolvedElectricalInterface;
use std::borrow::Cow;
use std::sync::Arc;

/// One original source selection and the actual ordered surface texture frames.
/// The live bulk material identities come from the receiving liquid samples,
/// not copies of the interface card's own surface declarations.
#[derive(Clone, Debug)]
pub struct FreeSurfaceElectricalInterfaceBinding {
    card: Arc<InterfaceSystemCard>,
    reference: Arc<ResolvedElectricalInterface>,
    texture_frames: [String; 2],
}
impl FreeSurfaceElectricalInterfaceBinding {
    /// Retain the original complete interface bundle and selected claims.
    /// Actual material/texture assignments are checked at every current epoch.
    pub fn try_new(
        card: Arc<InterfaceSystemCard>,
        reference: ResolvedElectricalInterface,
        texture_frames: [String; 2],
    ) -> Result<Self, Error> {
        let reference =
            ResolvedElectricalInterface::try_from_resolved(&card, reference.resolved().clone())?;
        Ok(Self {
            card,
            reference: Arc::new(reference),
            texture_frames,
        })
    }
}

/// Actual adjacent state used by this SAME bulk current evaluation.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceElectricalInterfaceSide<'a> {
    pub voxel: usize,
    pub coordinates: [usize; 3],
    pub specific_enthalpy_j_kg: f64,
    pub temperature_k: f64,
    /// Includes original bulk T/shear/transport/electrical source receipts.
    pub material: &'a ResolvedConductorState,
    /// Authored receiving texture frame, independent of the source declaration.
    pub texture_frame: &'a str,
}

/// Explicit instantaneous electrical state of an unchanged geometric face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FreeSurfaceElectricalConnectivity {
    /// The source supplies resistance; heat partition is independently declared.
    Closed { heat_fraction_to_first: f64 },
    /// No current, contact heat or motional/Lorentz transfer across this face.
    /// Original source support is still required; missing data is never Open.
    Open,
}

/// Physical query and connectivity supplied by the interface condition owner.
/// There is deliberately no independently copied resistance coefficient.
#[derive(Clone, Debug)]
pub struct FreeSurfaceElectricalInterfaceUpdate {
    pub query: QueryPoint,
    pub connectivity: FreeSurfaceElectricalConnectivity,
}

/// Ordered original source responses paired with one exact current epoch.
pub type FreeSurfaceElectricalInterfaceEpoch = SwitchingInterfaceEpoch<ResolvedElectricalInterface>;

type InterfaceLaw = dyn Fn(
        usize,
        &ElectricalContact,
        &ResolvedElectricalInterface,
        [FreeSurfaceElectricalInterfaceSide<'_>; 2],
    ) -> Result<FreeSurfaceElectricalInterfaceUpdate, Error>
    + Send
    + Sync;

#[derive(Clone)]
pub(super) struct MaterialInterfaces {
    bindings: Vec<FreeSurfaceElectricalInterfaceBinding>,
    resolve: Arc<InterfaceLaw>,
}
impl core::fmt::Debug for MaterialInterfaces {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("MaterialInterfaces")
            .field("bindings", &self.bindings)
            .finish_non_exhaustive()
    }
}

pub(super) struct ResolvedInterfaces<'a> {
    pub contacts: Cow<'a, [ElectricalContact]>,
    pub openings: Vec<[usize; 2]>,
    pub epoch: Option<Arc<FreeSurfaceElectricalInterfaceEpoch>>,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Evolve original material interfaces with the native current/heat/flow step.
    ///
    /// Bind exactly one original interface source per ordered authored contact.
    /// The callback receives the retained source selection and the SAME adjacent
    /// liquid samples used by bulk current. It supplies a complete physical
    /// QueryPoint and explicit Closed/Open state. The original material owner
    /// refreshes every property from the pinned source/claims; a newly available
    /// fallback claim cannot replace unsupported original data.
    ///
    /// Interface temperature, pressure, surface history and heat partition remain
    /// explicit caller models. No mean temperature or lattice-pressure EOS is
    /// guessed. Material/process/phase identities come from the actual liquid;
    /// the supplied texture frames are checked against the ordered interface
    /// source. Source support is checked even while electrically open.
    ///
    /// Resolved contacts/openings enter the original static or moving-current
    /// topology owner, for both electrode and isolated boundaries. Its ONE solve
    /// supplies ordered contact/bulk heat and reciprocal Lorentz force. An
    /// opening interrupts current loops without deleting liquid, thermal links,
    /// mass or enthalpy. Zero resistance remains perfect contact. Original
    /// authored coefficients are available as reference data but are not added
    /// to the resolved source resistance.
    ///
    /// The binding is retained automatically by all existing conduction, finite
    /// wall/gas, radiation, rheology and wetting steps. The callback must be pure
    /// with respect to physical history: candidate trials and endpoint support
    /// queries may repeat. Native cancellation, source, gauge, power/energy and
    /// split-motion gates publish the whole transaction or none of it.
    ///
    /// Geometry and the electrical core remain fixed/fully occupied. This adds
    /// no interface energy store, contact mechanics, induction or event-location
    /// solver. Connectivity is sampled on actual thermal substeps; discontinuous
    /// laws require declared time refinement. Replacing the electrical core
    /// requires rebinding its interfaces.
    pub fn with_material_electrical_interfaces(
        mut self,
        bindings: Vec<FreeSurfaceElectricalInterfaceBinding>,
        resolve: impl Fn(
            usize,
            &ElectricalContact,
            &ResolvedElectricalInterface,
            [FreeSurfaceElectricalInterfaceSide<'_>; 2],
        ) -> Result<FreeSurfaceElectricalInterfaceUpdate, Error>
        + Send
        + Sync
        + 'static,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        let core = self.electrical.as_mut().ok_or(Error::Invalid(
            "material interfaces require an existing electrical core",
        ))?;
        if bindings.len() != core.input.contacts.len() {
            return Err(Error::Invalid(
                "every authored electrical contact needs exactly one material binding",
            ));
        }
        // Safe ordered indexing precedes any source callback. The native current
        // owner still checks conducting participation and all solved topology.
        let mut faces = std::collections::BTreeSet::new();
        for contact in &core.input.contacts {
            poll(gate)?;
            let [a, b] = contact.cells;
            if a >= core.coordinates.len()
                || b >= core.coordinates.len()
                || !(0..6).any(|face| core.domain.neighbor(a, face) == Some(b))
                || !faces.insert((a.min(b), a.max(b)))
            {
                return Err(Error::Invalid(
                    "material electrical interfaces need unique adjacent core cells",
                ));
            }
        }
        poll(gate)?;
        core.interfaces = Some(MaterialInterfaces {
            bindings,
            resolve: Arc::new(resolve),
        });
        Ok(self)
    }
}

impl ElectricalCore {
    pub(super) fn resolve_interfaces<'a>(
        &'a self,
        time_s: f64,
        enthalpy: &[f64],
        temperature: &[f64],
        materials: &[ResolvedConductorState],
        gate: &CancelGate,
    ) -> Result<ResolvedInterfaces<'a>, Error> {
        let Some(law) = &self.interfaces else {
            return Ok(ResolvedInterfaces {
                contacts: Cow::Borrowed(&self.input.contacts),
                openings: Vec::new(),
                epoch: None,
            });
        };
        if enthalpy.len() != self.coordinates.len()
            || temperature.len() != self.coordinates.len()
            || materials.len() != self.coordinates.len()
        {
            return Err(Error::Invalid(
                "material electrical interfaces lost the actual bulk epoch",
            ));
        }
        let mut contacts = Vec::with_capacity(law.bindings.len());
        let mut openings = Vec::new();
        let mut samples = Vec::with_capacity(law.bindings.len());
        for (index, (&authored, binding)) in
            self.input.contacts.iter().zip(&law.bindings).enumerate()
        {
            poll(gate)?;
            let sides = [0, 1].map(|side| {
                let voxel = authored.cells[side];
                FreeSurfaceElectricalInterfaceSide {
                    voxel,
                    coordinates: self.coordinates[voxel],
                    specific_enthalpy_j_kg: enthalpy[voxel],
                    temperature_k: temperature[voxel],
                    material: &materials[voxel],
                    texture_frame: &binding.texture_frames[side],
                }
            });
            let update = (law.resolve)(index, &authored, &binding.reference, sides)?;
            poll(gate)?;
            let source = binding
                .reference
                .at_query_point(&binding.card, &update.query)?;
            source.require_surface_materials(
                sides[0].material.material().material(),
                sides[1].material.material().material(),
            )?;
            if source
                .surfaces()
                .iter()
                .zip(sides)
                .any(|(source, actual)| source.texture_frame != actual.texture_frame)
            {
                return Err(Error::Invalid(
                    "electrical interface source and actual ordered texture frames differ",
                ));
            }
            let response = match update.connectivity {
                FreeSurfaceElectricalConnectivity::Open => ElectricalInterfaceState::Open(source),
                FreeSurfaceElectricalConnectivity::Closed {
                    heat_fraction_to_first,
                } => {
                    let resistance = source.resistance_ohm_m2();
                    ElectricalInterfaceResponse::try_new(source, resistance, heat_fraction_to_first)
                        .map_err(|_| {
                            Error::Invalid(
                                "electrical interface heat partition must be finite and in [0,1]",
                            )
                        })?
                        .into()
                }
            };
            let sample = SwitchingInterfaceSample { authored, response };
            match sample.applied() {
                Some(contact) => contacts.push(contact),
                None => openings.push(authored.cells),
            }
            samples.push(sample);
        }
        poll(gate)?;
        Ok(ResolvedInterfaces {
            contacts: Cow::Owned(contacts),
            openings,
            epoch: Some(Arc::new(SwitchingInterfaceEpoch {
                time_s,
                contacts: samples,
            })),
        })
    }
}
