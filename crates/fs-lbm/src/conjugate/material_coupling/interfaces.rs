//! Source-bound electrical interfaces in the existing material feedback loop.
//!
//! Contact and bulk laws are re-evaluated at each proposed temperature, then
//! the original spatial current solver supplies ONE combined Joule source to
//! the thermal owner. Final contact receipts belong to the electrical field
//! evaluated at the returned temperature, not an earlier Picard iterate.

use fs_exec::CancelGate;
use fs_matdb::SurfaceSpec;
use fs_material::electrical_interface::ResolvedElectricalInterface;
use crate::conjugate::{
    ChtError, ElectricalContact, ElectricSetup, EnergySolution, ThermalSetup, Voxel, VoxelDomain,
};
use super::{MaterialCoupledSolution, MaterialCouplingConfig, invalid, poll, solve_with_electrical_owner};

/// One ordered source-supported contact and its immutable spatial assignment.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialContact {
    source: ResolvedElectricalInterface,
    face: ElectricalContact,
}

impl MaterialContact {
    /// Bind the actual receiving material/finish assignments to the source.
    /// Geometry adjacency and electrical participation are checked by the
    /// existing spatial solver. Heat partition is explicitly authored, not a
    /// material property inferred from the electrical resistance.
    pub fn try_new(
        source: ResolvedElectricalInterface,
        cells: [usize; 2],
        surfaces: [&SurfaceSpec; 2],
        heat_fraction_to_first: f64,
    ) -> Result<Self, ChtError> {
        if source.surfaces()[0] != *surfaces[0] || source.surfaces()[1] != *surfaces[1] {
            return Err(invalid("material.contact_surfaces", "source and receiving surfaces differ in material state, order or texture frame"));
        }
        if cells[0] == cells[1] || !heat_fraction_to_first.is_finite()
            || !(0.0..=1.0).contains(&heat_fraction_to_first)
        {
            return Err(invalid("material.contact_assignment", "needs distinct cells and a finite heat fraction in [0,1]"));
        }
        let face = ElectricalContact {
            cells, resistance_ohm_m2: source.resistance_ohm_m2(), heat_fraction_to_first,
        };
        Ok(Self { source, face })
    }
    /// Exact source query, original surface system and all property receipts.
    #[must_use]
    pub const fn source(&self) -> &ResolvedElectricalInterface { &self.source }
    /// Numerical lowering used by the existing current/heat solver.
    #[must_use]
    pub const fn face(&self) -> &ElectricalContact { &self.face }
}

/// Existing coupled fields plus the contacts used by the FINAL electrical solve.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialInterfaceSolution {
    /// Existing temperature/current/Joule fields and consistency defects.
    pub coupled: MaterialCoupledSolution,
    /// Final source-backed contacts in caller declaration order.
    pub contacts: Vec<MaterialContact>,
}

/// Couple source-supported electrical interfaces, bulk materials and heat.
///
/// `contacts(temperature)` resolves every assigned interface at the proposed
/// state. Interface temperature, pressure, history and other coordinates must
/// follow the caller's explicit model; no arithmetic mean is silently chosen.
/// The bulk law follows `ElectricSetup::solve_with_cell_resistivity`: declared
/// insulation is `None`, unsupported material data is an error, and conducting
/// fluid cells do not become mechanical solids.
///
/// The SAME existing feedback loop enforces both temperature and local Joule
/// consistency. Contact heat is already in `heated.power_w`; neither this
/// driver nor the thermal callback may add it again. The final source objects
/// are retained without re-querying or relabelling another material's receipt.
/// Fixed geometry, passive DC Ohmic interfaces; no contact mechanics, breakdown,
/// chemical evolution, thermal contact law or global convergence theorem.
///
/// # Errors
/// Invalid numerical budgets/fields, too many contact assignments, source,
/// electrical or thermal refusals, cancellation and exhausted coupling work.
/// No partial solution is returned. Captured callback state remains caller-owned.
#[allow(clippy::too_many_arguments)] // explicit independent laws, state and budgets
pub fn solve(
    domain: &VoxelDomain,
    electric: &ElectricSetup,
    thermal: &ThermalSetup,
    initial_temperature_k: &[f64],
    config: MaterialCouplingConfig,
    mut resistivity: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
    mut contacts: impl FnMut(&[f64]) -> Result<Vec<MaterialContact>, ChtError>,
    energy: impl FnMut(&ThermalSetup, &[f64]) -> Result<EnergySolution, ChtError>,
    gate: &CancelGate,
) -> Result<MaterialInterfaceSolution, ChtError> {
    let mut final_contacts = Vec::new();
    let coupled = solve_with_electrical_owner(
        domain, thermal, initial_temperature_k, config,
        |temperature| {
            poll(gate)?;
            let proposed = contacts(temperature)?;
            if proposed.len() > domain.cell_count().saturating_mul(3) {
                return Err(invalid("material.contacts", "contact count exceeds the number of Cartesian faces"));
            }
            let faces: Vec<_> = proposed.iter().map(|contact| *contact.face()).collect();
            let field = electric.solve_with_cell_contacts(domain, temperature, &faces, &mut resistivity, gate)?;
            final_contacts = proposed;
            Ok(field)
        }, energy, gate,
    )?;
    Ok(MaterialInterfaceSolution { coupled, contacts: final_contacts })
}


#[cfg(feature = "enthalpy-coupling")]
impl MaterialContact {
    /// Consume this exact source-supported surface assignment in the shared
    /// material/enthalpy contact lane. No property is re-queried or copied into
    /// an unrelated receipt. The expected ordered face is the current callback's
    /// authored topology; a stale or reversed assignment refuses before current.
    pub fn into_enthalpy_response(self, expected_cells: [usize; 2]) -> Result<
        crate::conjugate::enthalpy::transport::interfaces::ElectricalInterfaceResponse<Self>, ChtError,
    > {
        if self.face.cells != expected_cells {
            return Err(invalid("material.contact_assignment", "resolved source contact does not match the ordered enthalpy face"));
        }
        let face = self.face;
        crate::conjugate::enthalpy::transport::interfaces::ElectricalInterfaceResponse::try_new(
            self, face.resistance_ohm_m2, face.heat_fraction_to_first,
        ).map_err(|e| invalid("material.contact_response", &e.to_string()))
    }
}
