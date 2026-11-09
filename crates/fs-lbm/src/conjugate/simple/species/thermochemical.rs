//! Turbulent mass-frame mixing through the original finite gas inventories.
//!
//! The native face law supplies integrated species streams. The material owner
//! carries each original donor's partial enthalpy and reconstructs one common
//! temperature, pressure, caloric and transport endpoint per affected cell.

use super::{CancelGate, ChtError, FvFlow, MaterialAdvectionStep, VoxelDomain, poll};
use fs_material::thermochemical_gas::energy::{
    GasEnergyError, GasEnergySolve,
    parcel::{
        ClosedGasParcel,
        transfer::{
            snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
            species::network::{
                self, SpeciesNetworkCandidate, SpeciesNetworkError, SpeciesNetworkHeatTransfer,
                SpeciesNetworkSolve, SpeciesNetworkTransfer,
            },
        },
    },
};
use fs_math::dd::Dd;
use fs_qty::chemistry::SpeciesId;
use std::collections::BTreeMap;

mod adaptive;
mod heat;
pub use adaptive::TurbulentGasSpeciesAdvance;

/// One complete frozen turbulent-mixing snapshot and its original gas states.
#[derive(Clone, Debug, PartialEq)]
pub struct TurbulentGasSpeciesStep {
    /// Original source-bearing endpoints and species/partial-enthalpy balances.
    pub exchange: SpeciesNetworkCandidate,
    /// Entire positive physical interval sampled from the original state [s].
    pub duration_s: f64,
    /// Fluid-fluid faces with positive conductance, even at zero gradient.
    pub conducting_faces: usize,
    /// Largest absolute sum of the integrated species streams on one face [kg].
    pub maximum_face_mass_defect_kg: f64,
    /// Initial source-derived species partition against the accepted mass [kg].
    pub initial_carrier_partition_defect_kg: f64,
    /// Largest actual endpoint departure from the immutable carrier mass [kg].
    pub maximum_carrier_mass_defect_kg: f64,
    /// Sum of absolute endpoint departures from the immutable carrier mass [kg].
    pub absolute_carrier_mass_defect_kg: f64,
}

struct PreparedCell {
    /// Original typed basis conversion in the model's canonical species order.
    mass_fractions: Vec<f64>,
}

struct Face {
    cells: [u64; 2],
    conductance_kg_s: f64,
}

/// Frozen native geometry and turbulent coefficients for one physical window.
/// Source composition, temperature and partial enthalpy refresh at every probe.
struct PreparedTurbulentGasMixing<'a> {
    pub(super) domain: &'a VoxelDomain,
    pub(super) carrier: &'a MaterialAdvectionStep,
    faces: Vec<Face>,
    parcel_count: usize,
}

fn refuse(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.turbulent_gas_species",
        reason: reason.into(),
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn invalid(reason: &'static str) -> SpeciesNetworkError {
    SpeciesNetworkError::Invalid(reason)
}

fn check(cancelled: &mut dyn FnMut() -> bool) -> Result<(), SpeciesNetworkError> {
    if cancelled() {
        Err(SpeciesNetworkError::Cancelled)
    } else {
        Ok(())
    }
}

fn material_error(error: SpeciesNetworkError) -> ChtError {
    if error == SpeciesNetworkError::Cancelled {
        ChtError::Cancelled
    } else {
        refuse(&format!("original gas species exchange refused: {error}"))
    }
}

// Merge the ORIGINAL canonical identities. Exact absence means zero fraction;
// neither a fictitious trace species nor an index-based source is introduced.
fn visit_face(
    face: &Face,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    prepared: &BTreeMap<u64, PreparedCell>,
    duration_s: f64,
    mut visit: impl FnMut(&SpeciesId, f64) -> Result<(), SpeciesNetworkError>,
) -> Result<f64, SpeciesNetworkError> {
    let [left, right] = face.cells;
    let a = parcels[&left].model();
    let b = parcels[&right].model();
    let sources_a = a.thermodynamics().components();
    let sources_b = b.thermodynamics().components();
    let fractions_a = &prepared[&left].mass_fractions;
    let fractions_b = &prepared[&right].mass_fractions;
    let (mut i, mut j) = (0, 0);
    let mut net = Dd::ZERO;
    while i < sources_a.len() || j < sources_b.len() {
        let (species, fraction_a, fraction_b) = match (sources_a.get(i), sources_b.get(j)) {
            (Some(x), Some(y)) if x.species() == y.species() => {
                if x != y || a.transport_materials()[i] != b.transport_materials()[j] {
                    return Err(invalid(
                        "shared species has conflicting original caloric or transport sources",
                    ));
                }
                let row = (x.species(), fractions_a[i], fractions_b[j]);
                i += 1;
                j += 1;
                row
            }
            (Some(x), Some(y)) if x.species() < y.species() => {
                let row = (x.species(), fractions_a[i], 0.0);
                i += 1;
                row
            }
            (Some(x), None) => {
                let row = (x.species(), fractions_a[i], 0.0);
                i += 1;
                row
            }
            (_, Some(y)) => {
                let row = (y.species(), 0.0, fractions_b[j]);
                j += 1;
                row
            }
            (None, None) => break,
        };
        let gradient = fraction_a - fraction_b;
        if gradient == 0.0 {
            continue;
        }
        let mass = (face.conductance_kg_s * duration_s) * gradient;
        if !mass.is_finite() || mass == 0.0 {
            return Err(invalid(
                "nonzero turbulent species transfer is unrepresentable",
            ));
        }
        net = net + Dd::from_f64(mass);
        visit(species, mass)?;
    }
    let net = net.to_f64();
    if !net.is_finite() {
        return Err(invalid(
            "turbulent face mass-frame account is unrepresentable",
        ));
    }
    Ok(net)
}

impl FvFlow {
    /// Mix original ideal-gas material states with turbulent partial enthalpy.
    ///
    /// Parcel IDs are native flat voxel indices, with exactly one parcel for
    /// every fluid cell and none for solids. Each original parcel's mass must
    /// equal the immutable accepted carrier mass, and its volume must equal
    /// the native cell volume. Species and mass fractions come exclusively
    /// from that parcel's original thermochemical model and typed conversion.
    /// Common species on a conducting face must retain the same original
    /// caloric and transport sources. Source identities are never inferred.
    ///
    /// This flow's MOMENTUM eddy viscosity supplies `D_t = nu_t / Sc_t` with an
    /// explicitly positive Schmidt number, independently of heat-wall policy.
    /// The shared graded-cell resistance gives
    /// `G = A / (h_l/(2 rho_l D_l) + h_r/(2 rho_r D_r))` [kg/s]. Each stream is
    /// `dm_k = duration_s G (Y_lk - Y_rk)`, evaluated from ONE unchanged state.
    /// Solid, exterior and zero-diffusivity faces have no species flux.
    ///
    /// The original material snapshot owner combines all donor streams before
    /// rebuilding each affected cell once. Every stream carries its donor's
    /// actual partial enthalpy `h_k(T_donor)`. Its existing reference-invariant
    /// NASA caloric inversion conserves internal energy in the fixed volumes
    /// and returns the same updated source state for thermal, momentum and
    /// acoustic consumers. Formation energy is never added as a second heat.
    ///
    /// The returned source state owns the updated temperature, thermodynamic
    /// pressure, enthalpy and transport. The carrier supplies only mass and
    /// geometry admission here; its stored specific enthalpy is not a second
    /// caloric authority. SIMPLEC's gauge pressure is not thermodynamic pressure.
    /// Subsequent momentum/material coupling must consume these source endpoints.
    /// Matching the original flow geometry, mask and epoch remains the caller's
    /// obligation; dimensions and the accepted carrier mask are checked.
    ///
    /// This is a bounded, FIRST-ORDER explicit mixing component: at most the
    /// original 256 parcels and 4096 nonzero species streams. It is not the
    /// exponential/Strang isothermal kernel. Gross species withdrawals from
    /// every adjacent face share the original allowance; arriving material
    /// cannot replenish a donor during this interval. Excessive durations refuse
    /// without clipping. Refine time and space for accuracy. The mass tolerance
    /// additionally bounds source partition, each face's zero-sum mass defect,
    /// and the sum of endpoint carrier-mass departures. No fraction is repaired.
    ///
    /// Heat conduction, turbulent heat flux, molecular diffusion, chemical
    /// reactions and the next pressure/momentum evolution keep their own owners.
    /// This gradient-diffusion component supplies no independent entropy or
    /// turbulence-validation claim. Any refusal or cancellation publishes no
    /// partial result and leaves the flow, carrier and source parcels unchanged.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn diffuse_turbulent_gas_species(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        turbulent_schmidt: f64,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        controls: SpeciesNetworkSolve,
        gate: &CancelGate,
    ) -> Result<TurbulentGasSpeciesStep, ChtError> {
        if !positive(duration_s) {
            return Err(refuse("positive turbulent mixing duration required"));
        }
        PreparedTurbulentGasMixing::new(
            self,
            domain,
            carrier,
            parcels,
            turbulent_schmidt,
            endpoints,
            controls,
            gate,
        )?
        .step_candidate(parcels, duration_s, endpoints, controls, &[], &mut || {
            gate.is_requested()
        })
        .map_err(material_error)
    }
}

impl<'a> PreparedTurbulentGasMixing<'a> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        flow: &FvFlow,
        domain: &'a VoxelDomain,
        carrier: &'a MaterialAdvectionStep,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        turbulent_schmidt: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        controls: SpeciesNetworkSolve,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        let cells = domain.cell_count();
        if !positive(turbulent_schmidt)
            || !positive(controls.maximum_species_fraction)
            || controls.maximum_species_fraction >= 1.0
            || !positive(controls.mass_tolerance_kg)
            || !positive(controls.energy_tolerance_j)
            || parcels.is_empty()
            || parcels.len() > MAX_SNAPSHOT_PARCELS
            || parcels.len() != domain.fluid_count()
            || !parcels.keys().eq(endpoints.keys())
            || flow.field.dims() != domain.dims()
            || flow.eddy_viscosity_m2_s.len() != cells
        {
            return Err(refuse(
                "bounded gas/flow shape and positive mixing controls required",
            ));
        }
        carrier.validate_diffusion_carrier(domain, gate)?;
        // Check all source-to-carrier bindings and solid entries before preparing
        // any new workspace. Original source objects remain borrowed throughout.
        for (cell, &eddy) in flow.eddy_viscosity_m2_s.iter().enumerate() {
            if cell % 256 == 0 {
                poll(gate)?;
            }
            if !eddy.is_finite() || eddy < 0.0 || (!domain.is_fluid(cell) && eddy != 0.0) {
                return Err(refuse(
                    "eddy viscosity must be finite, nonnegative and zero in solids",
                ));
            }
        }
        for (&id, parcel) in parcels {
            poll(gate)?;
            let cell = usize::try_from(id).map_err(|_| refuse("parcel cell ID overflows"))?;
            if cell >= cells
                || !domain.is_fluid(cell)
                || parcel.mass_kg() != carrier.mass_kg[cell]
                || parcel.volume_m3() != domain.volume(cell)
            {
                return Err(refuse(
                    "source parcel must retain the exact accepted cell mass and volume",
                ));
            }
        }
        let mut diffusivities = BTreeMap::new();
        for &id in parcels.keys() {
            poll(gate)?;
            let eddy = flow.eddy_viscosity_m2_s[id as usize];
            let diffusivity = eddy / turbulent_schmidt;
            if !diffusivity.is_finite() || (eddy > 0.0 && diffusivity == 0.0) {
                return Err(refuse("positive turbulent diffusivity is unrepresentable"));
            }
            diffusivities.insert(id, diffusivity);
        }
        let mut faces = Vec::new();
        faces
            .try_reserve_exact(3 * parcels.len())
            .map_err(|_| refuse("unable to allocate turbulent gas faces"))?;
        for (&id, &diffusivity) in &diffusivities {
            poll(gate)?;
            if diffusivity == 0.0 {
                continue;
            }
            for side in [1, 3, 5] {
                let Some(next) = domain
                    .neighbor(id as usize, side)
                    .filter(|&next| domain.is_fluid(next))
                else {
                    continue;
                };
                let next_diffusivity = diffusivities[&(next as u64)];
                if next_diffusivity == 0.0 {
                    continue;
                }
                faces.push(Face {
                    cells: [id, next as u64],
                    conductance_kg_s: carrier.diffusion_face_conductance(
                        domain,
                        [id as usize, next],
                        side / 2,
                        [diffusivity, next_diffusivity],
                    )?,
                });
            }
        }
        poll(gate)?;
        Ok(Self {
            domain,
            carrier,
            faces,
            parcel_count: parcels.len(),
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn step_candidate(
        &self,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        duration_s: f64,
        endpoints: &BTreeMap<u64, GasEnergySolve>,
        controls: SpeciesNetworkSolve,
        heat: &[SpeciesNetworkHeatTransfer],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TurbulentGasSpeciesStep, SpeciesNetworkError> {
        check(cancelled)?;
        if !positive(duration_s) || heat.len() > MAX_SNAPSHOT_TRANSFERS {
            return Err(invalid(
                "positive mixing duration and bounded heat rows required",
            ));
        }
        self.audit_carrier(parcels, controls.mass_tolerance_kg, cancelled)?;
        let mut prepared = BTreeMap::new();
        let mut partition = Dd::ZERO;
        for (&id, parcel) in parcels {
            check(cancelled)?;
            let thermo = parcel.model().thermodynamics();
            let molar_masses: Vec<_> = thermo.components().iter().map(|s| s.molar_mass()).collect();
            let composition = thermo
                .mole_composition()
                .to_mass_fractions(&molar_masses)
                .map_err(|_| SpeciesNetworkError::Parcel {
                    id,
                    error: GasEnergyError::Invalid("original species basis conversion refused"),
                })?;
            let mut cell_partition = -Dd::from_f64(parcel.mass_kg());
            for &fraction in composition.fractions() {
                let mass = parcel.mass_kg() * fraction;
                if !positive(fraction) || fraction > 1.0 || !positive(mass) {
                    return Err(invalid(
                        "original positive species inventory is unrepresentable",
                    ));
                }
                cell_partition = cell_partition + Dd::from_f64(mass);
            }
            partition = partition + cell_partition.abs();
            prepared.insert(
                id,
                PreparedCell {
                    mass_fractions: composition.fractions().to_vec(),
                },
            );
        }
        let partition = partition.to_f64();
        if !partition.is_finite() || partition > controls.mass_tolerance_kg {
            return Err(invalid(
                "original source species do not partition the accepted carrier mass",
            ));
        }
        let mut transfer_count = 0usize;
        let mut maximum_face_mass_defect_kg: f64 = 0.0;
        for face in &self.faces {
            check(cancelled)?;
            let defect = visit_face(face, parcels, &prepared, duration_s, |_, _| {
                if transfer_count + heat.len() == MAX_SNAPSHOT_TRANSFERS {
                    return Err(invalid(
                        "turbulent species/heat streams exceed the original snapshot bound",
                    ));
                }
                transfer_count += 1;
                Ok(())
            })?
            .abs();
            if defect > controls.mass_tolerance_kg {
                return Err(invalid(
                    "turbulent face streams violate the zero-net-mass allowance",
                ));
            }
            maximum_face_mass_defect_kg = maximum_face_mass_defect_kg.max(defect);
        }
        // Count and admit the original streams before allocating their payload.
        // Face coefficients and typed compositions are already prepared once.
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(transfer_count)
            .map_err(|_| invalid("unable to allocate turbulent species streams"))?;
        for face in &self.faces {
            check(cancelled)?;
            visit_face(face, parcels, &prepared, duration_s, |species, mass| {
                let [left, right] = face.cells;
                let (donor, receiver) = if mass > 0.0 {
                    (left, right)
                } else {
                    (right, left)
                };
                transfers.push(SpeciesNetworkTransfer {
                    donor,
                    receiver,
                    species: species.clone(),
                    mass_kg: mass.abs(),
                });
                Ok(())
            })?;
        }
        let exchange = network::exchange_candidate_with_heat(
            parcels,
            &transfers,
            heat,
            endpoints,
            controls,
            &mut *cancelled,
        )?;
        let (maximum_carrier_mass_defect_kg, absolute_carrier_mass_defect_kg) =
            self.audit_carrier(exchange.parcels(), controls.mass_tolerance_kg, cancelled)?;
        check(cancelled)?;
        Ok(TurbulentGasSpeciesStep {
            exchange,
            duration_s,
            conducting_faces: self.faces.len(),
            maximum_face_mass_defect_kg,
            initial_carrier_partition_defect_kg: partition,
            maximum_carrier_mass_defect_kg,
            absolute_carrier_mass_defect_kg,
        })
    }

    /// Intermediate and final probes retain the SAME original fixed carrier.
    /// Existing representational mass defects are admitted, never reset to zero.
    fn audit_carrier(
        &self,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        mass_tolerance_kg: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(f64, f64), SpeciesNetworkError> {
        check(cancelled)?;
        if parcels.len() != self.parcel_count || !positive(mass_tolerance_kg) {
            return Err(invalid(
                "original gas parcel count and positive mass budget required",
            ));
        }
        let mut maximum_carrier_mass_defect_kg: f64 = 0.0;
        let mut carrier_defect = Dd::ZERO;
        for (&id, parcel) in parcels {
            check(cancelled)?;
            let cell = usize::try_from(id).map_err(|_| invalid("parcel cell ID overflows"))?;
            if cell >= self.domain.cell_count()
                || !self.domain.is_fluid(cell)
                || parcel.volume_m3() != self.domain.volume(cell)
            {
                return Err(invalid(
                    "gas probe changed its original cell identity or volume",
                ));
            }
            let defect = (parcel.mass_kg() - self.carrier.mass_kg[cell]).abs();
            maximum_carrier_mass_defect_kg = maximum_carrier_mass_defect_kg.max(defect);
            carrier_defect = carrier_defect + Dd::from_f64(defect);
        }
        let absolute_carrier_mass_defect_kg = carrier_defect.to_f64();
        if !absolute_carrier_mass_defect_kg.is_finite()
            || absolute_carrier_mass_defect_kg > mass_tolerance_kg
        {
            return Err(invalid(
                "turbulent gas endpoints violate the fixed carrier-mass allowance",
            ));
        }
        Ok((
            maximum_carrier_mass_defect_kg,
            absolute_carrier_mass_defect_kg,
        ))
    }
}
