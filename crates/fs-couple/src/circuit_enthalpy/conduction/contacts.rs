//! The native contact field, original interface sources and finite heat share
//! one thermal endpoint. No guessed interface temperature or heat split.

use super::*;
use fs_em::conduction::{ConductionContact, ConductionContactResponse};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, ResolvedElectricalInterface,
};
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_qty::{
    Pressure,
    semantic::{QuantityKind, ValueForm},
};

/// Explicit thermal ownership of one original surface contact. Contact indices
/// refer to the original conduction owner's `contacts()` in canonical trace
/// order, for either planar edges or matching tetrahedral faces.
#[derive(Clone, Copy, Debug)]
pub struct ConductionContactThermal<'a> {
    pub contact: usize,
    /// Original ordered interface card, bulk cards and thermal coordinates.
    /// The selected source claim stays pinned by the accepted field contact.
    pub point: ThermalInterfacePoint<'a>,
    /// Authored physical fraction of Joule heat entering surface A, in [0,1].
    /// Its remainder enters B. No effusivity or equal split is inferred.
    pub heat_fraction_to_first: f64,
}

/// A pressure coordinate owned by an admitted mechanical contact. The native
/// mechanics owner checks compressive traction and geometry; this shared
/// source owner checks the original coordinate, source pin and thermal state.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SolvedContactPressure<'a> {
    pub contact: usize,
    pub axis: &'a str,
    pub pressure_pa: f64,
}

/// Actual surface field and its heat destinations at the accepted endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductionContactHeat<R = ConductionContactResponse> {
    pub contact: usize,
    /// Derived from the two original adjacent volume cells' thermal routes, in
    /// source A/B order; these are not independent user-selected receivers.
    pub cells: [CellId; 2],
    pub heat_fraction_to_first: f64,
    pub response: R,
}

pub(super) struct PreparedContacts<'a, 'm> {
    bindings: Vec<&'a ConductionContactThermal<'m>>,
    cells: Vec<[CellId; 2]>,
    pressure_axes: BTreeMap<usize, String>,
}

// A refresh pins the selected claim, changing the selection receipt even at
// the same state. Compare its physical context and selected sources, retaining
// both original receipts rather than treating policy equality as physics.
fn same_interface(a: &ResolvedInterfaceStatePoint, b: &ResolvedInterfaceStatePoint) -> bool {
    a.card_identity() == b.card_identity()
        && a.surface_materials() == b.surface_materials()
        && a.axis_quantities() == b.axis_quantities()
        && a.query_point().len() == b.query_point().len()
        && a.query_point()
            .iter()
            .zip(b.query_point())
            .all(|((n, v), (m, w))| n == m && v.to_bits() == w.to_bits())
        && a.properties().len() == b.properties().len()
        && a.properties().iter().zip(b.properties()).all(|(p, q)| {
            p.requirement() == q.requirement()
                && p.source_axes() == q.source_axes()
                && p.answer().receipt.selected == q.answer().receipt.selected
                && p.value_si().to_bits() == q.value_si().to_bits()
        })
}

impl<'a, 'm> PreparedContacts<'a, 'm> {
    pub(super) fn new<C: PlanarConductorMaterial>(
        network: &ConductionThermalNetwork<'_>,
        problem: &PlanarConduction<C>,
        bindings: &'a [ConductionContactThermal<'m>],
        routes: &[CellId],
        accepted: &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error> {
        Self::new_for_contacts(
            network,
            problem.mesh(),
            problem.contacts(),
            problem.contact_liquid_mass_fraction_axes(),
            bindings,
            routes,
            accepted,
            cancel,
        )
    }

    // Linear and nonlinear current owners share original geometry/source/heat
    // admission without exposing a nonlinear law as a scalable frozen field.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new_for_contacts(
        network: &ConductionThermalNetwork<'_>,
        mesh: &fs_rep_mesh::TriComplex2,
        originals: &[ConductionContact],
        liquid_mass_fraction_axes: &[[Option<String>; 2]],
        bindings: &'a [ConductionContactThermal<'m>],
        routes: &[CellId],
        accepted: &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error> {
        check(cancel)?;
        if bindings.len() != originals.len() || liquid_mass_fraction_axes.len() != originals.len() {
            return Err(invalid(
                "every spatial electrical contact requires its source and heat partition",
            ));
        }
        if bindings.is_empty() {
            return Ok(Self {
                bindings: Vec::new(),
                cells: Vec::new(),
                pressure_axes: BTreeMap::new(),
            });
        }
        let trace = mesh
            .boundary_trace()
            .map_err(|_| invalid("spatial contact lost its original exterior traces"))?;
        let edge_key = |mut edge: [usize; 2]| {
            edge.sort_unstable();
            edge
        };
        let edges: BTreeMap<_, _> = trace
            .edges
            .iter()
            .map(|edge| {
                (
                    edge_key(edge.oriented_vertices.map(|v| v as usize)),
                    edge.source_face,
                )
            })
            .collect();
        let adjacent = originals
            .iter()
            .map(|contact| {
                check(cancel)?;
                let [first, second] = [contact.first_vertices, contact.second_vertices]
                    .map(|vertices| edges.get(&edge_key(vertices)).copied());
                match (first, second) {
                    (Some(first), Some(second)) => Ok([first, second]),
                    _ => Err(invalid("spatial contact has no original adjacent triangle")),
                }
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let sources = originals
            .iter()
            .map(|contact| &contact.source)
            .collect::<Vec<_>>();
        Self::new_for_sources(
            network,
            &sources,
            &adjacent,
            liquid_mass_fraction_axes,
            bindings,
            routes,
            accepted,
            &[],
            cancel,
        )
    }

    // Geometry owners supply their original adjacent volumes; source pinning,
    // caloric coordinates and two-sided heat ownership are identical in 2D/3D.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new_for_sources(
        network: &ConductionThermalNetwork<'_>,
        originals: &[&ResolvedElectricalInterface],
        adjacent_volumes: &[[usize; 2]],
        liquid_mass_fraction_axes: &[[Option<String>; 2]],
        bindings: &'a [ConductionContactThermal<'m>],
        routes: &[CellId],
        accepted: &MaterialNetworkSample,
        pressures: &[SolvedContactPressure<'_>],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error> {
        check(cancel)?;
        if bindings.len() != originals.len()
            || adjacent_volumes.len() != originals.len()
            || liquid_mass_fraction_axes.len() != originals.len()
        {
            return Err(invalid(
                "every spatial electrical contact requires its source, adjacent volumes and heat partition",
            ));
        }
        let mut ordered = vec![None; bindings.len()];
        for binding in bindings {
            check(cancel)?;
            if !binding.heat_fraction_to_first.is_finite()
                || !(0.0..=1.0).contains(&binding.heat_fraction_to_first)
                || binding.contact >= ordered.len()
                || ordered[binding.contact].replace(binding).is_some()
            {
                return Err(invalid(
                    "invalid, duplicate or foreign spatial contact heat binding",
                ));
            }
        }
        let bindings = ordered
            .into_iter()
            .map(|binding| {
                binding.ok_or_else(|| invalid("spatial contact thermal binding is absent"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut cells = Vec::with_capacity(bindings.len());
        for (adjacent, binding) in adjacent_volumes.iter().zip(&bindings) {
            check(cancel)?;
            if liquid_mass_fraction_axes[binding.contact]
                .each_ref()
                .map(|axis| axis.as_deref())
                != binding.point.liquid_fraction_axes
            {
                return Err(invalid(
                    "contact phase coordinates must retain the original native source ownership",
                ));
            }
            let mut receivers = [CellId(0); 2];
            for (side, face) in adjacent.iter().enumerate() {
                let cell = *routes.get(*face).ok_or_else(|| {
                    invalid("spatial contact adjacent volume lacks a heat receiver")
                })?;
                let bulk = network
                    .cells
                    .iter()
                    .find(|bulk| bulk.cell == cell)
                    .ok_or_else(|| {
                        invalid("spatial contact receiver lacks a bulk material binding")
                    })?;
                if bulk.point.material.content_hash()
                    != binding.point.materials[side].content_hash()
                {
                    return Err(invalid(
                        "spatial contact must retain both original bulk material cards",
                    ));
                }
                receivers[side] = cell;
            }
            cells.push(receivers);
        }
        let mut pressure_axes = BTreeMap::new();
        for pressure in pressures {
            check(cancel)?;
            let binding = bindings.get(pressure.contact).ok_or_else(|| {
                invalid("solved pressure requires an original electrical contact")
            })?;
            let axis = pressure.axis;
            if axis.trim().is_empty()
                || !binding.point.context.axes().contains_key(axis)
                || binding
                    .point
                    .context
                    .axis_quantities()
                    .get(axis)
                    .is_none_or(|quantity| {
                        quantity.dims() != Pressure::DIMS
                            || quantity.semantic_type().is_some_and(|semantic| {
                                semantic.kind() != QuantityKind::Pressure
                                    || semantic.form() != ValueForm::Static
                            })
                    })
                || binding.point.temperature_axes.contains(&axis)
                || binding.point.liquid_fraction_axes.contains(&Some(axis))
                || binding.point.mean_temperature_axis == Some(axis)
                || pressure_axes
                    .insert(pressure.contact, axis.to_owned())
                    .is_some()
            {
                return Err(invalid(
                    "solved contact pressure requires a unique distinct original static pressure coordinate in pascals",
                ));
            }
        }
        let prepared = Self {
            bindings,
            cells,
            pressure_axes,
        };
        let refreshed =
            prepared.refresh_sources_with_pressures(originals, accepted, pressures, cancel)?;
        if refreshed
            .iter()
            .zip(originals)
            .any(|(new, old)| !same_interface(new.resolved(), old.resolved()))
        {
            return Err(invalid(
                "accepted electrical contacts must match both actual thermal states and original source context",
            ));
        }
        Ok(prepared)
    }

    pub(super) fn refresh<C: PlanarConductorMaterial>(
        &self,
        problem: &PlanarConduction<C>,
        sample: &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<ConductionContact>, Error> {
        self.refresh_contacts(problem.contacts(), sample, cancel)
    }

    pub(super) fn refresh_contacts(
        &self,
        originals: &[ConductionContact],
        sample: &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<ConductionContact>, Error> {
        let originals_sources = originals
            .iter()
            .map(|contact| &contact.source)
            .collect::<Vec<_>>();
        let refreshed = self.refresh_sources(&originals_sources, sample, cancel)?;
        Ok(originals
            .iter()
            .zip(refreshed)
            .map(|(original, source)| ConductionContact {
                first_vertices: original.first_vertices,
                second_vertices: original.second_vertices,
                source,
            })
            .collect())
    }

    pub(super) fn refresh_sources(
        &self,
        originals: &[&ResolvedElectricalInterface],
        sample: &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<ResolvedElectricalInterface>, Error> {
        self.refresh_sources_with_pressures(originals, sample, &[], cancel)
    }

    pub(super) fn refresh_sources_with_pressures(
        &self,
        originals: &[&ResolvedElectricalInterface],
        sample: &MaterialNetworkSample,
        pressures: &[SolvedContactPressure<'_>],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<ResolvedElectricalInterface>, Error> {
        if originals.len() != self.bindings.len() {
            return Err(invalid(
                "spatial contact source coverage changed during feedback",
            ));
        }
        // A trial may update values, never the original ownership or coordinate.
        // Keeping this map beside the heat routes makes a missing/duplicate or
        // renamed mechanical coordinate fail before any source is published.
        let mut pressure_values = vec![None; self.bindings.len()];
        if pressures.len() != self.pressure_axes.len() {
            return Err(invalid(
                "solved contact pressure coverage changed during feedback",
            ));
        }
        for pressure in pressures {
            check(cancel)?;
            if !pressure.pressure_pa.is_finite()
                || pressure.pressure_pa < 0.0
                || self
                    .pressure_axes
                    .get(&pressure.contact)
                    .map(String::as_str)
                    != Some(pressure.axis)
                || pressure_values[pressure.contact]
                    .replace(pressure.pressure_pa)
                    .is_some()
            {
                return Err(invalid(
                    "solved pressure requires finite compression at each original contact coordinate",
                ));
            }
        }
        let mut contacts = Vec::with_capacity(self.bindings.len());
        for ((binding, cells), original) in self.bindings.iter().zip(&self.cells).zip(originals) {
            check(cancel)?;
            let phases = cells.map(|cell| {
                sample
                    .cells
                    .iter()
                    .find(|entry| entry.cell == cell)
                    .map(|entry| entry.phase)
            });
            let [Some(first), Some(second)] = phases else {
                return Err(invalid("spatial contact lost a receiving material state"));
            };
            let mut point = binding
                .point
                .at([first, second])
                .map_err(|error| material(MaterialNetworkError::Material(error)))?;
            if let Some(pressure_pa) = pressure_values[binding.contact] {
                let axis = &self.pressure_axes[&binding.contact];
                point = point
                    .with(axis, pressure_pa)
                    .map_err(|_| invalid("original mechanical pressure coordinate refused"))?;
            }
            let source = original
                .at_query_point(binding.point.interface, &point)
                .map_err(|error| {
                    material(MaterialNetworkError::Material(
                        ThermalMaterialError::Material(error),
                    ))
                })?;
            let property = source
                .resolved()
                .property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
                .ok_or_else(|| invalid("spatial electrical contact lost its source property"))?;
            if !property
                .source_axes()
                .iter()
                .all(|(axis, dims)| binding.point.is_response_axis(axis, *dims))
            {
                return Err(invalid(
                    "contact response requires owned thermal/phase or typed fixed external coordinates",
                ));
            }
            contacts.push(source);
        }
        Ok(contacts)
    }

    pub(super) fn heat<R: Clone>(
        &self,
        response: &[R],
    ) -> Result<Vec<ConductionContactHeat<R>>, Error> {
        if response.len() != self.bindings.len() {
            return Err(invalid(
                "spatial contact heat coverage changed during feedback",
            ));
        }
        Ok(self
            .bindings
            .iter()
            .zip(&self.cells)
            .zip(response)
            .map(|((binding, cells), response)| ConductionContactHeat {
                contact: binding.contact,
                cells: *cells,
                heat_fraction_to_first: binding.heat_fraction_to_first,
                response: response.clone(),
            })
            .collect())
    }
}

pub(super) fn deposit(
    contacts: &[ConductionContactHeat],
    powers: &mut BTreeMap<CellId, f64>,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), Error> {
    for contact in contacts {
        checkpoint()?;
        deposit_partition(
            contact.cells,
            contact.heat_fraction_to_first,
            contact.response.power_w,
            powers,
        )?;
    }
    Ok(())
}

pub(super) fn deposit_partition(
    cells: [CellId; 2],
    fraction: f64,
    power: f64,
    powers: &mut BTreeMap<CellId, f64>,
) -> Result<(), Error> {
    let first = finite(power * fraction)?;
    let second = finite(power - first)?;
    if power < 0.0
        || first < 0.0
        || second < 0.0
        || (power > 0.0 && ((fraction > 0.0 && first == 0.0) || (fraction < 1.0 && second == 0.0)))
    {
        return Err(invalid("spatial contact heat partition is unrepresentable"));
    }
    for (cell, power) in cells.into_iter().zip([first, second]) {
        let previous = powers.get(&cell).copied().unwrap_or(0.0);
        let combined = finite(previous + power)?;
        if power != 0.0 && previous == combined {
            return Err(invalid(
                "spatial contact heat disappeared in receiving cell power",
            ));
        }
        powers.insert(cell, combined);
    }
    Ok(())
}

impl ConductionThermalNetwork<'_> {
    /// Advance bulk current, surface voltage jumps, two-sided contact heat and
    /// every finite-cell enthalpy through the original nonlinear residual.
    /// Interface temperature/phase axes and fixed pressure/history come from
    /// each explicit `ThermalInterfacePoint`; original source choices stay pinned.
    /// Fraction axes must match the native `new_with_phase_contacts` ownership;
    /// their exact values come from the original caloric charts, even when
    /// temperatures remain constant throughout a latent-heat plateau.
    /// All bulk and surface Joule heat enters the receiving bodies exactly once.
    /// The accepted field and thermal history must be published together.
    ///
    /// Matching stationary traces, constant volume/density and zero interface
    /// storage retain the native owners' limitations. Heat fractions are fixed
    /// physical inputs for the interval; no spreading, opening or arcing law is
    /// supplied. First-order endpoint coupling needs time/mesh refinement.
    #[allow(clippy::too_many_arguments)]
    pub fn step_contact_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConduction<C>,
        drive: ConductionDrive<'_>,
        contacts: &[ConductionContactThermal<'_>],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionEnthalpyCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        self.admit_mesh_step(
            thermal,
            target_s,
            problem.mesh(),
            assembly,
            controls,
            initial_energy_increments_j,
        )?;
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive,
            contacts,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            |state| {
                prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |source| source,
            &cancel,
        )
    }
}
