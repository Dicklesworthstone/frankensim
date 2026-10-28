//! Stationary species balances on the existing prescribed upwind graph.
//!
//! Mass-balanced graphs need explicit inlet compositions. Acyclic sections use
//! direct mixing; inlet-fed recirculation uses the existing bounded dense solver.
//! Closed components with no inlet have no uniquely determined composition and
//! refuse. Species have no reaction, diffusion, condensation or storage rate.
//! The enthalpy observer uses these SAME directions and mass rates. This
//! stationary extension does not turn fixed thermal cells into a variable-
//! density transient flow solver.

use std::collections::{BTreeMap, BTreeSet};

use super::{Advected, MassBoundary, MassFlow, finite, invalid};
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation};

#[path = "species/recirculation.rs"]
mod recirculation;

/// Maximum unresolved cells in the dense recirculation solve, including any
/// downstream cells held back by a cycle. Acyclic graphs keep their 4096-cell cap.
pub const MAX_RECIRCULATION_CELLS: usize = 128;
/// Largest admitted mass-fraction error estimate from the scaled residual and
/// inverse-norm/roundoff check. This numerical allowance is not material data
/// uncertainty or a formal interval certificate.
pub const RECIRCULATION_FRACTION_TOLERANCE: f64 = 1e-10;

/// Frozen caloric-law composition in one explicit species/reference convention.
/// A common reference identifies the energy datum, not identical h(T) curves.
#[derive(Clone, Copy, Debug)]
pub struct CaloricComposition<'a> {
    /// Canonical species identities, including species absent from this mixture.
    pub species: &'a [String],
    /// Mass fractions in that exact order.
    pub mass_fractions: &'a [f64],
    /// Shared enthalpy reference including the constitutive convention.
    pub enthalpy_reference: &'a str,
}

/// Species carried by one already-declared inlet mass connection.
#[derive(Clone, Debug, PartialEq)]
pub struct InletSpecies {
    /// The MassBoundary::Inlet connection identity.
    pub boundary_id: u64,
    /// Authored mass fractions, required to sum to exact one.
    pub mass_fractions: Vec<f64>,
}

/// Accepted stationary compositions and the actual species-flux residual.
#[derive(Clone, Debug, PartialEq)]
pub struct StationarySpecies {
    species: Vec<String>,
    cells: BTreeMap<CellId, Vec<f64>>,
    inlets: BTreeMap<u64, Vec<f64>>,
    maximum_residual_kg_s: f64,
}

fn authored(fractions: &[f64], count: usize) -> Result<(), DepositionError> {
    if fractions.len() != count
        || fractions.iter().any(|x| !x.is_finite() || *x < 0.0 || *x > 1.0)
        || fractions.iter().sum::<f64>() != 1.0
    {
        return Err(invalid("stationary species need complete nonnegative unit mass fractions"));
    }
    Ok(())
}

/// Close only roundoff in a GENERATED simplex. Authored compositions never
/// use this operation. No concentration floor or finite mass source is added.
pub fn close_generated_simplex(fractions: &mut [f64]) -> Result<(), DepositionError> {
    if fractions.is_empty() || fractions.len() > 128
        || fractions.iter().any(|x| !x.is_finite() || *x < 0.0 || *x > 1.0)
    {
        return Err(invalid("invalid generated species simplex"));
    }
    let largest = fractions.iter().enumerate()
        .max_by(|(ia, a), (ib, b)| a.total_cmp(b).then_with(|| ib.cmp(ia)))
        .map(|(index, _)| index).expect("nonempty");
    for _ in 0..4 {
        let sum: f64 = fractions.iter().sum();
        if sum == 1.0 { return Ok(()); }
        let residual = 1.0 - sum;
        if !sum.is_finite() || residual.abs() > 32.0 * f64::EPSILON * fractions.len() as f64 {
            return Err(invalid("generated species simplex exceeds roundoff"));
        }
        fractions[largest] += residual;
        if !(0.0..=1.0).contains(&fractions[largest]) {
            return Err(invalid("generated species closure left the simplex"));
        }
    }
    // Adding a rounded total residual to the largest entry can alternate
    // between sums immediately below and above one. Close a component against
    // its actual prefix instead. A suffix is retained for a tiny final species
    // that cannot absorb a unit-scale rounding correction without distortion.
    // Only relative roundoff changes are allowed, and the actual ordered sum
    // remains the admission gate; no positive trace or exact zero is erased.
    let tolerance = 32.0 * f64::EPSILON * fractions.len() as f64;
    for index in (0..fractions.len()).rev() {
        let old = fractions[index];
        if old == 0.0 { continue; }
        let prefix: f64 = fractions[..index].iter().sum();
        let suffix: f64 = fractions[index + 1..].iter().sum();
        let closed = (1.0 - suffix) - prefix;
        if closed.is_finite() && closed > 0.0 && closed <= 1.0
            && (closed - old).abs() <= tolerance * old
        {
            fractions[index] = closed;
            if fractions.iter().sum::<f64>() == 1.0 { return Ok(()); }
            fractions[index] = old;
        }
    }
    Err(invalid("generated species simplex cannot represent unit sum"))
}

fn note_uniform(
    inputs: &mut BTreeMap<CellId, Option<Vec<f64>>>, cell: CellId, fractions: &[f64],
) {
    match inputs.entry(cell) {
        std::collections::btree_map::Entry::Vacant(entry) => {
            entry.insert(Some(fractions.to_vec()));
        }
        std::collections::btree_map::Entry::Occupied(mut entry) => {
            if entry.get().as_ref().is_some_and(|old| old.as_slice() != fractions) {
                entry.insert(None);
            }
        }
    }
}

impl StationarySpecies {
    fn solve(
        flows: &[MassFlow], boundaries: &[MassBoundary], species: Vec<String>,
        inlet_data: Vec<InletSpecies>, max_cells: usize, max_species: usize,
    ) -> Result<Self, DepositionError> {
        if species.is_empty() || species.len() > max_species.min(128)
            || inlet_data.len() > boundaries.len()
        {
            return Err(DepositionError::Budget("stationary species layout"));
        }
        let mut names = BTreeSet::new();
        if species.iter().any(|name| name.is_empty() || name.len() > 128 || !names.insert(name)) {
            return Err(invalid("stationary species identities must be bounded and unique"));
        }
        let count = species.len();
        let mut inlets = BTreeMap::new();
        for inlet in inlet_data {
            authored(&inlet.mass_fractions, count)?;
            if inlets.insert(inlet.boundary_id, inlet.mass_fractions).is_some() {
                return Err(invalid("duplicate stationary species inlet"));
            }
        }
        let expected: BTreeSet<_> = boundaries.iter().filter_map(|boundary| match boundary {
            MassBoundary::Inlet { id, .. } => Some(*id),
            MassBoundary::Outlet { .. } => None,
        }).collect();
        if expected != inlets.keys().copied().collect::<BTreeSet<_>>() {
            return Err(invalid("stationary species need exactly every declared inlet"));
        }
        let mut rates: BTreeMap<CellId, [f64; 2]> = BTreeMap::new();
        let mut edges: BTreeMap<CellId, Vec<(CellId, f64)>> = BTreeMap::new();
        let mut degree: BTreeMap<CellId, usize> = BTreeMap::new();
        let mut incoming_species: BTreeMap<CellId, Vec<f64>> = BTreeMap::new();
        let mut uniform_inputs: BTreeMap<CellId, Option<Vec<f64>>> = BTreeMap::new();
        for flow in flows {
            let Some((up, down, rate)) = flow.directed() else { continue; };
            let out = rates.entry(up).or_default();
            out[1] = finite(out[1] + rate)?;
            let incoming = rates.entry(down).or_default();
            incoming[0] = finite(incoming[0] + rate)?;
            edges.entry(up).or_default().push((down, rate));
            degree.entry(up).or_default();
            let n = degree.entry(down).or_default();
            *n = n.checked_add(1).ok_or(DepositionError::Arithmetic)?;
        }
        for boundary in boundaries {
            let (_, cell, rate) = boundary.identity();
            degree.entry(cell).or_default();
            let sums = rates.entry(cell).or_default();
            match boundary {
                MassBoundary::Inlet { id, .. } => {
                    sums[0] = finite(sums[0] + rate)?;
                    note_uniform(&mut uniform_inputs, cell, &inlets[id]);
                    let fluxes = incoming_species.entry(cell).or_insert_with(|| vec![0.0; count]);
                    for (flux, fraction) in fluxes.iter_mut().zip(&inlets[id]) {
                        *flux = finite(*flux + finite(rate * fraction)?)?;
                    }
                }
                MassBoundary::Outlet { .. } => { sums[1] = finite(sums[1] + rate)?; }
            }
        }
        if rates.is_empty() || rates.len() > max_cells.min(4096) {
            return Err(DepositionError::Budget("stationary species cells"));
        }
        for &[incoming, outgoing] in rates.values() {
            if incoming <= 0.0 || outgoing <= 0.0
                || (incoming - outgoing).abs() > 64.0 * f64::EPSILON * incoming.max(outgoing)
            {
                return Err(invalid("stationary species require balanced positive through-flow"));
            }
        }
        let mut ready: BTreeSet<_> = degree.iter().filter_map(|(&cell, &n)|
            if n == 0 { Some(cell) } else { None }).collect();
        let mut cells = BTreeMap::new();
        let mut maximum_residual_kg_s: f64 = 0.0;
        while let Some(cell) = ready.pop_first() {
            let fluxes = incoming_species.remove(&cell)
                .ok_or(invalid("stationary species component has no inlet"))?;
            let [incoming, outgoing] = rates[&cell];
            let fractions = if let Some(Some(exact)) = uniform_inputs.remove(&cell) {
                exact
            } else {
                let mut fractions: Vec<_> = fluxes.iter().map(|flux| flux / incoming).collect();
                close_generated_simplex(&mut fractions)?;
                fractions
            };
            for (&flux, &fraction) in fluxes.iter().zip(&fractions) {
                let residual = finite(flux - finite(outgoing * fraction)?)?.abs();
                if residual > 128.0 * f64::EPSILON * count as f64 * incoming.max(outgoing) {
                    return Err(invalid("stationary species flux does not close"));
                }
                maximum_residual_kg_s = maximum_residual_kg_s.max(residual);
            }
            if let Some(downstream) = edges.get(&cell) {
                for &(down, rate) in downstream {
                    note_uniform(&mut uniform_inputs, down, &fractions);
                    let fluxes = incoming_species.entry(down).or_insert_with(|| vec![0.0; count]);
                    for (flux, fraction) in fluxes.iter_mut().zip(&fractions) {
                        *flux = finite(*flux + finite(rate * fraction)?)?;
                    }
                    let n = degree.get_mut(&down).expect("registered downstream cell");
                    *n -= 1;
                    if *n == 0 { ready.insert(down); }
                }
            }
            cells.insert(cell, fractions);
        }
        if cells.len() != rates.len() {
            maximum_residual_kg_s = maximum_residual_kg_s.max(recirculation::complete(
                &rates, &edges, &incoming_species, uniform_inputs, &mut cells, count,
            )?);
        }
        Ok(Self { species, cells, inlets, maximum_residual_kg_s })
    }

    /// Canonical species order retained without renaming.
    pub fn species(&self) -> &[String] { &self.species }
    /// Frozen cell mass fractions obtained from the actual inlet/connection rates.
    pub fn cells(&self) -> &BTreeMap<CellId, Vec<f64>> { &self.cells }
    /// Actual largest absolute component mass-flow defect [kg/s].
    pub fn maximum_residual_kg_s(&self) -> f64 { self.maximum_residual_kg_s }

    pub(super) fn validate_laws(
        &self, thermal: &DepositionSimulation<'_>, boundaries: &[MassBoundary],
    ) -> Result<(), DepositionError> {
        let mut reference: Option<String> = None;
        let mut admit = |law: &dyn crate::enthalpy::ThermalLaw, fractions: &[f64]| {
            let descriptor = law.advected_composition()
                .ok_or(invalid("stationary species require composition-aware caloric laws"))?;
            if descriptor.species != self.species.as_slice() || descriptor.mass_fractions != fractions
                || descriptor.enthalpy_reference.is_empty()
            {
                return Err(invalid("stationary species do not match the receiving caloric law"));
            }
            if let Some(expected) = &reference {
                if descriptor.enthalpy_reference != expected.as_str() {
                    return Err(invalid("stationary species caloric laws have different enthalpy references"));
                }
            } else {
                reference = Some(descriptor.enthalpy_reference.to_owned());
            }
            Ok(())
        };
        for (&id, fractions) in &self.cells {
            let cell = thermal.cells().get(&id).ok_or(invalid("inactive stationary species cell"))?;
            admit(thermal.materials[&cell.material_id], fractions)?;
        }
        for boundary in boundaries {
            if let MassBoundary::Inlet { id, material_id, .. } = boundary {
                let law = thermal.materials.get(material_id)
                    .ok_or(invalid("unregistered stationary inlet caloric law"))?;
                admit(*law, &self.inlets[id])?;
            }
        }
        Ok(())
    }
}

impl<O> Advected<O> {
    /// Solve stationary species on this observer's existing mass-flow graph.
    /// Inlet-fed cycles use at most [`MAX_RECIRCULATION_CELLS`] unresolved cells;
    /// unreachable or ill-conditioned recirculation refuses without regularizing
    /// the flow or inventing a composition. The acyclic path is unchanged.
    /// Allows different immutable cell caloric laws only after their exact
    /// compositions and shared energy reference are checked on every advance.
    /// Inlet enthalpy remains the INLET law's enthalpy; it is never interpreted
    /// as the receiving mixture's temperature or a transferable temperature.
    pub fn with_stationary_species(
        mut self, species: Vec<String>, inlets: Vec<InletSpecies>,
        max_cells: usize, max_species: usize,
    ) -> Result<Self, DepositionError> {
        self.stationary_species = Some(StationarySpecies::solve(
            &self.flows, &self.boundaries, species, inlets, max_cells, max_species,
        )?);
        self.pending = None;
        self.boundary_pending = None;
        Ok(self)
    }

    /// The admitted species solution; absent on the unchanged original path.
    pub fn stationary_species(&self) -> Option<&StationarySpecies> {
        self.stationary_species.as_ref()
    }
}

#[cfg(test)]
#[path = "species/tests.rs"]
mod tests;
