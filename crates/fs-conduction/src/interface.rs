//! Finite-resistance contacts with separate temperature traces.
//!
//! [`ThermalInterfaces::new`] retains the original exact matching-P1 path.
//! [`ThermalInterfaces::with_nonmatching`] additionally binds explicitly named,
//! completely covered planar traces through common-refinement integration.
//! Both feed the same steady, nonlinear, transient and adjoint assemblers.
//! No perfect contact, face ownership, material data or gap closure is inferred.
mod matching;
pub mod nonmatching;
/// Source-backed temperature-dependent matching-contact advancement.
pub mod temperature;
/// Source-backed contact-area exchange on shared finite thermal inputs.
pub mod lumped_pair;
mod sensitivity;
pub use sensitivity::ContactResistanceGradient;

pub use matching::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    InterfaceFacePair, InterfaceFlux, InterfaceResistance, InterfaceSurface,
    ResistanceOrigin, ResistanceUncertainty, ResistanceValueOrigin,
    SeriesResistanceBudget, SeriesThermalResistance, SlabConductivityQuery, ThermalResistanceTerm,
};
pub use nonmatching::{NonmatchingOptions, NonmatchingSurface};

use std::collections::{BTreeMap, BTreeSet};
use fs_exec::Cx;
use fs_sparse::Coo;
use crate::{ConductionError, ConductionMesh, ThermalBoundary};

/// Complete declared contact operator. Matching contributions keep their
/// original integration/evaluation order. Nonmatching additions carry their
/// own checked geometry and explicit floating-point admission tolerances.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalInterfaces {
    matching: matching::ThermalInterfaces,
    // Some nonlinear nonmatching patches own even their exact subpatches.
    // The original matching binding remains the geometry/ownership witness.
    integration_matching: Option<matching::ThermalInterfaces>,
    nonmatching: Vec<nonmatching::Bound>,
}

impl ThermalInterfaces {
    /// Bind the original exact-coordinate matching interface model unchanged.
    pub fn new(mesh: &ConductionMesh, boundary: &ThermalBoundary,
        surfaces: Vec<InterfaceSurface>) -> Result<Self, ConductionError> {
        Ok(Self { matching: matching::ThermalInterfaces::new(mesh,boundary,surfaces)?,
            integration_matching: None, nonmatching: Vec::new() })
    }

    /// Bind matching and explicitly supplied nonmatching planar interfaces.
    /// Each declared face has one owner. Each nonmatching face must be fully
    /// covered within its explicit area tolerance, with no same-side overlap.
    /// All exact coincident candidates still require a declaration, as before.
    /// This does not search the whole mesh for UNDECLARED nonmatching contacts.
    ///
    /// Exact triangle pairs within a nonmatching patch are delegated to the
    /// original matching producer, not assembled twice. The remaining overlaps
    /// share the same name, original A-face resistance and signed heat report.
    ///
    /// # Errors
    /// All matching refusals, plus bounded planar-intersection/coverage,
    /// allocation-domain, shared ownership, or cancellation refusals.
    pub fn with_nonmatching(cx: &Cx<'_>, mesh: &ConductionMesh, boundary: &ThermalBoundary,
        mut matching: Vec<InterfaceSurface>, mut surfaces: Vec<NonmatchingSurface>) -> Result<Self, ConductionError> {
        cx.checkpoint().map_err(|_| ConductionError::Cancelled {stage:"nonmatching-contact",at:0})?;
        let candidates = Self::coincident_face_pairs(mesh)?;
        let exact: BTreeSet<_> = candidates.iter().map(|p|
            (p.side_a.min(p.side_b),p.side_a.max(p.side_b))).collect();
        let mut names = BTreeSet::new();
        let mut owned = BTreeSet::new();
        for surface in &matching {
            if !names.insert(surface.name().to_string()) {return Err(invalid("duplicate contact name"));}
            for pair in surface.face_pairs() {
                for slot in [pair.side_a,pair.side_b] {
                    if !owned.insert(slot) {return Err(invalid("a contact face has more than one owner"));}
                }
            }
        }
        surfaces.sort_by(|a,b|a.name.cmp(&b.name));
        let mut nonmatching = Vec::with_capacity(surfaces.len());
        for surface in surfaces {
            if !names.insert(surface.name.clone()) {return Err(invalid("duplicate contact name"));}
            for &slot in surface.side_a.iter().chain(&surface.side_b) {
                if !owned.insert(slot) {return Err(invalid("a contact face has more than one owner"));}
            }
            let bound = nonmatching::Bound::build(cx,mesh,boundary,&surface,&exact)?;
            let mut pairs = Vec::new();
            for pair in &candidates {
                if surface.side_a.binary_search(&pair.side_a).is_ok()
                    && surface.side_b.binary_search(&pair.side_b).is_ok() {
                    pairs.push(*pair);
                } else if surface.side_a.binary_search(&pair.side_b).is_ok()
                    && surface.side_b.binary_search(&pair.side_a).is_ok() {
                    pairs.push(InterfaceFacePair {side_a:pair.side_b,side_b:pair.side_a});
                }
            }
            if !pairs.is_empty() {
                matching.push(if surface.is_mapped() {
                    let resistances=pairs.iter().map(|pair|surface.resistance_for_face(pair.side_a).cloned()
                        .ok_or_else(||invalid("exact contact delegation lost its original A-face resistance")))
                        .collect::<Result<Vec<_>,_>>()?;
                    InterfaceSurface::new_mapped(surface.name,pairs,resistances)?
                } else {
                    let resistance=surface.resistance_for_face(pairs[0].side_a).cloned()
                        .ok_or_else(||invalid("exact contact delegation lost its original A-face resistance"))?;
                    InterfaceSurface::new(surface.name,pairs,resistance)?
                });
            }
            nonmatching.push(bound);
        }
        let matching = matching::ThermalInterfaces::new(mesh,boundary,matching)?;
        Ok(Self { matching, integration_matching: None, nonmatching })
    }

    fn matching_integration(&self) -> &matching::ThermalInterfaces {
        self.integration_matching.as_ref().unwrap_or(&self.matching)
    }

    fn refresh_matching(&mut self, surfaces: &[InterfaceSurface],
        quadrature_owned: &BTreeSet<String>) -> Result<(), ConductionError> {
        self.matching.update_surfaces(surfaces)?;
        self.integration_matching = (!quadrature_owned.is_empty())
            .then(|| self.matching.excluding(quadrature_owned));
        Ok(())
    }

    /// Exact matching candidates only. No proximity or nonmatching face search
    /// is silently introduced into legacy/native project admission.
    pub fn coincident_face_pairs(mesh:&ConductionMesh)->Result<Vec<InterfaceFacePair>,ConductionError> {
        matching::ThermalInterfaces::coincident_face_pairs(mesh)
    }
    pub(crate) fn require_no_undeclared(mesh:&ConductionMesh)->Result<(),ConductionError> {
        matching::ThermalInterfaces::require_no_undeclared(mesh)
    }
    pub(crate) fn validate_for(&self,mesh:&ConductionMesh,boundary:&ThermalBoundary)->Result<(),ConductionError> {
        self.matching.validate_for(mesh,boundary)?;
        for surface in &self.nonmatching {surface.validate_for(mesh,boundary)?;}
        Ok(())
    }
    pub(crate) fn assemble_into(&self,cx:&Cx<'_>,coo:&mut Coo)->Result<(),ConductionError> {
        self.matching_integration().assemble_into(cx,coo)?;
        for surface in &self.nonmatching {surface.assemble_into(cx,coo)?;}
        Ok(())
    }

    /// Number of named interfaces, not overlap triangles or delegated pieces.
    #[must_use]
    pub fn surface_count(&self)->usize {
        self.matching.surface_count()+self.nonmatching.iter()
            .filter(|s|self.matching.surface_is_mapped(&s.name).is_none()).count()
    }
    /// Whether a named interface carries a spatial resistance map.
    #[must_use]
    pub fn surface_is_mapped(&self,name:&str)->Option<bool> {
        if let Some(surface)=self.nonmatching.iter().find(|s|s.name==name) {Some(surface.is_mapped())}
        else {self.matching.surface_is_mapped(name)}
    }
    /// Matching: per-face-pair values in the original order. Nonmatching:
    /// resistance inherited from the original A face for each common-refinement
    /// integration triangle, including exactly matched delegated pieces.
    #[must_use]
    pub fn surface_face_resistances(&self,name:&str)->Option<Vec<f64>> {
        self.nonmatching.iter().find(|s|s.name==name).map(nonmatching::Bound::resistances)
            .or_else(||self.matching.surface_face_resistances(name))
    }
    /// Original source-backed values of a mapped nonmatching interface, keyed
    /// by canonical A-face slot. Uniform and matching-only surfaces return None.
    /// Each retained receipt is its actual source query; an aggregate heat
    /// report's representative receipt is not a spatially averaged material claim.
    #[must_use]
    pub fn nonmatching_face_resistance_sources(&self,name:&str)->Option<&[(usize,InterfaceResistance)]> {
        self.nonmatching.iter().find(|s|s.name==name).and_then(nonmatching::Bound::resistance_sources)
    }
    /// Whether a named surface uses explicitly admitted planar overlap traces.
    #[must_use]
    pub fn surface_is_nonmatching(&self,name:&str)->bool {
        self.nonmatching.iter().any(|s|s.name==name)
    }

    /// One signed A-to-B heat report per named interface. Matching-only calls
    /// return the original report directly, retaining legacy arithmetic order.
    /// A mapped nonmatching report retains the first canonical A-face receipt
    /// as a representative; inspect `nonmatching_face_resistance_sources` for
    /// every original local query. Conductance is the integrated local 1/R''.
    pub fn fluxes(&self,temperature:&[f64])->Result<Vec<InterfaceFlux>,ConductionError> {
        let original=self.matching_integration().fluxes(temperature)?;
        if self.nonmatching.is_empty(){return Ok(original);}
        let mut fluxes:BTreeMap<String,InterfaceFlux>=original.into_iter().map(|f|(f.interface.clone(),f)).collect();
        for surface in &self.nonmatching {
            if let Some(extra)=surface.flux(temperature)? {
                if let Some(total)=fluxes.get_mut(&extra.interface) {
                    total.area_m2=finite(total.area_m2+extra.area_m2)?;
                    total.conductance_w_per_k=finite(total.conductance_w_per_k+extra.conductance_w_per_k)?;
                    total.heat_rate_a_to_b_w=finite(total.heat_rate_a_to_b_w+extra.heat_rate_a_to_b_w)?;
                    total.mean_jump_k=finite(total.heat_rate_a_to_b_w/total.conductance_w_per_k)?;
                } else {fluxes.insert(extra.interface.clone(),extra);}
            }
            if surface.is_mapped() {
                if let Some(total)=fluxes.get_mut(&surface.name) {
                    total.receipt=surface.representative_resistance().receipt().clone();
                }
            }
        }
        Ok(fluxes.into_values().collect())
    }

    /// Contract a TOTAL coupled nodal-load adjoint with dK/dln(R'') for a
    /// nonmatching surface. Returns None for a matching-only/unknown name so
    /// existing matching sensitivity producers remain unchanged. Geometry,
    /// overlap topology and contact resistance law are held fixed.
    pub fn nonmatching_log_resistance_pullback(&self,cx:&Cx<'_>,name:&str,
        temperature:&[f64],nodal_load_adjoint:&[f64])->Result<Option<f64>,ConductionError> {
        self.nonmatching.iter().find(|s|s.name==name)
            .map(|s|s.log_resistance_pullback(cx,temperature,nodal_load_adjoint)).transpose()
    }
}
fn invalid(what:impl Into<String>)->ConductionError {
    ConductionError::Interface{interface:"<nonmatching-binding>".into(),what:what.into(),
        fix:"declare each contact face and name exactly once; preserve independent traces".into()}
}
fn finite(value:f64)->Result<f64,ConductionError> {
    if value.is_finite(){Ok(value)}else{Err(invalid("nonfinite combined contact flux"))}
}
