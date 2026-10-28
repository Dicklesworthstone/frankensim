//! One exterior boundary for a pressure/volume-flow network, exp(-i omega t).
//!
//! Unit outward opening flows become normal-velocity traces on the original
//! capped exterior mesh. Their work-conjugate pressure averages form the FULL
//! termination matrix. The original BEM solutions also produce microphone
//! fields, weighted by the flows solved with that matrix. Mutual loading and
//! coherent radiation therefore cannot accidentally describe different sources.
//!
//! This is a harmonic, rigid-wall coupling with one exterior atmosphere and
//! optional source-bound regional interior gases, not a causal time realization.
//! Openings retain one uniform-flow coordinate each. The caller
//! supplies physical caps and their correspondence to duct reference planes;
//! this does not infer that geometry or resolve interior evanescent modes.
//! Mesh/port-basis refinement and measured validation remain separate. Raw BEM
//! matrices are never symmetrized, diagonalized, fitted or repaired for passivity.

use super::*;
use fs_bem::{helmholtz::{self, Formulation, HelmholtzError, Medium, laplace::{self, PointField}},
    near_field::Geometry, panel3d::SpherePanels};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_duct::network::NetworkMedia;

/// One actual opening on a closed, outward-oriented, capped surface. All other
/// panels are rigid. Panels cannot belong to two different flow coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpeningPatch {
    /// Exact coupled-port label; patch order must equal the declared matrix order.
    pub port: String,
    /// Original surface panel indices, without new tessellation or area gains.
    pub panels: Vec<usize>,
}

/// Whole-call work and storage bounds, checked before allocating source traces.
#[derive(Clone, Copy, Debug)]
pub struct ExteriorNetworkBudget {
    /// Dense exterior panel cap, additionally bounded by the original BEM owner.
    pub max_panels: usize,
    /// Three source-panel arrays (trace, pressure, velocity) plus four receiver components per source.
    pub max_field_entries: usize,
    /// Conservative logical allowance n^3 + 3*p*n^2 + 4*p*n*r. Not a wall-time
    /// or quadrature cost certificate; BEM retains its own per-solve work limits.
    pub max_solve_terms: usize,
}

/// Refusals retain the original numerical owner. No partial response publishes.
#[derive(Debug)]
pub enum ExteriorNetworkError {
    /// Incomplete physical binding or exceeded declared work limit.
    Invalid(&'static str),
    /// Original duct propagation, boundary or Schur-solve refusal.
    Duct(DuctError),
    /// Original surface geometry, radiation or receiver refusal.
    Bem(HelmholtzError),
    /// Requested cancellation, without a partial response.
    Cancelled,
}
impl From<DuctError> for ExteriorNetworkError { fn from(e:DuctError)->Self { Self::Duct(e) } }
impl From<HelmholtzError> for ExteriorNetworkError { fn from(e:HelmholtzError)->Self { Self::Bem(e) } }
impl core::fmt::Display for ExteriorNetworkError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result { write!(f,"exterior network: {self:?}") }
}
impl std::error::Error for ExteriorNetworkError {}
fn bad(what:&'static str)->ExteriorNetworkError { ExteriorNetworkError::Invalid(what) }
fn poll(gate:&CancelGate)->Result<(),ExteriorNetworkError> {
    if gate.is_requested() { Err(ExteriorNetworkError::Cancelled) } else { Ok(()) }
}
fn finite(z:C64)->bool { z.re.is_finite() && z.im.is_finite() }

/// Immutable physical network, exterior geometry, opening basis and medium.
/// There is no cached frequency table or mutable source/receiver substitution.
pub struct RadiatingNetwork {
    driving: NetworkDrivingPoint,
    surface: SpherePanels,
    patches: Vec<OpeningPatch>,
    traces: Vec<Vec<C64>>,
    areas: Vec<f64>,
    gas: GasState,
    media: Option<NetworkMedia>,
    receivers: Vec<[f64;3]>,
    clearance_m: f64,
    source: ContentHash,
    boundary_hash: ContentHash,
    budget: ExteriorNetworkBudget,
}

/// One harmonic solution. Response arrays follow network node/port order;
/// outward flows follow `patches()`. Receiver fields follow `receivers()`.
#[derive(Debug)]
pub struct RadiatingNetworkSample {
    /// Original graph's unit-injected-flow transfer, not microphone pressure.
    pub response: NetworkDriveResponse,
    /// Actual source phasor, with owner-outward volume flow and explicit schema.
    pub inlet: HarmonicNetworkPort,
    /// Computed complete pressure/outward-flow matrix [Pa s/m^3].
    pub boundary: CoupledTermination,
    /// Actual solved opening flows [m^3/s], including the injected flow amplitude.
    pub outward_flows_m3_s: Vec<C64>,
    /// Actual coherent exterior pressure [Pa] and gradient [Pa/m]. No extra 1/r,
    /// delay or independent opening source is applied after the BEM observation.
    pub receiver_fields: Vec<PointField>,
    /// Raw mean power 1/2 Re sum conj(U_out)*p_opening [W]. A discrete surface
    /// work diagnostic, not an enclosed far-field or time-window energy claim.
    pub exterior_surface_power_w: f64,
}

impl RadiatingNetwork {
    /// Bind physical geometry once. The same coherent gas drives duct losses,
    /// acoustic propagation, the exterior load and all microphone observations.
    /// `source` identifies the caller's original geometry/material description;
    /// an additional boundary hash binds the exact mesh, trace and medium data.
    ///
    /// All receivers must be admitted exterior points of the original closed
    /// surface. The existing centroid pressure/gradient observation is retained;
    /// its clearance check is not a near-field quadrature-error certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn new(driving:NetworkDrivingPoint, surface:SpherePanels, patches:Vec<OpeningPatch>,
        gas:GasState, receivers:Vec<[f64;3]>, clearance_m:f64, source:ContentHash,
        budget:ExteriorNetworkBudget, gate:&CancelGate)->Result<Self,ExteriorNetworkError> {
        Self::new_inner(driving,surface,patches,gas,None,receivers,clearance_m,source,budget,gate)
    }

    /// Bind a complete gas record to every original network edge and retain a
    /// separate exterior record for BOTH radiation reaction and microphone
    /// pressure. All regions share static pressure; composition, temperature,
    /// transport and density may differ. Original thermochemical/pack identities
    /// remain available through [`Self::media`]. No mixing or mean flow is added.
    #[allow(clippy::too_many_arguments)]
    pub fn new_regional(driving:NetworkDrivingPoint, surface:SpherePanels, patches:Vec<OpeningPatch>,
        media:NetworkMedia, receivers:Vec<[f64;3]>, clearance_m:f64, source:ContentHash,
        budget:ExteriorNetworkBudget, gate:&CancelGate)->Result<Self,ExteriorNetworkError> {
        poll(gate)?;
        media.validate_for(driving.network())?;
        let gas=*media.exterior_gas();
        Self::new_inner(driving,surface,patches,gas,Some(media),receivers,clearance_m,source,budget,gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn new_inner(driving:NetworkDrivingPoint, surface:SpherePanels, patches:Vec<OpeningPatch>,
        gas:GasState, media:Option<NetworkMedia>, receivers:Vec<[f64;3]>, clearance_m:f64, source:ContentHash,
        budget:ExteriorNetworkBudget, gate:&CancelGate)->Result<Self,ExteriorNetworkError> {
        poll(gate)?;
        gas.validate_coherence().map_err(|_| bad("a coherent exterior gas state is required"))?;
        let n=surface.areas().len(); let p=patches.len(); let r=receivers.len();
        if n==0 || n>budget.max_panels || n>helmholtz::MAX_DENSE_PANELS
            || p==0 || p!=driving.coupled_labels().len() || p>fs_duct::network::MAX_NETWORK_PORTS
            || r==0 || r>fs_bem::near_field::MAX_RECEIVERS || source==ContentHash([0;32])
            || patches.iter().zip(driving.coupled_labels()).any(|(a,b)| &a.port!=b) {
            return Err(bad("nonempty bounded exterior, exact opening order, receivers and source are required"));
        }
        let entries=p.checked_mul(n).and_then(|v|v.checked_mul(3))
            .and_then(|v|p.checked_mul(r).and_then(|w|w.checked_mul(4)).and_then(|w|v.checked_add(w)))
            .ok_or_else(||bad("exterior field size overflow"))?;
        let terms=n.checked_mul(n).and_then(|n2|n2.checked_mul(n)
            .and_then(|v|p.checked_mul(n2).and_then(|w|w.checked_mul(3)).and_then(|w|v.checked_add(w))))
            .and_then(|v|p.checked_mul(n).and_then(|w|w.checked_mul(r)).and_then(|w|w.checked_mul(4))
                .and_then(|w|v.checked_add(w))).ok_or_else(||bad("exterior solve size overflow"))?;
        if entries>budget.max_field_entries || terms>budget.max_solve_terms {
            return Err(bad("complete exterior basis/solve exceeds its declared budget"));
        }
        // This original owner validates retained triangles, outward closed topology
        // and receiver separation. This is not a global self-intersection proof.
        Geometry::new(&surface,&receivers,clearance_m)?;
        let mut used=vec![false;n]; let mut traces=Vec::with_capacity(p); let mut areas=Vec::with_capacity(p);
        for patch in &patches {
            poll(gate)?;
            if patch.panels.is_empty() || patch.panels.len()>n { return Err(bad("opening needs actual surface panels")); }
            let mut area=0.0;
            for &panel in &patch.panels {
                if panel>=n || used[panel] { return Err(bad("opening panels must be valid and disjoint, including within a patch")); }
                used[panel]=true; area+=surface.areas()[panel];
            }
            if !area.is_finite() || area<=0.0 || !area.recip().is_finite() { return Err(bad("opening area cannot be normalized")); }
            let mut trace=vec![C64::ZERO;n];
            for &panel in &patch.panels { trace[panel]=C64::from_re(area.recip()); }
            areas.push(area); traces.push(trace);
        }
        let mut h=DomainHasher::new("fs-couple:network-exterior:v1"); h.update(source.as_bytes());
        for v in [n,p] { h.update(&(v as u64).to_le_bytes()); }
        for v in [gas.density,gas.sound_speed].into_iter()
            .chain(surface.centroids().iter().flatten().copied())
            .chain(surface.normals().iter().flatten().copied()).chain(surface.areas().iter().copied()) {
            h.update(&v.to_bits().to_le_bytes());
        }
        for triangle in surface.triangles().into_iter().flatten() {
            for v in triangle.iter().flatten() { h.update(&v.to_bits().to_le_bytes()); }
        }
        for (patch,trace) in patches.iter().zip(&traces) {
            h.update(&(patch.port.len() as u64).to_le_bytes()); h.update(patch.port.as_bytes());
            for v in trace { h.update(&v.re.to_bits().to_le_bytes()); }
        }
        let boundary_hash=h.finalize();
        poll(gate)?;
        Ok(Self{driving,surface,patches,traces,areas,gas,media,receivers,clearance_m,source,boundary_hash,budget})
    }
    /// Original duct geometry, scalar loads and coupled-port partition.
    #[must_use] pub fn driving_point(&self)->&NetworkDrivingPoint { &self.driving }
    /// Original capped exterior; all unselected panels are rigid.
    #[must_use] pub fn surface(&self)->&SpherePanels { &self.surface }
    /// Original opening labels and disjoint panel assignments.
    #[must_use] pub fn patches(&self)->&[OpeningPatch] { &self.patches }
    /// Physical patch areas used for both velocity and work projection.
    #[must_use] pub fn opening_areas_m2(&self)->&[f64] { &self.areas }
    /// Exterior medium; also the common interior gas for the uniform constructor.
    #[must_use] pub fn gas(&self)->GasState { self.gas }
    /// Complete original gas records on a regional assembly. A uniform assembly
    /// has no fabricated source records and returns `None`.
    #[must_use] pub fn media(&self)->Option<&NetworkMedia> { self.media.as_ref() }
    /// Original exterior pressure/gradient locations [m].
    #[must_use] pub fn receivers(&self)->&[[f64;3]] { &self.receivers }
    /// Caller-supplied description provenance, not a validation certificate.
    #[must_use] pub fn source(&self)->ContentHash { self.source }
    /// Exact exterior mesh, flow basis and acoustic medium identity.
    #[must_use] pub fn boundary_hash(&self)->ContentHash { self.boundary_hash }
    /// Original whole-call numerical limits.
    #[must_use] pub fn budget(&self)->ExteriorNetworkBudget { self.budget }

    /// Solve the entire exterior once for the complete unit-opening basis, then
    /// close the original duct-network Schur system with that SAME full matrix.
    /// All microphone fields reuse those same BEM solutions and solved flows.
    /// Cancellation is polled between original bounded synchronous operators;
    /// there is no new mid-LU cancellation or tighter latency guarantee.
    pub fn solve(&self,omega:f64,loss:LossModel,injected_volume_flow:C64,gate:&CancelGate)
        ->Result<RadiatingNetworkSample,ExteriorNetworkError> {
        poll(gate)?;
        if !omega.is_finite() || omega<=0.0 || !finite(injected_volume_flow) {
            return Err(bad("harmonic frequency must be positive and the physical source finite"));
        }
        let medium=Medium{density:self.gas.density,sound_speed:self.gas.sound_speed};
        let fields:Vec<_>=self.traces.iter().map(Vec::as_slice).collect();
        let solutions=helmholtz::solve_radiation_batch(&self.surface,omega/self.gas.sound_speed,
            medium,&fields,Formulation::BurtonMiller)?;
        poll(gate)?;
        let p=self.patches.len(); let mut z=vec![C64::ZERO;p*p];
        for (column,solution) in solutions.iter().enumerate() {
            for (row,patch) in self.patches.iter().enumerate() {
                let mut pressure=C64::ZERO;
                for &panel in &patch.panels {
                    pressure=pressure+solution.pressure[panel].scale(self.surface.areas()[panel]/self.areas[row]);
                }
                if !finite(pressure) { return Err(bad("projected opening pressure overflow")); }
                z[row*p+column]=pressure;
            }
        }
        let boundary=CoupledTermination{omega,labels:self.driving.coupled_labels().to_vec(),z,
            source_id:format!("fs-bem:network-openings:{}",self.boundary_hash)};
        // Complete Hermitian passivity is checked by the original duct owner.
        // A coarse active BEM matrix refuses; it is not repaired or diagonalized.
        let response=match &self.media {
            Some(media)=>self.driving.response_coupled_regional(media,omega,loss,&boundary)?,
            None=>self.driving.response_coupled(&self.gas,omega,loss,&boundary)?,
        };
        poll(gate)?;
        let labels:Vec<_>=self.driving.network().nodes.iter().filter_map(|n|match &n.boundary {
            NodeBoundary::Port{label}=>Some(label.as_str()), _=>None }).collect();
        let mut outward_flows_m3_s=Vec::with_capacity(p); let mut power=0.0;
        for patch in &self.patches {
            let index=labels.iter().position(|&label|label==patch.port).ok_or_else(||bad("opening label lost"))?;
            let flow=-response.port_inflows[index]*injected_volume_flow;
            let pressure=response.port_pressures[index]*injected_volume_flow;
            if !finite(flow) || !finite(pressure) { return Err(bad("physical opening phasor overflow")); }
            power+=0.5*(flow.conj()*pressure).re;
            outward_flows_m3_s.push(flow);
        }
        let mut receiver_fields=vec![PointField{pressure_pa:C64::ZERO,gradient_pa_per_m:[C64::ZERO;3]};self.receivers.len()];
        for (solution,&flow) in solutions.iter().zip(&outward_flows_m3_s) {
            poll(gate)?;
            let observations=laplace::observe_frequency(&self.surface,solution,&self.receivers,self.clearance_m)?;
            for (total,field) in receiver_fields.iter_mut().zip(observations) {
                total.pressure_pa=total.pressure_pa+field.pressure_pa*flow;
                for axis in 0..3 { total.gradient_pa_per_m[axis]=total.gradient_pa_per_m[axis]+field.gradient_pa_per_m[axis]*flow; }
            }
        }
        let pressure=response.impedance*injected_volume_flow;
        if !finite(pressure) || !power.is_finite() || receiver_fields.iter().any(|f|
            !finite(f.pressure_pa) || f.gradient_pa_per_m.iter().any(|&v|!finite(v))) {
            return Err(bad("combined physical observation overflow"));
        }
        poll(gate)?;
        Ok(RadiatingNetworkSample{inlet:HarmonicNetworkPort{omega,schema:self.driving.schema.clone(),
            pressure,outward_volume_flow:-injected_volume_flow}, response,boundary,outward_flows_m3_s,
            receiver_fields,exterior_surface_power_w:power})
    }
}

#[cfg(test)]
mod tests;
