//! Complete material valve images in the original mass-normalized solid basis.
//! The legacy card may omit zero-opening modes; this owner never does. Its
//! frozen K/C approximation and finite material support are admitted by the
//! same projection as the legacy path, with full damping blocks retained.

use super::{Approximation, Controls, invalid, poll, reduce_material};
use crate::linear3::visco::{TetViscoField, ViscoElements};
use crate::reduce::{ReduceError, ValveCardRequest};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::{CancelGate, Cx};
use fs_modal::ModePair;

/// Equilibrium-source comparators and the complete selected frozen image.
/// These are static responses, not broadband or omitted-mode error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalStaticCompliance {
    /// Mean opening per conjugate mean-opening force [m/N].
    pub exact_opening_force_m_per_n: f64,
    pub retained_opening_force_m_per_n: f64,
    /// Pressure-work volume per pascal [m3/Pa].
    pub exact_volume_m3_per_pa: f64,
    pub retained_volume_m3_per_pa: f64,
    /// Signed mean opening per pascal toward opening [m/Pa].
    pub exact_opening_m_per_pa: f64,
    pub retained_opening_m_per_pa: f64,
}

/// Immutable complete selected material image. The actual tetrahedral owner,
/// source laws, support map and all original M-orthonormal shapes stay attached.
/// A mode with zero mean opening can still carry pressure work or communicate
/// with another mode through C; it is never removed from this image.
///
/// Coordinates q have units m sqrt(kg), with M=I, K=diag(lambda), and the full
/// projected viscosity in approximation(). No independent-mode/reed conversion
/// is provided. The legacy MaterialValve keeps its coupled-block refusal.
pub struct MaterialModalValve {
    source: ViscoElements,
    surface: Option<super::MaterialValveSurface>,
    modes: Vec<ModePair>,
    approximation: Approximation,
    pressure_projection: Vec<f64>,
    opening_projection: Vec<f64>,
    nodes: usize,
    rest_opening_m: f64,
    width_m: f64,
    total_mass_kg: f64,
    source_id: String,
    identity: ContentHash,
    compliance: ModalStaticCompliance,
}

impl MaterialModalValve {
    pub fn source(&self) -> &ViscoElements {
        &self.source
    }
    /// Original exterior geometry; absent for unsupported higher-order traces.
    pub fn surface_geometry(&self) -> Option<&super::MaterialValveSurface> {
        self.surface.as_ref()
    }
    /// Original shapes in canonical free-DOF order; phi^T M phi = I.
    pub fn modes(&self) -> &[ModePair] {
        &self.modes
    }
    pub fn approximation(&self) -> &Approximation {
        &self.approximation
    }
    /// phi^T load_per_pa [m2/sqrt(kg)], positive toward opening.
    /// Its transpose is the swept-volume observation on modal velocity.
    pub fn pressure_projection_m2_per_sqrt_kg(&self) -> &[f64] {
        &self.pressure_projection
    }
    /// Mean opening per mass-normalized coordinate [1/sqrt(kg)].
    /// Exact zero entries remain present and distinct from pressure coupling.
    pub fn opening_projection_per_sqrt_kg(&self) -> &[f64] {
        &self.opening_projection
    }
    pub const fn rest_opening_m(&self) -> f64 {
        self.rest_opening_m
    }
    pub const fn width_m(&self) -> f64 {
        self.width_m
    }
    pub const fn total_mass_kg(&self) -> f64 {
        self.total_mass_kg
    }
    pub fn source_id(&self) -> &str {
        &self.source_id
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn compliance(&self) -> &ModalStaticCompliance {
        &self.compliance
    }
    pub const fn node_count(&self) -> usize {
        self.nodes
    }

    /// Project an arbitrary nodal force field with the SAME original shapes.
    /// Input N produces generalized force N/sqrt(kg). A dimensionless nodal
    /// force pattern per unit applied N produces the corresponding port map.
    /// Fixed components do no work. No mean-opening division is involved.
    pub fn project_nodal_force(
        &self,
        force_n: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<Vec<f64>, ReduceError> {
        poll_gate(gate)?;
        if force_n.len() != self.nodes || force_n.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid(
                "complete modal force needs one finite vector per source node",
            ));
        }
        let mut result = Vec::with_capacity(self.modes.len());
        for mode in &self.modes {
            poll_gate(gate)?;
            let value = mode
                .phi
                .iter()
                .zip(&self.source.equilibrium().free_dofs)
                .map(|(phi, &dof)| phi * force_n[dof / 3][dof % 3])
                .sum::<f64>();
            if !value.is_finite() {
                return Err(invalid("complete modal force projection overflow"));
            }
            result.push(value);
        }
        poll_gate(gate)?;
        Ok(result)
    }

    /// Expand original mass-normalized displacements or velocities into all
    /// source nodes. Fixed components remain zero. This is the exact transpose
    /// partner of project_nodal_force, not a separately normalized observer.
    pub fn expand_nodal(
        &self,
        coordinates: &[f64],
        gate: &CancelGate,
    ) -> Result<Vec<[f64; 3]>, ReduceError> {
        self.expand_nodal_cancellable(coordinates, || gate.is_requested())
    }

    /// Admit actual modal displacement at every original material quadrature
    /// point. Uses the same source shapes and B matrices as volume mechanics;
    /// no bound is inferred from modal amplitude or mean opening alone.
    pub fn validate_modal_displacement(
        &self,
        coordinates: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ReduceError> {
        let nodal = self.expand_nodal_cancellable(coordinates, &mut cancelled)?;
        let displacement: Vec<_> = nodal.into_iter().flatten().collect();
        self.source.validate_displacement(&displacement, cancelled)?;
        Ok(())
    }

    fn expand_nodal_cancellable(
        &self,
        coordinates: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, ReduceError> {
        let mut poll = || -> Result<(), ReduceError> {
            if cancelled() { Err(crate::linear3::TetElasticError::Cancelled.into()) }
            else { Ok(()) }
        };
        poll()?;
        if coordinates.len() != self.modes.len() || coordinates.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "complete modal observation needs all finite original coordinates",
            ));
        }
        let mut result = vec![[0.0; 3]; self.nodes];
        for (row, &dof) in self.source.equilibrium().free_dofs.iter().enumerate() {
            poll()?;
            let value = self
                .modes
                .iter()
                .zip(coordinates)
                .map(|(mode, q)| mode.phi[row] * q)
                .sum::<f64>();
            if !value.is_finite() {
                return Err(invalid("complete modal observation overflow"));
            }
            result[dof / 3][dof % 3] = value;
        }
        poll()?;
        Ok(result)
    }
}

fn poll_gate(gate: &CancelGate) -> Result<(), ReduceError> {
    if gate.is_requested() {
        Err(crate::linear3::TetElasticError::Cancelled.into())
    } else {
        Ok(())
    }
}

/// Reduce the SAME source material, mesh and pressure load with all selected
/// modes and full nonproportional viscosity. The explicit frozen-frequency
/// storage/loss dispersion checks, equilibrium compliance, source bands and
/// cancellation behavior remain mandatory. No authored diagonal loss is added.
///
/// Unlike the compatibility MaterialValve, this result must be consumed by a
/// coupled structural owner. At least one selected mode must move the declared
/// mean slit, but every zero-opening selected mode is also retained.
pub fn reduce_modal_pressure_loaded_valve(
    request: &ValveCardRequest<'_, TetViscoField<'_>>,
    nodal_load_per_pa_m2: &[[f64; 3]],
    controls: &Controls,
    cx: &Cx<'_>,
) -> Result<MaterialModalValve, ReduceError> {
    let (source, mut reduced) = reduce_material(request, nodal_load_per_pa_m2, controls, true, cx)?;
    let surface = super::MaterialValveSurface::from_request(request, controls.assembly, cx)?;
    let modes = reduced
        .complete_material_modes
        .take()
        .ok_or_else(|| invalid("complete material valve basis was not retained"))?;
    let approximation = reduced
        .material_loss
        .take()
        .ok_or_else(|| invalid("complete material valve viscosity was not retained"))?;
    let static_source = reduced
        .pressure_compliance
        .ok_or_else(|| invalid("complete material valve static response was not retained"))?;
    let axis = crate::reduce::unit(request.orifice.opening_axis)?;
    let free = &source.equilibrium().free_dofs;
    let mut pressure_projection = Vec::with_capacity(modes.len());
    let mut opening_projection = Vec::with_capacity(modes.len());
    let mut compliance = ModalStaticCompliance {
        exact_opening_force_m_per_n: static_source.exact_mean_force_m_per_n(),
        retained_opening_force_m_per_n: 0.0,
        exact_volume_m3_per_pa: static_source.exact_volume_m3_per_pa(),
        retained_volume_m3_per_pa: 0.0,
        exact_opening_m_per_pa: static_source.exact_opening_m_per_pa(),
        retained_opening_m_per_pa: 0.0,
    };
    for mode in &modes {
        poll(cx)?;
        let opening =
            crate::reduce::mean_face_translation(&mode.phi, free, request.orifice.face_nodes, axis);
        let pressure = crate::reduce::pressure_area(&mode.phi, free, nodal_load_per_pa_m2, 1.0);
        if !opening.is_finite() || !pressure.is_finite() {
            return Err(invalid(
                "complete material valve has a nonfinite physical projection",
            ));
        }
        compliance.retained_opening_force_m_per_n += opening * (opening / mode.lambda);
        compliance.retained_volume_m3_per_pa += pressure * (pressure / mode.lambda);
        compliance.retained_opening_m_per_pa += opening * (pressure / mode.lambda);
        opening_projection.push(opening);
        pressure_projection.push(pressure);
    }
    if [
        compliance.retained_opening_force_m_per_n,
        compliance.retained_volume_m3_per_pa,
        compliance.retained_opening_m_per_pa,
    ]
    .iter()
    .any(|v| !v.is_finite())
    {
        return Err(invalid("complete material valve static response overflow"));
    }
    let mass = reduced.card.total_mass_kg;
    if !mass.is_finite() || mass <= 0.0 {
        return Err(invalid(
            "complete material valve mass must be positive and finite",
        ));
    }
    // Identity includes every dark mode, coupled coefficient and signed port;
    // the compatibility card's smaller mode list cannot identify this image.
    let mut hash = DomainHasher::new("org.frankensim.fs-solid.complete-material-valve.v1");
    hash.update(reduced.card.identity.as_bytes());
    hash.update(&(request.nodes_m.len() as u64).to_le_bytes());
    hash.update(&(free.len() as u64).to_le_bytes());
    for &dof in free {
        hash.update(&(dof as u64).to_le_bytes());
    }
    hash.update(&(modes.len() as u64).to_le_bytes());
    for mode in &modes {
        hash.update(&mode.lambda.to_bits().to_le_bytes());
        hash.update(&(mode.phi.len() as u64).to_le_bytes());
        for value in &mode.phi {
            hash.update(&value.to_bits().to_le_bytes());
        }
    }
    for value in approximation
        .modal_viscosity_per_s
        .iter()
        .chain(&pressure_projection)
        .chain(&opening_projection)
    {
        hash.update(&value.to_bits().to_le_bytes());
    }
    for value in nodal_load_per_pa_m2.iter().flatten() {
        hash.update(&value.to_bits().to_le_bytes());
    }
    poll(cx)?;
    Ok(MaterialModalValve {
        source,
        surface,
        modes,
        approximation,
        pressure_projection,
        opening_projection,
        nodes: request.nodes_m.len(),
        rest_opening_m: reduced.card.rest_gap_m,
        width_m: reduced.card.width_m,
        total_mass_kg: mass,
        source_id: request.source_id.to_string(),
        identity: hash.finalize(),
        compliance,
    })
}
