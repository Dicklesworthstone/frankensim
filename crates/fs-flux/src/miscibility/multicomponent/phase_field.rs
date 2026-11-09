//! Conservative spatial demixing of original multicomponent liquids.
//!
//! The bulk free energy is the SAME homogeneous regular-solution law, including
//! its nonconvex branch. An independently authored positive constant kappa
//! [J/m] adds `E_e = kappa A/(4 d) sum_s (x_js-x_is)^2` on each graph face.
//! This species-symmetric convention reduces exactly to the binary coefficient
//! `kappa A/(2 d) (x_jA-x_iA)^2`. Kappa is not a surface tension [J/m^2].
//!
//! Every original component retains its own amount and carried low word. Full
//! Maxwell--Stefan transport acts on bulk plus gradient chemical potentials in
//! the fixed molar/volume frame, with equal fixed partial molar volumes. The
//! same complete original pair-diffusivity matrices and source species order
//! are retained; no dependent component or pairwise Fick replacement is used.
//!
//! The synchronous explicit step accepts only positive independent inventories,
//! closed species/capacity balances and descent of the actual TOTAL free energy.
//! It is not unconditionally stable: the square-gradient term introduces a
//! fourth-order spatial time-step restriction. A refused interval must be
//! refined by its caller. Geometry has the same orthogonal two-point or authored
//! compartment interpretation as the original graph.
//!
//! This owner deliberately never exposes its private zero-gradient network or
//! field. Nonconvex admission is confined to this positive-interface owner.
//! Prescribed-flow advection requires a separate consistent capillary work law.

use super::*;

mod solve;

/// Original multicomponent thermodynamics and molecular transport with a
/// positive species-symmetric interfacial coefficient on a fixed graph.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentPhaseFieldNetwork {
    bulk: MulticomponentMiscibilityNetwork,
    gradient_coefficient_j_m: f64,
    /// Canonical edge energy coefficient kappa A/(4 d) [J].
    edge_coefficients_j: Arc<[f64]>,
    workspace_bytes: usize,
}

/// All original independent component amounts in a closed isothermal phase
/// field. The private ordinary field supplies inventory storage only.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentPhaseField {
    network: MulticomponentPhaseFieldNetwork,
    bulk: MulticomponentMiscibilityField,
}

/// Complete-interval inventory and total-free-energy accounting.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentPhaseFieldStep {
    pub duration_s: f64,
    /// Canonical face then canonical species; positive means left to right.
    pub transferred_mol: Vec<Vec<f64>>,
    pub maxwell_stefan_evaluations: usize,
    pub maximum_gradient_residual_per_m: f64,
    pub maximum_molar_frame_residual_mol_m2_s: f64,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub bulk_free_energy_change_j: f64,
    pub interfacial_energy_change_j: f64,
    /// Original homogeneous free-energy change plus interface-energy change.
    pub free_energy_change_j: f64,
    /// Initial bulk-plus-gradient chemical dissipation [W]. This instantaneous
    /// diagnostic is separate from the complete finite-step energy inequality.
    pub dissipation_rate_w: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentPhaseFieldCandidate {
    field: MulticomponentPhaseField,
    step: MulticomponentPhaseFieldStep,
}

impl MulticomponentPhaseFieldCandidate {
    pub const fn field(&self) -> &MulticomponentPhaseField {
        &self.field
    }

    pub const fn step(&self) -> &MulticomponentPhaseFieldStep {
        &self.step
    }

    pub fn into_parts(self) -> (MulticomponentPhaseField, MulticomponentPhaseFieldStep) {
        (self.field, self.step)
    }
}

impl MulticomponentPhaseFieldNetwork {
    /// Retain original homogeneous thermodynamics, full molecular friction and
    /// graph geometry. A positive constant kappa regularizes nonconvex bulk
    /// laws; the ordinary zero-gradient constructor still requires convexity.
    ///
    /// The additional workspace bound covers retained edge coefficients,
    /// compensated cell gradients, reduced potentials and candidate wrappers.
    /// The original network's bound already covers simultaneous bulk candidates,
    /// independent inventory low words and the reusable original MS workspace.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        model: MulticomponentRegularSolutionV1,
        transport: MaxwellStefanV1,
        temperature: Temperature,
        concentration: AmountConcentration,
        cells: Vec<MulticomponentMiscibilityCell>,
        faces: Vec<MulticomponentMiscibilityFace>,
        diffusivity_matrices: Vec<Vec<f64>>,
        gradient_coefficient_j_m: f64,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if !positive(gradient_coefficient_j_m) {
            return Err(MulticomponentMiscibilityError::Invalid(
                "positive finite energetic multicomponent kappa required",
            ));
        }
        let rt = finite(R * temperature.value(), "phase-field RT")?;
        if !positive(rt) {
            return Err(MulticomponentMiscibilityError::Invalid(
                "positive representable phase-field RT required",
            ));
        }
        let extra = cells
            .len()
            .checked_mul(model.species().len())
            .and_then(|v| v.checked_mul(64))
            .and_then(|v| faces.len().checked_mul(8).and_then(|e| v.checked_add(e)))
            .and_then(|v| v.checked_add(3 * core::mem::size_of::<MulticomponentPhaseField>()))
            .ok_or(MulticomponentMiscibilityError::Invalid(
                "phase-field workspace size overflow",
            ))?;
        let bulk_budget = maximum_workspace_bytes.checked_sub(extra).ok_or(
            MulticomponentMiscibilityError::Workspace {
                required: extra,
                maximum: maximum_workspace_bytes,
            },
        )?;
        let bulk = MulticomponentMiscibilityNetwork::prepare(
            model,
            transport,
            temperature,
            concentration,
            cells,
            faces,
            diffusivity_matrices,
            bulk_budget,
            false,
            &mut cancelled,
        )
        .map_err(|error| match error {
            MulticomponentMiscibilityError::Workspace { required, .. } => {
                MulticomponentMiscibilityError::Workspace {
                    required: required.saturating_add(extra),
                    maximum: maximum_workspace_bytes,
                }
            }
            other => other,
        })?;
        let workspace_bytes = bulk.workspace_bytes().checked_add(extra).ok_or(
            MulticomponentMiscibilityError::Invalid("phase-field workspace size overflow"),
        )?;
        let mut edge_coefficients_j = Vec::with_capacity(bulk.faces().len());
        for face in bulk.faces() {
            poll(&mut cancelled)?;
            if face.area_m2 == 0.0 {
                edge_coefficients_j.push(0.0);
                continue;
            }
            // Reuse the original range-conscious product/quotient kernel:
            // kappa*A/(4*d) may be representable even when A/d alone is not.
            let coefficient = super::solve::pair_force(
                gradient_coefficient_j_m,
                face.area_m2,
                0.25,
                face.distance_m,
            )?;
            edge_coefficients_j.push(coefficient);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            bulk,
            gradient_coefficient_j_m,
            edge_coefficients_j: edge_coefficients_j.into(),
            workspace_bytes,
        })
    }

    pub fn model(&self) -> &MulticomponentRegularSolutionV1 {
        self.bulk.model()
    }

    pub fn transport(&self) -> &MaxwellStefanV1 {
        self.bulk.transport()
    }

    pub fn temperature(&self) -> Temperature {
        self.bulk.temperature()
    }

    pub fn concentration(&self) -> AmountConcentration {
        self.bulk.concentration()
    }

    pub fn cells(&self) -> &[MulticomponentMiscibilityCell] {
        self.bulk.cells()
    }

    pub fn faces(&self) -> &[MulticomponentMiscibilityFace] {
        self.bulk.faces()
    }

    pub fn capacities_mol(&self) -> &[f64] {
        self.bulk.capacities_mol()
    }

    pub fn species_count(&self) -> usize {
        self.bulk.species_count()
    }

    pub const fn gradient_coefficient_j_m(&self) -> f64 {
        self.gradient_coefficient_j_m
    }

    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Retain every independently offered original component. The underlying
    /// initializer performs its original positive-support and capacity checks,
    /// but its raw network/field never becomes a public zero-gradient owner.
    pub fn bind(
        &self,
        component_moles: Vec<Vec<f64>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentPhaseField> {
        poll(&mut cancelled)?;
        let bulk = self.bulk.bind(component_moles, &mut cancelled)?;
        let field = MulticomponentPhaseField {
            network: self.clone(),
            bulk,
        };
        // Qualify actual initial storage rather than allowing the first audit
        // to encounter an overflowing interfacial energy.
        field.interfacial_energy_j(&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(field)
    }
}

impl MulticomponentPhaseField {
    pub const fn network(&self) -> &MulticomponentPhaseFieldNetwork {
        &self.network
    }

    /// Published cell-major independent amounts; carried low words remain in
    /// conservative inventory and time-error accounts.
    pub fn component_moles(&self) -> &[f64] {
        self.bulk.component_moles()
    }

    pub fn mole_fractions(&self) -> &[f64] {
        self.bulk.mole_fractions()
    }

    pub fn cell_component_moles(&self, cell: usize) -> Option<&[f64]> {
        self.bulk.cell_component_moles(cell)
    }

    pub fn cell_mole_fractions(&self, cell: usize) -> Option<&[f64]> {
        self.bulk.cell_mole_fractions(cell)
    }

    pub const fn elapsed_s(&self) -> f64 {
        self.bulk.elapsed_s()
    }

    pub const fn workspace_bytes(&self) -> usize {
        self.network.workspace_bytes()
    }

    /// Actual nonnegative interface storage, summed in canonical face order.
    pub fn interfacial_energy_j(&self, mut cancelled: impl FnMut() -> bool) -> Result<f64> {
        let mut energy = Sum::default();
        for face in 0..self.network.faces().len() {
            poll(&mut cancelled)?;
            energy.add(solve::edge_energy_j(self, face)?)?;
        }
        poll(&mut cancelled)?;
        energy.value()
    }

    /// One synchronous explicit full-Maxwell--Stefan step. The existing solve
    /// control's free-energy tolerance applies to bulk PLUS interface energy.
    /// A refused candidate leaves this field and its original sources intact.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: MulticomponentMiscibilitySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentPhaseFieldCandidate> {
        solve::step(self, duration_s, controls, cancelled)
    }

    fn balances(
        &self,
        next: &Self,
        controls: MulticomponentMiscibilitySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64, f64, f64, f64)> {
        admit_controls(controls)?;
        let (components, capacity) = self.bulk.inventory_balances(
            &next.bulk,
            controls.amount_balance_tolerance_mol,
            cancelled,
        )?;
        let bulk = super::solve::free_energy_change(&self.bulk, &next.bulk, cancelled)?;
        let interfacial = solve::interfacial_energy_change_j(self, next, cancelled)?;
        let total = solve::admit_total_energy(bulk, interfacial, controls)?;
        Ok((components, capacity, bulk, interfacial, total))
    }
}

#[cfg(test)]
mod tests;
