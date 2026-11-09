//! Closed nonisothermal demixing with constant energetic square-gradient cost.
//!
//! An edge carries `E_ij = kappa A/(2 d) (x_i-x_j)^2`, allocated equally to
//! its endpoints. Each local energy is original homogeneous enthalpy plus its
//! half-edge allocation; entropy is the original homogeneous entropy. At fixed
//! local energy, the gradient contribution to the A-minus-B entropy potential
//! is `sum_j kappa A/(2 d) (1/T_i+1/T_j) (x_i-x_j)/N_i`. Neighbor temperatures
//! therefore enter the chemical force, and neighbors whose composition does
//! not change still recover temperature when their incident edge energy does.
//!
//! A face freezes temperature, source mobility and transported partial enthalpy
//! only for its convex-split composition solve. Ideal logarithms, negative chi
//! and the full graph gradient are implicit; positive chi is explicit. Original
//! calorics then recover all affected local temperatures, followed by implicit
//! Fourier exchange. This first-order split is not unconditionally entropy
//! stable: the complete original energy and entropy audits reject an inadmissible
//! interval. A smaller interval is required; the audits are never waived.
//!
//! This separate owner deliberately has no prescribed-flow advection API:
//! transporting interfacial energy requires its own mechanically consistent law.

use super::*;

/// Time-error control of composition, temperature and individual interface energies.
pub mod adaptive;
mod solve;

/// Accepted local temperatures and inventories with energetic interfaces.
///
/// The private bulk storage retains exactly the original graph and caloric law;
/// it is never exposed as a convex, zero-gradient transport field.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalPhaseField {
    bulk: BinaryNonisothermalField,
}

/// Conservative receipts for one complete canonical face sweep.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalPhaseFieldStep {
    pub duration_s: f64,
    /// Signed A transfer from canonical left to right [mol]; B is opposite.
    pub transferred_a_mol: Vec<f64>,
    /// Paired total local-energy transfer, including original partial enthalpy,
    /// authored reciprocal heat of transport and Fourier heat [J].
    pub transferred_energy_j: Vec<f64>,
    pub root_iterations: u64,
    pub caloric_iterations: u64,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    /// Original homogeneous enthalpy change plus interfacial energy change.
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
    pub homogeneous_enthalpy_change_j: f64,
    pub interfacial_energy_change_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryNonisothermalPhaseFieldCandidate {
    field: BinaryNonisothermalPhaseField,
    step: BinaryNonisothermalPhaseFieldStep,
}

impl BinaryNonisothermalPhaseFieldCandidate {
    pub const fn field(&self) -> &BinaryNonisothermalPhaseField {
        &self.field
    }
    pub const fn step(&self) -> &BinaryNonisothermalPhaseFieldStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BinaryNonisothermalPhaseField,
        BinaryNonisothermalPhaseFieldStep,
    ) {
        (self.field, self.step)
    }
}

impl BinaryMiscibilityNetwork {
    /// Initialize a closed local-temperature phase field from original liquid
    /// calorics, positive kappa and authored canonical thermal face data.
    ///
    /// The original chi(T), including excess heat capacity, may be nonconvex.
    /// Pure species endpoints are excluded because this local-temperature lane
    /// requires common positive species support. The initial temperatures are
    /// an initialization, not a conservative temperature event.
    pub fn bind_nonisothermal_phase_field(
        &self,
        initial: &BinaryMiscibilityState,
        caloric: BinaryLiquidCaloricModelV1,
        temperatures: Vec<Temperature>,
        faces: Vec<BinaryNonisothermalFace>,
        temperature_bracket: [Temperature; 2],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryNonisothermalPhaseField> {
        poll(&mut cancelled)?;
        self.admit(initial)?;
        let bracket = temperature_bracket.map(|t| t.value());
        let domain = caloric.temperature_domain().map(|t| t.value());
        if !positive(self.model().gradient_coefficient_j_m())
            || self.model().thermodynamics() != caloric.thermodynamics()
            || temperatures.len() != self.cells().len()
            || faces.len() != self.owner.faces.len()
            || faces.iter().any(|f| {
                !nonnegative(f.thermal_conductivity_w_m_k)
                    || !f.heat_of_transport_j_per_mol.is_finite()
            })
            || !positive(bracket[0])
            || !positive(bracket[1])
            || bracket[0] >= bracket[1]
            || bracket[0] < domain[0]
            || bracket[1] > domain[1]
        {
            return Err(BinaryNonisothermalError::Invalid(
                "original binary calorics, positive energetic kappa, supported thermal bracket and canonical thermal faces required",
            ));
        }
        // Bounds cover retained candidates, amount/temperature arrays, sparse
        // incident patches, receipt arrays and scalar caloric/root work. Source
        // material query storage remains owned by the original material source.
        let workspace_bytes = self
            .cells()
            .len()
            .checked_mul(4096)
            .and_then(|n| {
                self.owner
                    .faces
                    .len()
                    .checked_mul(1024)
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| n.checked_add(self.workspace_bytes()))
            .and_then(|n| n.checked_add(core::mem::size_of::<Owner>()))
            .ok_or(BinaryNonisothermalError::Invalid(
                "nonisothermal phase-field workspace size overflow",
            ))?;
        if workspace_bytes > maximum_workspace_bytes {
            return Err(BinaryNonisothermalError::Workspace {
                required: workspace_bytes,
                maximum: maximum_workspace_bytes,
            });
        }
        let field = BinaryNonisothermalPhaseField {
            bulk: BinaryNonisothermalField {
                owner: Arc::new(Owner {
                    network: self.clone(),
                    caloric,
                    faces: faces.into(),
                    temperature_bracket,
                    workspace_bytes,
                }),
                state: initial.clone(),
                temperatures,
            },
        };
        field.bulk.qualify(&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(field)
    }
}

impl BinaryNonisothermalPhaseField {
    pub fn network(&self) -> &BinaryMiscibilityNetwork {
        self.bulk.network()
    }
    pub fn caloric_model(&self) -> &BinaryLiquidCaloricModelV1 {
        self.bulk.caloric_model()
    }
    pub fn component_moles(&self) -> &[[f64; 2]] {
        self.bulk.component_moles()
    }
    pub fn mole_fractions(&self) -> &[[f64; 2]] {
        self.bulk.mole_fractions()
    }
    pub fn temperatures(&self) -> &[Temperature] {
        self.bulk.temperatures()
    }
    pub fn elapsed_s(&self) -> f64 {
        self.bulk.elapsed_s()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.bulk.workspace_bytes()
    }
    pub fn transport_faces(&self) -> &[BinaryNonisothermalFace] {
        self.bulk.transport_faces()
    }
    pub fn temperature_bracket(&self) -> [Temperature; 2] {
        self.bulk.temperature_bracket()
    }

    /// One bounded implicit composition/Fourier sweep. Publication is atomic
    /// after original species, capacity, total-energy and entropy checks.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: BinaryNonisothermalSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryNonisothermalPhaseFieldCandidate> {
        solve::step(self, duration_s, controls, cancelled)
    }

    fn balances(
        &self,
        next: &Self,
        controls: BinaryNonisothermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([f64; 2], f64, f64, f64, f64, f64)> {
        let mut amounts = [Sum::default(); 2];
        let mut homogeneous = Sum::default();
        let mut interfacial = Sum::default();
        let mut entropy = Sum::default();
        let mut capacity = 0.0_f64;
        for cell in 0..self.temperatures().len() {
            poll(cancelled)?;
            let mut residual = Sum::default();
            for s in 0..2 {
                amounts[s]
                    .add(next.bulk.state.amounts[cell][s] - self.bulk.state.amounts[cell][s])?;
                amounts[s].add(
                    next.bulk.state.amount_remainders[cell][s]
                        - self.bulk.state.amount_remainders[cell][s],
                )?;
                residual.add(next.bulk.state.amounts[cell][s])?;
                residual.add(next.bulk.state.amount_remainders[cell][s])?;
            }
            residual.add(-self.network().owner.capacities_mol[cell])?;
            capacity = capacity.max(residual.value()?.abs());
            let change = self.caloric_model().homogeneous_transport_change(
                &self.bulk.caloric_state(cell)?,
                &next.bulk.caloric_state(cell)?,
            )?;
            homogeneous.add(change.enthalpy_change_j)?;
            entropy.add(change.entropy_change_j_per_k)?;
        }
        for face in self.network().owner.faces.iter() {
            poll(cancelled)?;
            let [i, j] = face.physical.cells;
            let old = self.mole_fractions()[i][0] - self.mole_fractions()[j][0];
            let new = next.mole_fractions()[i][0] - next.mole_fractions()[j][0];
            interfacial.add(finite(
                (face.gradient_j * (new - old)) * (0.5 * new + 0.5 * old),
                "original interfacial energy change",
            )?)?;
        }
        let species = [amounts[0].value()?, amounts[1].value()?];
        if capacity > controls.amount_balance_tolerance_mol
            || species
                .iter()
                .any(|r| r.abs() > controls.amount_balance_tolerance_mol)
        {
            return Err(BinaryNonisothermalError::Species {
                residual_mol: species,
                tolerance_mol: controls.amount_balance_tolerance_mol,
            });
        }
        let homogeneous = homogeneous.value()?;
        let interfacial = interfacial.value()?;
        let mut total = Sum::default();
        total.add(homogeneous)?;
        total.add(interfacial)?;
        let total = total.value()?;
        if total.abs() > controls.energy_tolerance_j {
            return Err(BinaryNonisothermalError::Energy {
                residual_j: total,
                tolerance_j: controls.energy_tolerance_j,
            });
        }
        let entropy = entropy.value()?;
        if entropy < -controls.entropy_tolerance_j_per_k {
            return Err(BinaryNonisothermalError::Entropy {
                change_j_per_k: entropy,
                tolerance_j_per_k: controls.entropy_tolerance_j_per_k,
            });
        }
        Ok((species, capacity, total, entropy, homogeneous, interfacial))
    }
}

#[cfg(test)]
mod tests;
