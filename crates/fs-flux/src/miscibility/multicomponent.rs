//! Closed isothermal multicomponent liquids with full Maxwell–Stefan transport.
//!
//! Each component has its own conserved amount; none is reconstructed from a
//! carrier. A common fixed partial molar volume makes the zero-net-molar frame
//! preserve cell capacity. Every face uses the exact original chemical-potential
//! jump and one coupled friction solve, with an explicitly authored symmetric
//! molecular-diffusivity matrix. Arithmetic face compositions are a consistent
//! spatial approximation, not an exact finite-jump ideal Fick formula.
//!
//! This zero-interface-energy lane requires common positive species support
//! and a sufficient global convexity bound. The bound is inferred from the
//! original interaction spectrum at uniform composition, not merely from local
//! stability there. Nonconvex laws require an interfacial model and are refused.
//! One synchronous explicit step either preserves positive inventories and
//! decreases the actual original free energy within its allowance, or refuses.
//! There is no clipping, dependent-species projection or hidden subcycling.

use fs_math::det::ln;
use fs_qty::{AmountConcentration, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::solution::multicomponent::{
    MAX_REGULAR_SOLUTION_SPECIES_V1, MULTICOMPONENT_COMPOSITION_SUM_TOLERANCE_V1,
    MulticomponentRegularSolutionV1, MulticomponentSolutionErrorV1,
    stability::{MulticomponentStabilityControlV1, MulticomponentStabilityErrorV1},
};
use fs_thermochem::transport::maxwell_stefan::{
    MaxwellStefanErrorV1, MaxwellStefanMolarSolveV1, MaxwellStefanV1,
};
use std::sync::Arc;

mod solve;

/// Composition accuracy over complete physical intervals with bounded work.
pub mod adaptive;

/// Conservative prescribed volume flow for every original component.
pub mod advection;
/// Same-interval stirring and full molecular mixing, with coupled adaptation.
pub mod transport;

/// Local-temperature species and heat exchange using the original caloric law.
pub mod nonisothermal;

/// Nonconvex liquids with explicit species-symmetric interfacial energy.
pub mod phase_field;

/// Original liquid material properties retained through accepted endpoints.
#[cfg(any(feature = "continuum", feature = "experimental-miscibility-material"))]
pub mod material;

/// Fixed cell geometry; the common concentration sets its molar capacity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentMiscibilityCell {
    pub volume_m3: f64,
}

/// Canonical cell pair and an index into shared original molecular pair data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentMiscibilityFace {
    pub cells: [usize; 2],
    pub area_m2: f64,
    pub distance_m: f64,
    pub diffusivity_matrix: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MulticomponentMiscibilityError {
    Invalid(&'static str),
    Arithmetic(&'static str),
    Thermodynamics(MulticomponentSolutionErrorV1),
    Transport(MaxwellStefanErrorV1),
    Stability(MulticomponentStabilityErrorV1),
    Workspace {
        required: usize,
        maximum: usize,
    },
    /// A sufficient whole-simplex curvature lower bound could not be positive.
    Convexity {
        lower_bound: f64,
    },
    Positivity {
        cell: usize,
        species: usize,
        trial_mol: f64,
    },
    Conservation {
        component_residual_mol: Vec<f64>,
        maximum_capacity_residual_mol: f64,
        tolerance_mol: f64,
    },
    FreeEnergy {
        change_j: f64,
        tolerance_j: f64,
    },
    /// Prescribed flow must preserve the original closed cell volumes.
    VolumeDivergence {
        cell: usize,
        inflow_m3_s: f64,
        outflow_m3_s: f64,
        relative_residual: f64,
    },
    AdvectionWorkLimit {
        required_substeps: usize,
        maximum_substeps: usize,
    },
    Cancelled,
}
impl core::fmt::Display for MulticomponentMiscibilityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "multicomponent finite-volume miscibility: {self:?}")
    }
}
impl std::error::Error for MulticomponentMiscibilityError {}
impl From<MulticomponentSolutionErrorV1> for MulticomponentMiscibilityError {
    fn from(error: MulticomponentSolutionErrorV1) -> Self {
        if error == MulticomponentSolutionErrorV1::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermodynamics(error)
        }
    }
}
impl From<MaxwellStefanErrorV1> for MulticomponentMiscibilityError {
    fn from(error: MaxwellStefanErrorV1) -> Self {
        if error == MaxwellStefanErrorV1::Cancelled {
            Self::Cancelled
        } else {
            Self::Transport(error)
        }
    }
}
impl From<MulticomponentStabilityErrorV1> for MulticomponentMiscibilityError {
    fn from(error: MulticomponentStabilityErrorV1) -> Self {
        if error == MulticomponentStabilityErrorV1::Cancelled {
            Self::Cancelled
        } else {
            Self::Stability(error)
        }
    }
}
type Result<T> = core::result::Result<T, MulticomponentMiscibilityError>;

/// Original constitutive equation tolerances and complete-window allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MulticomponentMiscibilitySolve {
    pub maxwell_stefan: MaxwellStefanMolarSolveV1,
    pub amount_balance_tolerance_mol: f64,
    pub free_energy_tolerance_j: f64,
}
impl Default for MulticomponentMiscibilitySolve {
    fn default() -> Self {
        Self {
            maxwell_stefan: MaxwellStefanMolarSolveV1::default(),
            amount_balance_tolerance_mol: 1e-12,
            free_energy_tolerance_j: 1e-9,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Owner {
    model: MulticomponentRegularSolutionV1,
    transport: MaxwellStefanV1,
    temperature: Temperature,
    concentration: AmountConcentration,
    cells: Box<[MulticomponentMiscibilityCell]>,
    faces: Box<[MulticomponentMiscibilityFace]>,
    capacities_mol: Box<[f64]>,
    diffusivity_matrices: Box<[Vec<f64>]>,
    chi: Box<[f64]>,
    convexity_lower_bound: f64,
    workspace_bytes: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityNetwork {
    owner: Arc<Owner>,
}

/// Independent canonical component inventories in contiguous cell-major arrays.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityField {
    network: MulticomponentMiscibilityNetwork,
    amounts: Vec<f64>,
    amount_remainders: Vec<f64>,
    fractions: Vec<f64>,
    elapsed_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityStep {
    pub duration_s: f64,
    /// Canonical face then canonical species; positive means left to right.
    pub transferred_mol: Vec<Vec<f64>>,
    pub maxwell_stefan_evaluations: usize,
    pub maximum_gradient_residual_per_m: f64,
    pub maximum_molar_frame_residual_mol_m2_s: f64,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    pub free_energy_change_j: f64,
    /// Initial face chemical dissipation summed over the graph [W]. This is
    /// an instantaneous diagnostic, not a finite-step dissipation inequality.
    pub dissipation_rate_w: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentMiscibilityCandidate {
    field: MulticomponentMiscibilityField,
    step: MulticomponentMiscibilityStep,
}
impl MulticomponentMiscibilityCandidate {
    pub const fn field(&self) -> &MulticomponentMiscibilityField {
        &self.field
    }
    pub const fn step(&self) -> &MulticomponentMiscibilityStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        MulticomponentMiscibilityField,
        MulticomponentMiscibilityStep,
    ) {
        (self.field, self.step)
    }
}

impl MulticomponentMiscibilityNetwork {
    /// Bind original canonical species, a fixed isothermal law, common partial
    /// molar volume and explicit shared pair matrices. Each matrix is symmetric
    /// row-major n*n, with zero diagonal and positive off-diagonals. Faces are
    /// sorted by cell pair; their matrix index travels with their geometry.
    ///
    /// For a simplex tangent v, ideal curvature is at least 2||v||². At uniform
    /// composition the existing entropy-metric spectrum is I+(P chi P)/n, so
    /// `2+n*(lambda_min-1-uncertainty)>0` is a sufficient global convexity bound.
    /// It is conservative and includes numerical eigensolve uncertainty; it is
    /// not a necessary miscibility test or an interval arithmetic certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        model: MulticomponentRegularSolutionV1,
        transport: MaxwellStefanV1,
        temperature: Temperature,
        concentration: AmountConcentration,
        cells: Vec<MulticomponentMiscibilityCell>,
        faces: Vec<MulticomponentMiscibilityFace>,
        diffusivity_matrices: Vec<Vec<f64>>,
        maximum_workspace_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::prepare(
            model,
            transport,
            temperature,
            concentration,
            cells,
            faces,
            diffusivity_matrices,
            maximum_workspace_bytes,
            true,
            cancelled,
        )
    }

    /// Shared graph preparation. Only the positive-gradient-energy owner may
    /// request nonconvex preparation, and it never exposes this raw network or
    /// its raw fields. Public homogeneous construction retains its convexity gate.
    #[allow(clippy::too_many_arguments)]
    fn prepare(
        model: MulticomponentRegularSolutionV1,
        transport: MaxwellStefanV1,
        temperature: Temperature,
        concentration: AmountConcentration,
        cells: Vec<MulticomponentMiscibilityCell>,
        mut faces: Vec<MulticomponentMiscibilityFace>,
        diffusivity_matrices: Vec<Vec<f64>>,
        maximum_workspace_bytes: usize,
        require_convexity: bool,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let n = model.species().len();
        if model.species() != transport.species()
            || cells.is_empty()
            || cells.len() > super::MAX_MISCIBILITY_CELLS
            || faces.len() > super::MAX_MISCIBILITY_FACES
            || !positive(concentration.value())
            || cells.iter().any(|cell| !positive(cell.volume_m3))
        {
            return Err(MulticomponentMiscibilityError::Invalid(
                "same original species, positive concentration and finite cell volumes required",
            ));
        }
        let slots =
            cells
                .len()
                .checked_mul(n)
                .ok_or(MulticomponentMiscibilityError::Arithmetic(
                    "cell/species dimensions",
                ))?;
        let face_slots =
            faces
                .len()
                .checked_mul(n)
                .ok_or(MulticomponentMiscibilityError::Arithmetic(
                    "face/species dimensions",
                ))?;
        let matrix_slots = diffusivity_matrices.len().checked_mul(n * n).ok_or(
            MulticomponentMiscibilityError::Arithmetic("diffusivity table dimensions"),
        )?;
        let workspace_bytes = slots
            .checked_mul(512)
            .and_then(|v| face_slots.checked_mul(256).and_then(|w| v.checked_add(w)))
            .and_then(|v| matrix_slots.checked_mul(8).and_then(|w| v.checked_add(w)))
            .and_then(|v| (n * n).checked_mul(128).and_then(|w| v.checked_add(w)))
            .and_then(|v| cells.len().checked_mul(64).and_then(|w| v.checked_add(w)))
            .and_then(|v| faces.len().checked_mul(512).and_then(|w| v.checked_add(w)))
            .and_then(|v| v.checked_add(32 * 1024 + core::mem::size_of::<Owner>()))
            .ok_or(MulticomponentMiscibilityError::Arithmetic(
                "multicomponent workspace size",
            ))?;
        if workspace_bytes > maximum_workspace_bytes {
            return Err(MulticomponentMiscibilityError::Workspace {
                required: workspace_bytes,
                maximum: maximum_workspace_bytes,
            });
        }
        for matrix in &diffusivity_matrices {
            poll(&mut cancelled)?;
            if matrix.len() != n * n {
                return Err(MulticomponentMiscibilityError::Invalid(
                    "complete original pair diffusivity matrix required",
                ));
            }
            for i in 0..n {
                poll(&mut cancelled)?;
                if matrix[i * n + i] != 0.0 {
                    return Err(MulticomponentMiscibilityError::Invalid(
                        "zero molecular diffusivity diagonal required",
                    ));
                }
                for j in 0..i {
                    if !positive(matrix[i * n + j]) || matrix[i * n + j] != matrix[j * n + i] {
                        return Err(MulticomponentMiscibilityError::Invalid(
                            "positive exactly symmetric molecular pair diffusivities required",
                        ));
                    }
                }
            }
        }
        for face in &mut faces {
            poll(&mut cancelled)?;
            face.cells.sort();
            if face.cells[0] == face.cells[1]
                || face.cells[1] >= cells.len()
                || !nonnegative(face.area_m2)
                || !positive(face.distance_m)
                || face.diffusivity_matrix >= diffusivity_matrices.len()
            {
                return Err(MulticomponentMiscibilityError::Invalid(
                    "finite face geometry, distinct endpoints and original pair matrix required",
                ));
            }
        }
        faces.sort_by_key(|face| face.cells);
        if faces.windows(2).any(|pair| pair[0].cells == pair[1].cells) {
            return Err(MulticomponentMiscibilityError::Invalid(
                "duplicate physical face",
            ));
        }
        let uniform = vec![1.0 / n as f64; n];
        let state = model.evaluate(temperature, &uniform, &mut cancelled)?;
        let convexity_lower_bound = if n == 1 {
            2.0
        } else {
            let stability = state.local_stability(
                MulticomponentStabilityControlV1 {
                    max_sweeps: 64,
                    max_workspace_bytes: maximum_workspace_bytes,
                    ..MulticomponentStabilityControlV1::default()
                },
                &mut cancelled,
            )?;
            finite(
                2.0 + n as f64
                    * (stability.smallest_eigenvalue().ok_or(
                        MulticomponentMiscibilityError::Arithmetic("missing tangent spectrum"),
                    )? - 1.0
                        - stability.numerical_uncertainty()),
                "global simplex curvature lower bound",
            )?
        };
        if require_convexity && convexity_lower_bound <= 64.0 * f64::EPSILON * n as f64 {
            return Err(MulticomponentMiscibilityError::Convexity {
                lower_bound: convexity_lower_bound,
            });
        }
        let mut chi = Vec::with_capacity(n * n);
        for (&a, &b) in model
            .chi_constants()
            .iter()
            .zip(model.chi_inverse_temperature_k())
        {
            chi.push(finite(
                a + b / temperature.value(),
                "original isothermal interaction",
            )?);
        }
        let mut capacities = Vec::with_capacity(cells.len());
        for cell in &cells {
            poll(&mut cancelled)?;
            let capacity = finite(
                cell.volume_m3 * concentration.value(),
                "fixed molar capacity",
            )?;
            if !positive(capacity) {
                return Err(MulticomponentMiscibilityError::Arithmetic(
                    "positive molar capacity underflows",
                ));
            }
            capacities.push(capacity);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            owner: Arc::new(Owner {
                model,
                transport,
                temperature,
                concentration,
                cells: cells.into_boxed_slice(),
                faces: faces.into_boxed_slice(),
                capacities_mol: capacities.into_boxed_slice(),
                diffusivity_matrices: diffusivity_matrices.into_boxed_slice(),
                chi: chi.into_boxed_slice(),
                convexity_lower_bound,
                workspace_bytes,
            }),
        })
    }

    pub fn model(&self) -> &MulticomponentRegularSolutionV1 {
        &self.owner.model
    }
    pub fn transport(&self) -> &MaxwellStefanV1 {
        &self.owner.transport
    }
    pub fn temperature(&self) -> Temperature {
        self.owner.temperature
    }
    pub fn concentration(&self) -> AmountConcentration {
        self.owner.concentration
    }
    pub fn cells(&self) -> &[MulticomponentMiscibilityCell] {
        &self.owner.cells
    }
    pub fn faces(&self) -> &[MulticomponentMiscibilityFace] {
        &self.owner.faces
    }
    pub fn capacities_mol(&self) -> &[f64] {
        &self.owner.capacities_mol
    }
    pub fn species_count(&self) -> usize {
        self.model().species().len()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.owner.workspace_bytes
    }
    pub fn convexity_lower_bound(&self) -> f64 {
        self.owner.convexity_lower_bound
    }

    /// Retain every offered component independently. Only roundoff-sized
    /// capacity mismatch is admitted; amounts are never renormalized.
    pub fn bind(
        &self,
        component_moles: Vec<Vec<f64>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityField> {
        poll(&mut cancelled)?;
        let n = self.species_count();
        if component_moles.len() != self.cells().len()
            || component_moles.iter().any(|row| row.len() != n)
        {
            return Err(MulticomponentMiscibilityError::Invalid(
                "one independent amount per cell and canonical species required",
            ));
        }
        let mut amounts = Vec::with_capacity(component_moles.len() * n);
        let mut fractions = Vec::with_capacity(component_moles.len() * n);
        for (cell, row) in component_moles.into_iter().enumerate() {
            poll(&mut cancelled)?;
            let total = amount_sum(&row)?;
            let capacity = self.capacities_mol()[cell];
            if (total - capacity).abs() > MULTICOMPONENT_COMPOSITION_SUM_TOLERANCE_V1 * capacity {
                return Err(MulticomponentMiscibilityError::Invalid(
                    "offered inventory must match fixed molar capacity",
                ));
            }
            let x = self.fractions(&row, &mut cancelled)?;
            fractions.extend_from_slice(&x[..n]);
            amounts.extend(row);
        }
        poll(&mut cancelled)?;
        Ok(MulticomponentMiscibilityField {
            network: self.clone(),
            amount_remainders: vec![0.0; amounts.len()],
            amounts,
            fractions,
            elapsed_s: 0.0,
        })
    }

    fn fractions(
        &self,
        amounts: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; MAX_REGULAR_SOLUTION_SPECIES_V1]> {
        let total = amount_sum(amounts)?;
        let mut fractions = [0.0; MAX_REGULAR_SOLUTION_SPECIES_V1];
        for (s, &amount) in amounts.iter().enumerate() {
            fractions[s] = amount / total;
            if !positive(fractions[s]) {
                return Err(MulticomponentMiscibilityError::Arithmetic(
                    "positive component fraction is unrepresentable",
                ));
            }
        }
        let state =
            self.model()
                .evaluate(self.temperature(), &fractions[..amounts.len()], cancelled)?;
        fractions[..amounts.len()].copy_from_slice(state.mole_fractions());
        Ok(fractions)
    }
}

impl MulticomponentMiscibilityField {
    pub const fn network(&self) -> &MulticomponentMiscibilityNetwork {
        &self.network
    }
    /// Published cell-major amounts in canonical original species order.
    /// Thermodynamic and material states use these representable amounts;
    /// conservation and adaptive errors also retain their carried low words.
    pub fn component_moles(&self) -> &[f64] {
        &self.amounts
    }
    /// Cell-major, then canonical original species order.
    pub fn mole_fractions(&self) -> &[f64] {
        &self.fractions
    }
    pub fn cell_component_moles(&self, cell: usize) -> Option<&[f64]> {
        let n = self.network.species_count();
        let start = cell.checked_mul(n)?;
        self.amounts.get(start..start.checked_add(n)?)
    }
    pub fn cell_mole_fractions(&self, cell: usize) -> Option<&[f64]> {
        let n = self.network.species_count();
        let start = cell.checked_mul(n)?;
        self.fractions.get(start..start.checked_add(n)?)
    }
    pub const fn elapsed_s(&self) -> f64 {
        self.elapsed_s
    }
    pub fn workspace_bytes(&self) -> usize {
        self.network.workspace_bytes()
    }

    /// One synchronous explicit finite-volume step with bounded work. A large
    /// interval may refuse positivity or actual free-energy descent; refinement
    /// belongs to the separate adaptive owner, never to this physical step.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: MulticomponentMiscibilitySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentMiscibilityCandidate> {
        solve::step(self, duration_s, controls, cancelled)
    }

    fn balances(
        &self,
        next: &Self,
        controls: MulticomponentMiscibilitySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64, f64)> {
        let (components, capacity) =
            self.inventory_balances(next, controls.amount_balance_tolerance_mol, cancelled)?;
        let energy = solve::free_energy_change(self, next, cancelled)?;
        if energy > controls.free_energy_tolerance_j {
            return Err(MulticomponentMiscibilityError::FreeEnergy {
                change_j: energy,
                tolerance_j: controls.free_energy_tolerance_j,
            });
        }
        Ok((components, capacity, energy))
    }

    /// Shared independent-inventory account. Thermal and interfacial owners
    /// apply their complete energy audits without a homogeneous Gibbs gate.
    fn inventory_balances(
        &self,
        next: &Self,
        tolerance_mol: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64)> {
        let n = self.network.species_count();
        let mut components = vec![Sum::default(); n];
        let mut capacity = 0.0_f64;
        for cell in 0..self.network.cells().len() {
            poll(cancelled)?;
            let mut total = Sum::default();
            for (s, component) in components.iter_mut().enumerate() {
                let index = cell * n + s;
                component.add(next.amounts[index] - self.amounts[index])?;
                component.add(next.amount_remainders[index] - self.amount_remainders[index])?;
                total.add(next.amounts[index])?;
                total.add(next.amount_remainders[index])?;
            }
            total.add(-self.network.capacities_mol()[cell])?;
            capacity = capacity.max(total.value()?.abs());
        }
        let components = components
            .iter()
            .map(Sum::value)
            .collect::<Result<Vec<_>>>()?;
        if capacity > tolerance_mol || components.iter().any(|v| v.abs() > tolerance_mol) {
            return Err(MulticomponentMiscibilityError::Conservation {
                component_residual_mol: components,
                maximum_capacity_residual_mol: capacity,
                tolerance_mol,
            });
        }
        Ok((components, capacity))
    }
}

fn admit_controls(controls: MulticomponentMiscibilitySolve) -> Result<()> {
    let ms = controls.maxwell_stefan;
    if [
        controls.amount_balance_tolerance_mol,
        controls.free_energy_tolerance_j,
    ]
    .into_iter()
    .any(|v| !nonnegative(v))
        || !positive(ms.relative_tolerance)
        || ms.relative_tolerance > 1e-3
        || !positive(ms.absolute_gradient_tolerance_per_m)
        || !positive(ms.absolute_molar_flux_tolerance_mol_m2_s)
    {
        return Err(MulticomponentMiscibilityError::Invalid(
            "finite nonnegative physical and constitutive tolerances required",
        ));
    }
    Ok(())
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64, name: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(MulticomponentMiscibilityError::Arithmetic(name))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MulticomponentMiscibilityError::Cancelled)
    } else {
        Ok(())
    }
}
fn amount_sum(amounts: &[f64]) -> Result<f64> {
    let mut sum = Sum::default();
    for &amount in amounts {
        if !positive(amount) {
            return Err(MulticomponentMiscibilityError::Invalid(
                "common strictly positive species support required",
            ));
        }
        sum.add(amount)?;
    }
    let total = sum.value()?;
    if !positive(total) {
        return Err(MulticomponentMiscibilityError::Arithmetic(
            "positive total inventory",
        ));
    }
    Ok(total)
}

#[derive(Clone, Copy, Default)]
struct Sum {
    high: f64,
    low: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<()> {
        let next = finite(self.high + value, "compensated multicomponent sum")?;
        let error = if self.high.abs() >= value.abs() {
            (self.high - next) + value
        } else {
            (value - next) + self.high
        };
        self.low = finite(self.low + error, "compensated multicomponent remainder")?;
        self.high = next;
        Ok(())
    }
    fn value(&self) -> Result<f64> {
        finite(self.high + self.low, "compensated multicomponent value")
    }
}

/// Original deterministic ln(a/b), with exact reversal and close-value accuracy.
fn log_ratio(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let (lo, hi, sign) = if a > b { (b, a, 1.0) } else { (a, b, -1.0) };
    let delta = (hi - lo) / lo;
    let rounded = 1.0 + delta;
    sign * if rounded == 1.0 {
        delta
    } else if rounded.is_finite() {
        ln(rounded) * (delta / (rounded - 1.0))
    } else {
        ln(hi) - ln(lo)
    }
}

#[cfg(test)]
mod precision_tests;
#[cfg(test)]
mod tests;
