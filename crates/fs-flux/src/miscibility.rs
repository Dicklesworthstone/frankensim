//! Conservative finite-rate binary miscibility on a closed finite-volume graph.
//!
//! The original homogeneous regular-solution law supplies bulk chemical energy;
//! positive square-gradient energy regularizes spinodal decomposition. Cells
//! retain fixed volumes and the SAME site concentration, hence equal fixed
//! partial molar volumes. Each face exchanges A and B in opposite molar amounts.
//! This is zero net molar/volume flux. Unequal molar masses give a nonzero net
//! mass flux in this frame and changing local density; the barycentric mass
//! frame therefore differs from the fixed molar/volume frame used here.
//! Temperature and pressure are held externally; this is not caloric evolution,
//! advection, momentum, wetting/contact-line motion or a nucleation-rate model.
//!
//! Each physical face advances through one implicit convex-split exchange.
//! Logarithmic ideal terms, convex interactions and the FULL graph's gradient
//! energy are implicit; concave positive-chi interactions use the preceding
//! pair state. A monotone scalar solve retains both species without clipping or
//! a seeded trace. Canonical face sweeps are first order in time. Large steps
//! may be thermodynamically admissible but inaccurate: refine time and space.
//!
//! The convex-split construction follows the log/gradient-implicit principle of
//! Chen, Wang, Wang and Wise, JCP:X 3 (2019) 100031,
//! <https://doi.org/10.1016/j.jcpx.2019.100031>. The conservative pair splitting
//! here is a separate first-order graph discretization, not their global scheme.
//! For one solved exchange q, convexity gives delta F <= -q delta(mu_cs).
//! Every published candidate independently checks the whole graph's actual
//! free-energy change, that dissipation bound and both species balances.

use fs_thermochem::solution::BinarySolutionErrorV1;
use fs_thermochem::solution::gradient::BinarySquareGradientV1;
use std::sync::Arc;

#[path = "miscibility/solve.rs"]
mod solve;

pub mod adaptive;

/// Hard admission bounds keep every allocation and face sweep finite.
pub const MAX_MISCIBILITY_CELLS: usize = 65_536;
pub const MAX_MISCIBILITY_FACES: usize = 262_144;

/// One rigid finite volume. Its molar capacity is original c times volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryMiscibilityCell {
    pub volume_m3: f64,
}

/// One orthogonal two-point face or explicitly authored compartment connection.
/// The area/distance factor must be valid for the caller's spatial geometry;
/// arbitrary nonorthogonal meshes do not acquire a consistency certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryMiscibilityFace {
    pub cells: [usize; 2],
    pub area_m2: f64,
    pub distance_m: f64,
    /// Original binary Maxwell-Stefan diffusivity [m²/s], independent of chi.
    /// Explicit zero closes molecular transport, while gradient energy remains.
    pub diffusivity_m2_s: f64,
}

/// Scalar-root accuracy, actual conserved-inventory accuracy and finite work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryMiscibilitySolve {
    pub maximum_root_iterations: usize,
    pub absolute_transfer_tolerance_mol: f64,
    pub relative_transfer_tolerance: f64,
    /// Whole-step per-species and per-cell capacity allowance [mol].
    pub amount_balance_tolerance_mol: f64,
    /// Whole-step free-energy and convex-split dissipation allowance [J].
    /// This allowance is never multiplied by the number of faces.
    pub free_energy_tolerance_j: f64,
}

impl Default for BinaryMiscibilitySolve {
    fn default() -> Self {
        Self {
            maximum_root_iterations: 128,
            absolute_transfer_tolerance_mol: 1e-14,
            relative_transfer_tolerance: 1e-12,
            amount_balance_tolerance_mol: 1e-12,
            free_energy_tolerance_j: 1e-9,
        }
    }
}

/// No error mutates an accepted state or publishes a partially swept graph.
#[derive(Clone, Debug, PartialEq)]
pub enum BinaryMiscibilityError {
    Invalid(&'static str),
    Arithmetic(&'static str),
    Thermodynamics(BinarySolutionErrorV1),
    Workspace {
        required: usize,
        maximum: usize,
    },
    Root {
        face: usize,
        iterations: usize,
        residual_mol: f64,
    },
    Conservation {
        component_residual_mol: [f64; 2],
        capacity_residual_mol: f64,
    },
    FreeEnergy {
        change_j: f64,
        dissipation_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for BinaryMiscibilityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "binary finite-volume miscibility: {self:?}")
    }
}
impl std::error::Error for BinaryMiscibilityError {}
impl From<BinarySolutionErrorV1> for BinaryMiscibilityError {
    fn from(error: BinarySolutionErrorV1) -> Self {
        if error == BinarySolutionErrorV1::Cancelled {
            Self::Cancelled
        } else {
            Self::Thermodynamics(error)
        }
    }
}
type Result<T> = core::result::Result<T, BinaryMiscibilityError>;

#[derive(Clone, Debug, PartialEq)]
struct Face {
    physical: BinaryMiscibilityFace,
    gradient_j: f64,
    mobility_mol2_j_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
struct Owner {
    model: BinarySquareGradientV1,
    cells: Vec<BinaryMiscibilityCell>,
    capacities_mol: Vec<f64>,
    faces: Vec<Face>,
    adjacent_offsets: Vec<usize>,
    adjacent_faces: Vec<usize>,
    gradient_diagonal_j: Vec<f64>,
    workspace_bytes: usize,
}

/// Prepared immutable graph, shared with every accepted state.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityNetwork {
    owner: Arc<Owner>,
}

impl BinaryMiscibilityNetwork {
    /// Prepare canonical faces and CSR incidence once. Reordered or reversed
    /// input faces produce the same physical sweep. Duplicate connections refuse.
    /// The byte budget charges retained array payload and conservative transaction
    /// scratch, excluding allocator overhead and strings inside the pair model.
    /// On bounded-valence FV grids work is linear in face count times the bounded
    /// scalar iterations; general graph cache updates cost sum(degree²).
    pub fn new(
        model: BinarySquareGradientV1,
        cells: Vec<BinaryMiscibilityCell>,
        mut faces: Vec<BinaryMiscibilityFace>,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let n = cells.len();
        let m = faces.len();
        if n == 0 || n > MAX_MISCIBILITY_CELLS || m > MAX_MISCIBILITY_FACES {
            return Err(BinaryMiscibilityError::Invalid(
                "bounded nonempty finite-volume graph required",
            ));
        }
        let workspace_bytes = n
            .checked_mul(512)
            .and_then(|x| m.checked_mul(512).and_then(|y| x.checked_add(y)))
            .and_then(|x| x.checked_add(core::mem::size_of::<Owner>()))
            .ok_or(BinaryMiscibilityError::Arithmetic("workspace size"))?;
        if workspace_bytes > maximum_workspace_bytes {
            return Err(BinaryMiscibilityError::Workspace {
                required: workspace_bytes,
                maximum: maximum_workspace_bytes,
            });
        }
        let mut capacities_mol = Vec::with_capacity(n);
        for cell in &cells {
            poll(&mut cancelled)?;
            if !positive(cell.volume_m3) {
                return Err(BinaryMiscibilityError::Invalid(
                    "positive finite cell volumes required",
                ));
            }
            let capacity = cell.volume_m3 * model.molar_concentration_mol_m3();
            if !positive(capacity) {
                return Err(BinaryMiscibilityError::Arithmetic("cell molar capacity"));
            }
            capacities_mol.push(capacity);
        }
        for face in &mut faces {
            poll(&mut cancelled)?;
            face.cells.sort();
            if face.cells[0] == face.cells[1]
                || face.cells[1] >= n
                || !nonnegative(face.area_m2)
                || !positive(face.distance_m)
                || !nonnegative(face.diffusivity_m2_s)
            {
                return Err(BinaryMiscibilityError::Invalid(
                    "valid endpoints, nonnegative area/diffusivity and positive distance required",
                ));
            }
        }
        faces.sort_by_key(|face| face.cells);
        if faces.windows(2).any(|pair| pair[0].cells == pair[1].cells) {
            return Err(BinaryMiscibilityError::Invalid(
                "duplicate physical finite-volume face",
            ));
        }
        let mut prepared = Vec::with_capacity(m);
        let mut counts = vec![0usize; n];
        let mut diagonal = vec![Sum::default(); n];
        for physical in faces {
            poll(&mut cancelled)?;
            let geometry = finite(physical.area_m2 / physical.distance_m, "face geometry")?;
            if physical.area_m2 > 0.0 && geometry == 0.0 {
                return Err(BinaryMiscibilityError::Arithmetic(
                    "positive face geometry underflows",
                ));
            }
            let gradient_j = finite(
                model.gradient_coefficient_j_m() * geometry,
                "gradient face weight",
            )?;
            let mobility_mol2_j_s = finite(
                (model.molar_concentration_mol_m3() * physical.diffusivity_m2_s)
                    * (geometry / model.rt_j_per_mol()),
                "face mobility",
            )?;
            if (geometry > 0.0 && model.gradient_coefficient_j_m() > 0.0 && gradient_j == 0.0)
                || (geometry > 0.0 && physical.diffusivity_m2_s > 0.0 && mobility_mol2_j_s == 0.0)
            {
                return Err(BinaryMiscibilityError::Arithmetic(
                    "positive face coefficient underflows",
                ));
            }
            for &i in &physical.cells {
                counts[i] += 1;
                diagonal[i].add(gradient_j)?;
            }
            prepared.push(Face {
                physical,
                gradient_j,
                mobility_mol2_j_s,
            });
        }
        let mut adjacent_offsets = Vec::with_capacity(n + 1);
        adjacent_offsets.push(0);
        for count in counts {
            adjacent_offsets.push(adjacent_offsets.last().copied().unwrap() + count);
        }
        let mut cursor = adjacent_offsets[..n].to_vec();
        let mut adjacent_faces = vec![0usize; 2 * m];
        for (face, data) in prepared.iter().enumerate() {
            for &cell in &data.physical.cells {
                adjacent_faces[cursor[cell]] = face;
                cursor[cell] += 1;
            }
        }
        let gradient_diagonal_j = diagonal
            .into_iter()
            .map(|sum| sum.value())
            .collect::<Result<Vec<_>>>()?;
        poll(&mut cancelled)?;
        Ok(Self {
            owner: Arc::new(Owner {
                model,
                cells,
                capacities_mol,
                faces: prepared,
                adjacent_offsets,
                adjacent_faces,
                gradient_diagonal_j,
                workspace_bytes,
            }),
        })
    }

    pub fn model(&self) -> &BinarySquareGradientV1 {
        &self.owner.model
    }
    pub fn cells(&self) -> &[BinaryMiscibilityCell] {
        &self.owner.cells
    }
    pub fn faces(&self) -> impl ExactSizeIterator<Item = &BinaryMiscibilityFace> {
        self.owner.faces.iter().map(|face| &face.physical)
    }
    pub fn capacities_mol(&self) -> &[f64] {
        &self.owner.capacities_mol
    }
    pub fn workspace_bytes(&self) -> usize {
        self.owner.workspace_bytes
    }

    /// Initialize from independently represented fractions, retaining a real
    /// trace in either species order. Pure cells are admitted without a log floor.
    pub fn state_from_mole_fractions(
        &self,
        fractions: Vec<[f64; 2]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityState> {
        if fractions.len() != self.owner.cells.len() {
            return Err(BinaryMiscibilityError::Invalid(
                "one composition per finite volume required",
            ));
        }
        let mut amounts = Vec::with_capacity(fractions.len());
        for (fraction, &capacity) in fractions.into_iter().zip(&self.owner.capacities_mol) {
            poll(&mut cancelled)?;
            let x = self
                .model()
                .evaluate_composition(fraction)?
                .mole_fractions();
            let smaller = usize::from(x[1] < x[0]);
            let small = capacity * x[smaller];
            if x[smaller] > 0.0 && !positive(small) {
                return Err(BinaryMiscibilityError::Arithmetic(
                    "positive component inventory underflows",
                ));
            }
            let mut pair = [0.0; 2];
            pair[smaller] = small;
            pair[1 - smaller] = capacity - small;
            amounts.push(pair);
        }
        self.state_from_component_moles(amounts, cancelled)
    }

    /// Bind original paired inventories to these fixed cV capacities. No input
    /// amount is renormalized; a capacity mismatch beyond the original
    /// composition roundoff budget refuses.
    pub fn state_from_component_moles(
        &self,
        amounts: Vec<[f64; 2]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityState> {
        if amounts.len() != self.owner.cells.len() {
            return Err(BinaryMiscibilityError::Invalid(
                "one conserved pair per finite volume required",
            ));
        }
        let mut fractions = Vec::with_capacity(amounts.len());
        for (pair, &capacity) in amounts.iter().zip(&self.owner.capacities_mol) {
            poll(&mut cancelled)?;
            if pair.iter().any(|&n| !nonnegative(n)) {
                return Err(BinaryMiscibilityError::Invalid(
                    "finite nonnegative component amounts required",
                ));
            }
            fractions.push(self.fractions(*pair, capacity)?);
        }
        poll(&mut cancelled)?;
        Ok(BinaryMiscibilityState {
            owner: Arc::clone(&self.owner),
            amount_remainders: vec![[0.0; 2]; amounts.len()],
            amounts,
            fractions,
            elapsed_s: 0.0,
        })
    }

    /// Actual bulk plus square-gradient free energy [J], with no local flash.
    pub fn free_energy_j(
        &self,
        state: &BinaryMiscibilityState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64> {
        self.admit(state)?;
        let mut sum = Sum::default();
        for (&fraction, &capacity) in state.fractions.iter().zip(&self.owner.capacities_mol) {
            poll(&mut cancelled)?;
            sum.add(finite(
                capacity
                    * self
                        .model()
                        .evaluate_composition(fraction)?
                        .mixing_gibbs_j_per_mol(),
                "bulk free energy",
            )?)?;
        }
        for face in &self.owner.faces {
            poll(&mut cancelled)?;
            let [i, j] = face.physical.cells;
            let d = state.fractions[i][0] - state.fractions[j][0];
            sum.add(finite(
                (0.5 * face.gradient_j * d) * d,
                "gradient free energy",
            )?)?;
        }
        sum.value()
    }

    /// Physical chemical potentials of the current field, including every
    /// adjacent gradient term. A pure cell retains its divergent potential as
    /// None; implicit exchange does not evaluate an invented finite trace.
    pub fn exchange_potentials_j_per_mol(
        &self,
        state: &BinaryMiscibilityState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<Option<f64>>> {
        self.admit(state)?;
        let gradient = self.gradient(&state.fractions, &mut cancelled)?;
        state
            .fractions
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                poll(&mut cancelled)?;
                if x[0] == 0.0 || x[1] == 0.0 {
                    return Ok(None);
                }
                Ok(Some(finite(
                    self.model().exchange_potential_j_per_mol(x)?
                        + gradient[i].value()? / self.owner.capacities_mol[i],
                    "field chemical potential",
                )?))
            })
            .collect()
    }

    /// One atomic canonical face sweep. The scalar roots are implicit and
    /// positivity preserving; temporal accuracy remains first order. Every
    /// gradient contribution, both species and the whole-step energy inequality
    /// are checked before a new immutable state becomes visible.
    pub fn step_candidate(
        &self,
        initial: &BinaryMiscibilityState,
        duration_s: f64,
        controls: BinaryMiscibilitySolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityCandidate> {
        solve::step(self, initial, duration_s, controls, cancelled)
    }

    fn admit(&self, state: &BinaryMiscibilityState) -> Result<()> {
        if !Arc::ptr_eq(&self.owner, &state.owner) && self.owner != state.owner {
            return Err(BinaryMiscibilityError::Invalid(
                "state requires its original field model and finite-volume graph",
            ));
        }
        Ok(())
    }
    fn fractions(&self, amounts: [f64; 2], capacity: f64) -> Result<[f64; 2]> {
        let fractions = amounts.map(|amount| amount / capacity);
        if amounts
            .iter()
            .zip(fractions)
            .any(|(&amount, x)| amount > 0.0 && x == 0.0)
        {
            return Err(BinaryMiscibilityError::Arithmetic(
                "component lost in field normalization",
            ));
        }
        Ok(self
            .model()
            .evaluate_composition(fractions)?
            .mole_fractions())
    }
    fn gradient(
        &self,
        fractions: &[[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Sum>> {
        let mut result = vec![Sum::default(); fractions.len()];
        for face in &self.owner.faces {
            poll(cancelled)?;
            let [i, j] = face.physical.cells;
            let contribution = finite(
                face.gradient_j * (fractions[i][0] - fractions[j][0]),
                "gradient derivative",
            )?;
            result[i].add(contribution)?;
            result[j].add(-contribution)?;
        }
        Ok(result)
    }
    fn free_energy_change(
        &self,
        old: &BinaryMiscibilityState,
        new: &BinaryMiscibilityState,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        let mut sum = Sum::default();
        for i in 0..old.amounts.len() {
            poll(cancelled)?;
            sum.add(finite(
                self.owner.capacities_mol[i]
                    * self
                        .model()
                        .mixing_free_energy_change_j_per_mol(old.fractions[i], new.fractions[i])?,
                "bulk energy increment",
            )?)?;
        }
        for face in &self.owner.faces {
            poll(cancelled)?;
            let [i, j] = face.physical.cells;
            let a = old.fractions[i][0] - old.fractions[j][0];
            let b = new.fractions[i][0] - new.fractions[j][0];
            sum.add(finite(
                (face.gradient_j * (b - a)) * (0.5 * b + 0.5 * a),
                "gradient energy increment",
            )?)?;
        }
        sum.value()
    }
}

/// Conserved finite-volume inventories with exact original constitutive owner.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityState {
    owner: Arc<Owner>,
    amounts: Vec<[f64; 2]>,
    // Retain the low word of each accepted component exchange across steps.
    // Independently rounded debits and credits otherwise accumulate a false
    // capacity change, even when every individual transfer is conservative.
    amount_remainders: Vec<[f64; 2]>,
    fractions: Vec<[f64; 2]>,
    elapsed_s: f64,
}
impl BinaryMiscibilityState {
    pub fn component_moles(&self) -> &[[f64; 2]] {
        &self.amounts
    }
    pub fn mole_fractions(&self) -> &[[f64; 2]] {
        &self.fractions
    }
    pub const fn elapsed_s(&self) -> f64 {
        self.elapsed_s
    }
    pub fn model(&self) -> &BinarySquareGradientV1 {
        &self.owner.model
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityStep {
    pub duration_s: f64,
    /// Signed A transfer, with opposite B transfer, in canonical face order.
    pub transferred_a_mol: Vec<f64>,
    pub root_iterations: usize,
    pub maximum_root_residual_mol: f64,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub free_energy_before_j: f64,
    pub free_energy_after_j: f64,
    /// Direct increments, preserving small changes beside large bulk energies.
    pub free_energy_change_j: f64,
    /// Sum of accepted q times convex-split chemical-potential drop [J].
    pub convex_split_dissipation_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityCandidate {
    state: BinaryMiscibilityState,
    step: BinaryMiscibilityStep,
}
impl BinaryMiscibilityCandidate {
    pub const fn state(&self) -> &BinaryMiscibilityState {
        &self.state
    }
    pub const fn step(&self) -> &BinaryMiscibilityStep {
        &self.step
    }
    pub fn into_parts(self) -> (BinaryMiscibilityState, BinaryMiscibilityStep) {
        (self.state, self.step)
    }
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
        Err(BinaryMiscibilityError::Arithmetic(name))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(BinaryMiscibilityError::Cancelled)
    } else {
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Sum {
    main: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<()> {
        let next = finite(self.main + value, "compensated reduction")?;
        self.correction += if self.main.abs() >= value.abs() {
            (self.main - next) + value
        } else {
            (value - next) + self.main
        };
        self.main = next;
        Ok(())
    }
    fn value(self) -> Result<f64> {
        finite(
            self.main + self.correction,
            "compensated reduction correction",
        )
    }
}

#[cfg(test)]
#[path = "miscibility/tests.rs"]
mod tests;
