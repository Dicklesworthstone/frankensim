//! Finite storage in one-dimensional concentration films beside a membrane.
//!
//! Every film cell owns `area * width` of incompressible solvent and explicit
//! solute amounts. The SAME normal water flux crosses every interface: film
//! volumes stay fixed while the two finite bulk reservoirs lose/gain water.
//! Backward Euler balances every bulk and film solute together. At each water
//! trial, the original exponentially fitted transport coefficients form a
//! tridiagonal M-matrix, solved by the existing `fs_la` workspace in O(n) work.
//! The two terminal half cells and membrane use the original polarization law;
//! its reconstructed surface concentrations supply the osmotic driving force.
//!
//! This is an isothermal, neutral, ideal dilute film model with explicit constant
//! diffusivities and negligible film hydraulic loss. It does not infer crossflow,
//! turbulence, ionic fields, concentration-dependent diffusion, or heat storage.
//! Cell averages approximate center concentrations; spatial refinement and time
//! refinement are separate obligations. Each call accepts ONE first-order
//! implicit step, with bounded root work but no temporal error estimate or
//! uniqueness claim for permeating-solute roots. All inventories, including film
//! water, enter the balance report. Reported pressure work is not a closed heat
//! or chemical-energy balance. Failed/cancelled steps publish no state or time.

use super::{MAX_OSMOTIC_SPECIES, OsmoticCompartment, OsmoticNetworkError};
use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};
use fs_material::osmosis::polarization::{FilmPermeance, PolarizedMembraneTransport};
use fs_material::osmosis::{
    HydraulicPermeance, MembraneTransport, OsmosisError, ReflectionCoefficient, SolutePermeance,
};

/// Maximum finite cells on each side; all species share their geometry.
pub const MAX_TRANSIENT_FILM_CELLS: usize = 4096;
type Result<T> = core::result::Result<T, TransientOsmosisError>;

/// Original material/network refusals or a failed conservative linear solve.
#[derive(Debug)]
pub enum TransientOsmosisError {
    Network(OsmoticNetworkError),
    Linear(TridiagonalError),
    /// Residual of an original molar equation, before any clipping/correction.
    LinearResidual {
        species: usize,
        cell: usize,
        residual_mol: f64,
    },
    /// Total bulk-plus-film inventory changed beyond the requested allowance.
    InventoryBalance {
        species: usize,
        residual_mol: f64,
    },
}
impl core::fmt::Display for TransientOsmosisError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient concentration film: {self:?}")
    }
}
impl std::error::Error for TransientOsmosisError {}
impl From<OsmoticNetworkError> for TransientOsmosisError {
    fn from(value: OsmoticNetworkError) -> Self {
        Self::Network(value)
    }
}
impl From<OsmosisError> for TransientOsmosisError {
    fn from(value: OsmosisError) -> Self {
        Self::Network(value.into())
    }
}
impl From<TridiagonalError> for TransientOsmosisError {
    fn from(value: TridiagonalError) -> Self {
        if value == TridiagonalError::Cancelled {
            Self::Network(OsmoticNetworkError::Cancelled)
        } else {
            Self::Linear(value)
        }
    }
}
fn invalid(message: &'static str) -> TransientOsmosisError {
    OsmoticNetworkError::Invalid(message).into()
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("nonfinite film arithmetic"))
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(OsmoticNetworkError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn zeros(len: usize) -> Result<Box<[f64]>> {
    let mut data = Vec::new();
    data.try_reserve_exact(len)
        .map_err(|_| invalid("film workspace allocation failed"))?;
    data.resize(len, 0.0);
    Ok(data.into_boxed_slice())
}

/// Explicit uniform geometry, species diffusivities, and an initial profile.
/// Diffusivity and thickness are independent: a steady `D/thickness` coefficient
/// alone cannot specify either solvent storage or the diffusion time scale.
#[derive(Clone, Debug)]
pub struct TransientFilm {
    thickness_m: f64,
    cells: usize,
    diffusivity_m2_s: Box<[f64]>,
    initial_concentrations_mol_m3: Box<[f64]>,
}
impl TransientFilm {
    /// Concentrations are cell-major, with the original bulk species ordering
    /// within each cell. Cells on BOTH sides are ordered left to right.
    pub fn new(
        thickness_m: f64,
        cells: usize,
        diffusivity_m2_s: Vec<f64>,
        initial_concentrations_mol_m3: Vec<f64>,
    ) -> Result<Self> {
        if !positive(thickness_m)
            || cells == 0
            || cells > MAX_TRANSIENT_FILM_CELLS
            || diffusivity_m2_s.len() > MAX_OSMOTIC_SPECIES
            || initial_concentrations_mol_m3.len() != cells * diffusivity_m2_s.len()
            || diffusivity_m2_s.iter().any(|&d| !positive(d))
            || initial_concentrations_mol_m3
                .iter()
                .any(|&c| !nonnegative(c))
            || !positive(thickness_m / cells as f64)
        {
            return Err(invalid(
                "invalid finite film geometry, diffusivities or concentration layout",
            ));
        }
        Ok(Self {
            thickness_m,
            cells,
            diffusivity_m2_s: diffusivity_m2_s.into_boxed_slice(),
            initial_concentrations_mol_m3: initial_concentrations_mol_m3.into_boxed_slice(),
        })
    }
}

/// Numerical tolerances for ONE backward-Euler step, not adaptive time control.
#[derive(Clone, Copy, Debug)]
pub struct TransientStepControls {
    /// Relative solvent residual OR half-width of a sign-changing transfer bracket.
    pub relative_root_tolerance: f64,
    pub absolute_water_tolerance_m3: f64,
    pub relative_amount_tolerance: f64,
    pub absolute_amount_tolerance_mol: f64,
    /// Includes the zero-transfer trial and all bracketing/refused trials.
    pub maximum_root_evaluations: usize,
}
impl Default for TransientStepControls {
    fn default() -> Self {
        Self {
            relative_root_tolerance: 1e-9,
            absolute_water_tolerance_m3: 1e-18,
            relative_amount_tolerance: 1e-10,
            absolute_amount_tolerance_mol: 1e-18,
            maximum_root_evaluations: 128,
        }
    }
}

/// Accepted endpoint and the original bulk-plus-film balances, in SI units.
#[derive(Clone, Debug)]
pub struct TransientStepReport {
    pub interval_s: [f64; 2],
    pub solvent_transfer_m3: f64,
    /// Backward-Euler endpoint flux, positive left to right.
    pub volume_flux_m_s: f64,
    pub membrane_left_concentration_mol_m3: [f64; MAX_OSMOTIC_SPECIES],
    pub membrane_right_concentration_mol_m3: [f64; MAX_OSMOTIC_SPECIES],
    /// `(mean p_left - mean p_right) * transfer`, integrated exactly for the
    /// original linear bulk pressure/volume laws; may be negative in osmosis.
    pub pressure_work_into_solution_j: f64,
    pub water_balance_residual_m3: f64,
    pub solute_balance_residual_mol: [f64; MAX_OSMOTIC_SPECIES],
    pub species_count: usize,
    pub maximum_linear_residual_mol: f64,
    pub water_root_residual_m3: f64,
    pub root_evaluations: usize,
}

struct FilmGeometry {
    cells: usize,
    cell_volume_m3: f64,
    interior: MembraneTransport,
    outer_half_cell: MembraneTransport,
}
struct Workspace {
    trial_amounts: Box<[f64]>,
    lower: Box<[f64]>,
    diagonal: Box<[f64]>,
    upper: Box<[f64]>,
    rhs: Box<[f64]>,
    solution: Box<[f64]>,
    linear: TridiagonalWorkspace,
}
struct Trial {
    water: f64,
    residual: f64,
    faces: [[f64; MAX_OSMOTIC_SPECIES]; 2],
    balance: [f64; MAX_OSMOTIC_SPECIES],
    maximum_linear_residual: f64,
}

/// Two finite bulks and actual finite films. Prepared stepping allocates nothing.
/// Every species must be a neutral independent solute in a common ideal model;
/// nonideal solution owners are refused, rather than assigned scalar Fick laws.
pub struct TransientOsmoticPair {
    bulks: [OsmoticCompartment; 2],
    candidate_bulks: [OsmoticCompartment; 2],
    area_m2: f64,
    films: [Option<FilmGeometry>; 2],
    membrane: PolarizedMembraneTransport,
    film_amounts: Box<[f64]>,
    workspace: Workspace,
    species_count: usize,
    rows: usize,
    workspace_bytes: usize,
    time_s: f64,
}

// These are concentration-transfer operators only. Their zero hydraulic
// permeance is never sampled: the membrane supplies the SAME water flux to
// every interface; no artificial pressure/selectivity belongs to film cells.
fn film_operator(film: &TransientFilm, distance_m: f64) -> Result<MembraneTransport> {
    let solutes = film
        .diffusivity_m2_s
        .iter()
        .map(|&d| {
            let k = FilmPermeance::from_diffusion_layer(d, distance_m)?;
            Ok((
                ReflectionCoefficient::new(0.0)?,
                SolutePermeance::new(k.value())?,
            ))
        })
        .collect::<core::result::Result<Vec<_>, OsmosisError>>()?;
    Ok(MembraneTransport::new(
        HydraulicPermeance::new(0.0)?,
        solutes,
    )?)
}
fn prepare_film(film: &TransientFilm, area: f64) -> Result<FilmGeometry> {
    let dx = film.thickness_m / film.cells as f64;
    let volume = area * dx;
    if !positive(volume) {
        return Err(invalid("film cell water volume is not representable"));
    }
    Ok(FilmGeometry {
        cells: film.cells,
        cell_volume_m3: volume,
        interior: film_operator(film, dx)?,
        outer_half_cell: film_operator(film, 0.5 * dx)?,
    })
}

impl TransientOsmoticPair {
    /// Prepare all storage once. Initial film water/solute is supplied in
    /// addition to, never silently withdrawn from, the two bulk inventories.
    pub fn new(
        bulks: [OsmoticCompartment; 2],
        membrane_area_m2: f64,
        membrane: MembraneTransport,
        left_film: Option<TransientFilm>,
        right_film: Option<TransientFilm>,
        maximum_workspace_bytes: usize,
    ) -> Result<Self> {
        let count = membrane.solutes().len();
        let definitions = [left_film, right_film];
        let reference = bulks[0].solution();
        if !positive(membrane_area_m2)
            || count > MAX_OSMOTIC_SPECIES
            || bulks.iter().any(|b| {
                let s = b.solution();
                !s.is_ideal()
                    || s.solute_moles().len() != count
                    || s.temperature_k().to_bits() != reference.temperature_k().to_bits()
                    || s.solvent_molar_volume_m3_mol().to_bits()
                        != reference.solvent_molar_volume_m3_mol().to_bits()
            })
            || definitions
                .iter()
                .flatten()
                .any(|f| f.diffusivity_m2_s.len() != count)
        {
            return Err(invalid(
                "transient films require common neutral ideal species, solvent and temperature",
            ));
        }
        let film_cells: usize = definitions.iter().flatten().map(|f| f.cells).sum();
        let rows = film_cells + 2;
        // Conservative owned-heap bound: accepted/trial amounts, candidate
        // bulks, five matrix/RHS arrays, the original LU workspace, and all
        // boxed membrane/film coefficients. The initial supplied profile is
        // also covered while construction copies it into actual molar storage.
        let workspace_bytes = (3 * rows * count + 4 * count + 10 * rows) * size_of::<f64>()
            + rows * size_of::<bool>()
            + (5 * count) * size_of::<(ReflectionCoefficient, SolutePermeance)>()
            + 2 * count * size_of::<Option<FilmPermeance>>()
            + size_of::<Self>();
        if workspace_bytes > maximum_workspace_bytes {
            return Err(OsmoticNetworkError::Workspace {
                required_bytes: workspace_bytes,
                budget_bytes: maximum_workspace_bytes,
            }
            .into());
        }
        let half_films = |film: &Option<TransientFilm>| -> Result<Vec<Option<FilmPermeance>>> {
            match film {
                Some(f) => f
                    .diffusivity_m2_s
                    .iter()
                    .map(|&d| {
                        Ok(Some(FilmPermeance::from_diffusion_layer(
                            d,
                            0.5 * f.thickness_m / f.cells as f64,
                        )?))
                    })
                    .collect(),
                None => Ok(vec![None; count]),
            }
        };
        let membrane = PolarizedMembraneTransport::new(
            membrane,
            half_films(&definitions[0])?,
            half_films(&definitions[1])?,
        )?;
        let films = [
            definitions[0]
                .as_ref()
                .map(|f| prepare_film(f, membrane_area_m2))
                .transpose()?,
            definitions[1]
                .as_ref()
                .map(|f| prepare_film(f, membrane_area_m2))
                .transpose()?,
        ];
        let mut film_amounts = zeros(film_cells * count)?;
        let mut offset = 0;
        for side in 0..2 {
            if let (Some(definition), Some(geometry)) = (&definitions[side], &films[side]) {
                for &concentration in &definition.initial_concentrations_mol_m3 {
                    let amount = concentration * geometry.cell_volume_m3;
                    if !nonnegative(amount) || concentration > 0.0 && amount == 0.0 {
                        return Err(invalid("initial film amount is not representable"));
                    }
                    film_amounts[offset] = amount;
                    offset += 1;
                }
            }
        }
        let result = Self {
            candidate_bulks: bulks.clone(),
            bulks,
            area_m2: membrane_area_m2,
            films,
            membrane,
            film_amounts,
            workspace: Workspace {
                trial_amounts: zeros(rows * count)?,
                lower: zeros(rows - 1)?,
                diagonal: zeros(rows)?,
                upper: zeros(rows - 1)?,
                rhs: zeros(rows)?,
                solution: zeros(rows)?,
                linear: TridiagonalWorkspace::new(rows)?,
            },
            species_count: count,
            rows,
            workspace_bytes,
            time_s: 0.0,
        };
        finite(result.total_solvent_volume_m3())?;
        for species in 0..count {
            finite(result.total_amount(species))?;
        }
        Ok(result)
    }

    pub const fn bulks(&self) -> &[OsmoticCompartment; 2] {
        &self.bulks
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub const fn species_count(&self) -> usize {
        self.species_count
    }

    /// Side 0 is left, side 1 right; absent films have zero cells.
    pub fn film_cell_count(&self, side: usize) -> Option<usize> {
        self.films
            .get(side)
            .map(|f| f.as_ref().map_or(0, |f| f.cells))
    }
    /// Actual fixed water volume of one cell; None for an absent/invalid side.
    pub fn film_cell_volume_m3(&self, side: usize) -> Option<f64> {
        self.films.get(side)?.as_ref().map(|f| f.cell_volume_m3)
    }
    /// Actual cell inventory in the original species order. Both film profiles
    /// are indexed left to right, including the film on the right of the membrane.
    pub fn film_solute_moles(&self, side: usize, cell: usize) -> Option<&[f64]> {
        if cell >= self.films.get(side)?.as_ref()?.cells {
            return None;
        }
        let offset = if side == 0 { 0 } else { self.left_cells() };
        let begin = (offset + cell) * self.species_count;
        Some(&self.film_amounts[begin..begin + self.species_count])
    }
    pub fn total_solvent_volume_m3(&self) -> f64 {
        self.bulks
            .iter()
            .map(|b| b.solution().solvent_volume_m3())
            .sum::<f64>()
            + self
                .films
                .iter()
                .flatten()
                .map(|f| f.cells as f64 * f.cell_volume_m3)
                .sum::<f64>()
    }
    /// Include the two bulks AND every stored film cell in each species balance.
    pub fn total_solute_moles_into(&self, output: &mut [f64]) -> Result<()> {
        if output.len() != self.species_count {
            return Err(invalid("total amount output shape mismatch"));
        }
        for (i, amount) in output.iter_mut().enumerate() {
            *amount = self.total_amount(i);
        }
        Ok(())
    }
    fn left_cells(&self) -> usize {
        self.films[0].as_ref().map_or(0, |f| f.cells)
    }
    fn amount(&self, cell: usize, species: usize) -> f64 {
        if cell == 0 {
            self.bulks[0].solution().solute_moles()[species]
        } else if cell + 1 == self.rows {
            self.bulks[1].solution().solute_moles()[species]
        } else {
            self.film_amounts[(cell - 1) * self.species_count + species]
        }
    }
    fn total_amount(&self, species: usize) -> f64 {
        (0..self.rows).map(|row| self.amount(row, species)).sum()
    }
    fn volume(&self, cell: usize, water: f64) -> f64 {
        if cell == 0 {
            self.bulks[0].solution().solvent_volume_m3() - water
        } else if cell + 1 == self.rows {
            self.bulks[1].solution().solvent_volume_m3() + water
        } else if cell <= self.left_cells() {
            self.films[0]
                .as_ref()
                .expect("left film cell")
                .cell_volume_m3
        } else {
            self.films[1]
                .as_ref()
                .expect("right film cell")
                .cell_volume_m3
        }
    }

    fn evaluate(
        &mut self,
        area_time: f64,
        water: f64,
        controls: TransientStepControls,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        poll(cancelled)?;
        if *evaluations >= controls.maximum_root_evaluations {
            return Err(OsmoticNetworkError::NonlinearConvergence.into());
        }
        *evaluations += 1;
        let vl = self.volume(0, water);
        let vr = self.volume(self.rows - 1, water);
        if !positive(vl)
            || !positive(vr)
            || vl < self.bulks[0].minimum_volume_m3
            || vr < self.bulks[1].minimum_volume_m3
        {
            return Err(OsmoticNetworkError::Unbracketed.into());
        }
        let flux = finite(water / area_time)?;
        let middle = self.left_cells();
        let mut trial = Trial {
            water,
            residual: 0.0,
            faces: [[0.0; MAX_OSMOTIC_SPECIES]; 2],
            balance: [0.0; MAX_OSMOTIC_SPECIES],
            maximum_linear_residual: 0.0,
        };
        for species in 0..self.species_count {
            poll(cancelled)?;
            let membrane = self.membrane.coefficients(species, flux)?;
            for row in 0..self.rows {
                if row % 128 == 0 {
                    poll(cancelled)?;
                }
                self.workspace.diagonal[row] = 1.0;
                self.workspace.rhs[row] = self.amount(row, species);
            }
            for edge in 0..self.rows - 1 {
                if edge % 128 == 0 {
                    poll(cancelled)?;
                }
                let (kl, kr) = if edge == middle {
                    membrane.donor_coefficients()
                } else {
                    let side = usize::from(edge > middle);
                    let film = self.films[side].as_ref().expect("noncentral film edge");
                    let operator = if edge == 0 || edge + 2 == self.rows {
                        &film.outer_half_cell
                    } else {
                        &film.interior
                    };
                    operator.solute_stream_coefficients(species, flux)?
                };
                let a = area_time * kl / self.volume(edge, water);
                let b = area_time * kr / self.volume(edge + 1, water);
                if !nonnegative(a) || !nonnegative(b) {
                    return Err(invalid("film implicit coefficient overflow"));
                }
                // Solve for AMOUNTS. Each shared face contributes equal and
                // opposite rates; every column is strictly diagonally dominant.
                self.workspace.diagonal[edge] = finite(self.workspace.diagonal[edge] + a)?;
                self.workspace.diagonal[edge + 1] = finite(self.workspace.diagonal[edge + 1] + b)?;
                self.workspace.lower[edge] = -a;
                self.workspace.upper[edge] = -b;
            }
            let work = &mut self.workspace;
            work.linear.solve_into(
                &work.lower,
                &work.diagonal,
                &work.upper,
                &work.rhs,
                &mut work.solution,
                &mut *cancelled,
            )?;
            let mut before = 0.0;
            let mut after = 0.0;
            for row in 0..self.rows {
                if row % 128 == 0 {
                    poll(cancelled)?;
                }
                let value = work.solution[row];
                if !nonnegative(value) {
                    return Err(invalid("negative or nonfinite implicit film amount"));
                }
                let mut residual = work.diagonal[row].mul_add(value, -work.rhs[row]);
                if row > 0 {
                    residual = work.lower[row - 1].mul_add(work.solution[row - 1], residual);
                }
                if row + 1 < self.rows {
                    residual = work.upper[row].mul_add(work.solution[row + 1], residual);
                }
                let tolerance = controls.absolute_amount_tolerance_mol
                    + controls.relative_amount_tolerance * work.rhs[row].max(value);
                if !residual.is_finite() || residual.abs() > tolerance {
                    return Err(TransientOsmosisError::LinearResidual {
                        species,
                        cell: row,
                        residual_mol: residual,
                    });
                }
                trial.maximum_linear_residual = trial.maximum_linear_residual.max(residual.abs());
                before += work.rhs[row];
                after += value;
                work.trial_amounts[row * self.species_count + species] = value;
            }
            trial.balance[species] = after - before;
            if !trial.balance[species].is_finite()
                || trial.balance[species].abs()
                    > controls.absolute_amount_tolerance_mol
                        + controls.relative_amount_tolerance * before
            {
                return Err(TransientOsmosisError::InventoryBalance {
                    species,
                    residual_mol: trial.balance[species],
                });
            }
            let left = self.workspace.solution[middle] / self.volume(middle, water);
            let right = self.workspace.solution[middle + 1] / self.volume(middle + 1, water);
            (trial.faces[0][species], trial.faces[1][species]) =
                membrane.membrane_face_concentrations(left, right)?;
        }
        let predicted = self.membrane.membrane().water_flux(
            self.bulks[0].solution().temperature_k(),
            self.bulks[0].pressure.pressure_pa(vl)?,
            self.bulks[1].pressure.pressure_pa(vr)?,
            &trial.faces[0][..self.species_count],
            &trial.faces[1][..self.species_count],
        )?;
        trial.residual = finite(water - area_time * predicted.volume_flux_m_s)?;
        Ok(trial)
    }

    fn solve_water(
        &mut self,
        area_time: f64,
        controls: TransientStepControls,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial> {
        let converged = |trial: &Trial| {
            trial.residual.abs()
                <= controls.absolute_water_tolerance_m3
                    + controls.relative_root_tolerance * trial.water.abs()
        };
        let initial = self.evaluate(area_time, 0.0, controls, evaluations, cancelled)?;
        if converged(&initial) {
            return Ok(initial);
        }
        // Follow the initial flow direction from zero. Remote volume-floor
        // trials can overflow enrichment although a nearby physical root exists.
        let direction = if initial.residual < 0.0 { 1.0 } else { -1.0 };
        let donor = usize::from(direction < 0.0);
        let available = (self.bulks[donor].solution().solvent_volume_m3()
            - self.bulks[donor].minimum_volume_m3)
            .next_down();
        let mut distance = initial.residual.abs().min(0.125 * available);
        let mut previous = 0.0;
        let mut previous_residual = initial.residual;
        let mut refused = None;
        let (mut lo, mut hi, mut flo, mut fhi) = loop {
            if *evaluations >= controls.maximum_root_evaluations {
                return Err(OsmoticNetworkError::NonlinearConvergence.into());
            }
            if distance <= previous || distance > available {
                return Err(OsmoticNetworkError::Unbracketed.into());
            }
            match self.evaluate(
                area_time,
                direction * distance,
                controls,
                evaluations,
                cancelled,
            ) {
                Ok(trial) => {
                    if converged(&trial) {
                        return Ok(trial);
                    }
                    if direction * trial.residual >= 0.0 {
                        break if direction > 0.0 {
                            (previous, distance, previous_residual, trial.residual)
                        } else {
                            (-distance, -previous, trial.residual, previous_residual)
                        };
                    }
                    previous = distance;
                    previous_residual = trial.residual;
                    distance = refused.map_or_else(
                        || (2.0 * distance).min(available),
                        |limit| distance + 0.5 * (limit - distance),
                    );
                }
                Err(TransientOsmosisError::Network(OsmoticNetworkError::Cancelled)) => {
                    return Err(OsmoticNetworkError::Cancelled.into());
                }
                Err(_) => {
                    refused = Some(distance);
                    distance = previous + 0.5 * (distance - previous);
                }
            }
        };
        loop {
            let midpoint = lo + 0.5 * (hi - lo);
            // A steep osmotic force can amplify roundoff in the evaluated
            // residual. A narrow sign-changing bracket instead bounds the
            // transfer's location. Reevaluate it to stage the matching amounts.
            if hi - lo
                <= 2.0
                    * (controls.absolute_water_tolerance_m3
                        + controls.relative_root_tolerance * midpoint.abs())
            {
                return self.evaluate(area_time, midpoint, controls, evaluations, cancelled);
            }
            let fraction = -flo / (fhi - flo);
            let fraction = if fraction.is_finite() {
                fraction.clamp(0.1, 0.9)
            } else {
                0.5
            };
            let water = lo + fraction * (hi - lo);
            if water == lo || water == hi {
                return Err(OsmoticNetworkError::NonlinearConvergence.into());
            }
            let trial = self.evaluate(area_time, water, controls, evaluations, cancelled)?;
            if converged(&trial) {
                return Ok(trial);
            }
            if trial.residual < 0.0 {
                lo = water;
                flo = trial.residual;
            } else {
                hi = water;
                fhi = trial.residual;
            }
        }
    }

    /// Accept one fully implicit step. Water tolerances control its residual or
    /// sign-changing root bracket; amount tolerances check the original linear
    /// equations and total inventories. Neither estimates temporal truncation
    /// error. Refine dt separately; the actual water residual is always reported.
    pub fn step(
        &mut self,
        duration_s: f64,
        controls: TransientStepControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransientStepReport> {
        poll(&mut cancelled)?;
        if !positive(duration_s)
            || !positive(controls.relative_root_tolerance)
            || controls.relative_root_tolerance >= 1.0
            || !positive(controls.absolute_water_tolerance_m3)
            || !positive(controls.relative_amount_tolerance)
            || controls.relative_amount_tolerance >= 1.0
            || !positive(controls.absolute_amount_tolerance_mol)
            || !(2..=4096).contains(&controls.maximum_root_evaluations)
        {
            return Err(invalid(
                "invalid transient film duration or root/balance controls",
            ));
        }
        let end = finite(self.time_s + duration_s)?;
        let area_time = self.area_m2 * duration_s;
        if end <= self.time_s || !positive(area_time) {
            return Err(OsmoticNetworkError::StepUnderflow.into());
        }
        let before_water = self.total_solvent_volume_m3();
        let mut evaluations = 0;
        let trial = self.solve_water(area_time, controls, &mut evaluations, &mut cancelled)?;
        let vl = self.volume(0, trial.water);
        let vr = self.volume(self.rows - 1, trial.water);
        let count = self.species_count;
        self.candidate_bulks[0]
            .solution
            .replace_inventories(vl, &self.workspace.trial_amounts[..count])?;
        self.candidate_bulks[1]
            .solution
            .replace_inventories(vr, &self.workspace.trial_amounts[(self.rows - 1) * count..])?;
        let pressure_before = self.bulks[0].pressure_pa()? - self.bulks[1].pressure_pa()?;
        let pressure_after =
            self.candidate_bulks[0].pressure_pa()? - self.candidate_bulks[1].pressure_pa()?;
        let work = finite((0.5 * pressure_before + 0.5 * pressure_after) * trial.water)?;
        let film_water: f64 = self
            .films
            .iter()
            .flatten()
            .map(|f| f.cell_volume_m3 * f.cells as f64)
            .sum();
        let water_residual = finite((vl + vr + film_water) - before_water)?;
        if water_residual.abs()
            > controls.absolute_water_tolerance_m3 + 32.0 * f64::EPSILON * before_water
        {
            return Err(invalid("bulk-plus-film water balance failed"));
        }
        let report = TransientStepReport {
            interval_s: [self.time_s, end],
            solvent_transfer_m3: trial.water,
            volume_flux_m_s: trial.water / area_time,
            membrane_left_concentration_mol_m3: trial.faces[0],
            membrane_right_concentration_mol_m3: trial.faces[1],
            pressure_work_into_solution_j: work,
            water_balance_residual_m3: water_residual,
            solute_balance_residual_mol: trial.balance,
            species_count: count,
            maximum_linear_residual_mol: trial.maximum_linear_residual,
            water_root_residual_m3: trial.residual,
            root_evaluations: evaluations,
        };
        poll(&mut cancelled)?;
        // The sole publication boundary: all source, root, linear and balance
        // refusals have completed. Prepared copies cannot fail or allocate.
        self.film_amounts
            .copy_from_slice(&self.workspace.trial_amounts[count..(self.rows - 1) * count]);
        core::mem::swap(&mut self.bulks, &mut self.candidate_bulks);
        self.time_s = end;
        Ok(report)
    }
}

#[cfg(test)]
#[path = "transient/tests.rs"]
mod tests;
