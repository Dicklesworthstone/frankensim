//! History-dependent gravity sag of a slender, already spanning deposited road.
//!
//! Reuses `deposition::LayeredStrip` for every axial cell's material histories,
//! thermal eigenstrain, implicit section compliance and energy audit. Only the
//! statically determinate beam kinematics/load transfer live here. Each cell has
//! constant curvature and axial strain; integrating twice gives a C1 piecewise
//! quadratic displacement. The EXACT cell-average equilibrium moment is the
//! work-conjugate load, not the midpoint moment. Thus section work equals the
//! integral of the declared piecewise-constant dead load against displacement.
//!
//! Local frame: x along the horizontal reference span, z upward. Supports are
//! a left pin/right roller or left clamp/free right end. The ordinary path has zero axial resultant; `drawing` adds explicitly declared
//! tensile drawing force or explicit slack with small-deflection geometric stiffening.
//! No compression/buckling, inertia, shear deformation,
//! adhesion, fracture or molten-strand flow is inferred. Explicit frictionless
//! point supports are available separately through `supports`.
//! A complete initial seed must be supplied. `growth` extends a clamped seed
//! at its current tip without activating future material. Complete layers may be born
//! stress-free on the current shape only with an explicit perfect-bond assumption.
//! This is bounded issue #1 R5 functionality, not whole-print qualification.

use std::collections::{BTreeMap, BTreeSet};
use crate::deposition::{self, AxialControl, BendingControl, Control, LayerGeometry, LayeredStrip};
use fs_exec::Cx;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::thermovisco::{Energy, ThermoMaxwell};
use fs_material::thermovisco::source::ResolvedThermoMaxwell;

/// Later-pass footprint capture on the actual deformed surface.
pub mod contact;
/// Solved unilateral support reactions on the same history-dependent beam.
pub mod supports;
/// Equal-and-opposite contact forces between deformable neighbouring roads.
pub mod neighbours;
/// Stress-free cell birth at the advancing tip of a clamped road.
pub mod growth;
/// Declared drawing tension/slack with beam-column geometric stiffness.
pub mod drawing;

/// Fixed transverse restraints; both allow free axial thermal contraction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Support {
    /// w(0)=w(L)=0; end rotations are free, no end moments.
    SimplySupported,
    /// w(0)=w'(0)=0; right end is traction-free.
    CantileverLeft,
}

/// Retained-state limits, including TOTAL Maxwell branch states across cells.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Budget {
    pub max_cells: usize,
    pub max_layers_per_cell: usize,
    pub max_branch_states: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct Config {
    /// Strain, rotation, force/moment and whole-beam joule gates. The energy
    /// allowance is divided between sections before their individual updates.
    pub section: deposition::Config,
    /// Actual maximum |w| allowed, metres; not a measured discretization bound.
    pub maximum_deflection_m: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Invalid(&'static str),
    Budget,
    Cancelled,
    Section { cell: u64, error: deposition::Error },
    Regime(&'static str),
    Energy(f64),
}
impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "deposited bridge refused: {self:?}")
    }
}
impl std::error::Error for Error {}

/// One axial reference interval and its physical section history.
#[derive(Debug, Clone, PartialEq)]
pub struct Cell<'a> {
    id: u64,
    start_m: f64,
    section: LayeredStrip<'a>,
    mass_per_length_kg_m: f64,
    // Accepted left-end reference-line displacement, rotation and axial motion.
    w_m: f64,
    slope: f64,
    u_m: f64,
}
impl<'a> Cell<'a> {
    pub fn id(&self) -> u64 { self.id }
    pub fn start_m(&self) -> f64 { self.start_m }
    pub fn length_m(&self) -> f64 { self.section.length_m() }
    pub fn section(&self) -> &LayeredStrip<'a> { &self.section }
    pub fn mass_kg(&self) -> f64 { self.mass_per_length_kg_m * self.length_m() }
    fn height(&self, s: f64) -> f64 {
        self.w_m + s * (self.slope + 0.5 * s * self.section.curvature_per_m())
    }
    fn integral_height(&self) -> f64 {
        let h = self.length_m();
        h * (self.w_m + 0.5 * h * self.slope + h*h*self.section.curvature_per_m()/6.0)
    }
}

/// Immutable accepted shape and all material histories. Clone is a checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct Bridge<'a> {
    cells: Vec<Cell<'a>>,
    length_m: f64,
    time_s: f64,
    support: Support,
    budget: Budget,
    branch_states: usize,
}

#[derive(Debug, Clone)]
pub struct Report {
    /// Upward vertical reactions at left and right, N (signed for added loads).
    pub vertical_reactions_n: [f64; 2],
    /// Counteracting clamp reaction, N m; zero for pin/roller support.
    pub clamp_reaction_nm: f64,
    /// Work-conjugate cell-average moments, same order as `Bridge::cells`.
    pub moments_nm: Vec<f64>,
    pub maximum_sag_m: f64,
    pub maximum_absolute_deflection_m: f64,
    pub maximum_absolute_rotation_rad: f64,
    pub material_point_updates: usize,
    /// Extensive material accounts from the actual section updates, J.
    pub energy_j: Energy,
    /// Integral of endpoint dead load against actual displacement increment, J.
    /// Endpoint quadrature is not exact time-integrated actuator work.
    pub load_work_j: f64,
    pub load_work_residual_j: f64,
}

impl<'a> Bridge<'a> {
    /// Ordered `(cell_id, reference_length_m)` intervals tile the complete span.
    /// Empty sections carry neither weight nor support until a layer is born.
    pub fn new(cells: &[(u64, f64)], support: Support, budget: Budget) -> Result<Self, Error> {
        if cells.is_empty() || cells.len() > budget.max_cells || budget.max_layers_per_cell == 0 {
            return Err(Error::Budget);
        }
        let mut ids = BTreeSet::new();
        let mut output = Vec::with_capacity(cells.len());
        let mut length_m = 0.0;
        for &(id, length) in cells {
            positive(length)?;
            if !ids.insert(id) { return Err(Error::Invalid("duplicate axial cell ID")); }
            let section = LayeredStrip::new(length, deposition::Budget {
                max_layers: budget.max_layers_per_cell, max_branch_states: budget.max_branch_states,
            }).map_err(|error| Error::Section { cell: id, error })?;
            output.push(Cell { id, start_m: length_m, section, mass_per_length_kg_m: 0.0,
                w_m: 0.0, slope: 0.0, u_m: 0.0 });
            let end = finite(length_m + length)?;
            if end <= length_m { return Err(Error::Invalid("unrepresentable axial interval")); }
            length_m = end;
        }
        Ok(Self { cells: output, length_m, time_s: 0.0, support, budget, branch_states: 0 })
    }
    pub fn cells(&self) -> &[Cell<'a>] { &self.cells }
    pub fn length_m(&self) -> f64 { self.length_m }
    pub fn time_s(&self) -> f64 { self.time_s }
    pub fn mass_kg(&self) -> f64 { self.cells.iter().map(Cell::mass_kg).sum() }
    pub fn stored_energy_j(&self) -> f64 {
        self.cells.iter().map(|cell| cell.section.stored_energy_j()).sum()
    }

    /// Deposit a COMPLETE rectangular road on the existing section geometry.
    /// Perfect bonding is explicitly selected by this call, never inferred from
    /// thermal contact. New points use the current total strain and curvature as
    /// their stress-free birth reference; old points and shape are unchanged.
    /// One birth temperature per axial cell is required. All lengths/areas stay
    /// reference quantities; this does not solve path laying or bead coalescence.
    pub fn deposit_bonded_layer(&self, layer: u64, geometry: LayerGeometry,
        density_kg_m3: f64, law: &'a ThermoMaxwell, temperatures_k: &BTreeMap<u64, f64>,
        cancelled: &mut impl FnMut() -> bool) -> Result<Self, Error> {
        poll(cancelled)?;
        positive(density_kg_m3)?;
        if temperatures_k.len() != self.cells.len() { return Err(Error::Invalid("one birth temperature per cell required")); }
        let extra = self.cells.len().checked_mul(2).and_then(|n| n.checked_mul(law.branch_count()))
            .ok_or(Error::Budget)?;
        let branches = self.branch_states.checked_add(extra).ok_or(Error::Budget)?;
        if branches > self.budget.max_branch_states { return Err(Error::Budget); }
        let mass_per_length = finite(density_kg_m3 * geometry.area_m2)?;
        positive(mass_per_length)?;
        let lo = finite(geometry.centroid_z_m - 0.5 * geometry.thickness_m)?;
        let hi = finite(geometry.centroid_z_m + 0.5 * geometry.thickness_m)?;
        if hi <= lo { return Err(Error::Invalid("nonpositive/unrepresentable layer thickness")); }
        // A beam theory requires slender geometry; this is an admission screen,
        // not evidence that shear deformation is negligible for a given material.
        let width = geometry.area_m2 / geometry.thickness_m;
        positive(width)?;
        let mut bottom = lo;
        let mut top = hi;
        let mut touches = self.cells[0].section.layers().is_empty();
        for previous in self.cells[0].section.layers().values() {
            let a = previous.geometry.centroid_z_m - 0.5 * previous.geometry.thickness_m;
            let b = previous.geometry.centroid_z_m + 0.5 * previous.geometry.thickness_m;
            if lo < b && hi > a { return Err(Error::Invalid("overlapping layer cross sections")); }
            touches |= lo == b || hi == a;
            bottom = bottom.min(a); top = top.max(b);
        }
        if !touches { return Err(Error::Invalid("perfectly bonded layers must share a section face")); }
        if (top-bottom).max(width) > self.length_m/10.0 {
            return Err(Error::Regime("beam span must be at least ten times section depth and width")); }
        let mut next = self.clone();
        for cell in &mut next.cells {
            poll(cancelled)?;
            let temperature = *temperatures_k.get(&cell.id).ok_or(Error::Invalid("missing birth temperature"))?;
            cell.section = cell.section.deposit(layer, geometry, law, temperature)
                .map_err(|error| Error::Section { cell: cell.id, error })?;
            cell.mass_per_length_kg_m = finite(cell.mass_per_length_kg_m + mass_per_length)?;
        }
        finite(next.mass_kg())?;
        next.branch_states = branches;
        poll(cancelled)?;
        Ok(next)
    }

    /// Deposit a complete bonded layer whose viscoelastic law AND reference
    /// mass density come from one source-compatible material bundle.
    ///
    /// This is the source-bound counterpart of `deposit_bonded_layer`: geometry
    /// and birth temperatures remain specimen inputs, while density is taken from
    /// the same-card caloric bundle accepted by `ResolvedThermoMaxwell`. A material
    /// swap therefore changes inertia/weight and constitutive history together
    /// instead of permitting an unrelated hand-authored density to survive.
    ///
    /// The caloric chart is used here only as source authority for reference
    /// density; no heat equation or phase evolution is performed by this method.
    pub fn deposit_bonded_material_layer(&self, layer: u64, geometry: LayerGeometry,
        source: &'a ResolvedThermoMaxwell, caloric: &ResolvedHeatCapacityCurve,
        temperatures_k: &BTreeMap<u64, f64>, cancelled: &mut impl FnMut() -> bool)
        -> Result<Self, Error> {
        poll(cancelled)?;
        let density = source.reference_density_kg_m3(caloric)
            .map_err(|_| Error::Invalid("mechanical law and reference density must share one admitted material source"))?;
        let next = self.deposit_bonded_layer(layer, geometry, density, source.law(),
            temperatures_k, cancelled)?;
        poll(cancelled)?;
        Ok(next)
    }

    /// Physical coordinates (x,z) of a reference point. z is measured from the
    /// section datum; x is the REFERENCE coordinate, not an inverse world query.
    /// Uses the same strain/curvature as stress: X=x+u-z*theta, Z=z+w.
    pub fn deformed_position_m(&self, x: f64, z: f64) -> Result<[f64; 2], Error> {
        let (cell, s) = self.locate(x)?;
        finite(z)?;
        let theta = cell.slope + cell.section.curvature_per_m()*s;
        Ok([finite(x + cell.u_m + cell.section.axial_strain()*s - z*theta)?,
            finite(z + cell.height(s))?])
    }
    fn locate(&self, x: f64) -> Result<(&Cell<'a>, f64), Error> {
        if !x.is_finite() || x < 0.0 || x > self.length_m || self.cells[0].section.layers().is_empty() {
            return Err(Error::Invalid("query needs a born span and an in-range reference position")); }
        let index = self.cells.partition_point(|cell| cell.start_m <= x).saturating_sub(1);
        let cell = &self.cells[index];
        Ok((cell, x-cell.start_m))
    }

    /// One quasistatic, backward-Euler history interval. Temperatures are keyed
    /// by `(axial_cell_id, layer_id)`. `additional_downward_n_m` has one signed
    /// dead line load per axial cell. Gravity acts on every born reference mass.
    /// No born material, guessed temperature, or partial successful history can
    /// leak from a refused call. Work is linear in cells for a fixed section size;
    /// no dense global matrix. Refine spatial cells and thermal/time steps separately.
    pub fn advance_controlled(&self, target_s: f64, temperatures_k: &BTreeMap<(u64, u64), f64>,
        gravity_m_s2: f64, additional_downward_n_m: &[f64], config: Config,
        cancelled: &mut impl FnMut() -> bool) -> Result<(Self, Report), Error> {
        self.advance_with_point_forces(target_s, temperatures_k, gravity_m_s2,
            additional_downward_n_m, &[], config, cancelled)
    }

    // Point forces are upward N at REFERENCE x. They come from the checked
    // support solve, and enter both the moment diagram and displacement work.
    fn advance_with_point_forces(&self, target_s: f64, temperatures_k: &BTreeMap<(u64, u64), f64>,
        gravity_m_s2: f64, additional_downward_n_m: &[f64], point_forces: &[(f64, f64)],
        config: Config, cancelled: &mut impl FnMut() -> bool) -> Result<(Self, Report), Error> {
        poll(cancelled)?;
        positive(finite(target_s-self.time_s)?)?;
        if !gravity_m_s2.is_finite() || gravity_m_s2 < 0.0 { return Err(Error::Invalid("gravity must be finite and nonnegative")); }
        positive(config.maximum_deflection_m)?;
        let count = self.cells.len();
        let layers = self.cells[0].section.layers().len();
        let expected = count.checked_mul(layers).ok_or(Error::Budget)?;
        if layers == 0 || temperatures_k.len() != expected || additional_downward_n_m.len() != count {
            return Err(Error::Invalid("born span, complete temperatures and one line load per cell required")); }
        let load = self.load_field(gravity_m_s2, additional_downward_n_m, point_forces, cancelled)?;
        let mut next = self.clone();
        let mut energy = Energy::default();
        let mut updates = 0usize;
        let mut section_config = config.section;
        section_config.energy_tolerance_j /= count as f64;
        for (cell, &moment) in next.cells.iter_mut().zip(&load.moments) {
            poll(cancelled)?;
            let mut temperatures = BTreeMap::new();
            for &layer in cell.section.layers().keys() {
                let t = *temperatures_k.get(&(cell.id, layer)).ok_or(Error::Invalid("missing section temperature"))?;
                temperatures.insert(layer, t);
            }
            let control = Control { axial: AxialControl::ForceN(0.0), bending: BendingControl::MomentNm(moment) };
            let (section, report) = cell.section.advance_to(target_s, &temperatures, control, section_config,
                || if cancelled() { Err(deposition::Error::Cancelled) } else { Ok(()) })
                .map_err(|error| if error == deposition::Error::Cancelled { Error::Cancelled }
                    else { Error::Section { cell: cell.id, error } })?;
            cell.section = section;
            add_energy(&mut energy, report.energy_j)?;
            updates = updates.checked_add(report.material_point_updates).ok_or(Error::Budget)?;
        }
        // Integrate the section curvatures exactly on their piecewise-constant
        // representation. The pin/roller integration constant enforces w(L)=0.
        let mut slope = 0.0;
        if self.support == Support::SimplySupported {
            for cell in &next.cells {
                poll(cancelled)?;
                let h = cell.length_m();
                slope = finite(slope - cell.section.curvature_per_m()*h
                    *((self.length_m-cell.start_m-0.5*h)/self.length_m))?;
            }
        }
        let mut w = 0.0;
        let mut u = 0.0;
        let mut max_abs: f64 = 0.0;
        let mut sag: f64 = 0.0;
        let mut rotation: f64 = 0.0;
        let mut work = 0.0;
        for ((cell, old), &q) in next.cells.iter_mut().zip(&self.cells).zip(&load.distributed) {
            poll(cancelled)?;
            let h = cell.length_m();
            let kappa = cell.section.curvature_per_m();
            cell.w_m = w; cell.slope = slope; cell.u_m = u;
            let end = finite(cell.height(h))?;
            for value in [w, end] { max_abs = max_abs.max(value.abs()); sag = sag.max(-value); }
            if kappa != 0.0 {
                let critical = -slope/kappa;
                if critical > 0.0 && critical < h {
                    let value = finite(cell.height(critical))?;
                    max_abs = max_abs.max(value.abs()); sag = sag.max(-value);
                }
            }
            rotation = rotation.max(slope.abs());
            slope = finite(slope+kappa*h)?;
            rotation = rotation.max(slope.abs());
            u = finite(u+cell.section.axial_strain()*h)?;
            w = end;
            work = finite(work-q*(cell.integral_height()-old.integral_height()))?;
        }
        for &(x, force) in point_forces {
            poll(cancelled)?;
            let (old, s_old) = self.locate(x)?;
            let (new, s_new) = next.locate(x)?;
            work = finite(work + force*(new.height(s_new)-old.height(s_old)))?;
        }
        if max_abs > config.maximum_deflection_m || rotation > config.section.max_rotation_rad {
            return Err(Error::Regime("integrated shape exceeds the declared deflection/rotation range")); }
        next.time_s = target_s;
        let work_defect = finite(work-energy.mechanical_work)?;
        let storage = finite(next.stored_energy_j()-self.stored_energy_j())?;
        let balance = finite(storage+energy.dissipated+energy.numerical_dissipation
            -work-energy.thermal_strain_work-energy.modulus_work)?;
        for defect in [work_defect, balance, finite(storage-energy.stored_change)?] {
            if defect.abs() > config.section.energy_tolerance_j { return Err(Error::Energy(defect)); }
        }
        energy.residual = balance;
        poll(cancelled)?;
        Ok((next, Report { vertical_reactions_n: load.reactions, clamp_reaction_nm: load.clamp,
            moments_nm: load.moments, maximum_sag_m: sag, maximum_absolute_deflection_m: max_abs,
            maximum_absolute_rotation_rad: rotation, material_point_updates: updates,
            energy_j: energy, load_work_j: work, load_work_residual_j: work_defect }))
    }

    fn load_field(&self, gravity: f64, extra: &[f64], point_forces: &[(f64, f64)],
        cancelled: &mut impl FnMut() -> bool) -> Result<LoadField, Error> {
        if !gravity.is_finite() || gravity < 0.0 || extra.len() != self.cells.len() {
            return Err(Error::Invalid("invalid gravity or incomplete cell loads"));
        }
        let mut distributed = Vec::with_capacity(self.cells.len());
        let mut weight = 0.0;
        let mut first_moment = 0.0;
        for (cell, &added) in self.cells.iter().zip(extra) {
            poll(cancelled)?;
            let q = finite(gravity*cell.mass_per_length_kg_m + finite(added)?)?;
            let force = finite(q*cell.length_m())?;
            weight = finite(weight+force)?;
            first_moment = finite(first_moment+force*(cell.start_m+0.5*cell.length_m()))?;
            distributed.push(q);
        }
        for &(x, force) in point_forces {
            poll(cancelled)?;
            if !x.is_finite() || !(0.0..=self.length_m).contains(&x) {
                return Err(Error::Invalid("point force outside reference span"));
            }
            weight = finite(weight-finite(force)?)?;
            first_moment = finite(first_moment-force*x)?;
        }
        let right = match self.support { Support::SimplySupported => first_moment/self.length_m, Support::CantileverLeft => 0.0 };
        let left = finite(weight-right)?;
        let clamp = match self.support { Support::SimplySupported => 0.0, Support::CantileverLeft => first_moment };
        let mut passed_force = 0.0;
        let mut passed_moment = 0.0;
        let mut moments = Vec::with_capacity(self.cells.len());
        for (cell, &q) in self.cells.iter().zip(&distributed) {
            poll(cancelled)?;
            let h = cell.length_m();
            let x = cell.start_m;
            let mut moment = finite(-clamp+(left-passed_force)*(x+0.5*h)+passed_moment-q*h*h/6.0)?;
            for &(position, force) in point_forces {
                // Integral of the positive-part ramp / cell length. This
                // factorization avoids subtracting two nearly equal squares.
                let a = (x-position).max(0.0);
                let b = (x+h-position).max(0.0);
                moment = finite(moment+force*(0.5*((b-a)/h)*(a+b)))?;
            }
            moments.push(moment);
            passed_force = finite(passed_force+q*h)?;
            passed_moment = finite(passed_moment+q*h*(x+0.5*h))?;
        }
        Ok(LoadField { distributed, moments, reactions: [left, right], clamp })
    }

    pub fn advance_to(&self, cx: &Cx<'_>, target_s: f64, temperatures_k: &BTreeMap<(u64, u64), f64>,
        gravity_m_s2: f64, additional_downward_n_m: &[f64], config: Config) -> Result<(Self, Report), Error> {
        self.advance_controlled(target_s, temperatures_k, gravity_m_s2, additional_downward_n_m,
            config, &mut || cx.checkpoint().is_err())
    }
}

struct LoadField {
    distributed: Vec<f64>,
    moments: Vec<f64>,
    reactions: [f64; 2],
    clamp: f64,
}

fn add_energy(sum: &mut Energy, part: Energy) -> Result<(), Error> {
    for (target, value) in [(&mut sum.mechanical_work, part.mechanical_work),
        (&mut sum.thermal_strain_work, part.thermal_strain_work), (&mut sum.modulus_work, part.modulus_work),
        (&mut sum.stored_change, part.stored_change), (&mut sum.dissipated, part.dissipated),
        (&mut sum.numerical_dissipation, part.numerical_dissipation), (&mut sum.residual, part.residual)] {
        *target = finite(*target+value)?;
    }
    Ok(())
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() { Ok(value) } else { Err(Error::Invalid("nonfinite beam arithmetic")) }
}
fn positive(value: f64) -> Result<(), Error> {
    if value.is_finite() && value > 0.0 { Ok(()) } else { Err(Error::Invalid("positive finite quantity required")) }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}

#[cfg(test)]
mod tests;
