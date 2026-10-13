//! History-dependent, perfectly bonded deposited-strip mechanics (bounded R2).
//!
//! Uses the existing fiber-section convention `epsilon(z)=epsilon0-z*kappa`,
//! `N=sum(A*sigma)`, `M=-sum(A*z*sigma)`. Each rectangular layer has two Gauss
//! points, so axial/bending stiffness and quadratic stored energy are integrated
//! exactly for the represented uniform-layer law. `fs-material::thermovisco`
//! owns temperature response, branch history, consistent tangent and work audit.
//!
//! This is an axially UNIFORM, uniaxial Euler-Bernoulli section over a declared
//! reference length. Layers span that entire length and are explicitly assumed
//! perfectly bonded; a thermal contact does not establish that assumption.
//! Geometry, slenderness and nonoverlapping layer volumes are caller obligations.
//! No plate/torsion/shear response, Poisson contraction, inertia, gravity, sag,
//! adhesive detachment prediction or geometric feedback to conduction is claimed.
//! Force/moment controls prescribe quasi-static endpoints; release takes an
//! explicit positive duration, not an inferred impact or fracture event.
//! Resolve a prescribed unloading ramp when its work matters: an instantaneous
//! force jump can leave finite algorithmic dissipation even as dt approaches zero.

use std::collections::BTreeMap;
use fs_material::thermovisco::{Energy, Point, ThermoMaxwell};

/// Rectangle dimensions in the strip's local section frame (x along the strip).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LayerGeometry {
    /// Cross-sectional area, m2; volume is area times reference strip length.
    pub area_m2: f64,
    /// Centroid z coordinate from the section reference axis, m.
    pub centroid_z_m: f64,
    /// Thickness in z, m. Width is area/thickness.
    pub thickness_m: f64,
}
impl LayerGeometry {
    fn points(self) -> [f64; 2] {
        let offset = 0.5*self.thickness_m*0.577_350_269_189_625_8;
        [self.centroid_z_m-offset, self.centroid_z_m+offset]
    }
}

/// Axial displacement/strain restraint or applied tensile force.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AxialControl { Strain(f64), ForceN(f64) }
/// Bending restraint or applied work-conjugate section moment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BendingControl { CurvaturePerM(f64), MomentNm(f64) }
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Control { pub axial: AxialControl, pub bending: BendingControl }
impl Control {
    /// Fixed reference length and zero curvature. Reactions are outputs.
    pub const HELD: Self = Self { axial: AxialControl::Strain(0.0), bending: BendingControl::CurvaturePerM(0.0) };
    /// No axial force or bending moment. This does not predict when adhesion fails.
    pub const RELEASED: Self = Self { axial: AxialControl::ForceN(0.0), bending: BendingControl::MomentNm(0.0) };
}

/// Explicit retained-state ceilings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Budget { pub max_layers: usize, pub max_branch_states: usize }
/// Numerical gates and a bounded small-deformation geometry regime.
#[derive(Debug, Clone, Copy)]
pub struct Config {
    pub force_tolerance_n: f64,
    pub moment_tolerance_nm: f64,
    pub energy_tolerance_j: f64,
    /// Maximum extreme-fiber absolute strain, no greater than 0.05.
    pub max_strain: f64,
    /// Maximum absolute end rotation, no greater than 0.1 rad.
    pub max_rotation_rad: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Error {
    Invalid(&'static str),
    Budget,
    Cancelled,
    Material(fs_material::thermovisco::Error),
    Regime(&'static str),
    Equilibrium { force_n: f64, moment_nm: f64 },
    Energy(f64),
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "deposited strip refusal: {self:?}") }
}
impl std::error::Error for Error {}
impl From<fs_material::thermovisco::Error> for Error {
    fn from(value: fs_material::thermovisco::Error) -> Self { Self::Material(value) }
}

/// An accepted layer and its two independent material-point histories.
#[derive(Debug, Clone, PartialEq)]
pub struct Layer<'a> {
    pub geometry: LayerGeometry,
    pub birth_time_s: f64,
    points: [Point<'a>; 2],
}
impl<'a> Layer<'a> {
    #[must_use]
    pub fn material_points(&self) -> &[Point<'a>; 2] { &self.points }
}

/// Accepted thermal-history-dependent section; clone is an in-memory checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct LayeredStrip<'a> {
    length_m: f64,
    time_s: f64,
    strain: f64,
    curvature_per_m: f64,
    layers: BTreeMap<u64, Layer<'a>>,
    budget: Budget,
}

/// Actual endpoint resultants, reactions and extensive work in joules.
#[derive(Debug, Clone, Copy, Default)]
pub struct Report {
    pub axial_force_n: f64,
    pub bending_moment_nm: f64,
    pub axial_reaction_n: f64,
    pub bending_reaction_nm: f64,
    pub force_residual_n: f64,
    pub moment_residual_nm: f64,
    pub material_point_updates: usize,
    /// Every Energy member is multiplied by reference volume here: units J.
    pub energy_j: Energy,
}

/// Affine endpoint curvature at zero axial force and fixed temperature/duration.
/// This is the derivative of the same implicit section solve, not a fitted EI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BendingResponse {
    pub zero_moment_curvature_per_m: f64,
    /// d(curvature)/d(moment), in 1/(N m2), strictly positive.
    pub compliance_per_n_m2: f64,
}

impl<'a> LayeredStrip<'a> {
    pub fn new(length_m: f64, budget: Budget) -> Result<Self, Error> {
        positive(length_m)?;
        Ok(Self { length_m, time_s: 0.0, strain: 0.0, curvature_per_m: 0.0,
            layers: BTreeMap::new(), budget })
    }
    #[must_use]
    pub fn length_m(&self) -> f64 { self.length_m }
    #[must_use]
    pub fn time_s(&self) -> f64 { self.time_s }
    #[must_use]
    pub fn axial_strain(&self) -> f64 { self.strain }
    #[must_use]
    pub fn curvature_per_m(&self) -> f64 { self.curvature_per_m }
    #[must_use]
    pub fn layers(&self) -> &BTreeMap<u64, Layer<'a>> { &self.layers }
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.layers.values().map(|layer| 0.5*layer.geometry.area_m2*self.length_m
            * layer.points.iter().map(Point::stored_energy_j_m3).sum::<f64>()).sum()
    }

    /// Condense axial freedom without aging history or testing a fictitious
    /// unit-load material state. Final strain/rate admission belongs to advance_to.
    pub fn bending_response(&self, target_s: f64, temperatures_k: &BTreeMap<u64, f64>,
        mut checkpoint: impl FnMut() -> Result<(), Error>) -> Result<BendingResponse, Error> {
        checkpoint()?;
        if self.layers.is_empty() { return Err(Error::Invalid("unborn section has no bending compliance")); }
        let (scale, k, history) = self.linearization(target_s-self.time_s, temperatures_k, &mut checkpoint)?;
        let free = solve(k, [-history[0], -history[1]], [None, None])?;
        let unit = solve(k, [0.0, 1.0], [None, None])?;
        let compliance = finite(unit[1]/scale/scale)?;
        positive(compliance)?;
        checkpoint()?;
        Ok(BendingResponse { zero_moment_curvature_per_m: finite(free[1]/scale)?,
            compliance_per_n_m2: compliance })
    }

    fn linearization(&self, dt: f64, temperatures_k: &BTreeMap<u64, f64>,
        checkpoint: &mut impl FnMut() -> Result<(), Error>)
        -> Result<(f64, [[f64; 2]; 2], [f64; 2]), Error> {
        positive(dt)?;
        if temperatures_k.len() != self.layers.len() || self.layers.keys().any(|id| !temperatures_k.contains_key(id)) {
            return Err(Error::Invalid("one endpoint temperature per active layer required"));
        }
        let scale = self.layers.values().map(|l| l.geometry.centroid_z_m.abs()+0.5*l.geometry.thickness_m)
            .fold(0.0_f64, f64::max);
        let scale = if scale > 0.0 { scale } else { 1.0 };
        let mut k = [[0.0; 2]; 2]; let mut history = [0.0; 2];
        for (&id, layer) in &self.layers {
            checkpoint()?;
            for (point, z) in layer.points.iter().zip(layer.geometry.points()) {
                let affine = point.linearize(temperatures_k[&id], dt)?;
                let b = [1.0, -z/scale]; let area = 0.5*layer.geometry.area_m2;
                for i in 0..2 {
                    history[i] = finite(history[i]+area*b[i]*affine.offset_pa)?;
                    for j in 0..2 { k[i][j] = finite(k[i][j]+area*b[i]*b[j]*affine.tangent_pa)?; }
                }
            }
        }
        Ok((scale, k, history))
    }

    /// Birth at the current accepted time and CURRENT section strain/curvature.
    /// Both Gauss points start stress free; no prior layer history is reset.
    /// Stable IDs normally correspond to the thermal deposition cell IDs.
    pub fn deposit(&self, id: u64, geometry: LayerGeometry, law: &'a ThermoMaxwell,
        temperature_k: f64) -> Result<Self, Error> {
        if self.layers.contains_key(&id) { return Err(Error::Invalid("layer ID already active")); }
        let existing: usize = self.layers.values().map(|l| 2*l.points[0].law().branch_count()).sum();
        if self.layers.len() >= self.budget.max_layers
            || 2*law.branch_count() > self.budget.max_branch_states.saturating_sub(existing) {
            return Err(Error::Budget);
        }
        positive(geometry.area_m2)?; positive(geometry.thickness_m)?;
        finite(geometry.centroid_z_m)?;
        positive(finite(geometry.area_m2*self.length_m)?)?;
        positive(finite(geometry.area_m2/geometry.thickness_m)?)?;
        let z = geometry.points();
        if !z[0].is_finite() || !z[1].is_finite() || z[0] == z[1] {
            return Err(Error::Invalid("layer Gauss-point positions are not representable"));
        }
        let layer = Layer { geometry, birth_time_s: self.time_s, points: [
            law.birth(temperature_k, self.strain-z[0]*self.curvature_per_m)?,
            law.birth(temperature_k, self.strain-z[1]*self.curvature_per_m)?,
        ] };
        let mut next = self.clone(); next.layers.insert(id, layer); Ok(next)
    }

    /// Deformed position of a reference point (x,z) in the local strip plane, m.
    /// Rigid translation and rotation at x=0 are fixed by the coordinate convention.
    /// Callers retain the world-frame transform; no rendering-only deformation is
    /// substituted for the strain/curvature used by the actual material update.
    pub fn deformed_position_m(&self, x_m: f64, z_m: f64) -> Result<[f64; 2], Error> {
        if !x_m.is_finite() || !(0.0..=self.length_m).contains(&x_m) || !z_m.is_finite() {
            return Err(Error::Invalid("position must be finite and x within reference length"));
        }
        Ok([finite(x_m*(1.0+self.strain-z_m*self.curvature_per_m))?,
            finite(z_m+0.5*self.curvature_per_m*x_m*x_m)?])
    }

    /// One physical backward-Euler mechanics interval, using endpoint cell
    /// temperatures. Exactly one temperature per active layer is required.
    /// No accepted state is mutated by trials, refusals or cancellation.
    /// `checkpoint` can be bound directly to an executor Cx cancellation gate.
    pub fn advance_to(&self, target_s: f64, temperatures_k: &BTreeMap<u64, f64>,
        control: Control, config: Config, mut checkpoint: impl FnMut() -> Result<(), Error>)
        -> Result<(Self, Report), Error> {
        checkpoint()?;
        let dt = target_s-self.time_s; positive(dt)?;
        positive(config.force_tolerance_n)?; positive(config.moment_tolerance_nm)?;
        positive(config.energy_tolerance_j)?; positive(config.max_strain)?; positive(config.max_rotation_rad)?;
        if config.max_strain > 0.05 || config.max_rotation_rad > 0.1 {
            return Err(Error::Regime("small-strain/rotation ceiling cannot exceed 0.05 / 0.1 rad"));
        }
        let (scale, k, history) = self.linearization(dt, temperatures_k, &mut checkpoint)?;
        let (prescribed0, force) = match control.axial {
            AxialControl::Strain(e) => (Some(finite(e)?), 0.0),
            AxialControl::ForceN(n) => (None, finite(n)?),
        };
        let (prescribed1, moment) = match control.bending {
            BendingControl::CurvaturePerM(c) => (Some(finite(c*scale)?), 0.0),
            BendingControl::MomentNm(m) => (None, finite(m)?),
        };
        let rhs = [finite(force-history[0])?, finite(moment/scale-history[1])?];
        let x = if self.layers.is_empty() {
            if force != 0.0 || moment != 0.0 { return Err(Error::Invalid("unborn material cannot carry loads")); }
            [prescribed0.unwrap_or(self.strain), prescribed1.unwrap_or(self.curvature_per_m*scale)]
        } else {
            solve(k, rhs, [prescribed0, prescribed1])?
        };
        let strain = x[0]; let curvature = finite(x[1]/scale)?;
        if !strain.is_finite() || strain.abs()+curvature.abs()*scale > config.max_strain
            || curvature.abs()*self.length_m > config.max_rotation_rad {
            return Err(Error::Regime("predicted geometry leaves bounded strain/rotation regime"));
        }
        let mut next = self.clone(); let mut report = Report::default();
        for (&id, layer) in &self.layers {
            checkpoint()?;
            let area = 0.5*layer.geometry.area_m2; let volume = area*self.length_m;
            let tolerance = config.energy_tolerance_j/(2.0*self.layers.len() as f64)/volume;
            let mut updated = layer.clone();
            for (q, z) in layer.geometry.points().into_iter().enumerate() {
                let step = layer.points[q].advance(temperatures_k[&id], strain-z*curvature, dt, tolerance)?;
                report.axial_force_n += area*step.stress_pa;
                report.bending_moment_nm -= area*z*step.stress_pa;
                add_energy(&mut report.energy_j, step.energy, volume);
                report.material_point_updates += 1; updated.points[q] = step.point;
            }
            next.layers.insert(id, updated);
        }
        next.time_s = target_s; next.strain = strain; next.curvature_per_m = curvature;
        if prescribed0.is_some() { report.axial_reaction_n = report.axial_force_n; }
        else { report.force_residual_n = finite(report.axial_force_n-force)?; }
        if prescribed1.is_some() { report.bending_reaction_nm = report.bending_moment_nm; }
        else { report.moment_residual_nm = finite(report.bending_moment_nm-moment)?; }
        finite(report.axial_force_n)?; finite(report.bending_moment_nm)?;
        if report.force_residual_n.abs() > config.force_tolerance_n
            || report.moment_residual_nm.abs() > config.moment_tolerance_nm {
            return Err(Error::Equilibrium { force_n: report.force_residual_n, moment_nm: report.moment_residual_nm });
        }
        // Independent global storage and work-conjugacy checks, not just a sum
        // of constitutive defects that could cancel or omit a material point.
        let work = self.length_m*(report.axial_force_n*(strain-self.strain)
            +report.bending_moment_nm*(curvature-self.curvature_per_m));
        let storage = finite(next.stored_energy_j()-self.stored_energy_j())?;
        let e = &mut report.energy_j;
        let defect = finite(storage+e.dissipated+e.numerical_dissipation
            -e.mechanical_work-e.thermal_strain_work-e.modulus_work)?;
        if defect.abs() > config.energy_tolerance_j
            || finite(work-e.mechanical_work)?.abs() > config.energy_tolerance_j
            || finite(storage-e.stored_change)?.abs() > config.energy_tolerance_j {
            return Err(Error::Energy(defect));
        }
        e.residual = defect;
        checkpoint()?;
        Ok((next, report))
    }
}

fn solve(k: [[f64; 2]; 2], rhs: [f64; 2], prescribed: [Option<f64>; 2]) -> Result<[f64; 2], Error> {
    let result = match prescribed {
        [Some(a), Some(b)] => [a, b],
        [None, Some(b)] => { positive(k[0][0])?; [(rhs[0]-k[0][1]*b)/k[0][0], b] },
        [Some(a), None] => { positive(k[1][1])?; [a, (rhs[1]-k[1][0]*a)/k[1][1]] },
        [None, None] => {
            positive(k[0][0])?;
            let ratio = k[0][1]/k[0][0]; let schur = k[1][1]-ratio*k[0][1];
            if !schur.is_finite() || schur <= 64.0*f64::EPSILON*k[1][1].abs() {
                return Err(Error::Invalid("section bending stiffness is singular/ill-conditioned; choose a nearby reference axis"));
            }
            let b = (rhs[1]-ratio*rhs[0])/schur;
            [(rhs[0]-k[0][1]*b)/k[0][0], b]
        }
    };
    finite(result[0])?; finite(result[1])?; Ok(result)
}
fn add_energy(out: &mut Energy, e: Energy, volume: f64) {
    out.mechanical_work += volume*e.mechanical_work;
    out.thermal_strain_work += volume*e.thermal_strain_work;
    out.modulus_work += volume*e.modulus_work;
    out.stored_change += volume*e.stored_change;
    out.dissipated += volume*e.dissipated;
    out.numerical_dissipation += volume*e.numerical_dissipation;
    out.residual += volume*e.residual;
}
fn finite(v: f64) -> Result<f64, Error> {
    if v.is_finite() { Ok(v) } else { Err(Error::Invalid("nonfinite strip arithmetic")) }
}
fn positive(v: f64) -> Result<(), Error> {
    if v.is_finite() && v > 0.0 { Ok(()) } else { Err(Error::Invalid("positive finite SI quantity required")) }
}

#[cfg(test)]
#[path = "deposition/tests.rs"]
mod tests;
