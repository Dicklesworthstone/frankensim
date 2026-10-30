//! Finite source gas and the ACTUAL inertial, thermally responsive tetrahedral wall.
//!
//! Reference pressure faces supply one native load footprint b. The same b
//! defines affine cavity volume V = V_ref + b dot (u - u_ref). Native Newmark
//! work is therefore (p_start+p_end)/2 times swept volume, less exterior-pressure
//! work, with additional fixed reference loads accounted separately.
//!
//! Two conservative-coordinate residuals close gas energy and geometric volume.
//! Each trial uses the original frozen gas caloric inversion, original spatial
//! heat substeps and original reciprocal thermoelastic/Newmark candidate. There
//! is no replacement piston, reduced modal basis, mechanical law or gas EOS.
//!
//! Gas composition/mass, reference faces and film geometry remain fixed. The
//! gas is well mixed. Its coupled ENDPOINT temperature is held during the
//! original explicit wall substeps; that thermal quadrature is first order.
//! Pressure uses native endpoint work, not the energy-midpoint gas law. Refine
//! the complete interval for time accuracy. Small-deformation/fully-solid source
//! restrictions are those of the original wall owner; no follower-load,
//! finite-strain volume, chemistry, nozzle flow or acoustic propagation is implied.

mod heat;

use super::*;
use crate::solid_enthalpy::gas::SolidGasContact;
use fs_material::thermal::source::ThermalMaterialError;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_solid::loads3::{SurfaceLoad3, VolumeLoad3};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress};
use std::cell::Cell;
use std::convert::Infallible;

pub use heat::CavityContactHeat;
use heat::{Film, FilmObserver};

#[derive(Debug)]
pub enum CavityError<E = Infallible> {
    Invalid(&'static str),
    Gas(GasEnergyError),
    Binding(SolidEnthalpyError),
    Solid(ThermoelasticFeedbackError),
    Thermal(DepositionError),
    Layer(ThermalMaterialError),
    Observer(E),
    Newton(NewtonError),
    Interrupted(SolverRunProgress),
    NotConverged(NewtonReport),
    TrialBudget,
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for CavityError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite gas/deformable material cavity: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CavityError<E> {}
impl<E> From<DepositionError> for CavityError<E> {
    fn from(value: DepositionError) -> Self { Self::Thermal(value) }
}
impl<E> From<GasEnergyError> for CavityError<E> {
    fn from(value: GasEnergyError) -> Self { Self::Gas(value) }
}
impl<E> From<SolidEnthalpyError> for CavityError<E> {
    fn from(value: SolidEnthalpyError) -> Self { Self::Binding(value) }
}
impl<E> From<ThermalMaterialError> for CavityError<E> {
    fn from(value: ThermalMaterialError) -> Self { Self::Layer(value) }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), CavityError<E>> {
    cx.checkpoint().map_err(|_| CavityError::Cancelled)
}
fn number<E>(value: f64) -> Result<f64, CavityError<E>> {
    if value.is_finite() { Ok(value) }
    else { Err(CavityError::Invalid("unrepresentable cavity energy, volume or force")) }
}
fn feedback<E>(error: ThermoelasticFeedbackError<CavityError<E>>) -> CavityError<E> {
    use ThermoelasticFeedbackError as F;
    let original = match error {
        F::Observer(error) => return error,
        F::Invalid(message) => F::Invalid(message),
        F::Source(error) => F::Source(error),
        F::Solid(error) => F::Solid(error),
        F::Thermal(error) => F::Thermal(error),
        F::NotConverged { iterations, heat_defect_j, combined_residual_j } =>
            F::NotConverged { iterations, heat_defect_j, combined_residual_j },
        F::Cancelled => F::Cancelled,
    };
    CavityError::Solid(original)
}
fn difference_error<E>(error: FiniteDifferenceError<CavityError<E>>) -> CavityError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => CavityError::Invalid(message),
    }
}
fn dot_motion<E>(
    forces: &[(usize, f64)], before: &[[f64; 3]], after: &[[f64; 3]],
) -> Result<f64, CavityError<E>> {
    if before.len() != after.len() { return Err(CavityError::Invalid("native nodal dimensions differ")); }
    let mut result = 0.0;
    for &(dof, coefficient) in forces {
        let a = before.get(dof / 3).ok_or(CavityError::Invalid("invalid pressure/load footprint"))?;
        let b = after.get(dof / 3).ok_or(CavityError::Invalid("invalid pressure/load footprint"))?;
        result = number(result + coefficient * (b[dof % 3] - a[dof % 3]))?;
    }
    Ok(result)
}

/// Numerical controls plus independent absolute physical admission gates.
#[derive(Clone, Copy, Debug)]
pub struct CavitySolve {
    pub gas: GasEnergySolve,
    pub solid: ThermoelasticSolve,
    pub thermal: StepConfig,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Finite-difference increment in the two J-scaled coordinates.
    pub difference_step_j: f64,
    pub newton: NewtonKrylovConfig,
    /// Positive supported physical volume interval, including the accepted start.
    pub volume_interval_m3: [f64; 2],
    pub volume_tolerance_m3: f64,
    /// Bounds gas, native pressure-work, interface and enclosing energy defects.
    pub energy_tolerance_j: f64,
}

/// Retained gas and structural history. Loads, volume reference and material
/// films are private so subsequent intervals cannot replace one accepted owner.
pub struct ThermoelasticCavity<'s, 'm> {
    solid: DynamicThermoelasticSolid<'s, 'm>,
    gas: ClosedGasParcel,
    faces: Vec<usize>,
    footprint_m2: Vec<(usize, f64)>,
    reference_loads_n: Vec<(usize, f64)>,
    exterior_pressure_pa: f64,
    reference_volume_m3: f64,
    reference_displacement_m: Vec<[f64; 3]>,
    films: Vec<Film>,
}

/// Physical account of one accepted interval, in original reference geometry.
pub struct CavityBalance {
    pub solid: DynamicThermoelasticReport,
    pub contacts: Vec<CavityContactHeat>,
    pub gas_energy_change_j: f64,
    /// Actual equal/opposite film heat; positive enters the material cells.
    pub gas_to_solid_heat_j: f64,
    pub external_gas_heat_j: f64,
    /// Original solid external heat AFTER removing internal gas-film exchange.
    pub external_solid_heat_j: f64,
    pub gas_pressure_work_j: f64,
    pub exterior_pressure_work_j: f64,
    pub reference_load_work_j: f64,
    pub swept_volume_m3: f64,
    pub geometric_volume_m3: f64,
    /// Stored gas volume minus the reference-anchored affine geometric volume.
    pub volume_residual_m3: f64,
    pub gas_energy_residual_j: f64,
    /// Native applied work minus gas/exterior/reference-load work.
    pub pressure_work_residual_j: f64,
    pub interface_heat_residual_j: f64,
    /// Delta(U_gas+H_cells+U_thermoelastic+K) + D + W_exterior
    /// - Q_gas - Q_other - W_reference - incoming cell energy.
    pub combined_residual_j: f64,
}

pub struct CavityCandidate<'s, 'm, 't, O> {
    pub thermal: DepositionSimulation<'t>,
    pub cavity: ThermoelasticCavity<'s, 'm>,
    /// The caller's original observer history, advanced only on this candidate.
    pub observer: O,
    pub balance: CavityBalance,
    pub nonlinear: NewtonReport,
    pub trials: usize,
}

struct Trial<'s, 'm, 't, O> {
    thermal: DepositionSimulation<'t>,
    solid: DynamicThermoelasticSolid<'s, 'm>,
    gas: ClosedGasParcel,
    observer: O,
    balance: CavityBalance,
}

impl<'s, 'm> ThermoelasticSolid<'s, 'm> {
    /// Launch a cavity from the actual loaded equilibrium and nodal velocity.
    ///
    /// Supply every moving cavity-facing reference surface. Native positive
    /// pressure acts inward on the solid; its identical footprint is positive
    /// cavity volume work. Omitted cavity surfaces must remain stationary. The
    /// caller owns cavity topology and the positive gas volume at this reference.
    ///
    /// The initial equilibrium load must be exactly this source gas pressure
    /// minus the fixed exterior pressure on these faces, plus the supplied
    /// fixed reference nodal loads. Later pressures are generated internally.
    /// Reference loads may represent prestress/external loading; their actual
    /// displacement work remains external. Supports remain stationary.
    ///
    /// Film patches retain their original material cells/claims and exposed-area
    /// budget. Empty patches declare an adiabatic cavity, not missing data.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_cavity(
        self, thermal: &DepositionSimulation<'_>, gas: ClosedGasParcel,
        faces: &[usize], exterior_pressure_pa: f64, reference_loads_n: &[(usize, f64)],
        contacts: &[SolidGasContact], velocity_m_s: &[[f64; 3]],
        transient: TransientConfig3, cx: &Cx<'_>,
    ) -> Result<ThermoelasticCavity<'s, 'm>, CavityError> {
        poll(cx)?;
        if faces.is_empty() || !exterior_pressure_pa.is_finite() || exterior_pressure_pa < 0.0 {
            return Err(CavityError::Invalid("nonempty native pressure faces and nonnegative exterior pressure required"));
        }
        self.admit_thermal(thermal, true, cx).map_err(CavityError::Solid)?;
        let source = self.source;
        let unit = source.geometry().assemble(
            &[SurfaceLoad3::Pressure { faces, pascals: 1.0 }], VolumeLoad3::None, &[], cx,
        ).map_err(SolidEnthalpyError::from)?;
        let footprint_m2 = unit.nodal_forces_n().to_vec();
        if footprint_m2.iter().all(|(_, value)| *value == 0.0) {
            return Err(CavityError::Invalid("nonzero native pressure footprint required"));
        }
        let loads = source.geometry().assemble(
            &[SurfaceLoad3::Pressure { faces, pascals: gas.state().state().pressure - exterior_pressure_pa }],
            VolumeLoad3::None, reference_loads_n, cx,
        ).map_err(SolidEnthalpyError::from)?;
        let films = heat::prepare(source, thermal, &gas, contacts, cx)?;
        let solid = self.begin_transient(&loads, velocity_m_s, transient, cx).map_err(CavityError::Solid)?;
        let reference_displacement_m = solid.nodal_motion(cx).map_err(CavityError::Solid)?[0].clone();
        poll(cx)?;
        Ok(ThermoelasticCavity {
            reference_volume_m3: gas.volume_m3(), gas, solid, faces: faces.to_vec(),
            footprint_m2, reference_loads_n: reference_loads_n.to_vec(), exterior_pressure_pa,
            reference_displacement_m, films,
        })
    }
}

impl<'s, 'm> ThermoelasticCavity<'s, 'm> {
    pub fn gas(&self) -> &ClosedGasParcel { &self.gas }
    pub fn solid(&self) -> &DynamicThermoelasticSolid<'s, 'm> { &self.solid }
    pub fn time_s(&self) -> f64 { self.solid.time_s() }

    /// Advance gas pressure, real material motion and reciprocal spatial heat.
    ///
    /// The two guesses are gas Delta U and p_start times Delta V, both joules.
    /// Each nonlinear/Jacobian probe starts from the SAME accepted gas, motion,
    /// thermal and observer histories. No trial changes a native clock or leaks
    /// its forces, temperatures or material source receipts into another trial.
    ///
    /// Pass the original material-heat/network observer when internal transport
    /// is source driven. All its source, conductance, reservoir, internal-transfer,
    /// interval, event and timestep hooks are retained. This adapter supplies
    /// only gas films; it never duplicates the observer's heat mechanisms.
    ///
    /// Volume is anchored to the initial accepted displacement, not accumulated
    /// from approximate roots. Gas energy uses original-law increments, and the
    /// final native pressure-work and total-energy gates are independent of the
    /// Newton residual norm. Energy closure is not a time-accuracy certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self, thermal: &DepositionSimulation<'t>, observer: &O, target_s: f64,
        external_gas_heat_j: f64, controls: CavitySolve, initial_coordinates_j: [f64; 2],
        cx: &Cx<'_>,
    ) -> Result<CavityCandidate<'s, 'm, 't, O>, CavityError<O::Error>> {
        poll(cx)?;
        let [vmin, vmax] = controls.volume_interval_m3;
        let p0 = self.gas.state().state().pressure;
        let tol = controls.energy_tolerance_j;
        if !target_s.is_finite() || target_s <= self.time_s() || !external_gas_heat_j.is_finite()
            || !vmin.is_finite() || vmin <= 0.0 || !vmax.is_finite() || vmax < vmin
            || self.gas.volume_m3() < vmin || self.gas.volume_m3() > vmax
            || !controls.volume_tolerance_m3.is_finite() || controls.volume_tolerance_m3 <= 0.0
            || !tol.is_finite() || tol <= 0.0 || controls.solid.energy_tolerance_j > tol
            || controls.maximum_trials == 0 || controls.maximum_newton_iterations == 0
            || !p0.is_finite() || p0 <= 0.0
            || self.gas.mass_kg() * controls.gas.energy_tolerance_j_kg > tol
        {
            return Err(CavityError::Invalid("positive physical interval, bounded volume/work and compatible native energy tolerances required"));
        }
        self.solid.admit::<CavityError<O::Error>>(thermal, true, cx).map_err(feedback)?;
        let before = self.solid.nodal_motion(cx).map_err(CavityError::Solid)?;
        let before_geometric = number(self.reference_volume_m3
            + dot_motion(&self.footprint_m2, &self.reference_displacement_m, &before[0])?)?;
        if (before_geometric - self.gas.volume_m3()).abs() > controls.volume_tolerance_m3 {
            return Err(CavityError::Invalid("accepted gas and native wall volume no longer agree"));
        }
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            poll(cx)?;
            if x.len() != 2 || !x.iter().all(|value| value.is_finite()) {
                return Err(CavityError::Invalid("two finite conservative cavity coordinates required"));
            }
            if trials.get() >= controls.maximum_trials { return Err(CavityError::TrialBudget); }
            trials.set(trials.get() + 1);
            let increment = number(x[1] / p0)?;
            let volume = number(self.gas.volume_m3() + increment)?;
            if x[1] != 0.0 && (increment == 0.0 || volume == self.gas.volume_m3()) {
                return Err(CavityError::Invalid("gas volume coordinate is not representable"));
            }
            if volume < vmin || volume > vmax {
                return Err(CavityError::Invalid("trial gas volume left its supported interval"));
            }
            let (gas, gas_energy_change_j) = self.gas.at_volume_internal_energy_increment(
                volume, x[0], controls.gas, || cx.checkpoint().is_err(),
            )?;
            let p1 = gas.state().state().pressure;
            let loads = self.solid.source.geometry().assemble(
                &[SurfaceLoad3::Pressure { faces: &self.faces, pascals: p1 - self.exterior_pressure_pa }],
                VolumeLoad3::None, &self.reference_loads_n, cx,
            ).map_err(SolidEnthalpyError::from)?;
            let film = FilmObserver::new(self, gas.state(), observer.clone());
            let candidate = self.solid.advance(
                thermal, &film, target_s, &loads, controls.thermal, controls.solid, cx,
            ).map_err(feedback)?;
            let (thermal_next, solid, film, report) = candidate.into_parts();
            let after = solid.nodal_motion(cx).map_err(CavityError::Solid)?;
            let swept = dot_motion(&self.footprint_m2, &before[0], &after[0])?;
            let geometric = number(self.reference_volume_m3
                + dot_motion(&self.footprint_m2, &self.reference_displacement_m, &after[0])?)?;
            if geometric < vmin || geometric > vmax {
                return Err(CavityError::Invalid("native wall left the positive supported cavity volume interval"));
            }
            let gas_work = number((0.5 * p0 + 0.5 * p1) * swept)?;
            let exterior_work = number(self.exterior_pressure_pa * swept)?;
            let reference_work = dot_motion(&self.reference_loads_n, &before[0], &after[0])?;
            let pressure_defect = number(report.mechanical_work_j - gas_work + exterior_work - reference_work)?;
            let external_solid = number(report.external_heat_j - film.heat_j)?;
            let gas_defect = number(gas_energy_change_j + gas_work - external_gas_heat_j + film.heat_j)?;
            let total = number(gas_energy_change_j + report.thermal_audit.stored_energy_change_j
                + report.internal_energy_change_j + report.kinetic_energy_change_j + report.damping_export_j
                + exterior_work - reference_work - report.thermal_audit.incoming_energy_j
                - external_solid - external_gas_heat_j)?;
            let interface_defect = number(film.samples.iter().map(|row| row.heat_into_solid_j).sum::<f64>() - film.heat_j)?;
            let balance = CavityBalance {
                solid: report, contacts: film.samples, gas_energy_change_j, gas_to_solid_heat_j: film.heat_j,
                external_gas_heat_j, external_solid_heat_j: external_solid, gas_pressure_work_j: gas_work,
                exterior_pressure_work_j: exterior_work, reference_load_work_j: reference_work,
                swept_volume_m3: swept, geometric_volume_m3: geometric,
                volume_residual_m3: number(volume - geometric)?, gas_energy_residual_j: gas_defect,
                pressure_work_residual_j: pressure_defect, interface_heat_residual_j: interface_defect,
                combined_residual_j: total,
            };
            poll(cx)?;
            Ok(Trial { thermal: thermal_next, solid, gas, observer: film.inner, balance })
        };
        let problem = FallibleFiniteDifference::new(2, controls.difference_step_j, |x: &[f64]| {
            let candidate = evaluate(x)?;
            Ok::<_, CavityError<O::Error>>(vec![
                candidate.balance.gas_energy_residual_j,
                number(p0 * candidate.balance.volume_residual_m3)?,
            ])
        }).map_err(difference_error)?;
        let construction = NewtonKrylovState::new(&problem, initial_coordinates_j.to_vec(), controls.newton);
        if let Some(error) = problem.take_failure() { return Err(difference_error(error)); }
        let mut solver = construction.map_err(CavityError::Newton)?;
        let result = solver.run_controlled(&problem, controls.maximum_newton_iterations,
            &mut || cx.checkpoint().is_err() || problem.has_failed());
        if let Some(error) = problem.take_failure() { return Err(difference_error(error)); }
        poll(cx)?;
        if result.progress != SolverRunProgress::Complete {
            return Err(CavityError::Interrupted(result.progress));
        }
        if !result.report.converged { return Err(CavityError::NotConverged(result.report)); }
        let trial = evaluate(&solver.x)?;
        let b = &trial.balance;
        if b.volume_residual_m3.abs() > controls.volume_tolerance_m3
            || [b.gas_energy_residual_j, b.pressure_work_residual_j, b.interface_heat_residual_j,
                b.combined_residual_j].iter().any(|value| value.abs() > tol)
        {
            return Err(CavityError::Invalid("gas, swept volume, pressure work or total energy did not close"));
        }
        poll(cx)?;
        Ok(CavityCandidate {
            thermal: trial.thermal,
            cavity: ThermoelasticCavity {
                solid: trial.solid, gas: trial.gas, faces: self.faces.clone(),
                footprint_m2: self.footprint_m2.clone(),
                reference_loads_n: self.reference_loads_n.clone(),
                exterior_pressure_pa: self.exterior_pressure_pa,
                reference_volume_m3: self.reference_volume_m3,
                reference_displacement_m: self.reference_displacement_m.clone(),
                films: self.films.clone(),
            },
            observer: trial.observer, balance: trial.balance,
            nonlinear: result.report, trials: trials.get(),
        })
    }
}
