//! Closed 3D winding terminals on the original harmonic A-phi field.
//!
//! Profiles, Whitney source pairings and omitted end-flux factors come from
//! `SpatialWindingBasis`. Source-bound stranded resistance is integrated by the
//! same binder as transient spatial windings. Resolved eddy-current cells stay
//! disjoint. Only the at-most-16-terminal impedance is dense: every unit field
//! and actual drive uses the native sparse harmonic owner and its full physical
//! acceptance gates. Cached unit cochains seed repeated drives without replacing
//! their actual-amplitude current, flux, heat and mechanical-power checks.
//!
//! Peak phasors use exp(-i omega t), so V = R I - i omega (F^T A + L_end I).
//! Prescribed motion in separate eddy regions can make Z nonreciprocal or active;
//! no symmetry or passivity gate is imposed on that terminal matrix. Strands
//! remain stationary, source-bound and unresolved internally. Skin/proximity
//! effects in the separate resolved conductors are retained. No strand-scale
//! eddies, end resistance, material refresh, saturation, or circuit graph is
//! inferred. Omitted end storage must not duplicate flux already in the mesh.

use super::super::{
    SpatialInduction, geometry,
    nonlinear::winding::{SpatialNonlinearWindingError, SpatialWindingConductor},
    winding_kernel::{SourceBinding, bind_profiles},
};
use super::{SpatialHarmonicInductionControls, SpatialHarmonicInductionResponse};
use crate::{
    planar::eddy::harmonic::circuit::{HarmonicCircuitError, terminal_solve},
    spatial::{
        SpatialMagneticError, SpatialMagneticMaterial,
        vector_potential::winding::{SpatialWindingBasis, SpatialWindingError},
    },
    winding::end::EndWindingError,
};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_math::c64::C64;
use fs_solver::SolveReport;

/// Native field allowances and bounded terminal preparation/drive work.
#[derive(Clone, Copy, Debug)]
pub struct SpatialHarmonicWindingControls {
    pub field: SpatialHarmonicInductionControls,
    /// Number of original physical windings, from 1 through 16.
    pub max_windings: usize,
    /// Additional profiles, unit cochains/reports, copper and terminal payload
    /// in real scalar slots. Native field/Krylov budgets apply separately.
    pub maximum_response_scalars: usize,
    /// True scaled residual of the native small complex terminal solve.
    pub relative_terminal_tolerance: f64,
    /// Absolute defect of each ACTUAL reconstructed voltage under voltage drive.
    pub voltage_tolerance_v: f64,
    /// Absolute real/reactive terminal, spatial-copper and field pairing defects.
    pub power_tolerance_w: f64,
}

#[derive(Debug)]
pub enum SpatialHarmonicWindingError {
    Invalid(&'static str),
    Field(SpatialMagneticError),
    Basis(SpatialWindingError),
    Sources(SpatialNonlinearWindingError),
    End(EndWindingError),
    Terminal(HarmonicCircuitError),
    Balance(&'static str),
    Cancelled,
}
impl core::fmt::Display for SpatialHarmonicWindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial harmonic windings: {self:?}")
    }
}
impl std::error::Error for SpatialHarmonicWindingError {}
impl From<SpatialMagneticError> for SpatialHarmonicWindingError {
    fn from(error: SpatialMagneticError) -> Self {
        match error {
            SpatialMagneticError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<SpatialWindingError> for SpatialHarmonicWindingError {
    fn from(error: SpatialWindingError) -> Self {
        match error {
            SpatialWindingError::Magnetic(SpatialMagneticError::Cancelled) => Self::Cancelled,
            other => Self::Basis(other),
        }
    }
}
impl From<SpatialNonlinearWindingError> for SpatialHarmonicWindingError {
    fn from(error: SpatialNonlinearWindingError) -> Self {
        match error {
            SpatialNonlinearWindingError::Cancelled => Self::Cancelled,
            other => Self::Sources(other),
        }
    }
}
impl From<EndWindingError> for SpatialHarmonicWindingError {
    fn from(error: EndWindingError) -> Self {
        Self::End(error)
    }
}
impl From<HarmonicCircuitError> for SpatialHarmonicWindingError {
    fn from(error: HarmonicCircuitError) -> Self {
        match error {
            HarmonicCircuitError::Cancelled => Self::Cancelled,
            other => Self::Terminal(other),
        }
    }
}
type Result<T> = core::result::Result<T, SpatialHarmonicWindingError>;
fn invalid(message: &'static str) -> SpatialHarmonicWindingError {
    SpatialHarmonicWindingError::Invalid(message)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable harmonic winding arithmetic"))
    }
}
fn phasor(value: C64) -> Result<C64> {
    finite(value.re)?;
    finite(value.im)?;
    Ok(value)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialHarmonicWindingError::Cancelled)
    } else {
        Ok(())
    }
}

/// Shared original field, closed winding profiles and fixed conductor sources.
pub struct SpatialHarmonicWindings<'s, 'b, 'm, C: SpatialConductorMaterial = ResolvedConductorState>
{
    field: &'s SpatialInduction<'m, C>,
    basis: &'b SpatialWindingBasis<'m>,
    profiles: Vec<Vec<[f64; 3]>>,
    sources: SourceBinding,
    controls: SpatialHarmonicWindingControls,
}

/// Frequency-specific impedance and cached sparse-field unit solutions.
pub struct SpatialHarmonicWindingFrequency<
    'p,
    's,
    'b,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    prepared: &'p SpatialHarmonicWindings<'s, 'b, 'm, C>,
    omega_rad_s: f64,
    velocity_m_s: Vec<[f64; 3]>,
    impedance_ohm: Vec<C64>,
    solutions: Vec<Vec<f64>>,
    reports: Vec<SolveReport>,
}

/// Original stranded bundle's intrinsic complex fields and cycle-mean heat.
#[derive(Clone, Debug)]
pub struct SpatialHarmonicWindingCell {
    pub source: SpatialWindingConductor,
    pub volume_m3: f64,
    pub conductor_volume_m3: f64,
    /// Within the actual filled conductor, not averaged over the whole cell.
    pub current_density_a_m2: [C64; 3],
    pub electric_field_v_m: [C64; 3],
    pub mean_power_w: f64,
}

/// Actual phase-sensitive native field and its physical terminal/heat account.
pub struct SpatialHarmonicWindingResponse<
    's,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    pub field: SpatialHarmonicInductionResponse<'s, 'm, C>,
    pub currents_a: Vec<C64>,
    pub terminal_voltages_v: Vec<C64>,
    pub winding_cells: Vec<SpatialHarmonicWindingCell>,
    pub mean_winding_power_w: f64,
    pub mean_end_energy_j: f64,
    /// Resolved active-field plus explicitly omitted end-region storage.
    pub mean_magnetic_energy_j: f64,
    /// sum(conj(I) V)/2; inductive reactive power is negative in this convention.
    pub supplied_power_va: C64,
    pub power_residual_w: f64,
    pub reactive_power_residual_var: f64,
    /// Actual voltage minus commanded voltage; zero for prescribed-current drives.
    pub maximum_terminal_voltage_defect_v: f64,
    /// Present only for prescribed-voltage drives; the small native LU's residual.
    pub terminal_solver_relative_residual: Option<f64>,
}
impl<C: SpatialConductorMaterial> SpatialHarmonicWindingResponse<'_, '_, C> {
    /// Heat only: prescribed-motion power is kept in the original field account.
    pub fn mean_dissipated_power_w(&self) -> f64 {
        self.mean_winding_power_w + self.field.balance().mean_ohmic_power_w
    }
}

impl<'s, 'b, 'm, C: SpatialConductorMaterial> SpatialHarmonicWindings<'s, 'b, 'm, C> {
    /// Bind the exact same magnetic owner and AC boundary traces as the basis.
    /// Every active profile cell needs its original scalar conductor source and
    /// explicit strand fill. Multiple insulated bundles can share a cell only
    /// within total fill one; they cannot overlap resolved eddy conductivity.
    pub fn new(
        field: &'s SpatialInduction<'m, C>,
        basis: &'b SpatialWindingBasis<'m>,
        sources: Vec<SpatialWindingConductor>,
        controls: SpatialHarmonicWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let count = basis.winding_count();
        let mesh = field.problem().mesh();
        let unknowns = field.prepared.mass.nrows();
        let payload = mesh
            .tets
            .len()
            .checked_mul(3)
            .and_then(|v| unknowns.checked_mul(2).and_then(|u| v.checked_add(u)))
            .and_then(|v| {
                controls
                    .field
                    .maximum_cycles
                    .checked_mul(8)
                    .and_then(|h| h.checked_add(8))
                    .and_then(|h| v.checked_add(h))
            })
            .and_then(|v| v.checked_mul(count))
            .and_then(|v| {
                mesh.tets
                    .len()
                    .checked_mul(6)
                    .and_then(|u| v.checked_add(u))
            })
            .and_then(|v| unknowns.checked_mul(2).and_then(|u| v.checked_add(u)))
            .and_then(|v| {
                mesh.edges
                    .len()
                    .checked_mul(2)
                    .and_then(|u| v.checked_add(u))
            })
            .and_then(|v| {
                mesh.vertex_count
                    .checked_mul(3)
                    .and_then(|u| v.checked_add(u))
            })
            .and_then(|v| sources.len().checked_mul(32).and_then(|u| v.checked_add(u)))
            .and_then(|v| {
                count
                    .checked_mul(count)
                    .and_then(|u| u.checked_mul(8))
                    .and_then(|u| v.checked_add(u))
            });
        if count == 0
            || count > controls.max_windings
            || controls.max_windings > 16
            || payload.is_none_or(|v| v > controls.maximum_response_scalars)
            || !positive(controls.relative_terminal_tolerance)
            || controls.relative_terminal_tolerance >= 1.0
            || !positive(controls.voltage_tolerance_v)
            || !positive(controls.power_tolerance_w)
            || !std::ptr::eq(field.problem(), basis.problem())
            || basis
                .unit_field(0)
                .is_none_or(|f| f.prescribed_edges() != field.prepared.fixed_edges.as_slice())
        {
            return Err(invalid(
                "bounded winding payload, positive allowances and identical original field/traces required",
            ));
        }
        let mut profiles = Vec::with_capacity(count);
        for winding in 0..count {
            poll(cancelled)?;
            profiles.push(
                basis
                    .unit_field(winding)
                    .ok_or(invalid("missing original winding profile"))?
                    .impressed_current_density_a_m2()
                    .to_vec(),
            );
        }
        let sources = bind_profiles(
            geometry(field.problem()),
            &profiles,
            |cell| field.conductors()[cell].is_some(),
            sources,
            None,
            |cell, source| match &field.problem().materials()[cell] {
                SpatialMagneticMaterial::Linear(material) if source.material() == material => {
                    Ok(())
                }
                _ => Err(SpatialNonlinearWindingError::Invalid(
                    "stranded resistance requires its exact original magnetic material bundle",
                )),
            },
            cancelled,
        )?;
        Ok(Self {
            field,
            basis,
            profiles,
            sources,
            controls,
        })
    }

    pub fn basis(&self) -> &'b SpatialWindingBasis<'m> {
        self.basis
    }
    pub fn field(&self) -> &'s SpatialInduction<'m, C> {
        self.field
    }
    pub fn winding_conductors(&self) -> &[SpatialWindingConductor] {
        &self.sources.conductors
    }
    pub fn resistance_ohm(&self) -> &[f64] {
        &self.sources.resistance_ohm
    }

    /// Prepare Z at one frequency and prescribed steady velocity. Only eddy
    /// regions may move; every vertex incident on a stranded cell must be still.
    /// Unit solves retain full Ampere/charge/flux and real/reactive-power gates.
    pub fn at_frequency(
        &self,
        omega_rad_s: f64,
        velocity_m_s: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicWindingFrequency<'_, 's, 'b, 'm, C>> {
        poll(cancelled)?;
        if !positive(omega_rad_s)
            || velocity_m_s.len() != self.field.problem().mesh().vertex_count
            || velocity_m_s.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "positive frequency and complete finite velocity required",
            ));
        }
        for source in &self.sources.conductors {
            poll(cancelled)?;
            if self.field.problem().mesh().tets[source.cell]
                .iter()
                .any(|&v| velocity_m_s[v as usize] != [0.0; 3])
            {
                return Err(invalid(
                    "fixed stranded profiles require stationary winding vertices",
                ));
            }
        }
        let count = self.basis.winding_count();
        let mut impedance_ohm = vec![C64::ZERO; count * count];
        let mut solutions = Vec::with_capacity(count);
        let mut reports = Vec::with_capacity(count);
        for source in 0..count {
            poll(cancelled)?;
            let current: Vec<_> = self.profiles[source]
                .iter()
                .map(|j| j.map(C64::from_re))
                .collect();
            let field = self.field.solve_harmonic(
                &current,
                velocity_m_s,
                omega_rad_s,
                self.controls.field,
                cancelled,
            )?;
            let flux = pair(self.basis, field.potential_wb(), cancelled)?;
            let mut unit = vec![0.0; count];
            unit[source] = 1.0;
            let end_flux = self.basis.end_inductance_model().flux_linkages_wb(&unit)?;
            for receiver in 0..count {
                let resistance = if receiver == source {
                    self.sources.resistance_ohm[receiver]
                } else {
                    0.0
                };
                impedance_ohm[receiver * count + source] = phasor(C64::new(
                    resistance + omega_rad_s * flux[receiver].im,
                    -omega_rad_s * finite(flux[receiver].re + end_flux[receiver])?,
                ))?;
            }
            // All full-field gates have passed. Retain only the seed and report;
            // unused per-unit B/J/cell recovery payload drops here.
            solutions.push(field.solution);
            reports.push(field.report);
        }
        poll(cancelled)?;
        Ok(SpatialHarmonicWindingFrequency {
            prepared: self,
            omega_rad_s,
            velocity_m_s: velocity_m_s.to_vec(),
            impedance_ohm,
            solutions,
            reports,
        })
    }
}

fn pair(
    basis: &SpatialWindingBasis<'_>,
    potential: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>> {
    let real: Vec<_> = potential.iter().map(|v| v.re).collect();
    let imaginary: Vec<_> = potential.iter().map(|v| v.im).collect();
    let real = basis.pair_potential_wb(&real, cancelled)?;
    let imaginary = basis.pair_potential_wb(&imaginary, cancelled)?;
    real.into_iter()
        .zip(imaginary)
        .map(|(r, i)| phasor(C64::new(r, i)))
        .collect()
}

fn power(voltages: &[C64], currents: &[C64], cancelled: &mut impl FnMut() -> bool) -> Result<C64> {
    let mut sum = [0.0_f64; 2];
    let mut correction = [0.0; 2];
    for (&voltage, &current) in voltages.iter().zip(currents) {
        poll(cancelled)?;
        let value = phasor((current.conj() * voltage).scale(0.5))?;
        for (axis, term) in [value.re, value.im].into_iter().enumerate() {
            let next = finite(sum[axis] + term)?;
            correction[axis] = finite(
                correction[axis]
                    + if sum[axis].abs() >= term.abs() {
                        (sum[axis] - next) + term
                    } else {
                        (term - next) + sum[axis]
                    },
            )?;
            sum[axis] = next;
        }
    }
    phasor(C64::new(sum[0] + correction[0], sum[1] + correction[1]))
}

impl<'s, 'm, C: SpatialConductorMaterial> SpatialHarmonicWindingFrequency<'_, 's, '_, 'm, C> {
    pub fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    /// Row-major terminal Z, preserving moving-field nonreciprocity. Inductive
    /// imaginary parts are negative under exp(-i omega t).
    pub fn impedance_ohm(&self) -> &[C64] {
        &self.impedance_ohm
    }
    pub fn unit_reports(&self) -> &[SolveReport] {
        &self.reports
    }

    /// Drive every original winding with a prescribed peak current phasor.
    /// Cached unit cochains seed ONE native actual-amplitude solve. A small
    /// seed residual can avoid Krylov work; full physical recovery and gates
    /// always run. End voltage/energy use original factors, never Z's expanded
    /// end block, so cancelling shared end flux cannot erase resolved physics.
    pub fn current_response(
        &self,
        currents_a: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicWindingResponse<'s, 'm, C>> {
        poll(cancelled)?;
        let prepared = self.prepared;
        let count = prepared.basis.winding_count();
        if currents_a.len() != count
            || currents_a
                .iter()
                .any(|i| !i.re.is_finite() || !i.im.is_finite())
        {
            return Err(invalid(
                "one finite peak current per original winding required",
            ));
        }
        let dimension = self.solutions[0].len() / 2;
        let mut seed = vec![0.0; 2 * dimension];
        let mut current = vec![[C64::ZERO; 3]; prepared.field.problem().mesh().tets.len()];
        for (w, &amplitude) in currents_a.iter().enumerate() {
            for (cell, profile) in prepared.profiles[w].iter().enumerate() {
                poll(cancelled)?;
                for axis in 0..3 {
                    current[cell][axis] =
                        phasor(current[cell][axis] + amplitude.scale(profile[axis]))?;
                }
            }
            for coordinate in 0..dimension {
                poll(cancelled)?;
                let value = C64::new(
                    self.solutions[w][coordinate],
                    self.solutions[w][dimension + coordinate],
                ) * amplitude;
                seed[coordinate] = finite(seed[coordinate] + value.re)?;
                seed[dimension + coordinate] = finite(seed[dimension + coordinate] + value.im)?;
            }
        }
        let field = prepared.field.solve_harmonic_seeded(
            &current,
            &self.velocity_m_s,
            self.omega_rad_s,
            prepared.controls.field,
            Some(&seed),
            cancelled,
        )?;
        let flux = pair(prepared.basis, field.potential_wb(), cancelled)?;
        let real: Vec<_> = currents_a.iter().map(|i| i.re).collect();
        let imaginary: Vec<_> = currents_a.iter().map(|i| i.im).collect();
        let end = prepared.basis.end_inductance_model();
        let end_real = end.flux_linkages_wb(&real)?;
        let end_imaginary = end.flux_linkages_wb(&imaginary)?;
        let mean_end_energy_j =
            finite(0.5 * end.stored_energy_j(&real)? + 0.5 * end.stored_energy_j(&imaginary)?)?;
        let mut field_voltages = Vec::with_capacity(count);
        let mut terminal_voltages_v = Vec::with_capacity(count);
        for w in 0..count {
            let active = phasor(C64::new(
                self.omega_rad_s * flux[w].im,
                -self.omega_rad_s * flux[w].re,
            ))?;
            field_voltages.push(active);
            terminal_voltages_v.push(phasor(
                active
                    + currents_a[w].scale(prepared.sources.resistance_ohm[w])
                    + C64::new(
                        self.omega_rad_s * end_imaginary[w],
                        -self.omega_rad_s * end_real[w],
                    ),
            )?);
        }
        let mut winding_cells = Vec::with_capacity(prepared.sources.conductors.len());
        let mut mean_winding_power_w = 0.0;
        let geometry = geometry(prepared.field.problem());
        for (source, &resistance) in prepared
            .sources
            .conductors
            .iter()
            .zip(&prepared.sources.cell_resistance_ohm)
        {
            poll(cancelled)?;
            let volume_m3 = geometry.geometry.geometry.vol_signed[source.cell].abs();
            let conductor_volume_m3 = finite(volume_m3 * source.fill_fraction)?;
            if !positive(conductor_volume_m3) {
                return Err(invalid("stranded conductor volume is not representable"));
            }
            let current = currents_a[source.winding];
            let mut current_density_a_m2 = [C64::ZERO; 3];
            let mut electric_field_v_m = [C64::ZERO; 3];
            for axis in 0..3 {
                current_density_a_m2[axis] = phasor(current.scale(
                    prepared.profiles[source.winding][source.cell][axis] / source.fill_fraction,
                ))?;
                electric_field_v_m[axis] =
                    phasor(current_density_a_m2[axis].scale(source.conductor.resistivity_ohm_m()))?;
            }
            let mean_power_w = finite(
                0.5 * (resistance * current.re) * current.re
                    + 0.5 * (resistance * current.im) * current.im,
            )?;
            mean_winding_power_w = finite(mean_winding_power_w + mean_power_w)?;
            winding_cells.push(SpatialHarmonicWindingCell {
                source: source.clone(),
                volume_m3,
                conductor_volume_m3,
                current_density_a_m2,
                electric_field_v_m,
                mean_power_w,
            });
        }
        let supplied_power_va = power(&terminal_voltages_v, currents_a, cancelled)?;
        // Check the original active source pairing independently, so large end
        // storage cannot mask a wrong winding orientation or field voltage.
        let active_power = power(&field_voltages, currents_a, cancelled)?;
        let native = field.balance();
        let active_defect = phasor(
            active_power
                - C64::new(
                    native.mean_supplied_power_w,
                    native.reactive_supplied_power_var,
                ),
        )?;
        let power_residual_w = finite(
            mean_winding_power_w + native.mean_ohmic_power_w + native.mean_mechanical_power_w
                - supplied_power_va.re,
        )?;
        let reactive_power_residual_var = finite(
            (supplied_power_va.im - native.reactive_supplied_power_var)
                + 2.0 * self.omega_rad_s * mean_end_energy_j
                + native.reactive_power_residual_var,
        )?;
        if active_defect.re.abs() > prepared.controls.power_tolerance_w
            || active_defect.im.abs() > prepared.controls.power_tolerance_w
            || power_residual_w.abs() > prepared.controls.power_tolerance_w
            || reactive_power_residual_var.abs() > prepared.controls.power_tolerance_w
        {
            return Err(SpatialHarmonicWindingError::Balance(
                "stranded copper, original field, prescribed motion and end power do not close",
            ));
        }
        let mean_magnetic_energy_j = finite(native.mean_magnetic_energy_j + mean_end_energy_j)?;
        poll(cancelled)?;
        Ok(SpatialHarmonicWindingResponse {
            field,
            currents_a: currents_a.to_vec(),
            terminal_voltages_v,
            winding_cells,
            mean_winding_power_w,
            mean_end_energy_j,
            mean_magnetic_energy_j,
            supplied_power_va,
            power_residual_w,
            reactive_power_residual_var,
            maximum_terminal_voltage_defect_v: 0.0,
            terminal_solver_relative_residual: None,
        })
    }

    /// Drive every winding with prescribed peak terminal voltage, using the
    /// existing scaled complex LU on the small terminal matrix. Singular or
    /// unresolved terminal systems refuse; no diagonal leakage is invented.
    /// After that solve, recover the actual native field and require EVERY
    /// reconstructed voltage to meet the original absolute voltage allowance.
    pub fn voltage_response(
        &self,
        voltages_v: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicWindingResponse<'s, 'm, C>> {
        poll(cancelled)?;
        if voltages_v.len() != self.prepared.basis.winding_count()
            || voltages_v
                .iter()
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(invalid(
                "one finite peak voltage per original winding required",
            ));
        }
        let (currents, residual) = terminal_solve(
            &self.impedance_ohm,
            voltages_v,
            self.prepared.controls.relative_terminal_tolerance,
            cancelled,
        )?;
        let mut response = self.current_response(&currents, cancelled)?;
        let mut maximum = 0.0_f64;
        for (&actual, &target) in response.terminal_voltages_v.iter().zip(voltages_v) {
            poll(cancelled)?;
            maximum = maximum.max(finite(phasor(actual - target)?.abs())?);
        }
        if maximum > self.prepared.controls.voltage_tolerance_v {
            return Err(SpatialHarmonicWindingError::Balance(
                "actual harmonic winding voltage misses its commanded terminal value",
            ));
        }
        response.maximum_terminal_voltage_defect_v = maximum;
        response.terminal_solver_relative_residual = Some(residual);
        Ok(response)
    }
}
