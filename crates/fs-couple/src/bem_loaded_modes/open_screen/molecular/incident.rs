//! Sound-driven open-shell motion using the original two-sided work pairing.
use super::*;
use fs_bem::helmholtz::bc::IncidentField;
use fs_bem::reciprocal_radiation::open_screen::{
    MolecularScreenIllumination, MolecularScreenIlluminationLimits,
    PreparedMolecularScreenIlluminationReceivers,
};

/// Complete boundary illumination and mechanical-batch limits.
#[derive(Debug, Clone, Copy)]
pub struct MolecularOpenScreenIncidentOptions {
    pub boundary: MolecularScreenIlluminationLimits,
    pub maximum_load_cases: usize,
    /// cases*(112*modes+1024)+48*modes^2+320*modes+4096 logical bytes.
    /// Includes force columns, motion and total-reaction results, conversion
    /// headers and mechanical scratch. Source/boundary/receiver rows and
    /// allocator metadata have their separate limits.
    pub maximum_response_bytes: usize,
    /// 32*modes^3+(320*cases+256)*modes^2, excluding boundary/source sampling.
    pub maximum_work: u64,
    pub physical_residual_tolerance: f64,
    pub power_relative_tolerance: f64,
}

impl MolecularOpenScreenIncidentOptions {
    fn mechanical(self) -> OpenScreenHarmonicOptions {
        OpenScreenHarmonicOptions {
            boundary: self.boundary.boundary.ports,
            maximum_load_cases: self.maximum_load_cases,
            maximum_response_bytes: self.maximum_response_bytes,
            maximum_work: self.maximum_work,
            physical_residual_tolerance: self.physical_residual_tolerance,
            power_relative_tolerance: self.power_relative_tolerance,
        }
    }
}

/// Cycle-mean powers [W] for peak phasors. Sound can transfer energy INTO the
/// shell, so total_exterior_w may be negative and is never clipped.
#[derive(Debug, Clone, Copy)]
pub struct MolecularOpenScreenIncidentPower {
    pub external_mechanical_input_w: f64,
    /// Blocked-shell force applied to the actual loaded shell velocity.
    pub blocked_acoustic_input_w: f64,
    pub structural_viscosity_w: f64,
    /// Work of the outgoing motion-only field; includes molecular gas heat.
    pub motion_exterior_w: f64,
    /// Direct blocked pressure-jump work against actual motion. This exchange
    /// term retains incident/motion interference; no difference of large
    /// independently computed radiation powers is used to recover it.
    pub blocked_pressure_work_w: f64,
    /// Direct work of the TOTAL jump against actual structural velocity.
    /// It is not outgoing flux, scattered power or molecular absorption alone.
    pub total_exterior_w: f64,
    /// External mechanical input - structural heat - total exterior work.
    pub balance_w: f64,
    pub apparent_power_scale_w: f64,
    pub relative_balance: f64,
}

#[derive(Debug, Clone)]
pub struct MolecularOpenScreenIncidentResponse {
    /// The actual sound/force-driven motion. This member's reaction and power
    /// concern motion radiation and the effective external+blocked force.
    pub motion: MolecularOpenScreenHarmonicResponse,
    pub blocked_acoustic_force: Vec<C64>,
    /// Force ON the shell from the total jump: blocked force - Z*velocity.
    pub total_fluid_reaction: Vec<C64>,
    pub power: MolecularOpenScreenIncidentPower,
}

/// All source columns and blocked fields share one boundary LU; all physical
/// load cases share one mechanical LU. Original incident sources remain borrowed.
pub struct MolecularOpenScreenIncidentBatch<'p, 's, 'i> {
    radiation: MolecularScreenIllumination<'p, 's, 'i>,
    responses: Vec<MolecularOpenScreenIncidentResponse>,
}

pub struct PreparedMolecularOpenScreenIncidentReceivers<'b, 'p, 's, 'i> {
    receiver: PreparedMolecularScreenIlluminationReceivers<'b, 'p, 's, 'i>,
    responses: &'b [MolecularOpenScreenIncidentResponse],
}

impl<'s> PreparedOpenScreenModes<'s> {
    /// Illuminate the original shell and solve its resulting loaded motion.
    /// One complete external mechanical-force vector is required per incident
    /// case; an explicit zero vector selects purely sound-driven motion.
    /// Incoming pressure is continuous across a zero-thickness screen, while
    /// blocked scattering supplies the pressure jump that drives the shell.
    ///
    /// The original P1 boundary/source maps, molecular material, full K/C and
    /// all mutual entries are retained. This stationary harmonic owner does
    /// not infer an absorption cross-section or a causal time-domain response.
    pub fn harmonic_molecular_incident_response_batch<'p, 'i>(
        &'p self,
        exterior: &MolecularExterior,
        omega: f64,
        incidents: &[IncidentField<'i>],
        mechanical_forces: &[&[C64]],
        options: MolecularOpenScreenIncidentOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularOpenScreenIncidentBatch<'p, 's, 'i>, BemModesError> {
        poll(&mut cancelled)?;
        let mechanical = options.mechanical();
        self.source
            .admit_harmonic_batch(omega, mechanical_forces, mechanical)?;
        let m = self.source.mode_count() as u128;
        let c = incidents.len() as u128;
        let bytes = c * (112 * m + 1024) + 48 * m * m + 320 * m + 4096;
        let work = 32 * m * m * m + (320 * c + 256) * m * m;
        let gas = exterior.source().gas();
        if incidents.len() != mechanical_forces.len()
            || bytes > options.maximum_response_bytes as u128
            || work > options.maximum_work as u128
            || gas.density.to_bits() != self.source.medium.density.to_bits()
            || gas.sound_speed.to_bits() != self.source.medium.sound_speed.to_bits()
        {
            return Err(bad(
                "molecular screen illumination needs matched gas, paired complete forces and admitted response resources",
            ));
        }
        let rows: Vec<_> = self.source.weights.iter().map(Vec::as_slice).collect();
        let radiation = self
            .screen
            .solve_molecular_illumination(
                exterior,
                omega,
                &rows,
                incidents,
                options.boundary,
                &mut cancelled,
            )
            .map_err(BemModesError::MolecularScreen)?;
        let mut effective = Vec::with_capacity(incidents.len());
        for (external, blocked) in mechanical_forces.iter().zip(radiation.blocked_forces()) {
            poll(&mut cancelled)?;
            let force: Vec<_> = external.iter().zip(blocked).map(|(&a, &b)| a + b).collect();
            if force.iter().any(|&z| !finite(z)) {
                return Err(bad("molecular screen combined mechanical force overflows"));
            }
            effective.push(force);
        }
        let refs: Vec<_> = effective.iter().map(Vec::as_slice).collect();
        let solved = self.source.solve_harmonic_batch(
            omega,
            &refs,
            mechanical,
            radiation.motion().impedance(),
            radiation.maximum_backward_residual(),
            |velocity| radiation.motion().resisting_efforts(velocity),
            &mut cancelled,
        )?;
        let mut responses = Vec::with_capacity(incidents.len());
        for (case, response) in solved.into_iter().enumerate() {
            poll(&mut cancelled)?;
            let motion = MolecularOpenScreenHarmonicResponse::from_mechanical(response);
            let blocked = &radiation.blocked_forces()[case];
            let total = radiation
                .total_fluid_reaction(case, &motion.modal_velocity, &mut cancelled)
                .map_err(BemModesError::ScreenBoundary)?;
            let power = self.source.illuminated_work(
                mechanical_forces[case],
                blocked,
                &motion.modal_velocity,
                &motion.fluid_reaction,
                &total,
                options.power_relative_tolerance,
                &mut cancelled,
            )?;
            responses.push(MolecularOpenScreenIncidentResponse {
                motion,
                blocked_acoustic_force: blocked.clone(),
                total_fluid_reaction: total,
                power,
            });
        }
        poll(&mut cancelled)?;
        Ok(MolecularOpenScreenIncidentBatch {
            radiation,
            responses,
        })
    }
}

#[derive(Default, Clone, Copy)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }
    fn finish(self) -> f64 {
        self.value + self.correction
    }
}

impl OpenScreenLoadedModes {
    #[allow(clippy::too_many_arguments)]
    fn illuminated_work(
        &self,
        external: &[C64],
        blocked: &[C64],
        velocity: &[C64],
        motion_reaction: &[C64],
        total_reaction: &[C64],
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularOpenScreenIncidentPower, BemModesError> {
        let amplitude = external
            .iter()
            .chain(blocked)
            .chain(velocity)
            .chain(motion_reaction)
            .chain(total_reaction)
            .map(|z| z.abs())
            .fold(0.0, f64::max);
        let amplitude = if amplitude == 0.0 { 1.0 } else { amplitude };
        let scale = |value: C64| -> Result<C64, BemModesError> {
            let result = divide(value, amplitude);
            if !finite(result)
                || (value.re != 0.0 && result.re == 0.0)
                || (value.im != 0.0 && result.im == 0.0)
            {
                Err(bad(
                    "illuminated screen work normalization is unrepresentable",
                ))
            } else {
                Ok(result)
            }
        };
        let mut sums = [Sum::default(); 6]; // external, blocked, heat, motion, total, envelope
        let m = self.mode_count();
        for i in 0..m {
            poll(cancelled)?;
            let v = scale(velocity[i])?;
            let (f, b, r, t) = (
                scale(external[i])?,
                scale(blocked[i])?,
                scale(motion_reaction[i])?,
                scale(total_reaction[i])?,
            );
            sums[0].add(0.5 * (f * v.conj()).re);
            sums[1].add(0.5 * (b * v.conj()).re);
            sums[3].add(-0.5 * (r * v.conj()).re);
            sums[4].add(-0.5 * (t * v.conj()).re);
            sums[5].add(0.5 * (f.abs() + b.abs() + r.abs() + t.abs()) * v.abs());
            for j in 0..m {
                let u = scale(velocity[j])?;
                sums[2].add(0.5 * self.viscosity[i * m + j] * (v.conj() * u).re);
                sums[5].add(0.5 * self.viscosity[i * m + j].abs() * v.abs() * u.abs());
            }
        }
        let s = sums.map(Sum::finish);
        let balance = s[0] - s[2] - s[4];
        let relative = if s[5] == 0.0 {
            balance.abs()
        } else {
            balance.abs() / s[5]
        };
        if s.iter().any(|v| !v.is_finite()) || !relative.is_finite() {
            return Err(bad("illuminated screen signed work overflows"));
        }
        if relative > tolerance {
            return Err(BemModesError::PowerBalance {
                relative,
                tolerance,
            });
        }
        let restore = |v: f64| -> Result<f64, BemModesError> {
            let value = (v * amplitude) * amplitude;
            if !value.is_finite() || (v != 0.0 && value == 0.0) {
                Err(bad("illuminated screen physical work is unrepresentable"))
            } else {
                Ok(value)
            }
        };
        Ok(MolecularOpenScreenIncidentPower {
            external_mechanical_input_w: restore(s[0])?,
            blocked_acoustic_input_w: restore(s[1])?,
            structural_viscosity_w: restore(s[2])?,
            motion_exterior_w: restore(s[3])?,
            blocked_pressure_work_w: restore(-s[1])?,
            total_exterior_w: restore(s[4])?,
            balance_w: restore(balance)?,
            apparent_power_scale_w: restore(s[5])?,
            relative_balance: relative,
        })
    }
}

impl<'p, 's, 'i> MolecularOpenScreenIncidentBatch<'p, 's, 'i> {
    pub fn responses(&self) -> &[MolecularOpenScreenIncidentResponse] {
        &self.responses
    }
    pub fn radiation(&self) -> &MolecularScreenIllumination<'p, 's, 'i> {
        &self.radiation
    }
    pub fn total_pressure_jump(
        &self,
        case: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, BemModesError> {
        poll(&mut cancelled)?;
        let response = self
            .responses
            .get(case)
            .ok_or_else(|| bad("illuminated screen load-case index out of range"))?;
        self.radiation
            .total_pressure_jump(case, &response.motion.modal_velocity, cancelled)
            .map_err(BemModesError::ScreenBoundary)
    }
    pub fn prepare_receivers<'b>(
        &'b self,
        points_m: &[[f64; 3]],
        separation_m: f64,
        options: near_field::Options,
        maximum_row_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularOpenScreenIncidentReceivers<'b, 'p, 's, 'i>, BemModesError> {
        poll(&mut cancelled)?;
        let receiver = self
            .radiation
            .prepare_receivers(
                points_m,
                separation_m,
                options,
                maximum_row_bytes,
                &mut cancelled,
            )
            .map_err(BemModesError::ScreenBoundary)?;
        Ok(PreparedMolecularOpenScreenIncidentReceivers {
            receiver,
            responses: &self.responses,
        })
    }
}

impl PreparedMolecularOpenScreenIncidentReceivers<'_, '_, '_, '_> {
    pub fn points_m(&self) -> &[[f64; 3]] {
        self.receiver.points_m()
    }
    pub fn kernel_evaluations(&self) -> usize {
        self.receiver.kernel_evaluations()
    }
    /// Original incident + blocked scattering + motion radiation, with one
    /// coherent outgoing-jump integration and analytic incident added once.
    pub fn evaluate(
        &self,
        case: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, BemModesError> {
        poll(&mut cancelled)?;
        let response = self
            .responses
            .get(case)
            .ok_or_else(|| bad("illuminated screen load-case index out of range"))?;
        self.receiver
            .evaluate(case, &response.motion.modal_velocity, cancelled)
            .map_err(BemModesError::ScreenBoundary)
    }
}

#[cfg(test)]
mod tests;
