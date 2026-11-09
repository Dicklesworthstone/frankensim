//! Complete original material, interconnection, forcing and clock directions.

use super::{
    ForcedOrbitSensitivity, ForcedSensitivityDirection, ForcedSensitivityError, PeriodicPortInput,
    PhsForcedSensitivityError, PhsOrbitError, PreparedPhsForcedSensitivity, bad,
};
use crate::phs_orbit::sideband::PhsStorageParameterAction;
use fs_orbit::ForcedParameterPoint;
use fs_phs::{
    PhsError, PhsStructureParameterDirection, PortHamiltonian, PreparedPhsParameterDirection,
};
use std::cell::RefCell;

/// One complete physical direction at the unchanged original state and phase.
/// None holds the corresponding original law fixed. Matrix derivatives are
/// signed directions, not newly fitted passive models. The original physical
/// input is retained by the preparation and cannot be replaced by this query.
pub struct PhsForcedParameterDirection<'a> {
    pub structure: PhsStructureParameterDirection<'a>,
    /// Complete partial_p grad H at fixed state, including material memory.
    pub storage: Option<&'a PhsStorageParameterAction<'a>>,
    /// Fourier-coefficient derivative at fixed phase; same base clock, ports
    /// and admitted bandwidth. A phase change is authored through these complex
    /// coefficients. The separate clock derivative does not change coefficients.
    pub input: Option<&'a PeriodicPortInput>,
}

impl PreparedPhsForcedSensitivity<'_> {
    /// Differentiate simultaneous J/R/G, storage, Fourier input and prescribed
    /// frequency changes. The original pHS equation has the port s*I, so its
    /// complete fixed-phase residual column is
    /// `i*n*Omega_p*X_n - Fourier[(J-R)h_p+(J_p-R_p)e+G_p*u+G*u_p]`.
    /// A Fourier source has no additional clock term at FIXED PHASE: adding its
    /// fixed-time derivative here would count the carrier change twice.
    ///
    /// All original coordinates and nonlinear/memory couplings remain. The
    /// returned waveform coefficients are fixed-phase derivatives; a fixed-time
    /// observation also differentiates its advancing phase. The source stays
    /// unchanged, and the existing coupled factor and two AFT grids are reused.
    pub fn solve_parameter(
        &mut self,
        chart: ForcedSensitivityDirection,
        direction: PhsForcedParameterDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ForcedOrbitSensitivity, PhsForcedSensitivityError> {
        if cancelled() {
            return Err(ForcedSensitivityError::Cancelled.into());
        }
        let n = self.system.state_dim();
        with_parameter_direction(
            self.system,
            &self.input,
            self.prepared.harmonics(),
            direction,
            |action| {
                self.prepared.solve_split_direction(
                    chart,
                    action,
                    |point, out| {
                        ode_clock_column(n, chart.angular_frequency_derivative, point, out)
                    },
                    &mut cancelled,
                )
            },
        )
    }
}

pub(in crate::phs_orbit) fn ode_clock_column(
    n: usize,
    omega_derivative: f64,
    point: ForcedParameterPoint<'_>,
    out: &mut [f64],
) -> bool {
    out.fill(0.0);
    for harmonic in 1..=point.harmonics {
        for state in 0..n {
            let re = n + (harmonic - 1) * 2 * n + state;
            let im = re + n;
            out[re] = -(harmonic as f64) * omega_derivative * point.coefficients[im];
            out[im] = harmonic as f64 * omega_derivative * point.coefficients[re];
        }
    }
    out.iter().all(|v| v.is_finite())
}

pub(super) fn with_parameter_direction<T>(
    system: &PortHamiltonian,
    input: &PeriodicPortInput,
    harmonics: usize,
    direction: PhsForcedParameterDirection<'_>,
    run: impl FnOnce(
        &mut dyn FnMut(f64, &[f64], &mut [f64]) -> bool,
    ) -> Result<T, ForcedSensitivityError>,
) -> Result<T, PhsForcedSensitivityError> {
    with_parameter_samples(system, input, harmonics, direction, |samples| {
        run(&mut |time, state, out| samples.borrow_mut().vector_field(time, state, out))
    })
}

// One original physical sample owner serves both the residual column and the
// direct objective partial. The numeric consumer invokes its two callbacks
// sequentially; sharing this scratch avoids duplicate source assembly/buffers.
pub(super) fn with_parameter_samples<'s, 'd, T>(
    system: &'s PortHamiltonian,
    input: &'s PeriodicPortInput,
    harmonics: usize,
    direction: PhsForcedParameterDirection<'d>,
    run: impl FnOnce(&RefCell<ParameterSamples<'s, 'd>>) -> Result<T, ForcedSensitivityError>,
) -> Result<T, PhsForcedSensitivityError> {
    if let Some(derivative) = direction.input {
        if derivative.omega_rad_s() != input.omega_rad_s()
            || derivative.port_count() != system.port_dim()
            || derivative.highest_nonzero_harmonic() > harmonics
        {
            return Err(bad(
                "physical input direction requires original clock, ports and retained bandwidth",
            )
            .into());
        }
    }
    let physical = system.prepare_parameter_direction(direction.structure)?;
    let (n, m) = (system.state_dim(), system.port_dim());
    let samples = RefCell::new(ParameterSamples {
        system,
        input,
        physical,
        direction,
        effort: vec![f64::NAN; n],
        effort_partial: vec![0.0; n],
        ports: vec![f64::NAN; m],
        ports_partial: vec![0.0; m],
        source_refused: false,
        physical_error: None,
    });
    let result = run(&samples);
    let samples = samples.into_inner();
    if samples.source_refused {
        return Err(PhsOrbitError::SourceRefused.into());
    }
    if let Some(error) = samples.physical_error {
        return Err(error.into());
    }
    result.map_err(Into::into)
}

pub(super) struct ParameterSamples<'s, 'd> {
    system: &'s PortHamiltonian,
    input: &'s PeriodicPortInput,
    physical: PreparedPhsParameterDirection<'s, 'd>,
    direction: PhsForcedParameterDirection<'d>,
    effort: Vec<f64>,
    effort_partial: Vec<f64>,
    ports: Vec<f64>,
    ports_partial: Vec<f64>,
    source_refused: bool,
    physical_error: Option<PhsError>,
}

impl ParameterSamples<'_, '_> {
    fn load(&mut self, time: f64, state: &[f64]) -> bool {
        self.effort.fill(f64::NAN);
        self.system.storage().gradient(state, &mut self.effort);
        if self.effort.iter().any(|x| !x.is_finite())
            || !self.input.sample_into(time, &mut self.ports)
        {
            self.source_refused = true;
            return false;
        }
        if let Some(storage) = self.direction.storage {
            self.effort_partial.fill(f64::NAN);
            if !storage(state, &mut self.effort_partial)
                || self.effort_partial.iter().any(|v| !v.is_finite())
            {
                return false;
            }
        }
        if let Some(input) = self.direction.input {
            if !input.sample_into(time, &mut self.ports_partial) {
                return false;
            }
        }
        true
    }

    pub(super) fn vector_field(&mut self, time: f64, state: &[f64], out: &mut [f64]) -> bool {
        if !self.load(time, state) {
            return false;
        }
        match self.physical.vector_field_from_effort_into(
            &self.effort,
            &self.effort_partial,
            &self.ports,
            &self.ports_partial,
            out,
        ) {
            Ok(()) => true,
            Err(error) => {
                self.physical_error = Some(error);
                false
            }
        }
    }

    pub(super) fn objective_partial(
        &mut self,
        time: f64,
        state: &[f64],
        objective: super::power::PhsPowerObjective,
        omega_derivative: f64,
        out: &mut f64,
    ) -> bool {
        if !self.load(time, state) {
            return false;
        }
        let power = if objective.supplied() {
            self.physical.supplied_power_from_effort(
                &self.effort,
                &self.effort_partial,
                &self.ports,
                &self.ports_partial,
            )
        } else {
            self.physical
                .resistance_power_from_effort(&self.effort, &self.effort_partial)
        };
        let direct = match power {
            Ok(value) => value,
            Err(error) => {
                self.physical_error = Some(error);
                return false;
            }
        };
        let original = if objective.per_cycle() && omega_derivative != 0.0 {
            let power = if objective.supplied() {
                self.system
                    .supplied_power_from_effort(&self.effort, &self.ports)
            } else {
                self.system.resistance_power_from_effort(&self.effort)
            };
            match power {
                Ok(value) => value,
                Err(error) => {
                    self.physical_error = Some(error);
                    return false;
                }
            }
        } else {
            0.0
        };
        let Some(value) =
            objective.direct_partial(self.input.omega_rad_s(), omega_derivative, direct, original)
        else {
            return false;
        };
        *out = value;
        true
    }
}
