//! Complete original material, interconnection, forcing and clock directions.

use super::{
    ForcedOrbitSensitivity, ForcedSensitivityDirection, ForcedSensitivityError, PeriodicPortInput,
    PhsForcedSensitivityError, PhsOrbitError, PreparedPhsForcedSensitivity, bad,
};
use crate::phs_orbit::sideband::PhsStorageParameterAction;
use fs_orbit::ForcedParameterPoint;
use fs_phs::{PhsStructureParameterDirection, PortHamiltonian};

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

pub(super) fn ode_clock_column(
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
    let mut effort = vec![f64::NAN; n];
    let mut effort_partial = vec![0.0; n];
    let mut ports = vec![f64::NAN; m];
    let mut ports_partial = vec![0.0; m];
    let mut source_refused = false;
    let mut physical_error = None;
    let result = run(&mut |time, state, out| {
        effort.fill(f64::NAN);
        system.storage().gradient(state, &mut effort);
        if effort.iter().any(|x| !x.is_finite()) || !input.sample_into(time, &mut ports) {
            source_refused = true;
            return false;
        }
        if let Some(storage) = direction.storage {
            effort_partial.fill(f64::NAN);
            if !storage(state, &mut effort_partial) || effort_partial.iter().any(|v| !v.is_finite())
            {
                return false;
            }
        }
        if let Some(input) = direction.input {
            if !input.sample_into(time, &mut ports_partial) {
                return false;
            }
        }
        match physical.vector_field_from_effort_into(
            &effort,
            &effort_partial,
            &ports,
            &ports_partial,
            out,
        ) {
            Ok(()) => true,
            Err(error) => {
                physical_error = Some(error);
                false
            }
        }
    });
    if source_refused {
        return Err(PhsOrbitError::SourceRefused.into());
    }
    if let Some(error) = physical_error {
        return Err(error.into());
    }
    result.map_err(Into::into)
}
