//! Original Neumann actions and admitted traces for simultaneous physical solves.
use super::*;

impl IterativeRadiation<'_> {
    /// Apply the ORIGINAL velocity-to-right-side map B*(i*omega*rho*v).
    /// This shares all triangle, directional, self and near corrections with
    /// `start`; no acoustic solve or impedance column is constructed.
    /// Coupled callers must prepare a bare Neumann owner, without Robin data.
    pub fn coupled_velocity_image(
        &self,
        velocity: &[C64],
        gate: &CancelGate,
    ) -> Result<Vec<C64>, IterativeError> {
        self.admit_coupled_neumann(gate)?;
        self.boundary_rhs(velocity, gate)
    }

    /// Existing native A-block right preconditioner, for a larger simultaneous
    /// solve. It changes solver coordinates only, never the physical operator.
    pub fn coupled_pressure_precondition(
        &self,
        pressure: &[C64],
        gate: &CancelGate,
    ) -> Result<Vec<C64>, IterativeError> {
        self.admit_coupled_neumann(gate)?;
        if pressure.len() != self.surface.areas().len() || pressure.iter().any(|&v| !finite(v)) {
            return Err(IterativeError::Invalid {
                what: "one finite coupled pressure coordinate per panel required",
            });
        }
        let answer = self.inverse(pressure);
        poll(gate)?;
        if answer.iter().any(|&v| !finite(v)) {
            return Err(IterativeError::Invalid {
                what: "coupled pressure preconditioner produced a nonfinite value",
            });
        }
        Ok(answer)
    }

    /// Admit independently solved PHYSICAL Neumann traces against the same
    /// original FMM equation, then produce the native observer-bound field.
    /// The reported backward defect is ||A p-B q||/(||A p||+||B q||).
    /// It is neither a relative RHS residual nor a spatial/FMM error bound.
    /// The preparation's unchanged tolerance gates this additional check.
    pub fn finish_coupled_traces(
        &self,
        pressure: Vec<C64>,
        velocity: Vec<C64>,
        gate: &CancelGate,
    ) -> Result<(RadiationSolution, f64), IterativeError> {
        self.admit_coupled_neumann(gate)?;
        let a = self.apply(&pressure, gate)?;
        let b = self.boundary_rhs(&velocity, gate)?;
        let scale = a
            .iter()
            .chain(&b)
            .fold(0.0_f64, |s, z| s.max(z.re.abs()).max(z.im.abs()));
        let defect = if scale == 0.0 {
            0.0
        } else {
            let left = a
                .iter()
                .flat_map(|z| [z.re / scale, z.im / scale])
                .collect::<Vec<_>>();
            let right = b
                .iter()
                .flat_map(|z| [z.re / scale, z.im / scale])
                .collect::<Vec<_>>();
            let difference = left
                .iter()
                .zip(&right)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>();
            fs_solver::norm2(&difference) / (fs_solver::norm2(&left) + fs_solver::norm2(&right))
        };
        if !defect.is_finite() || defect > self.options.tol {
            return Err(IterativeError::Invalid {
                what: "coupled physical Neumann trace residual exceeds the original tolerance",
            });
        }
        let (radiated_power, radiated_power_roundoff_interval, power_diagnostics) =
            evaluate_radiated_power(&pressure, &velocity, self.surface.areas(), self.medium);
        poll(gate)?;
        Ok((
            RadiationSolution {
                pressure,
                velocity,
                k: self.k,
                medium: self.medium,
                surface_fingerprint: surface_fingerprint(self.surface),
                panels_per_wavelength: self.ppw,
                radiated_power,
                radiated_power_roundoff_interval,
                power_diagnostics,
                condition_lower_bound: 1.0,
                dense_cap_utilization: self.surface.areas().len() as f64 / MAX_DENSE_PANELS as f64,
            },
            defect,
        ))
    }

    /// Original FMM preparation resources with counters supplied by a composed
    /// solve. One composed matvec may apply both A and B; these counters do not
    /// claim to count individual FMM potentials, near kernels or factor solves.
    pub fn coupled_resource_receipt(
        &self,
        matvecs: usize,
        iterations: usize,
        solve_wall_s: f64,
    ) -> Result<ResourceReceipt, IterativeError> {
        if !solve_wall_s.is_finite() || solve_wall_s < 0.0 {
            return Err(IterativeError::Invalid {
                what: "finite nonnegative coupled solve duration required",
            });
        }
        Ok(self.receipt_for(solve_wall_s, matvecs, iterations))
    }

    fn admit_coupled_neumann(&self, gate: &CancelGate) -> Result<(), IterativeError> {
        poll(gate)?;
        if self.admittance.is_some() {
            return Err(IterativeError::Invalid {
                what: "simultaneous Neumann actions require an owner without Robin admittance",
            });
        }
        Ok(())
    }
}
