//! Bounded local Robin coordinates, including exact pressure-release patches.
//!
//! With w=rho*c*v, r=(Z-rho*c)/(Z+rho*c), a=(1+r)/2 and b=(1-r)/2,
//! the passive wall is b*p+a*w=0. Unlike beta=rho*c/Z, these coefficients
//! remain finite at a mechanical resonance (Z=0, r=-1). This is the same
//! original CBIE/HBIE pair A*p=i*k*B*w, not a different spatial solver.
//!
//! For prescribed wall data h, write p=a*x+conj(b)*h/s and
//! w=-b*x+conj(a)*h/s, s=|a|^2+|b|^2 >= 1/2. Two in-place column
//! transforms give (A*a+i*k*B*b)*x=(-A*conj(b)+i*k*B*conj(a))*h/s.
//! Thus neither endpoint requires an infinite coefficient, a clipped law,
//! extra dense matrix storage, nor an artificial wall resistance.

use super::*;

enum ReflectionFactor<'a> {
    Rigid(PreparedBoundary<'a>),
    Local {
        drive: Vec<C64>,
        lu: LuComplex,
        condition_lower_bound: f64,
    },
}

/// One dense local boundary factorization in bounded reflection coordinates.
///
/// Coefficients are normal-incidence pressure reflection in exp(-i omega t).
/// `+1` is rigid and `-1` is pressure release. Passive laws require |r|<=1
/// within arithmetic roundoff; no coefficient is clipped. A prescribed drive
/// retains the convention v=v_drive-Y*p wherever finite Y exists.
///
/// This owner reuses the existing closed-surface geometry, CBIE/Burton-Miller
/// matrices, dense work ceiling and receiver traces. It does not supply an
/// iterative, nonlocal, moving-wall or nonlinear boundary model. Preparation
/// has the original dense LU cancellation granularity.
pub struct PreparedReflectionBoundary<'a> {
    surface: &'a SpherePanels,
    k: f64,
    medium: Medium,
    pressure_coefficient: Vec<C64>,
    flow_coefficient: Vec<C64>,
    inverse_norm: Vec<f64>,
    declared_active: bool,
    panels_per_wavelength: f64,
    factor: ReflectionFactor<'a>,
}

impl<'a> PreparedReflectionBoundary<'a> {
    /// Admit and copy one coefficient per original panel, then factor once.
    pub fn new(
        surface: &'a SpherePanels,
        k: f64,
        medium: Medium,
        reflection: &[C64],
        declared_active: bool,
        formulation: Formulation,
    ) -> Result<Self, HelmholtzError> {
        let ppw = admit_resolution(surface, k)?;
        incident::validate_frequency(k, medium)?;
        validate_closed(surface)?;
        let n = surface.areas().len();
        if reflection.len() != n
            || reflection
                .iter()
                .any(|&r| !finite(r) || (!declared_active && r.abs() > 1.0 + 64.0 * f64::EPSILON))
        {
            return Err(bad(
                "one finite, declared-passive reflection per boundary panel required",
            ));
        }
        let pressure_coefficient: Vec<_> = reflection
            .iter()
            .map(|&r| (C64::ONE + r).scale(0.5))
            .collect();
        let flow_coefficient: Vec<_> = reflection
            .iter()
            .map(|&r| (C64::ONE - r).scale(0.5))
            .collect();
        let inverse_norm: Vec<_> = pressure_coefficient
            .iter()
            .zip(&flow_coefficient)
            .map(|(&a, &b)| 1.0 / (a.norm_sq() + b.norm_sq()))
            .collect();
        if inverse_norm.iter().any(|&s| !s.is_finite() || s <= 0.0) {
            return Err(bad("unrepresentable local reflection normalization"));
        }
        let factor = if reflection.iter().all(|&r| r == C64::ONE) {
            ReflectionFactor::Rigid(PreparedBoundary::new(
                surface,
                k,
                medium,
                Admittance::Rigid,
                formulation,
            )?)
        } else {
            let (mut matrix, mut drive) =
                assemble_dense_with_static(surface, k, alpha_for(formulation, k), None);
            let ik = C64::new(0.0, k);
            for i in 0..n {
                for j in 0..n {
                    let index = i * n + j;
                    let a = matrix[index];
                    let ikb = ik * drive[index];
                    let p = pressure_coefficient[j];
                    let f = flow_coefficient[j];
                    matrix[index] = a * p + ikb * f;
                    drive[index] = (ikb * p.conj() - a * f.conj()).scale(inverse_norm[j]);
                }
            }
            if matrix.iter().chain(&drive).any(|&z| !finite(z)) {
                return Err(bad("nonfinite reflection boundary operator"));
            }
            let lu = lu_complex(&matrix, n).map_err(|_| HelmholtzError::Singular)?;
            let condition_lower_bound = condition_lower_bound(&matrix, n, &lu);
            ReflectionFactor::Local {
                drive,
                lu,
                condition_lower_bound,
            }
        };
        Ok(Self {
            surface,
            k,
            medium,
            pressure_coefficient,
            flow_coefficient,
            inverse_norm,
            declared_active,
            panels_per_wavelength: ppw,
            factor,
        })
    }

    /// Prescribed outward drive velocity [m/s] through the prepared wall.
    /// A drive on a pressure-release patch has zero pressure work, as the
    /// same finite-impedance equations require in their Z->0 limit.
    pub fn solve(&self, drive_velocity: &[C64]) -> Result<BoundarySolution, HelmholtzError> {
        if drive_velocity.len() != self.surface.areas().len()
            || drive_velocity.iter().any(|&v| !finite(v))
        {
            return Err(bad("reflection boundary needs one finite drive per panel"));
        }
        if let ReflectionFactor::Rigid(owner) = &self.factor {
            let mut result = owner.solve(drive_velocity)?;
            result.declared_active = self.declared_active;
            return Ok(result);
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let h: Vec<_> = drive_velocity
            .iter()
            .zip(&self.pressure_coefficient)
            .map(|(&v, &a)| a * v.scale(rho_c))
            .collect();
        let (radiation, mut wall_amplitude) = self.trace(&h)?;
        for i in 0..wall_amplitude.len() {
            wall_amplitude[i] = wall_amplitude[i]
                + self.flow_coefficient[i].conj()
                    * drive_velocity[i].scale(rho_c * self.inverse_norm[i]);
        }
        let absorbed_power_w = self.wall_work(&wall_amplitude)?;
        let driving_power_w: f64 = radiation
            .pressure
            .iter()
            .zip(drive_velocity)
            .zip(self.surface.areas())
            .map(|((&p, &v), &area)| 0.5 * (p * v.conj()).re * area)
            .sum();
        if !driving_power_w.is_finite() {
            return Err(bad("nonfinite reflection boundary drive work"));
        }
        Ok(BoundarySolution {
            radiation,
            absorbed_power_w,
            driving_power_w,
            declared_active: self.declared_active,
        })
    }

    /// Scatter a coherently sampled incident field, including BOTH its
    /// pressure and velocity. The outgoing trace is separate from the total
    /// wall trace; absorption uses the latter, including interference.
    pub fn scatter(&self, incident: &IncidentTraces) -> Result<ScatteringSolution, HelmholtzError> {
        incident.validate_binding(self.surface, self.k, self.medium)?;
        if let ReflectionFactor::Rigid(owner) = &self.factor {
            let mut result = owner.scatter(incident)?;
            result.declared_active = self.declared_active;
            return Ok(result);
        }
        let rho_c = self.medium.density * self.medium.sound_speed;
        let h: Vec<_> = (0..self.surface.areas().len())
            .map(|i| {
                -(self.flow_coefficient[i] * incident.pressure_pa()[i]
                    + self.pressure_coefficient[i] * incident.normal_velocity_m_s()[i].scale(rho_c))
            })
            .collect();
        let (scattered, mut wall_amplitude) = self.trace(&h)?;
        for i in 0..wall_amplitude.len() {
            wall_amplitude[i] = wall_amplitude[i]
                + (self.pressure_coefficient[i].conj() * incident.pressure_pa()[i]
                    - self.flow_coefficient[i].conj()
                        * incident.normal_velocity_m_s()[i].scale(rho_c))
                .scale(self.inverse_norm[i]);
        }
        let absorbed_power_w = self.wall_work(&wall_amplitude)?;
        Ok(ScatteringSolution {
            scattered,
            incident_pressure_pa: incident.pressure_pa().to_vec(),
            incident_normal_velocity_m_s: incident.normal_velocity_m_s().to_vec(),
            absorbed_power_w,
            declared_active: self.declared_active,
        })
    }

    fn trace(&self, h: &[C64]) -> Result<(RadiationSolution, Vec<C64>), HelmholtzError> {
        let ReflectionFactor::Local {
            drive,
            lu,
            condition_lower_bound,
        } = &self.factor
        else {
            return Err(bad("rigid reflection must use the original Neumann owner"));
        };
        let n = h.len();
        let mut x = vec![C64::ZERO; n];
        for i in 0..n {
            let mut sum = ComplexSum::default();
            for j in 0..n {
                sum.add(drive[i * n + j] * h[j]);
            }
            x[i] = sum.finish();
        }
        lu.solve(&mut x);
        let rho_c = self.medium.density * self.medium.sound_speed;
        let mut pressure = Vec::with_capacity(n);
        let mut velocity = Vec::with_capacity(n);
        for i in 0..n {
            let p = self.pressure_coefficient[i];
            let f = self.flow_coefficient[i];
            let particular = h[i].scale(self.inverse_norm[i]);
            pressure.push(p * x[i] + f.conj() * particular);
            velocity.push((-f * x[i] + p.conj() * particular).scale(1.0 / rho_c));
        }
        if pressure.iter().chain(&velocity).any(|&z| !finite(z)) {
            return Err(bad("nonfinite reflection boundary trace"));
        }
        let (radiated_power, radiated_power_roundoff_interval, power_diagnostics) =
            evaluate_radiated_power(&pressure, &velocity, self.surface.areas(), self.medium);
        Ok((
            RadiationSolution {
                pressure,
                velocity,
                k: self.k,
                medium: self.medium,
                surface_fingerprint: surface_fingerprint(self.surface),
                panels_per_wavelength: self.panels_per_wavelength,
                radiated_power,
                radiated_power_roundoff_interval,
                power_diagnostics,
                condition_lower_bound: *condition_lower_bound,
                dense_cap_utilization: n as f64 / MAX_DENSE_PANELS as f64,
            },
            x,
        ))
    }

    fn wall_work(&self, wall_amplitude: &[C64]) -> Result<f64, HelmholtzError> {
        let mut work = ComplexSum::default();
        let mut apparent = 0.0;
        let inverse_rho_c = 1.0 / (self.medium.density * self.medium.sound_speed);
        for (i, (&amplitude, &area)) in wall_amplitude.iter().zip(self.surface.areas()).enumerate()
        {
            // Evaluate the actual passive wall branch, not the subtraction
            // of nearly cancelling incident/scattered velocities. Rigid and
            // pressure-release patches then have exactly zero absorption.
            let p = self.pressure_coefficient[i] * amplitude;
            let v = (self.flow_coefficient[i] * amplitude).scale(inverse_rho_c);
            work.add((p * v.conj()).scale(0.5 * area));
            apparent += 0.5 * p.abs() * v.abs() * area;
        }
        let work = work.finish().re;
        if !work.is_finite()
            || !apparent.is_finite()
            || (!self.declared_active && work < -1e-10 * apparent.max(f64::MIN_POSITIVE))
        {
            return Err(bad("reflection boundary failed passive-work check"));
        }
        Ok(work)
    }
}
