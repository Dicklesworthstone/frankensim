//! Driven response of the SAME retained pressure/displacement/material pencil.
//! One assembly and LU per real frequency, shared across all admitted loads.
use super::*;

/// Cold solve/output limits, additional to the model's boundary and state caps.
#[derive(Debug, Clone, Copy)]
pub struct BemHarmonicOptions {
    /// Nonempty force batch, at most 256 complete structural vectors.
    pub max_load_cases: usize,
    /// Logical response payload ceiling: cases*(32*full_dimension+512) bytes.
    /// Solver scratch remains under the model's separately declared byte cap.
    pub max_response_bytes: usize,
    /// Conservative LU, back-substitution and audit units:
    /// 32*dimension^3 + 128*cases*dimension^2, excluding BEM assembly.
    pub max_work: u64,
    /// Relative closure of input=Kelvin+Maxwell+surface work, normalized by
    /// the sum of absolute complex-work magnitudes, not tiny real input alone.
    /// This is a numerical balance gate, not a continuum passivity certificate.
    pub power_relative_tolerance: f64,
}

/// Cycle-average powers for PEAK complex amplitudes under exp(-i omega t).
/// Raw signed surface work is retained; a negative BEM radiation artifact is
/// never clipped or repaired, even when the discrete equations balance.
#[derive(Debug, Clone, Copy)]
pub struct BemHarmonicPower {
    pub input_w: f64,
    pub kelvin_w: f64,
    pub relaxation_w: f64,
    pub radiated_w: f64,
    pub balance_w: f64,
    pub apparent_power_scale_w: f64,
    pub relative_balance: f64,
}

/// Actual force-scaled physical response, NOT eigenvector-normalized values.
/// Source modal coordinates, skin panel order and material branch order stay
/// unchanged. No finite-distance receiver or acoustic storage is inferred.
#[derive(Debug, Clone)]
pub struct BemHarmonicResponse {
    pub angular_frequency_rad_s: f64,
    pub modal_displacement: Vec<C64>,
    pub modal_velocity: Vec<C64>,
    pub surface_pressure_pa: Vec<C64>,
    pub surface_normal_velocity_m_s: Vec<C64>,
    pub elastic_memory_sqrt_j: Vec<C64>,
    pub mechanical_residual: f64,
    pub boundary_residual: f64,
    pub memory_residual: f64,
    pub power: BemHarmonicPower,
}

impl BemLoadedModes {
    /// Solve every structural force against one ORIGINAL full-state operator.
    /// Fixed two-sided balancing changes only the linear solve. Returned states
    /// must pass the same unscaled componentwise row checks as resonance poles,
    /// now including their physical RHS, and an independent work/heat balance.
    ///
    /// All forcing columns share one BEM assembly and one complex LU. No acoustic
    /// or material state elimination, fitting, diagonal damping replacement or
    /// frequency interpolation. An exact lossless pole may be singular and is
    /// refused, not broadened by artificial resistance. Failed batches return
    /// no partial responses; the immutable model remains reusable.
    pub fn harmonic_response_batch(
        &self, omega: f64, forces: &[&[C64]], options: &BemHarmonicOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<BemHarmonicResponse>, BemModesError> {
        poll(&mut cancelled)?;
        let (m, p, n) = (self.mode_count(), self.panel_count(), self.dimension());
        let cases = forces.len();
        if !omega.is_finite() || omega <= 0.0 || cases == 0
            || options.max_load_cases == 0 || options.max_load_cases > 256
            || cases > options.max_load_cases
            || !options.power_relative_tolerance.is_finite()
            || options.power_relative_tolerance <= 0.0 || options.power_relative_tolerance >= 1.0
            || forces.iter().any(|f| f.len() != m || f.iter().any(|&z| !finite(z)))
        {
            return Err(bad("harmonic response needs a positive frequency, complete finite force batch and valid limits"));
        }
        let bytes = cases as u128*(32*n as u128+512);
        let work = 32*(n as u128).pow(3)+128*cases as u128*(n as u128).pow(2);
        if bytes > options.max_response_bytes as u128 || work > options.max_work as u128 {
            return Err(bad("harmonic response exceeds its declared output or solve-work budget"));
        }
        let raw = self.assemble(C64::from_re(omega), &mut 0, &mut cancelled)?;
        let mut balanced = raw.value.matrix.clone();
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                let at = i*n+j;
                balanced[at] = balanced[at].scale(self.column_scales[j]/self.row_scales[i]);
                if !finite(balanced[at]) {
                    return Err(bad("balanced response operator is not representable"));
                }
            }
        }
        let lu = lu_complex_with_cancel(&balanced, n, &mut cancelled)
            .map_err(BemModesError::ResponseSolve)?;
        drop(balanced);
        let mut results = Vec::with_capacity(cases);
        for &force in forces {
            poll(&mut cancelled)?;
            let mut state = vec![C64::ZERO; n];
            for i in 0..m { state[i] = divide(force[i], self.row_scales[i]); }
            lu.solve_with_cancel(&mut state, &mut cancelled).map_err(BemModesError::ResponseSolve)?;
            for (x,&scale) in state.iter_mut().zip(&self.column_scales) { *x = x.scale(scale); }
            if state.iter().any(|&z| !finite(z)) {
                return Err(bad("physical harmonic response overflow"));
            }
            let residuals = self.check_physical_rows(&raw, &state, Some(force), &mut cancelled)?;
            let velocity: Vec<_> = state[..m].iter().map(|&x| C64::new(0.0,-omega)*x).collect();
            let mut normal_velocity = vec![C64::ZERO; p];
            for (i,out) in normal_velocity.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                for j in 0..m { *out = *out+velocity[j].scale(self.weights[j][i]); }
            }
            if velocity.iter().chain(&normal_velocity).any(|&z| !finite(z)) {
                return Err(bad("physical harmonic velocity overflow"));
            }
            let power = self.response_power(force, &velocity, &state[m..m+p],
                &normal_velocity, &state[m+p..], options.power_relative_tolerance, &mut cancelled)?;
            let memories = state.split_off(m+p);
            let pressure = state.split_off(m);
            results.push(BemHarmonicResponse {
                angular_frequency_rad_s: omega,
                modal_displacement: state,
                modal_velocity: velocity,
                surface_pressure_pa: pressure,
                surface_normal_velocity_m_s: normal_velocity,
                elastic_memory_sqrt_j: memories,
                mechanical_residual: residuals[0], boundary_residual: residuals[1], memory_residual: residuals[2],
                power,
            });
        }
        poll(&mut cancelled)?;
        Ok(results)
    }

    fn response_power(
        &self, force: &[C64], velocity: &[C64], pressure: &[C64], normal: &[C64], memory: &[C64],
        tolerance: f64, cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<BemHarmonicPower, BemModesError> {
        // One common amplitude rescaling preserves all power products. This
        // avoids squaring an extreme physical drive before computing its audit.
        // It is not a comparison or norm between unlike physical quantities.
        let amplitude = force.iter().map(|z| z.abs()).fold(0.0, f64::max);
        let amplitude = if amplitude == 0.0 { 1.0 } else { amplitude };
        let normalized = |x: &[C64]| -> Result<Vec<C64>, BemModesError> {
            let mut result = Vec::with_capacity(x.len());
            for &z in x {
                let v = divide(z, amplitude);
                if !finite(v) || (z.re != 0.0 && v.re == 0.0) || (z.im != 0.0 && v.im == 0.0) {
                    return Err(bad("harmonic work normalization is not representable"));
                }
                result.push(v);
            }
            Ok(result)
        };
        let f = normalized(force)?;
        let v = normalized(velocity)?;
        let p = normalized(pressure)?;
        let u = normalized(normal)?;
        let e = normalized(memory)?;
        let mut sums = [0.0; 5];
        let mut corrections = [0.0; 5];
        let mut add = |channel: usize, value: f64| {
            let adjusted = value-corrections[channel];
            let next = sums[channel]+adjusted;
            corrections[channel] = (next-sums[channel])-adjusted;
            sums[channel] = next;
        };
        let m = self.mode_count();
        for i in 0..m {
            poll(cancelled)?;
            add(0, 0.5*(f[i]*v[i].conj()).re);
            add(4, 0.5*f[i].abs()*v[i].abs());
            for j in 0..m {
                let c = self.viscosity[i*m+j];
                add(1, 0.5*(v[i].conj()*v[j]).re*c);
                add(4, 0.5*v[i].abs()*v[j].abs()*c.abs());
            }
        }
        for (&x,&tau) in e.iter().zip(&self.memory.times) {
            poll(cancelled)?;
            let heat = 0.5*x.norm_sq()/tau;
            add(2, heat);
            add(4, heat);
        }
        for ((&p,&u),&area) in p.iter().zip(&u).zip(self.surface.areas()) {
            poll(cancelled)?;
            add(3, 0.5*(p*u.conj()).re*area);
            add(4, 0.5*p.abs()*u.abs()*area);
        }
        let balance = sums[0]-sums[1]-sums[2]-sums[3];
        let relative = if sums[4] == 0.0 { balance.abs() } else { balance.abs()/sums[4] };
        if sums.iter().any(|x| !x.is_finite()) || !relative.is_finite() {
            return Err(bad("harmonic power accounting overflow"));
        }
        if relative > tolerance {
            return Err(BemModesError::PowerBalance { relative, tolerance });
        }
        let restore = |x: f64| -> Result<f64, BemModesError> {
            let value = (x*amplitude)*amplitude;
            if !value.is_finite() || (x != 0.0 && value == 0.0) {
                return Err(bad("physical harmonic power is not representable"));
            }
            Ok(value)
        };
        Ok(BemHarmonicPower {
            input_w: restore(sums[0])?, kelvin_w: restore(sums[1])?,
            relaxation_w: restore(sums[2])?, radiated_w: restore(sums[3])?,
            balance_w: restore(balance)?, apparent_power_scale_w: restore(sums[4])?,
            relative_balance: relative,
        })
    }
}

#[cfg(test)]
#[path = "harmonic_tests.rs"]
mod tests;
