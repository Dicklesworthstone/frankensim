//! Exact elimination of the original independent wall histories at s=-i omega.
//!
//! For each original arm, y=sqrt(g r) v/(s+r), so its resisting effort is
//! g s/(s+r) v. With the SAME signed wall channels V,T and Bv=K^+ V^T,
//! the remaining canonical equations are
//!
//! (s I + Bv Dv V) q - M^-1 p = 0,
//! K q + (s I + T^T Dt T M^-1) p = B U.
//!
//! Dv includes the original Poiseuille drag. Every eliminated history is
//! reconstructed after the solve; the original source then evaluates heat.
//! This is a Schur complement, not a new attenuation law or a reduced basis.

use super::*;
use fs_la::eigen_complex::{LuComplex, lu_complex_with_cancel};
use fs_phs::harmonic::{HarmonicBudget, HarmonicError, MAX_HARMONIC_DIMENSION};

/// Original cylindrical modal source with frequency-domain wall elimination.
/// The dense dimension is only 2*coordinates; wall-pole/channel storage is
/// linear. Source geometry, all acoustic modes and all diffusion arms remain.
#[derive(Debug)]
pub struct ModalViscothermalHarmonic {
    source: ModalGalerkinSource,
    viscous: Vec<f64>,
    thermal: Vec<f64>,
    input: Vec<f64>,
    budget: HarmonicBudget,
}

/// One factorization reused for different coherent physical port flows.
pub struct ModalViscothermalFactor<'a> {
    model: &'a ModalViscothermalHarmonic,
    omega: f64,
    matrix: Vec<C64>,
    row_scales: Vec<f64>,
    column_scales: Vec<f64>,
    lu: LuComplex,
}

/// Peak exp(-i omega t) response in the original source coordinates.
#[derive(Debug, Clone)]
pub struct ModalViscothermalResponse {
    /// Complete original [q, canonical p, viscous histories, thermal histories].
    pub state: Vec<C64>,
    /// Physical nodal modal pressures M^-1 p [Pa].
    pub pressure_pa: Vec<C64>,
    /// Pressure conjugate to each supplied inward modal flow [Pa].
    pub port_pressure_pa: Vec<C64>,
    /// Cycle-mean physical acoustic plus wall-memory storage [J].
    pub stored_energy_j: f64,
    /// Cycle-mean signed work supplied by all inward flows [W].
    pub supplied_power_w: f64,
    /// Independently evaluated original viscous and thermal resistor heat [W].
    pub heat: ModalViscothermalHeat,
    /// Largest componentwise condensed-equation or reconstructed-memory residual.
    pub relative_residual: f64,
}

fn fail(what: &'static str) -> HarmonicError {
    HarmonicError::Invalid { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HarmonicError> {
    if cancelled() {
        Err(HarmonicError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn real_row(row: &[f64], values: &[C64]) -> C64 {
    row.iter()
        .zip(values)
        .fold(C64::ZERO, |sum, (&a, &v)| sum + v.scale(a))
}
fn excess(terms: &[ZkDiffusiveTerm], s: C64) -> C64 {
    terms.iter().fold(C64::ZERO, |sum, term| {
        sum + (s / (s + C64::from_re(term.rate_per_s))).scale(term.gain_per_s)
    })
}

impl ModalGalerkin {
    /// Prepare the original modal viscothermal physics for harmonic queries
    /// without building dense matrices over the independent wall histories.
    ///
    /// Uses exactly the geometry, gas, signed wall projections, kinetic port
    /// lift, positive Bessel arms and thin-layer admission of
    /// [`Self::with_viscothermal`]. Every acoustic coordinate is retained, and
    /// each response reconstructs all original memory coordinates and heat.
    /// The source's spatial discretization is unchanged. This steady harmonic
    /// solve does not advance transients or relax the source's physical domain.
    ///
    /// `budget.max_dimension` bounds the 2*coordinates acoustic solve, rather
    /// than acoustic-plus-wall states. Byte/work bounds include all retained
    /// channels, reconstructed states, and factor/solve scratch. A frequency
    /// outside the admitted wall band refuses. Memory count may exceed the
    /// exact time propagator's cap without dropping poles or modes.
    pub fn into_viscothermal_harmonic(
        self,
        resolution: ModalViscothermalResolution,
        budget: HarmonicBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalViscothermalHarmonic, ModalTimeError> {
        let (model, admission) = prepare_law(&self.source, resolution, &mut cancelled)?;
        let n = self.source.coordinates();
        let dimension = 2 * n;
        let ports = self.source.port_count();
        let channels = self.source.elements + self.source.nodes.len();
        let d = dimension as u128;
        let states = admission.states as u128;
        let poles = resolution.poles as u128;
        let bytes = 192 * d * d
            + 64 * states
            + 64 * d * (ports + channels + 1) as u128
            + 96 * poles * poles;
        let work = 80 * d * d * d
            + 32 * d * d * (ports + channels + 1) as u128
            + 64 * states
            + 80 * poles * poles * poles;
        if dimension > budget.max_dimension
            || dimension > MAX_HARMONIC_DIMENSION
            || bytes > budget.max_dense_bytes as u128
            || work > budget.max_dense_work as u128
            || !budget.relative_residual_tolerance.is_finite()
            || !(0.0..1.0).contains(&budget.relative_residual_tolerance)
            || budget.relative_residual_tolerance == 0.0
        {
            return Err(invalid(
                "condensed modal wall solve exceeds declared resources or accuracy domain",
            ));
        }
        let (source, system) = self.into_parts();
        let (_, resistance, input) = system.structure();
        if system.state_dim() != dimension || resistance.iter().any(|&v| v != 0.0) {
            return Err(invalid(
                "condensed wall solve requires the complete original lossless operator",
            ));
        }
        let input = input.to_vec();
        let signs: Vec<_> = source
            .radial_eigenvalues
            .iter()
            .map(|&gamma| {
                if bessel_j0_certified(gamma) > 0.0 {
                    1.0
                } else {
                    -1.0
                }
            })
            .collect();
        let (velocity, pressure) = wall_channels(
            &source.mass,
            source.modes.0,
            source.nodes.len(),
            source.elements,
            source.profile.length(),
            source.profile.inlet_radius(),
            source.gas.density,
            &signs,
        )?;
        let lifted = lift_velocity(
            &source.stiffness,
            source.modes.0,
            source.end,
            &velocity,
            &mut cancelled,
        )?;
        let mut viscous = vec![0.0; n * n];
        let mut thermal = vec![0.0; n * n];
        for (b, v) in lifted.iter().zip(&velocity) {
            cancelled_at(&mut cancelled, 0)?;
            for i in 0..n {
                for j in 0..n {
                    viscous[i * n + j] += b[i] * v[j];
                }
            }
        }
        for t in &pressure {
            cancelled_at(&mut cancelled, 0)?;
            let physical: Vec<f64> = (0..n)
                .map(|j| (0..n).map(|k| t[k] * source.inverse_mass[k * n + j]).sum())
                .collect();
            for i in 0..n {
                for j in 0..n {
                    thermal[i * n + j] += t[i] * physical[j];
                }
            }
        }
        if viscous.iter().chain(&thermal).any(|v| !v.is_finite()) {
            return Err(invalid(
                "condensed modal wall projection is unrepresentable",
            ));
        }
        let mut source = source;
        source.viscothermal = Some(ViscothermalSource {
            model,
            admission,
            velocity,
            pressure,
            states: admission.states,
        });
        cancelled_at(&mut cancelled, source.elements)?;
        Ok(ModalViscothermalHarmonic {
            source,
            viscous,
            thermal,
            input,
            budget,
        })
    }
}

impl ModalViscothermalHarmonic {
    /// Original source, including every wall arm, signed channel and admission.
    pub const fn source(&self) -> &ModalGalerkinSource {
        &self.source
    }

    /// Dense solve size. The complete reconstructed state can be much larger.
    pub fn solve_dimension(&self) -> usize {
        2 * self.source.coordinates()
    }

    /// Factor within the source's admitted physical band, with the existing
    /// cancellable complex LU owner and independent row/column equilibration.
    pub fn factor(
        &self,
        omega: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalViscothermalFactor<'_>, HarmonicError> {
        poll(&mut cancelled)?;
        let wall = self
            .source
            .viscothermal
            .as_ref()
            .expect("prepared wall source");
        if !omega.is_finite()
            || !(wall.model.band_rad_s()[0]..=wall.model.band_rad_s()[1]).contains(&omega)
        {
            return Err(fail(
                "condensed modal frequency is outside the original wall band",
            ));
        }
        let s = C64::new(0.0, -omega);
        let dv =
            C64::from_re(wall.model.series_drag_per_s()) + excess(wall.model.series_terms(), s);
        let dt = excess(wall.model.shunt_terms(), s);
        let n = self.source.coordinates();
        let d = 2 * n;
        let mut matrix = vec![C64::ZERO; d * d];
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                matrix[i * d + j] = dv.scale(self.viscous[i * n + j]);
                matrix[i * d + n + j] = C64::from_re(-self.source.inverse_mass[i * n + j]);
                matrix[(n + i) * d + j] = C64::from_re(self.source.stiffness[i * n + j]);
                matrix[(n + i) * d + n + j] = dt.scale(self.thermal[i * n + j]);
            }
            matrix[i * d + i] = matrix[i * d + i] + s;
            matrix[(n + i) * d + n + i] = matrix[(n + i) * d + n + i] + s;
        }
        if matrix.iter().any(|&v| !finite(v)) {
            return Err(fail("condensed modal operator overflow"));
        }
        let mut scaled = matrix.clone();
        let mut row_scales = vec![0.0; d];
        let mut column_scales = vec![0.0; d];
        for (i, row) in scaled.chunks_exact_mut(d).enumerate() {
            poll(&mut cancelled)?;
            let scale = row.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            if scale == 0.0 {
                return Err(HarmonicError::Singular { index: i });
            }
            if !scale.is_finite() {
                return Err(fail("condensed modal row scale overflow"));
            }
            row_scales[i] = scale;
            for v in row {
                *v = C64::new(v.re / scale, v.im / scale);
            }
        }
        for j in 0..d {
            poll(&mut cancelled)?;
            let scale = (0..d)
                .map(|i| scaled[i * d + j].abs())
                .fold(0.0_f64, f64::max);
            if scale == 0.0 {
                return Err(HarmonicError::Singular { index: j });
            }
            column_scales[j] = scale;
            for i in 0..d {
                let value = scaled[i * d + j];
                scaled[i * d + j] = C64::new(value.re / scale, value.im / scale);
            }
        }
        let lu = lu_complex_with_cancel(&scaled, d, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(ModalViscothermalFactor {
            model: self,
            omega,
            matrix,
            row_scales,
            column_scales,
            lu,
        })
    }
}

impl ModalViscothermalFactor<'_> {
    /// Solve coherent modal volume flows INTO every physical face port.
    /// Reuses the factor and restores every original wall memory analytically;
    /// source heat is evaluated from those memories, never an energy deficit.
    pub fn response(
        &self,
        input: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ModalViscothermalResponse, HarmonicError> {
        poll(&mut cancelled)?;
        let source = &self.model.source;
        let n = source.coordinates();
        let d = 2 * n;
        let ports = source.port_count();
        if input.len() != ports || input.iter().any(|&v| !finite(v)) {
            return Err(fail(
                "one finite inward modal flow per original port is required",
            ));
        }
        let rhs: Vec<_> = self
            .model
            .input
            .chunks_exact(ports)
            .map(|row| real_row(row, input))
            .collect();
        let mut state = self.solve(&rhs, &mut cancelled)?;
        let mut correction = vec![C64::ZERO; d];
        let mut relative = 0.0_f64;
        for pass in 0..=4 {
            relative = 0.0;
            for i in 0..d {
                poll(&mut cancelled)?;
                let mut residual = -rhs[i];
                let mut denominator = rhs[i].abs();
                for (j, &value) in state.iter().enumerate() {
                    let a = self.matrix[i * d + j];
                    residual = residual + a * value;
                    denominator += a.abs() * value.abs();
                }
                if !finite(residual) || !denominator.is_finite() {
                    return Err(fail("condensed modal equation residual overflow"));
                }
                relative = relative.max(if denominator == 0.0 {
                    0.0
                } else {
                    residual.abs() / denominator
                });
                correction[i] = -residual;
            }
            if relative <= self.model.budget.relative_residual_tolerance || pass == 4 {
                break;
            }
            let delta = self.solve(&correction, &mut cancelled)?;
            for (x, dx) in state.iter_mut().zip(delta) {
                *x = *x + dx;
            }
        }
        let pressure_pa: Vec<_> = source
            .inverse_mass
            .chunks_exact(n)
            .map(|row| real_row(row, &state[n..]))
            .collect();
        let wall = source.viscothermal.as_ref().expect("prepared wall source");
        let s = C64::new(0.0, -self.omega);
        state.reserve(source.state_dimension() - d);
        for (channels, field, terms) in [
            (
                &wall.velocity,
                state[..n].to_vec(),
                wall.model.series_terms(),
            ),
            (
                &wall.pressure,
                pressure_pa.clone(),
                wall.model.shunt_terms(),
            ),
        ] {
            let values: Vec<_> = channels.iter().map(|row| real_row(row, &field)).collect();
            for term in terms {
                poll(&mut cancelled)?;
                let coupling = det::sqrt(term.gain_per_s) * det::sqrt(term.rate_per_s);
                let a = s + C64::from_re(term.rate_per_s);
                for value in &values {
                    let b = value.scale(coupling);
                    let memory = b / a;
                    let residual = a * memory - b;
                    let scale = a.abs() * memory.abs() + b.abs();
                    if !finite(memory) || !finite(residual) || !scale.is_finite() {
                        return Err(fail("original modal wall history is unrepresentable"));
                    }
                    if scale != 0.0 {
                        relative = relative.max(residual.abs() / scale);
                    }
                    state.push(memory);
                }
            }
        }
        if relative > self.model.budget.relative_residual_tolerance {
            return Err(HarmonicError::Residual {
                relative,
                tolerance: self.model.budget.relative_residual_tolerance,
            });
        }
        let heat = source
            .viscothermal_mean_heat_w(&state)
            .map_err(|_| fail("original modal wall heat source refused reconstructed histories"))?;
        let mut effort = vec![C64::ZERO; d];
        for i in 0..n {
            effort[i] = real_row(&source.stiffness[i * n..(i + 1) * n], &state[..n]);
            effort[n + i] = pressure_pa[i];
        }
        let mut stored_energy_j = (0..d)
            .map(|i| (state[i].conj() * effort[i]).re)
            .sum::<f64>();
        stored_energy_j += state[d..].iter().map(|v| v.norm_sq()).sum::<f64>();
        stored_energy_j *= 0.25;
        let port_pressure_pa: Vec<_> = (0..ports)
            .map(|p| {
                (0..d).fold(C64::ZERO, |sum, i| {
                    sum + effort[i].scale(self.model.input[i * ports + p])
                })
            })
            .collect();
        let supplied_power_w = 0.5
            * port_pressure_pa
                .iter()
                .zip(input)
                .map(|(&p, &u)| (p * u.conj()).re)
                .sum::<f64>();
        if !stored_energy_j.is_finite()
            || stored_energy_j < 0.0
            || !supplied_power_w.is_finite()
            || pressure_pa
                .iter()
                .chain(&port_pressure_pa)
                .chain(&state)
                .any(|&v| !finite(v))
        {
            return Err(fail(
                "original modal pressure, work or storage is unrepresentable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(ModalViscothermalResponse {
            state,
            pressure_pa,
            port_pressure_pa,
            stored_energy_j,
            supplied_power_w,
            heat,
            relative_residual: relative,
        })
    }

    fn solve(
        &self,
        rhs: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HarmonicError> {
        let mut result: Vec<_> = rhs
            .iter()
            .zip(&self.row_scales)
            .map(|(&v, &scale)| C64::new(v.re / scale, v.im / scale))
            .collect();
        if result.iter().any(|&v| !finite(v)) {
            return Err(fail("modal flow right-hand side overflow"));
        }
        self.lu.solve_with_cancel(&mut result, &mut *cancelled)?;
        for (v, &scale) in result.iter_mut().zip(&self.column_scales) {
            *v = C64::new(v.re / scale, v.im / scale);
        }
        if result.iter().any(|&v| !finite(v)) {
            return Err(fail("condensed modal solution overflow"));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
