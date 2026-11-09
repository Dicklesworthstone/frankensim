//! Exact affine loading of the original constrained Kelvin/Maxwell history.
//!
//! The existing preparation supplies orthonormal energy modes with
//! `z_i' = -r_i z_i + h_i . u'`. Each original Prony ramp supplies its endpoint,
//! mean and independent nonnegative integral `r_i * integral(z_i² dt)`.
//! Reconstructed ORIGINAL tensor stresses supply work and both load tangents.
//! The residual Kelvin Gram supplies the remaining positive heat directly.

use super::*;
use crate::visco::PronyRamp;
use crate::visco_tensor::WEIGHTS;

/// Accepted affine in-plane loading; all energies are per original volume.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorKelvinPlaneStressRampStep {
    /// Full physical endpoint stress [Pa], with zero transverse traction.
    pub stress: Voigt,
    /// Full physical endpoint strain, retaining normal relaxation.
    pub strain: Voigt,
    /// Integrated strain change before addition to the old rounded strain.
    pub strain_increment: Voigt,
    /// Compatible endpoint strain rate [1/s].
    pub strain_rate: Voigt,
    /// Exact interval-mean physical stress [Pa], including in-plane reactions.
    pub mean_stress: Voigt,
    /// Original equilibrium and Maxwell spring storage [J/m³].
    pub stored_energy: f64,
    /// Independently evaluated local quadratic storage change [J/m³].
    pub stored_energy_change: f64,
    /// Signed original stress work [J/m³], including all imposed strain ports.
    pub work: f64,
    /// Positive original Kelvin plus Maxwell heat integral [J/m³].
    pub dissipated: f64,
    /// Work minus storage change minus physical heat; never charged as heat.
    pub balance_residual: f64,
    /// Conservative full-path component bounds in physical tensor-shear order.
    /// These are numerical sufficient bounds, not outward-rounded certificates.
    pub strain_bounds: [Voigt; 2],
}

/// Exact derivatives with respect to the three endpoint engineering strains.
///
/// Columns are `[epsilon_aa, epsilon_bb, gamma_ab]` in the prepared material
/// plane; rows use the full physical `[xx, yy, zz, xy, yz, zx]` tensor-shear
/// convention. Initial history, coefficients and duration are held fixed.
/// These are small source maps for implicit loading, not parameter derivatives
/// or an admission of a finite load outside the original material domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorKelvinPlaneStressRampTangent {
    /// Full endpoint strain derivative.
    pub strain: [[f64; 3]; 6],
    /// Full endpoint stress derivative [Pa].
    pub stress: [[f64; 3]; 6],
    /// Full interval-mean stress derivative [Pa].
    pub mean_stress: [[f64; 3]; 6],
    /// Compatible endpoint strain-rate derivative [1/s].
    pub strain_rate: [[f64; 3]; 6],
}

impl PreparedTensorKelvinPlaneStress {
    fn ramp_coefficients(&self, duration: f64) -> Result<Vec<PronyRamp>> {
        if !duration.is_finite() || duration <= 0.0 {
            return Err(invalid(
                "Kelvin plane-stress ramps require positive finite duration",
            ));
        }
        self.rates
            .iter()
            .map(|&rate| {
                let r = duration * rate;
                if !r.is_finite() || r <= 0.0 {
                    return Err(invalid(
                        "unrepresentable Kelvin plane-stress ramp relaxation number",
                    ));
                }
                Ok(PronyRamp::new(r))
            })
            .collect()
    }

    // Recover the stress from original spring coordinates, independently of
    // prepared modal input/output identities. `rate` is a full Mandel rate.
    fn ramp_stress(&self, strain: &Voigt, q: &[Voigt], rate: &[Dd; 6]) -> Result<Voigt> {
        let (spring, mut scale) = self.ramp_spring(strain, q)?;
        let viscosity = self.model.kv_pa_s().expect("admitted Kelvin tensor");
        let tolerance = 8192.0 * f64::EPSILON * (self.state_coordinates() + 6) as f64;
        let mut stress = [0.0; 6];
        for i in 0..6 {
            let mut value = spring[i];
            for j in 0..6 {
                let term = Dd::from_f64(viscosity[i][j]) * rate[j];
                value = value + term;
                scale[i] += term.abs().to_f64();
            }
            stress[i] = number(value)?.to_f64();
            if self.outside.contains(&i) {
                if !scale[i].is_finite() || stress[i].abs() > tolerance * scale[i] {
                    return Err(invalid(
                        "original Kelvin ramp plane-stress constraint is unresolved",
                    ));
                }
                stress[i] = 0.0;
            }
        }
        Ok(stress)
    }

    fn ramp_spring(&self, strain: &Voigt, q: &[Voigt]) -> Result<([Dd; 6], Voigt)> {
        let eq = apply(&self.model.equilibrium_gram, strain);
        let mut spring = [Dd::ZERO; 6];
        let mut scale = [0.0; 6];
        for i in 0..6 {
            for row in 0..6 {
                let term =
                    Dd::from_f64(self.model.equilibrium_gram[row][i]) * Dd::from_f64(eq[row]);
                spring[i] = spring[i] + term;
                scale[i] += term.abs().to_f64();
            }
            for (j, branch) in self.model.branches().iter().enumerate() {
                let gram = branch.resolved_gram_rows(&self.model.factors[j]);
                for row in 0..6 {
                    let term = Dd::from_f64(gram[row][i]) * Dd::from_f64(q[j][row]);
                    spring[i] = spring[i] + term;
                    scale[i] += term.abs().to_f64();
                }
            }
            number(spring[i])?;
            if !scale[i].is_finite() {
                return Err(invalid("unrepresentable original Kelvin spring stress"));
            }
        }
        Ok((spring, scale))
    }

    fn ramp_endpoint_rate(
        &self,
        strain: &Voigt,
        q: &[Voigt],
        mut rate: [Dd; 6],
    ) -> Result<[Dd; 6]> {
        let (spring, _) = self.ramp_spring(strain, q)?;
        let free = if let Some(constraint) = &self.constraint {
            constraint.endpoint_rate(self, &spring, &rate, q)?
        } else {
            let viscosity = self.model.kv_pa_s().expect("admitted Kelvin tensor");
            let rhs = self.outside.map(|axis| {
                let mut value = -spring[axis];
                for &j in &self.inside {
                    value = value - Dd::from_f64(viscosity[axis][j]) * rate[j];
                }
                value
            });
            let root = self
                .viscosity
                .as_ref()
                .expect("definite transverse Kelvin factor");
            root.solve(root.solve(rhs, true)?, false)?
        };
        for i in 0..3 {
            rate[self.outside[i]] = number(free[i])?;
        }
        Ok(rate)
    }

    // State increments are lifted into the ORIGINAL normal strains and branch
    // coordinates; semidefinite modes already live in the algebraic nullspace.
    fn ramp_state_increment(&self, amplitudes: &[f64]) -> Result<Vec<Dd>> {
        let n = self.state_coordinates();
        let m = self.rates.len();
        let mut result = vec![Dd::ZERO; n];
        for row in 0..n {
            let mut sum = Sum::default();
            for mode in 0..m {
                sum.add([self.modes[row * m + mode], amplitudes[mode], 1.0, 1.0]);
            }
            result[row] = number(sum.retained_value())?;
        }
        Ok(result)
    }

    /// Advance a hold, reversal or arbitrary affine engineering strain ramp.
    ///
    /// The ramp starts at the accepted in-plane strain; the supplied three
    /// values are its endpoint. Duration must be positive: this API introduces
    /// no instantaneous strain jump. The original transverse rate enforces
    /// TOTAL zero traction, including viscosity, throughout the path. Positive
    /// transverse viscosity permits a traction/rate switch without a memory
    /// reset; algebraic null directions must already satisfy spring traction.
    ///
    /// This reuses the native constrained modes and scalar exact-ramp owner.
    /// Each step is O(n²), n=3+6*branches, with no new matrix factorization and
    /// storage bounded by the admitted preparation. Full original history and
    /// source strain limits survive. Refusal/cancellation changes no field.
    /// Coefficients and temperature are fixed; an arbitrary ramp provides no
    /// frequency-band calibration claim for a law fitted on a finite band.
    #[allow(clippy::too_many_lines)]
    pub fn advance(
        &self,
        history: &mut TensorPronyHistory,
        engineering: [f64; 3],
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorKelvinPlaneStressRampStep> {
        poll(&mut cancelled)?;
        if !Arc::ptr_eq(&self.model, &history.model) || engineering.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "Kelvin ramp requires its original history and finite engineering strains",
            ));
        }
        let ramps = self.ramp_coefficients(duration_s)?;
        if let Some(constraint) = &self.constraint {
            constraint.admit_history(self, history)?;
        }
        let n = self.state_coordinates();
        let m = self.rates.len();
        let old_u =
            core::array::from_fn::<_, 3, _>(|j| history.strain[self.inside[j]] / ENGINEERING[j]);
        let target_inside = core::array::from_fn::<_, 3, _>(|j| engineering[j] * ENGINEERING[j]);
        // Subtract in the ORIGINAL state coordinates before changing units.
        // A representable held endpoint must not acquire a spurious viscous
        // impulse from an engineering-coordinate round trip at tiny duration.
        let inside_delta = core::array::from_fn::<_, 3, _>(|j| {
            Dd::from_f64(target_inside[j]) - Dd::from_f64(history.strain[self.inside[j]])
        });
        let du = core::array::from_fn::<_, 3, _>(|j| {
            (inside_delta[j] / Dd::from_f64(ENGINEERING[j])).to_f64()
        });
        if du.iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable Kelvin strain increment"));
        }
        let mut lag = self.outside.map(|i| Dd::from_f64(history.strain[i]));
        for i in 0..3 {
            for j in 0..3 {
                lag[i] =
                    lag[i] - Dd::from_f64(self.map[j][self.outside[i]]) * Dd::from_f64(old_u[j]);
            }
        }
        let mut state = vec![Dd::ZERO; n];
        state[..3].copy_from_slice(&self.equilibrium.apply(lag)?);
        for row in 3..n {
            state[row] = Dd::from_f64(history.q[(row - 3) / 6][(row - 3) % 6]);
        }
        let mut increments = zeros(m)?;
        let mut mean_increments = zeros(m)?;
        let mut heat = Sum::default();
        let mut lower = [0.0; 6];
        let mut upper = [0.0; 6];
        let mut mapped: [Dd; 6] = core::array::from_fn(|i| {
            (0..3).fold(Dd::ZERO, |s, j| {
                s + Dd::from_f64(self.map[j][i]) * Dd::from_f64(du[j])
            })
        });
        for j in 0..3 {
            mapped[self.inside[j]] = inside_delta[j];
        }
        for i in 0..6 {
            lower[i] = mapped[i].to_f64().min(0.0);
            upper[i] = mapped[i].to_f64().max(0.0);
        }
        for mode in 0..m {
            poll(&mut cancelled)?;
            let mut initial = Sum::default();
            for row in 0..n {
                add(&mut initial, state[row], self.modes[row * m + mode]);
            }
            let initial = number(initial.retained_value())?.to_f64();
            let drive = number((0..3).fold(Dd::ZERO, |s, j| {
                s + Dd::from_f64(self.forcing[mode][j]) * Dd::from_f64(du[j])
            }))?
            .to_f64();
            let (_, _, loss) = ramps[mode].apply(initial, drive);
            increments[mode] = ramps[mode].spring_increment(initial, drive);
            let c = ramps[mode].apply(0.0, 1.0).1;
            // 1 - phi1(r) = r*phi2(r), using the native small-r series.
            mean_increments[mode] = c * drive - ((duration_s * self.rates[mode]) * c) * initial;
            if ![loss, increments[mode], mean_increments[mode]]
                .iter()
                .all(|v| v.is_finite())
                || loss < 0.0
            {
                return Err(invalid("unrepresentable Kelvin ramp mode or physical heat"));
            }
            heat.add([loss, 1.0, 1.0, 1.0]);
            // Every driven scalar mode is monotone during an affine ramp.
            // Summing its endpoint extrema bounds the whole physical path.
            for i in 0..6 {
                let excursion = self.strain_modes[mode][i] * increments[mode];
                lower[i] += excursion.min(0.0);
                upper[i] += excursion.max(0.0);
            }
        }
        let bounds = [
            core::array::from_fn(|i| history.strain[i] + lower[i]),
            core::array::from_fn(|i| history.strain[i] + upper[i]),
        ];
        let maximum = core::array::from_fn(|i| bounds[0][i].abs().max(bounds[1][i].abs()));
        self.model.validate_mandel_strain(&maximum).map_err(|_| {
            invalid("full Kelvin ramp strain envelope exceeds the original material domain")
        })?;
        let change = self.ramp_state_increment(&increments)?;
        let mean_change = self.ramp_state_increment(&mean_increments)?;
        let normal = self
            .equilibrium
            .solve(core::array::from_fn(|i| change[i]), false)?;
        let normal_mean = self
            .equilibrium
            .solve(core::array::from_fn(|i| mean_change[i]), false)?;
        let mut delta = mapped;
        let mut mean_delta = mapped.map(|v| Dd::from_f64(0.5) * v);
        for i in 0..3 {
            delta[self.outside[i]] = delta[self.outside[i]] + normal[i];
            mean_delta[self.outside[i]] = mean_delta[self.outside[i]] + normal_mean[i];
        }
        let mut target = [0.0; 6];
        let mut mean_strain = [0.0; 6];
        for i in 0..6 {
            target[i] = number(Dd::from_f64(history.strain[i]) + delta[i])?.to_f64();
            mean_strain[i] = number(Dd::from_f64(history.strain[i]) + mean_delta[i])?.to_f64();
        }
        for j in 0..3 {
            target[self.inside[j]] = target_inside[j];
        }
        self.model
            .validate_mandel_strain(&target)
            .map_err(|_| invalid("Kelvin ramp endpoint exceeds the original strain domain"))?;
        let mut q = history.q.clone();
        let mut mean_q = q.clone();
        for row in 3..n {
            q[(row - 3) / 6][(row - 3) % 6] = number(state[row] + change[row])?.to_f64();
            mean_q[(row - 3) / 6][(row - 3) % 6] = number(state[row] + mean_change[row])?.to_f64();
        }
        let mean_rate = delta.map(|v| v / Dd::from_f64(duration_s));
        let mut known_rate = [Dd::ZERO; 6];
        for j in 0..3 {
            known_rate[self.inside[j]] = inside_delta[j] / Dd::from_f64(duration_s);
        }
        let rate = self.ramp_endpoint_rate(&target, &q, known_rate)?;
        let stress = self.ramp_stress(&target, &q, &rate)?;
        let mean_stress = self.ramp_stress(&mean_strain, &mean_q, &mean_rate)?;
        let work = number((0..3).fold(Dd::ZERO, |s, j| {
            s + Dd::from_f64(mean_stress[self.inside[j]]) * delta[self.inside[j]]
        }))?;
        for row in &self.direct_gram {
            let channel = number((0..3).fold(Dd::ZERO, |s, j| {
                s + Dd::from_f64(row[j]) * Dd::from_f64(du[j])
            }))?;
            add(&mut heat, channel * channel / Dd::from_f64(duration_s), 1.0);
        }
        let eq0 = apply(&self.model.equilibrium_gram, &history.strain);
        let eq1 = apply(&self.model.equilibrium_gram, &target);
        let eq_delta = self
            .model
            .equilibrium_gram
            .map(|row| (0..6).fold(Dd::ZERO, |s, i| s + Dd::from_f64(row[i]) * delta[i]));
        let mut storage_change = Dd::ZERO;
        let mut exchange = 0.0;
        let mut stored = 0.5 * dot(&eq1, &eq1);
        for i in 0..n + 3 {
            let (old, increment) = if i < 6 {
                (eq0[i], eq_delta[i])
            } else {
                (history.q[(i - 6) / 6][(i - 6) % 6], change[i - 3])
            };
            let term = (Dd::from_f64(old) + Dd::from_f64(0.5) * increment) * increment;
            storage_change = storage_change + term;
            exchange += term.abs().to_f64();
        }
        for memory in &q {
            stored += 0.5 * dot(memory, memory);
        }
        let heat = number(heat.retained_value())?;
        let residual = number(work - storage_change - heat)?.to_f64();
        let total_work = number(Dd::from_f64(history.work) + work)?.to_f64();
        let total_heat = number(Dd::from_f64(history.heat) + heat)?.to_f64();
        let tolerance = 8192.0 * f64::EPSILON * (n + 6) as f64;
        if q.iter().any(|v| !finite(v))
            || !stored.is_finite()
            || !exchange.is_finite()
            || stored < 0.0
            || heat.lt(Dd::ZERO)
            || residual.abs() > tolerance * (exchange + work.abs().to_f64() + heat.to_f64())
        {
            return Err(invalid(
                "independent Kelvin ramp work/storage/heat balance is unresolved",
            ));
        }
        let result = TensorKelvinPlaneStressRampStep {
            stress: tensor(&stress),
            strain: tensor(&target),
            strain_increment: tensor(&delta.map(Dd::to_f64)),
            strain_rate: tensor(&rate.map(Dd::to_f64)),
            mean_stress: tensor(&mean_stress),
            stored_energy: stored,
            stored_energy_change: storage_change.to_f64(),
            work: work.to_f64(),
            dissipated: heat.to_f64(),
            balance_residual: residual,
            strain_bounds: bounds.map(|v| tensor(&v)),
        };
        poll(&mut cancelled)?;
        history.strain = target;
        history.q = q;
        history.stress = stress;
        history.stored = stored;
        history.work = total_work;
        history.heat = total_heat;
        Ok(result)
    }

    /// Consistent endpoint and interval-mean load maps for an affine ramp.
    /// The three analytical columns reuse the prepared modes and original
    /// tensor stress/rate reconstruction. No perturbed history is advanced.
    /// The same positive duration and representability refusals apply.
    pub fn ramp_tangent(
        &self,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorKelvinPlaneStressRampTangent> {
        poll(&mut cancelled)?;
        let ramps = self.ramp_coefficients(duration_s)?;
        let mut answer = TensorKelvinPlaneStressRampTangent {
            strain: [[0.0; 3]; 6],
            stress: [[0.0; 3]; 6],
            mean_stress: [[0.0; 3]; 6],
            strain_rate: [[0.0; 3]; 6],
        };
        for col in 0..3 {
            poll(&mut cancelled)?;
            let increments: Vec<_> = ramps
                .iter()
                .enumerate()
                .map(|(mode, ramp)| ramp.apply(0.0, self.forcing[mode][col]).0)
                .collect();
            let means: Vec<_> = ramps
                .iter()
                .enumerate()
                .map(|(mode, ramp)| ramp.apply(0.0, self.forcing[mode][col]).1)
                .collect();
            let change = self.ramp_state_increment(&increments)?;
            let mean_change = self.ramp_state_increment(&means)?;
            let free = self
                .equilibrium
                .solve(core::array::from_fn(|i| change[i]), false)?;
            let mean_free = self
                .equilibrium
                .solve(core::array::from_fn(|i| mean_change[i]), false)?;
            let mut strain = self.map[col].map(Dd::from_f64);
            let mut mean_strain = strain.map(|v| Dd::from_f64(0.5) * v);
            for i in 0..3 {
                strain[self.outside[i]] = strain[self.outside[i]] + free[i];
                mean_strain[self.outside[i]] = mean_strain[self.outside[i]] + mean_free[i];
            }
            let q: Vec<Voigt> = (0..self.model.branches().len())
                .map(|j| core::array::from_fn(|i| change[3 + 6 * j + i].to_f64()))
                .collect();
            let mean_q: Vec<Voigt> = (0..q.len())
                .map(|j| core::array::from_fn(|i| mean_change[3 + 6 * j + i].to_f64()))
                .collect();
            let mut known_rate = [Dd::ZERO; 6];
            known_rate[self.inside[col]] =
                Dd::from_f64(ENGINEERING[col]) / Dd::from_f64(duration_s);
            let rate = self.ramp_endpoint_rate(&strain.map(Dd::to_f64), &q, known_rate)?;
            let stress = self.ramp_stress(&strain.map(Dd::to_f64), &q, &rate)?;
            let mean_stress = self.ramp_stress(
                &mean_strain.map(Dd::to_f64),
                &mean_q,
                &strain.map(|v| v / Dd::from_f64(duration_s)),
            )?;
            for i in 0..6 {
                answer.strain[i][col] = number(strain[i])?.to_f64() / WEIGHTS[i];
                answer.stress[i][col] = stress[i] / WEIGHTS[i];
                answer.mean_stress[i][col] = mean_stress[i] / WEIGHTS[i];
                answer.strain_rate[i][col] = number(rate[i])?.to_f64() / WEIGHTS[i];
            }
        }
        poll(&mut cancelled)?;
        Ok(answer)
    }
}
