//! Incremental motion about a held static deflection of the original SOS law.
//!
//! For the source potential `V`, a supplied displacement `a` requires the
//! constant generalized force `f0 = grad V(a)`. With `q` measured FROM that
//! equilibrium, the mechanical potential is `V(a+q)-V(a)-f0.q`. It contains
//! the original quartic terms AND the cubic terms responsible for asymmetric
//! motion and even harmonics. This is a coordinate translation under a dead
//! load, not a curvature model, a new shell theory, or an equilibrium solver.
//!
//! Each channel is expanded using `b=E a`, `s=a.E a`, `t=b.q`, `u=q.E q`:
//! `delta V = c (s*u/2 + t*t + t*u + u*u/4)`. Neither nearby source energies
//! nor nearby source forces are subtracted. All preparation is cold; warm
//! energy, effort, tangent and midpoint-remainder evaluations allocate nothing.

use fs_math::dd::Dd;
use fs_phs::{ContinuousStorageTangent, QuadraticStorage, Storage, linear_exact::QuadraticPart};

use crate::{NlModalError, SosModalStorage};

struct ReferenceChannel {
    product: Vec<f64>,
    form: f64,
    diagonal: bool,
}

/// The unchanged nonlinear modal law translated to a caller-supplied static
/// deflection with nonnegative tangent curvature under its held force.
///
/// State remains interleaved `[q_0, p_0, q_1, p_1, ...]`, with unit modal mass;
/// only displacement origin changes. Drive supplied to the assembled pHS is
/// ADDITIONAL to [`Self::preload_forces`]. Applying that force a second time
/// would count the preload twice. Damping and ports retain the source basis.
///
/// The complete local stiffness is admitted by the existing pHS PSD owner.
/// Negative tangent curvature refuses; no diagonalization, clipping or stiffness
/// floor changes the source. A semidefinite Hessian does NOT establish local
/// stability: cubic terms can destabilize a null direction. Neither stability,
/// global minimality nor validity of the source geometry is certified here.
pub struct PreloadedModalStorage {
    source: SosModalStorage,
    reference: Vec<f64>,
    preload: Vec<f64>,
    stiffness: QuadraticStorage,
    channels: Vec<ReferenceChannel>,
}

fn invalid(what: &'static str) -> NlModalError {
    NlModalError::Parameter { what }
}

// Twice the exact midpoint's displacement from the represented midpoint.
// Differences avoid overflowing a+b; error-free sums retain a small endpoint
// even when it lies far below the other endpoint's floating-point spacing.
fn midpoint_offset_twice(a: f64, b: f64, midpoint: f64) -> Dd {
    (Dd::from_f64(a) - Dd::from_f64(midpoint)) + (Dd::from_f64(b) - Dd::from_f64(midpoint))
}

impl PreloadedModalStorage {
    /// Translate the source to a declared equilibrium, retaining all modes.
    ///
    /// This does not solve for a deflection or infer a pressure projection.
    /// The caller must establish that the returned held generalized force is
    /// supplied by the intended physical supports/load. The source matrices
    /// must be exactly symmetric: their conservative scalar potential, effort
    /// and tangent then describe the same translated law without a repair.
    ///
    /// Assembly costs O(channels * modes^2), retains the full coupled local
    /// stiffness, and is followed by the existing pHS symmetric-PSD check. Warm
    /// evaluation costs O(modes) per diagonal stress channel and O(modes^2)
    /// per dense channel, with no state-dependent preparation.
    ///
    /// # Errors
    /// Invalid source/reference, nonfinite or unrepresentable stiffness/load,
    /// asymmetric stress matrices, or negative local tangent curvature.
    pub fn new(
        source: SosModalStorage,
        equilibrium_displacements: &[f64],
    ) -> Result<Self, NlModalError> {
        source.validate()?;
        let n = source.omegas.len();
        if equilibrium_displacements.len() != n
            || equilibrium_displacements.iter().any(|q| !q.is_finite())
        {
            return Err(invalid(
                "preload requires one finite displacement per source mode",
            ));
        }
        let mut stiffness = vec![0.0; n * n];
        let mut preload = vec![0.0; n];
        for i in 0..n {
            let squared = source.omegas[i] * source.omegas[i];
            stiffness[i * n + i] = squared;
            preload[i] = squared * equilibrium_displacements[i];
        }
        let mut channels = Vec::with_capacity(source.channels.len());
        for channel in &source.channels {
            let e = &channel.coupling;
            if (0..n).any(|i| (0..i).any(|j| e[i * n + j] != e[j * n + i])) {
                return Err(invalid(
                    "preload requires exactly symmetric conservative stress channels",
                ));
            }
            let diagonal = (0..n).all(|i| (0..n).all(|j| i == j || e[i * n + j] == 0.0));
            if channel.coefficient == 0.0 {
                channels.push(ReferenceChannel {
                    product: vec![0.0; n],
                    form: 0.0,
                    diagonal,
                });
                continue;
            }
            let product: Vec<f64> = (0..n)
                .map(|i| {
                    if diagonal {
                        e[i * n + i] * equilibrium_displacements[i]
                    } else {
                        (0..n)
                            .map(|j| e[i * n + j] * equilibrium_displacements[j])
                            .sum()
                    }
                })
                .collect();
            let form: f64 = equilibrium_displacements
                .iter()
                .zip(&product)
                .map(|(a, b)| a * b)
                .sum();
            if !form.is_finite() || product.iter().any(|v| !v.is_finite()) {
                return Err(invalid("preload stress contraction is not representable"));
            }
            for i in 0..n {
                preload[i] += channel.coefficient * form * product[i];
                for j in 0..=i {
                    let entry =
                        channel.coefficient * (form * e[i * n + j] + 2.0 * product[i] * product[j]);
                    stiffness[i * n + j] += entry;
                    if i != j {
                        stiffness[j * n + i] += entry;
                    }
                }
            }
            channels.push(ReferenceChannel {
                product,
                form,
                diagonal,
            });
        }
        if stiffness.iter().chain(&preload).any(|v| !v.is_finite()) {
            return Err(invalid(
                "preload force or local stiffness is not representable",
            ));
        }
        let stiffness = QuadraticStorage::new(stiffness, n).map_err(NlModalError::Phs)?;
        Ok(Self {
            source,
            reference: equilibrium_displacements.to_vec(),
            preload,
            stiffness,
            channels,
        })
    }

    /// Original, immutable modal law and unchanged source frequencies/basis.
    pub fn source(&self) -> &SosModalStorage {
        &self.source
    }

    /// Number of original mass-normalized displacement coordinates.
    pub fn mode_count(&self) -> usize {
        self.source.omegas.len()
    }

    /// Static displacement in the source's original modal coordinates.
    pub fn equilibrium_displacements(&self) -> &[f64] {
        &self.reference
    }

    /// Constant generalized force that must be held to sustain this reference.
    /// Incremental pHS drive excludes this already-accounted dead load.
    pub fn preload_forces(&self) -> &[f64] {
        &self.preload
    }

    /// Complete row-major local tangent `d^2 V / dq^2` at the reference.
    /// Signed off-diagonal entries are retained in the original basis.
    pub fn local_stiffness(&self) -> &[f64] {
        self.stiffness.matrix()
    }

    fn valid_state(&self, state: &[f64]) -> bool {
        state.len() == 2 * self.mode_count() && state.iter().all(|v| v.is_finite())
    }

    fn row(&self, channel: usize, row: usize, state: &[f64]) -> f64 {
        let n = self.mode_count();
        let e = &self.source.channels[channel].coupling;
        if self.channels[channel].diagonal {
            e[row * n + row] * state[2 * row]
        } else {
            (0..n).map(|j| e[row * n + j] * state[2 * j]).sum()
        }
    }

    fn forms(&self, channel: usize, state: &[f64]) -> (f64, f64) {
        let mut t = 0.0;
        let mut u = 0.0;
        for i in 0..self.mode_count() {
            t += self.channels[channel].product[i] * state[2 * i];
            u += state[2 * i] * self.row(channel, i, state);
        }
        (t, u)
    }

    /// Exact continuous effort derivative in the incremental state basis.
    /// This is the Hessian of the full translated potential, including every
    /// cubic/quartic interaction, not just [`Self::local_stiffness`].
    ///
    /// # Errors
    /// Incomplete/nonfinite state or direction, wrong output size, or arithmetic
    /// overflow. Shape/input refusal leaves output unchanged; arithmetic refusal
    /// leaves its contents unspecified.
    pub fn hessian_action(
        &self,
        state: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> Result<(), NlModalError> {
        if !self.valid_state(state) || !self.valid_state(direction) || out.len() != state.len() {
            return Err(invalid(
                "preloaded tangent requires complete finite state/direction/output",
            ));
        }
        let n = self.mode_count();
        for i in 0..n {
            out[2 * i] = self.source.omegas[i] * self.source.omegas[i] * direction[2 * i];
            out[2 * i + 1] = direction[2 * i + 1];
        }
        for (j, channel) in self.source.channels.iter().enumerate() {
            if channel.coefficient == 0.0 {
                continue;
            }
            let (t, u) = self.forms(j, state);
            let mut bv = 0.0;
            let mut qev = 0.0;
            for i in 0..n {
                bv += self.channels[j].product[i] * direction[2 * i];
                qev += direction[2 * i] * self.row(j, i, state);
            }
            let form = self.channels[j].form + 2.0 * t + u;
            for i in 0..n {
                out[2 * i] += channel.coefficient
                    * (form * self.row(j, i, direction)
                        + 2.0 * (bv + qev) * (self.channels[j].product[i] + self.row(j, i, state)));
            }
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(invalid("preloaded tangent arithmetic is not representable"));
        }
        Ok(())
    }
}

impl ContinuousStorageTangent for PreloadedModalStorage {
    fn apply(&self, state: &[f64], direction: &[f64], out: &mut [f64]) -> bool {
        self.hessian_action(state, direction, out).is_ok()
    }
}

impl Storage for PreloadedModalStorage {
    fn continuous_tangent(&self) -> Option<&dyn ContinuousStorageTangent> {
        Some(self)
    }

    fn hamiltonian(&self, state: &[f64]) -> f64 {
        if !self.valid_state(state) {
            return f64::NAN;
        }
        let mut energy = 0.0;
        for (i, omega) in self.source.omegas.iter().enumerate() {
            energy += f64::midpoint(
                state[2 * i + 1] * state[2 * i + 1],
                omega * omega * state[2 * i] * state[2 * i],
            );
        }
        for (j, channel) in self.source.channels.iter().enumerate() {
            if channel.coefficient == 0.0 {
                continue;
            }
            let (t, u) = self.forms(j, state);
            let c = channel.coefficient;
            energy +=
                0.5 * c * self.channels[j].form * u + c * t * t + c * t * u + 0.25 * c * u * u;
        }
        energy
    }

    fn gradient(&self, state: &[f64], out: &mut [f64]) {
        if !self.valid_state(state) || out.len() != state.len() {
            out.fill(f64::NAN);
            return;
        }
        let n = self.mode_count();
        for i in 0..n {
            out[2 * i] = self.source.omegas[i] * self.source.omegas[i] * state[2 * i];
            out[2 * i + 1] = state[2 * i + 1];
        }
        for (j, channel) in self.source.channels.iter().enumerate() {
            if channel.coefficient == 0.0 {
                continue;
            }
            let (t, u) = self.forms(j, state);
            let changed_form = 2.0 * t + u;
            for i in 0..n {
                out[2 * i] += channel.coefficient
                    * (changed_form * self.channels[j].product[i]
                        + (self.channels[j].form + changed_form) * self.row(j, i, state));
            }
        }
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        if !self.valid_state(a) || !self.valid_state(b) || !self.valid_state(mid) {
            return Some(f64::NAN);
        }
        let n = self.mode_count();
        let mut remainder = Dd::ZERO;
        let half = Dd::from_f64(0.5);
        let quarter = Dd::from_f64(0.25);
        for (i, omega) in self.source.omegas.iter().enumerate() {
            let dq = Dd::from_f64(b[2 * i]) - Dd::from_f64(a[2 * i]);
            let dp = Dd::from_f64(b[2 * i + 1]) - Dd::from_f64(a[2 * i + 1]);
            let oq = midpoint_offset_twice(a[2 * i], b[2 * i], mid[2 * i]);
            let op = midpoint_offset_twice(a[2 * i + 1], b[2 * i + 1], mid[2 * i + 1]);
            remainder = remainder + half * (Dd::from_f64(omega * omega) * dq * oq + dp * op);
        }
        for (k, channel) in self.source.channels.iter().enumerate() {
            if channel.coefficient == 0.0 {
                continue;
            }
            // d=b-a, o=(a-mid)+(b-mid), z=reference+mid mathematically.
            // D=d.E.d, O=o.E.o, C=d.E.o, M=z.E.d, N=z.E.o, S=z.E.z.
            // The complete remainder is c*((M+C/2)*(N+(D+O)/4)+S*C/2).
            // This includes both the ideal-midpoint cubic remainder and the
            // change of effort from the supplied midpoint to the exact mean.
            // Form z's contractions from the stored reference coefficients;
            // neither z nor a difference of source energies is constructed.
            let mut d_ed = Dd::ZERO;
            let mut o_eo = Dd::ZERO;
            let mut d_eo = Dd::ZERO;
            let mut z_ed = Dd::ZERO;
            let mut z_eo = Dd::ZERO;
            let mut z_ez = Dd::from_f64(self.channels[k].form);
            for i in 0..n {
                let di = Dd::from_f64(b[2 * i]) - Dd::from_f64(a[2 * i]);
                let oi = midpoint_offset_twice(a[2 * i], b[2 * i], mid[2 * i]);
                let mi = Dd::from_f64(mid[2 * i]);
                let (ed, eo, em) = if self.channels[k].diagonal {
                    let e = Dd::from_f64(channel.coupling[i * n + i]);
                    (e * di, e * oi, e * mi)
                } else {
                    let mut ed = Dd::ZERO;
                    let mut eo = Dd::ZERO;
                    let mut em = Dd::ZERO;
                    for j in 0..n {
                        let e = Dd::from_f64(channel.coupling[i * n + j]);
                        ed = ed + e * (Dd::from_f64(b[2 * j]) - Dd::from_f64(a[2 * j]));
                        eo = eo + e * midpoint_offset_twice(a[2 * j], b[2 * j], mid[2 * j]);
                        em = em + e * Dd::from_f64(mid[2 * j]);
                    }
                    (ed, eo, em)
                };
                let reference_product = Dd::from_f64(self.channels[k].product[i]);
                d_ed = d_ed + di * ed;
                o_eo = o_eo + oi * eo;
                d_eo = d_eo + di * eo;
                z_ed = z_ed + reference_product * di + mi * ed;
                z_eo = z_eo + reference_product * oi + mi * eo;
                z_ez = z_ez + Dd::from_f64(2.0) * reference_product * mi + mi * em;
            }
            remainder = remainder
                + Dd::from_f64(channel.coefficient)
                    * ((z_ed + half * d_eo) * (z_eo + quarter * (d_ed + o_eo))
                        + half * z_ez * d_eo);
        }
        Some(remainder.to_f64())
    }

    fn quadratic_part(&self) -> Option<QuadraticPart> {
        let n = self.mode_count();
        let dimension = 2 * n;
        let mut matrix = vec![0.0; dimension * dimension];
        for i in 0..n {
            matrix[(2 * i + 1) * dimension + 2 * i + 1] = 1.0;
            for j in 0..n {
                matrix[2 * i * dimension + 2 * j] = self.local_stiffness()[i * n + j];
            }
        }
        Some(QuadraticPart {
            matrix,
            exact: self.source.channels.iter().all(|channel| {
                channel.coefficient == 0.0 || channel.coupling.iter().all(|entry| *entry == 0.0)
            }),
        })
    }
}
