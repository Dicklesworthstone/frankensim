//! Explicit magnetic storage outside the resolved winding field.
//!
//! An independent self-inductance and signed shared-flux modes describe
//! L = diag(L_self) + sum(ell * n * n^T), with nonnegative coefficients.
//! The caller supplies physical end-region data; these terms must not duplicate
//! flux already resolved by the field. No leakage, inverse or regularization is
//! inferred. Dependent and zero modes are valid passive storage.
//!
//! Flux, energy and finite-step work use the original modes, never a quadratic
//! form expanded into cancelling matrix entries. The matrix is cached only for
//! inspection. At most 16 windings and 16 shared modes bound every operation.

/// The small terminal block is independent of field mesh size.
pub const MAX_END_WINDINGS: usize = 16;
/// Maximum independently authored shared end-region flux modes.
pub const MAX_END_MODES: usize = 16;

/// One shared magnetic path outside the resolved field domain.
#[derive(Clone, Debug, PartialEq)]
pub struct EndWindingMode {
    /// Nonnegative reference inductance [H]; zero explicitly neglects this mode.
    pub inductance_h: f64,
    /// Signed dimensionless winding participation in original terminal order.
    /// Reversing one winding reverses its entry, including its mutual polarity.
    pub turns: Vec<f64>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EndWindingError {
    Invalid(&'static str),
    Arithmetic,
}
impl core::fmt::Display for EndWindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "end-winding storage: {self:?}")
    }
}
impl std::error::Error for EndWindingError {}
type Result<T> = core::result::Result<T, EndWindingError>;

/// Fixed passive end storage, shared by field extraction and circuit owners.
#[derive(Clone, Debug, PartialEq)]
pub struct EndWindingInductance {
    diagonal_h: Vec<f64>,
    modes: Vec<EndWindingMode>,
    matrix_h: Vec<f64>,
}

/// Quadratic-storage quantities for a finite current increment.
/// Endpoint energies are observations; stored_change_j is evaluated directly
/// from modal midpoint currents and increments, not by subtracting them.
#[derive(Clone, Debug)]
pub struct EndWindingStep {
    pub flux_increment_wb: Vec<f64>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    /// Right-end current paired with the finite flux increment [J].
    pub endpoint_work_j: f64,
    /// Backward-Euler loss: half of delta_I^T L_end delta_I.
    /// This is numerical dissipation, never physical conductor heat.
    pub numerical_dissipation_j: f64,
}

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EndWindingError::Arithmetic)
    }
}

// Multiplying the smallest and largest magnitudes first avoids a gratuitous
// intermediate overflow/underflow when the complete three-factor product fits.
fn product(a: f64, b: f64, c: f64) -> Result<f64> {
    if a == 0.0 || b == 0.0 || c == 0.0 {
        return Ok(0.0);
    }
    let mut values = [a, b, c];
    values.sort_by(|x, y| x.abs().total_cmp(&y.abs()));
    finite(finite(values[0] * values[2])? * values[1])
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<()> {
        let next = finite(self.value + term)?;
        let lost = if self.value.abs() >= term.abs() {
            (self.value - next) + term
        } else {
            (term - next) + self.value
        };
        self.correction = finite(self.correction + lost)?;
        self.value = next;
        Ok(())
    }
    fn finish(self) -> Result<f64> {
        finite(self.value + self.correction)
    }
}
fn projection(turns: &[f64], currents: &[f64]) -> Result<f64> {
    let mut sum = Sum::default();
    for (&n, &i) in turns.iter().zip(currents) {
        if n == 0.0 {
            continue;
        }
        let term = finite(n * i)?;
        sum.add(term)?;
        sum.add(finite(n.mul_add(i, -term))?)?;
    }
    sum.finish()
}

#[derive(Default)]
struct Work {
    before: Sum,
    after: Sum,
    change: Sum,
    endpoint: Sum,
    numerical: Sum,
}
impl Work {
    fn add(&mut self, l: f64, before: f64, after: f64, increment: f64) -> Result<()> {
        self.before.add(product(l, 0.5 * before, before)?)?;
        self.after.add(product(l, 0.5 * after, after)?)?;
        self.change
            .add(product(l, f64::midpoint(before, after), increment)?)?;
        self.endpoint.add(product(l, after, increment)?)?;
        self.numerical
            .add(product(l, 0.5 * increment, increment)?)?;
        Ok(())
    }
}

impl EndWindingInductance {
    /// Admit authored self-inductances and shared flux modes once.
    /// Positive-semidefinite storage follows from the factors themselves;
    /// singular modes need neither inversion nor invented series leakage.
    pub fn new(diagonal_h: Vec<f64>, modes: Vec<EndWindingMode>) -> Result<Self> {
        let count = diagonal_h.len();
        if count == 0 || count > MAX_END_WINDINGS || modes.len() > MAX_END_MODES {
            return Err(EndWindingError::Invalid(
                "bounded nonempty winding block required",
            ));
        }
        if diagonal_h.iter().any(|l| !l.is_finite() || *l < 0.0)
            || modes.iter().any(|mode| {
                !mode.inductance_h.is_finite()
                    || mode.inductance_h < 0.0
                    || mode.turns.len() != count
                    || mode.turns.iter().any(|n| !n.is_finite())
            })
        {
            return Err(EndWindingError::Invalid(
                "finite nonnegative coefficients and complete signed modes required",
            ));
        }
        let mut matrix_h = vec![0.0; count * count];
        for row in 0..count {
            for column in row..count {
                let mut sum = Sum::default();
                if row == column {
                    sum.add(diagonal_h[row])?;
                }
                for mode in &modes {
                    sum.add(product(
                        mode.inductance_h,
                        mode.turns[row],
                        mode.turns[column],
                    )?)?;
                }
                let value = sum.finish()?;
                matrix_h[row * count + column] = value;
                matrix_h[column * count + row] = value;
            }
        }
        Ok(Self {
            diagonal_h,
            modes,
            matrix_h,
        })
    }

    pub fn winding_count(&self) -> usize {
        self.diagonal_h.len()
    }
    /// Independent diagonal contribution, before shared-mode self terms.
    pub fn diagonal_h(&self) -> &[f64] {
        &self.diagonal_h
    }
    pub fn modes(&self) -> &[EndWindingMode] {
        &self.modes
    }
    /// Full symmetric row-major matrix [H], including shared-mode self terms.
    /// Evolution uses the original factors to preserve cancelling flux modes.
    pub fn matrix_h(&self) -> &[f64] {
        &self.matrix_h
    }

    fn admit_currents(&self, currents: &[f64]) -> Result<()> {
        if currents.len() != self.winding_count() || currents.iter().any(|i| !i.is_finite()) {
            return Err(EndWindingError::Invalid(
                "complete finite terminal currents required",
            ));
        }
        Ok(())
    }

    pub fn flux_linkages_wb(&self, currents: &[f64]) -> Result<Vec<f64>> {
        self.admit_currents(currents)?;
        let mut flux = [Sum::default(); MAX_END_WINDINGS];
        for (w, (&l, &i)) in self.diagonal_h.iter().zip(currents).enumerate() {
            flux[w].add(product(l, i, 1.0)?)?;
        }
        for mode in &self.modes {
            if mode.inductance_h == 0.0 {
                continue;
            }
            let current = projection(&mode.turns, currents)?;
            for (w, &n) in mode.turns.iter().enumerate() {
                flux[w].add(product(mode.inductance_h, n, current)?)?;
            }
        }
        flux[..self.winding_count()]
            .iter()
            .map(|s| s.finish())
            .collect()
    }

    pub fn stored_energy_j(&self, currents: &[f64]) -> Result<f64> {
        self.admit_currents(currents)?;
        let mut energy = Sum::default();
        for (&l, &i) in self.diagonal_h.iter().zip(currents) {
            energy.add(product(l, 0.5 * i, i)?)?;
        }
        for mode in &self.modes {
            if mode.inductance_h == 0.0 {
                continue;
            }
            let current = projection(&mode.turns, currents)?;
            energy.add(product(mode.inductance_h, 0.5 * current, current)?)?;
        }
        energy.finish()
    }

    /// Evaluate flux and work from the SAME finite current increment.
    /// Modal increments are formed before energy evaluation, preserving small
    /// changes on a large magnetic bias. Perfectly cancelling shared modes
    /// remain zero even when an expanded quadratic form would overflow.
    pub fn step(&self, before: &[f64], after: &[f64]) -> Result<EndWindingStep> {
        self.admit_currents(before)?;
        self.admit_currents(after)?;
        let count = self.winding_count();
        let mut increments = [0.0; MAX_END_WINDINGS];
        for w in 0..count {
            increments[w] = after[w] - before[w];
        }
        let mut flux = [Sum::default(); MAX_END_WINDINGS];
        let mut work = Work::default();
        for (w, &l) in self.diagonal_h.iter().enumerate() {
            if l == 0.0 {
                continue;
            }
            let increment = finite(increments[w])?;
            flux[w].add(product(l, increment, 1.0)?)?;
            work.add(l, before[w], after[w], increment)?;
        }
        for mode in &self.modes {
            if mode.inductance_h == 0.0 {
                continue;
            }
            let old = projection(&mode.turns, before)?;
            let new = projection(&mode.turns, after)?;
            let finite_support = mode
                .turns
                .iter()
                .zip(&increments)
                .all(|(&n, di)| n == 0.0 || di.is_finite());
            let increment = if finite_support {
                projection(&mode.turns, &increments[..count])?
            } else {
                // An individual current difference may overflow although its
                // shared mode remains small. Cancel half-increments first,
                // rather than subtracting rounded, strongly biased projections.
                let mut half = [0.0; MAX_END_WINDINGS];
                for w in 0..count {
                    if mode.turns[w] != 0.0 {
                        half[w] = f64::midpoint(-before[w], after[w]);
                        if half[w] == 0.0 && increments[w] != 0.0 {
                            return Err(EndWindingError::Arithmetic);
                        }
                    }
                }
                finite(2.0 * projection(&mode.turns, &half[..count])?)?
            };
            for (w, &n) in mode.turns.iter().enumerate() {
                flux[w].add(product(mode.inductance_h, n, increment)?)?;
            }
            work.add(mode.inductance_h, old, new, increment)?;
        }
        Ok(EndWindingStep {
            flux_increment_wb: flux[..count]
                .iter()
                .map(|s| s.finish())
                .collect::<Result<_>>()?,
            stored_before_j: work.before.finish()?,
            stored_after_j: work.after.finish()?,
            stored_change_j: work.change.finish()?,
            endpoint_work_j: work.endpoint.finish()?,
            numerical_dissipation_j: work.numerical.finish()?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 2e-14 * a.abs().max(b.abs()).max(1.0),
            "{a} != {b}"
        );
    }

    #[test]
    fn signed_shared_modes_match_independent_flux_and_storage() {
        let model = EndWindingInductance::new(
            vec![0.02, 0.03],
            vec![
                EndWindingMode {
                    inductance_h: 0.2,
                    turns: vec![1.0, -2.0],
                },
                EndWindingMode {
                    inductance_h: 0.1,
                    turns: vec![3.0, 1.0],
                },
            ],
        )
        .unwrap();
        for (&actual, expected) in model.matrix_h().iter().zip([1.12, -0.1, -0.1, 0.93]) {
            near(actual, expected);
        }
        let current = [2.0, -1.0];
        let flux = model.flux_linkages_wb(&current).unwrap();
        near(flux[0], 2.34);
        near(flux[1], -1.13);
        near(model.stored_energy_j(&current).unwrap(), 2.905);
        let step = model.step(&[0.3, 0.7], &current).unwrap();
        let di = [1.7, -1.7];
        near(step.flux_increment_wb[0], 1.12 * di[0] - 0.1 * di[1]);
        near(step.flux_increment_wb[1], -0.1 * di[0] + 0.93 * di[1]);
        near(
            step.stored_after_j - step.stored_before_j,
            step.stored_change_j,
        );
        near(
            step.stored_change_j + step.numerical_dissipation_j,
            step.endpoint_work_j,
        );
        near(
            step.endpoint_work_j,
            current[0] * step.flux_increment_wb[0] + current[1] * step.flux_increment_wb[1],
        );
    }

    #[test]
    fn null_modes_survive_large_opposing_currents_without_expanded_energy() {
        let model = EndWindingInductance::new(
            vec![0.0, 0.0],
            vec![EndWindingMode {
                inductance_h: 1.0,
                turns: vec![1.0, -1.0],
            }],
        )
        .unwrap();
        let currents = [f64::MAX, f64::MAX];
        assert_eq!(model.stored_energy_j(&currents).unwrap(), 0.0);
        assert_eq!(model.flux_linkages_wb(&currents).unwrap(), [0.0, 0.0]);
        let step = model.step(&currents, &[-f64::MAX, -f64::MAX]).unwrap();
        assert_eq!(step.flux_increment_wb, [0.0, 0.0]);
        assert_eq!(step.stored_change_j, 0.0);
        assert_eq!(step.endpoint_work_j, 0.0);
        assert_eq!(step.numerical_dissipation_j, 0.0);
    }

    #[test]
    fn modal_projection_retains_small_component_between_cancelling_currents() {
        let model = EndWindingInductance::new(
            vec![0.0; 3],
            vec![EndWindingMode {
                inductance_h: 1.0,
                turns: vec![1.0, 1.0, -1.0],
            }],
        )
        .unwrap();
        let currents = [1e16, 1.0, 1e16];
        assert_eq!(model.stored_energy_j(&currents).unwrap(), 0.5);
        assert_eq!(model.flux_linkages_wb(&currents).unwrap(), [1.0, 1.0, -1.0]);
        let biased = EndWindingInductance::new(vec![1.0], vec![]).unwrap();
        let step = biased.step(&[1e12], &[1e12 + 1.0]).unwrap();
        assert_eq!(step.stored_change_j, 1e12 + 0.5);
        assert_eq!(step.endpoint_work_j, 1e12 + 1.0);
        assert_eq!(step.numerical_dissipation_j, 0.5);
        assert_ne!(
            step.stored_after_j - step.stored_before_j,
            step.stored_change_j
        );
    }

    #[test]
    fn overflowed_current_reversal_preserves_small_increment_on_large_modal_bias() {
        let model = EndWindingInductance::new(
            vec![0.0; 4],
            vec![EndWindingMode {
                inductance_h: 1.0,
                turns: vec![1.0, 1.0, -1.0, 1.0],
            }],
        )
        .unwrap();
        let step = model
            .step(
                &[f64::MAX, 1e16, f64::MAX, 1.0],
                &[-f64::MAX, 1e16, -f64::MAX, 2.0],
            )
            .unwrap();
        assert_eq!(step.flux_increment_wb, [1.0, 1.0, -1.0, 1.0]);
        assert_eq!(step.numerical_dissipation_j, 0.5);
        near(step.stored_change_j, 1e16 + 1.5);
        near(step.endpoint_work_j, 1e16 + 2.0);
    }

    #[test]
    fn winding_polarity_transforms_flux_but_preserves_energy_and_loss() {
        let make = |sign| {
            EndWindingInductance::new(
                vec![0.02, 0.03],
                vec![EndWindingMode {
                    inductance_h: 0.4,
                    turns: vec![1.0, -2.0 * sign],
                }],
            )
            .unwrap()
        };
        let (a, b) = (make(1.0), make(-1.0));
        let (old, new) = ([0.2, -0.5], [0.7, 0.3]);
        let sa = a.step(&old, &new).unwrap();
        let sb = b.step(&[old[0], -old[1]], &[new[0], -new[1]]).unwrap();
        near(sa.flux_increment_wb[0], sb.flux_increment_wb[0]);
        near(sa.flux_increment_wb[1], -sb.flux_increment_wb[1]);
        near(sa.stored_change_j, sb.stored_change_j);
        near(sa.numerical_dissipation_j, sb.numerical_dissipation_j);
        near(sa.endpoint_work_j, sb.endpoint_work_j);
        near(
            a.stored_energy_j(&new).unwrap(),
            b.stored_energy_j(&[new[0], -new[1]]).unwrap(),
        );
    }

    #[test]
    fn admission_refuses_invalid_or_unrepresentable_models_and_currents() {
        assert!(EndWindingInductance::new(vec![], vec![]).is_err());
        assert!(EndWindingInductance::new(vec![0.0; MAX_END_WINDINGS + 1], vec![]).is_err());
        let mode = EndWindingMode {
            inductance_h: 1.0,
            turns: vec![1.0],
        };
        assert!(
            EndWindingInductance::new(vec![0.0], vec![mode.clone(); MAX_END_MODES + 1]).is_err()
        );
        assert!(EndWindingInductance::new(vec![-1.0], vec![]).is_err());
        assert!(EndWindingInductance::new(vec![f64::NAN], vec![]).is_err());
        for invalid in [
            EndWindingMode {
                inductance_h: -1.0,
                turns: vec![1.0],
            },
            EndWindingMode {
                inductance_h: f64::INFINITY,
                turns: vec![1.0],
            },
            EndWindingMode {
                inductance_h: 1.0,
                turns: vec![1.0, 2.0],
            },
            EndWindingMode {
                inductance_h: 1.0,
                turns: vec![f64::NAN],
            },
            EndWindingMode {
                inductance_h: f64::MAX,
                turns: vec![2.0],
            },
        ] {
            assert!(EndWindingInductance::new(vec![0.0], vec![invalid]).is_err());
        }
        let model = EndWindingInductance::new(
            vec![0.0],
            vec![EndWindingMode {
                inductance_h: 0.0,
                turns: vec![f64::MAX],
            }],
        )
        .unwrap();
        assert_eq!(model.stored_energy_j(&[f64::MAX]).unwrap(), 0.0);
        assert!(model.flux_linkages_wb(&[]).is_err());
        assert!(model.step(&[0.0], &[f64::NAN]).is_err());
    }
}
