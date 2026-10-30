//! Small joint constraint algebra; the complete mechanical operator stays matrix-free.
use super::*;

pub(crate) struct Contacts {
    pub columns: Vec<Vec<f64>>,
    pub norms: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    norm2: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub norm1: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub phases: [LinearNetworkPhase; MAX_LINEAR_NETWORK_CONTACTS],
    pub static_capacity: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub kinetic_capacity: [f64; MAX_LINEAR_NETWORK_CONTACTS],
    pub velocity_coeff: [[f64; ORDER + 1]; MAX_LINEAR_NETWORK_CONTACTS],
    pub reaction_coeff: [[f64; ORDER + 1]; MAX_LINEAR_NETWORK_CONTACTS],
    gram: [f64; MAX_LINEAR_NETWORK_CONTACTS * MAX_LINEAR_NETWORK_CONTACTS],
    lower: [f64; MAX_LINEAR_NETWORK_CONTACTS * MAX_LINEAR_NETWORK_CONTACTS],
    active: [usize; MAX_LINEAR_NETWORK_CONTACTS],
    active_count: usize,
    minimum_pivot: f64,
    active_mask: u16,
    projection_bound: f64,
    roundoff_gain: f64,
}

impl Contacts {
    pub fn new(columns: Vec<Vec<f64>>, minimum_pivot: f64) -> Result<Self, Error> {
        let count = columns.len();
        if count == 0
            || count > MAX_LINEAR_NETWORK_CONTACTS
            || !minimum_pivot.is_finite()
            || minimum_pivot <= 0.0
            || minimum_pivot >= 1.0
        {
            return Err(Error::Invalid(
                "bounded independent junction columns and normalized pivot",
            ));
        }
        let n = columns[0].len();
        if n == 0
            || count > n
            || columns
                .iter()
                .any(|c| c.len() != n || c.iter().any(|x| !x.is_finite()))
        {
            return Err(Error::Invalid("complete independent junction work columns"));
        }
        let mut norms = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        let mut norm2 = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        let mut norm1 = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        for i in 0..count {
            norm2[i] = dot(&columns[i], &columns[i])?;
            norms[i] = finite(det::sqrt(norm2[i]))?;
            if norms[i] <= 0.0 {
                return Err(Error::Invalid("nonzero junction work column"));
            }
            norm1[i] = finite(columns[i].iter().map(|v| v.abs()).sum())?;
        }
        let mut gram = [0.0; MAX_LINEAR_NETWORK_CONTACTS * MAX_LINEAR_NETWORK_CONTACTS];
        for i in 0..count {
            for j in 0..=i {
                let value = finite(
                    columns[i]
                        .iter()
                        .zip(&columns[j])
                        .map(|(a, b)| (a / norms[i]) * (b / norms[j]))
                        .sum(),
                )?;
                gram[i * MAX_LINEAR_NETWORK_CONTACTS + j] = value;
                gram[j * MAX_LINEAR_NETWORK_CONTACTS + i] = value;
            }
        }
        let mut result = Self {
            columns,
            norms,
            norm2,
            norm1,
            gram,
            phases: [LinearNetworkPhase::Sticking; MAX_LINEAR_NETWORK_CONTACTS],
            static_capacity: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            kinetic_capacity: [0.0; MAX_LINEAR_NETWORK_CONTACTS],
            velocity_coeff: [[0.0; ORDER + 1]; MAX_LINEAR_NETWORK_CONTACTS],
            reaction_coeff: [[0.0; ORDER + 1]; MAX_LINEAR_NETWORK_CONTACTS],
            lower: [0.0; MAX_LINEAR_NETWORK_CONTACTS * MAX_LINEAR_NETWORK_CONTACTS],
            active: [0; MAX_LINEAR_NETWORK_CONTACTS],
            active_count: 0,
            minimum_pivot,
            active_mask: u16::MAX,
            projection_bound: 1.0,
            roundoff_gain: 1.0,
        };
        // Full rank makes every subsequently selected principal constraint block
        // positive definite. No regularization changes physical contact forces.
        result.factor()?;
        result.phases.fill(LinearNetworkPhase::Fresh);
        Ok(result)
    }

    pub fn count(&self) -> usize {
        self.columns.len()
    }

    pub fn factor(&mut self) -> Result<(), Error> {
        let mask = (0..self.count())
            .filter(|i| self.phases[*i] == LinearNetworkPhase::Sticking)
            .fold(0u16, |mask, i| mask | (1u16 << i));
        if self.active_mask == mask {
            return Ok(());
        }
        // A failed factor must not leave a previous cache key on partial data.
        self.active_mask = u16::MAX;
        self.active_count = 0;
        self.lower.fill(0.0);
        for i in 0..self.count() {
            if self.phases[i] == LinearNetworkPhase::Sticking {
                self.active[self.active_count] = i;
                self.active_count += 1;
            }
        }
        for i in 0..self.active_count {
            for j in 0..=i {
                let mut value =
                    self.gram[self.active[i] * MAX_LINEAR_NETWORK_CONTACTS + self.active[j]];
                for k in 0..j {
                    value -= self.lower[i * MAX_LINEAR_NETWORK_CONTACTS + k]
                        * self.lower[j * MAX_LINEAR_NETWORK_CONTACTS + k];
                }
                value = finite(value)?;
                self.lower[i * MAX_LINEAR_NETWORK_CONTACTS + j] = if i == j {
                    if value <= self.minimum_pivot {
                        return Err(Error::Invalid(
                            "junction columns are dependent or unresolved at the declared pivot",
                        ));
                    }
                    det::sqrt(value)
                } else {
                    value / self.lower[j * MAX_LINEAR_NETWORK_CONTACTS + j]
                };
            }
        }
        self.projection_bound = self.compute_projection_bound();
        self.roundoff_gain = 1.0;
        if self.active_count > 1 {
            // Roundoff in the normalized Gram solve can be amplified even though
            // the exact mass projector has norm one. Carry this conditioning in
            // trajectory/work uncertainty rather than hiding it in a pivot test.
            let mut rows = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
            for j in 0..self.active_count {
                let mut rhs = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
                rhs[j] = 1.0;
                self.solve(&mut rhs)?;
                for i in 0..self.active_count {
                    rows[i] = finite(rows[i] + rhs[i].abs())?;
                }
            }
            self.roundoff_gain = finite(
                1.0 + self.active_count as f64
                    * rows[..self.active_count]
                        .iter()
                        .copied()
                        .fold(0.0, f64::max),
            )?;
        }
        self.active_mask = mask;
        Ok(())
    }

    fn solve(&self, rhs: &mut [f64; MAX_LINEAR_NETWORK_CONTACTS]) -> Result<(), Error> {
        for i in 0..self.active_count {
            let mut value = rhs[i];
            for j in 0..i {
                value -= self.lower[i * MAX_LINEAR_NETWORK_CONTACTS + j] * rhs[j];
            }
            rhs[i] = finite(value / self.lower[i * MAX_LINEAR_NETWORK_CONTACTS + i])?;
        }
        for i in (0..self.active_count).rev() {
            let mut value = rhs[i];
            for j in i + 1..self.active_count {
                value -= self.lower[j * MAX_LINEAR_NETWORK_CONTACTS + i] * rhs[j];
            }
            rhs[i] = finite(value / self.lower[i * MAX_LINEAR_NETWORK_CONTACTS + i])?;
        }
        Ok(())
    }

    pub const fn projection_norm(&self) -> f64 {
        self.projection_bound
    }
    pub const fn roundoff_gain(&self) -> f64 {
        self.roundoff_gain
    }

    fn compute_projection_bound(&self) -> f64 {
        if self.active_count == 0 {
            1.0
        } else if self.active_count == 1 {
            let i = self.active[0];
            let norm2 = self.norm2[i];
            self.columns[i]
                .iter()
                .map(|b| (1.0 - b * b / norm2).abs() + b.abs() * (self.norm1[i] - b.abs()) / norm2)
                .fold(0.0, f64::max)
        } else {
            // An orthogonal projector has spectral norm <=1, hence its infinity
            // norm is <=sqrt(n). This avoids forming a dense n-by-n operator.
            det::sqrt(self.columns[0].len() as f64)
        }
    }

    /// Apply every held sliding effort, then solve ALL static reactions jointly.
    pub fn reactions(
        &self,
        acceleration: &mut [f64],
        constant: bool,
    ) -> Result<[f64; MAX_LINEAR_NETWORK_CONTACTS], Error> {
        let mut result = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        if constant {
            for i in 0..self.count() {
                if let LinearNetworkPhase::Sliding { direction } = self.phases[i] {
                    result[i] = finite(direction as f64 * self.kinetic_capacity[i])?;
                    for (a, b) in acceleration.iter_mut().zip(&self.columns[i]) {
                        *a = finite(*a - b * result[i])?;
                    }
                }
            }
        }
        if self.active_count == 1 {
            // Retain the exact scalar projection, including exactly stationary
            // states, without a needless square-root / Cholesky round trip.
            let i = self.active[0];
            result[i] = finite(dot(&self.columns[i], acceleration)? / self.norm2[i])?;
            for (a, b) in acceleration.iter_mut().zip(&self.columns[i]) {
                *a = finite(*a - b * result[i])?;
            }
            return Ok(result);
        }
        let mut rhs = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        for (j, &i) in self.active[..self.active_count].iter().enumerate() {
            rhs[j] = finite(dot(&self.columns[i], acceleration)? / self.norms[i])?;
        }
        self.solve(&mut rhs)?;
        for (j, &i) in self.active[..self.active_count].iter().enumerate() {
            result[i] = finite(rhs[j] / self.norms[i])?;
        }
        for (j, &i) in self.active[..self.active_count].iter().enumerate() {
            for (a, b) in acceleration.iter_mut().zip(&self.columns[i]) {
                *a = finite(*a - (b / self.norms[i]) * rhs[j])?;
            }
        }
        if self.active_count == acceleration.len() {
            // Full-rank constraints spanning the complete mass-normalized space
            // leave no acceleration; do not manufacture roundoff vibration.
            acceleration.fill(0.0);
        }
        Ok(result)
    }

    pub fn coefficient(
        &mut self,
        velocity: &[f64],
        acceleration: &mut [f64],
        order: usize,
    ) -> Result<(), Error> {
        let reaction = self.reactions(acceleration, order == 0)?;
        for i in 0..self.count() {
            self.velocity_coeff[i][order] = dot(&self.columns[i], velocity)?;
            self.reaction_coeff[i][order] = reaction[i];
        }
        Ok(())
    }

    /// Remove only the numerical violation of the active constraints. The
    /// caller accounts for the complete correction in its work uncertainty.
    pub fn project(
        &self,
        value: &mut [f64],
        reference: Option<&[f64]>,
    ) -> Result<(f64, f64), Error> {
        if self.active_count == 1 {
            let i = self.active[0];
            let drift = finite(
                self.columns[i]
                    .iter()
                    .enumerate()
                    .map(|(k, b)| b * (value[k] - reference.map_or(0.0, |r| r[k])))
                    .sum(),
            )?;
            let scale = finite(drift / self.norm2[i])?;
            let mut maximum: f64 = 0.0;
            for (v, b) in value.iter_mut().zip(&self.columns[i]) {
                let correction = finite(b * scale)?;
                *v = finite(*v - correction)?;
                maximum = maximum.max(correction.abs());
            }
            return Ok((maximum, finite(scale.abs() * self.norms[i])?));
        }
        let mut rhs = [0.0; MAX_LINEAR_NETWORK_CONTACTS];
        for (j, &i) in self.active[..self.active_count].iter().enumerate() {
            rhs[j] = finite(
                self.columns[i]
                    .iter()
                    .enumerate()
                    .map(|(k, b)| {
                        (b / self.norms[i]) * (value[k] - reference.map_or(0.0, |r| r[k]))
                    })
                    .sum(),
            )?;
        }
        self.solve(&mut rhs)?;
        let mut maximum: f64 = 0.0;
        let mut squared: f64 = 0.0;
        for k in 0..value.len() {
            let correction = if self.active_count == value.len() {
                finite(value[k] - reference.map_or(0.0, |r| r[k]))?
            } else {
                finite(
                    self.active[..self.active_count]
                        .iter()
                        .enumerate()
                        .map(|(j, &i)| self.columns[i][k] / self.norms[i] * rhs[j])
                        .sum(),
                )?
            };
            value[k] = finite(value[k] - correction)?;
            maximum = maximum.max(correction.abs());
            squared = finite(squared + correction * correction)?;
        }
        Ok((maximum, det::sqrt(squared)))
    }
}
