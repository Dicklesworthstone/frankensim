//! Linear reversible immobile storage on the original aqueous transport owner.
//!
//! Per bulk transverse area, b_eq=B*c and db/dt=k*(B*c-b). Backward Euler
//! eliminates b locally: alpha=k*dt/(1+k*dt), b1=(1-alpha)*b0+alpha*B*c1.
//! The concentration system gains alpha*B on its storage diagonal and alpha*b0
//! on its RHS. Its water storage, Darcy advection and theta*D/tau diffusion
//! coefficients remain the original physical quantities.

use super::transport::{
    self, DiluteSoluteControls, DiluteSoluteError, DiluteSoluteReport, DiluteSoluteStep,
    ImmobileStorage, Result, finite, integrated_flux, nonnegative, poll, sum,
};

/// Explicit independent dilute linear sorption law for one cell and species.
///
/// Capacity B [m] is dry solid mass per transverse area [kg/m²] times the
/// linear distribution coefficient Kd [m³/kg]. It stays fixed as water changes.
/// The caller qualifies that linear isotherm and, when used, the relaxation
/// rate at the actual material/solute state. These are not water volumes or
/// maximum site loadings. No ion competition, nonlinear saturation, chemical
/// reaction, heat of sorption or hydraulic feedback is supplied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LinearSoluteSorption {
    /// Instantaneous equilibrium. Original b0 must already equal B*c0 within
    /// the explicit absolute equation tolerance; no initial projection occurs.
    Equilibrium { capacity_m: f64 },
    /// Finite linear-driving-force exchange. A zero rate leaves the original
    /// bound pool isolated, including when its initial loading is nonzero.
    Kinetic {
        capacity_m: f64,
        rate_s_inverse: f64,
    },
}

/// Original transport accounts plus separately retained aqueous/bound storage.
#[derive(Clone, Debug, PartialEq)]
pub struct SorbingDiluteSoluteReport<const N: usize> {
    /// Amount totals/changes/residuals here count aqueous PLUS sorbed solute.
    /// Face integrals still contain aqueous advection and diffusion only.
    pub total: DiluteSoluteReport<N>,
    pub aqueous_amount_before_mol_m2: [f64; N],
    pub aqueous_amount_after_mol_m2: [f64; N],
    pub aqueous_amount_change_mol_m2: [f64; N],
    pub sorbed_amount_before_mol_m2: [f64; N],
    pub sorbed_amount_after_mol_m2: [f64; N],
    pub sorbed_amount_change_mol_m2: [f64; N],
    /// Defect of the represented bound change against the locally eliminated
    /// backward-Euler exchange, or equilibrium storage equation [mol/m²].
    pub maximum_sorption_equation_defect_mol_m2: [f64; N],
}

#[derive(Clone, Debug, PartialEq)]
pub struct SorbingDiluteSoluteStep<const N: usize> {
    pub aqueous_amount_mol_m2: Vec<[f64; N]>,
    pub sorbed_amount_mol_m2: Vec<[f64; N]>,
    pub concentration_mol_m3: Vec<[f64; N]>,
    /// Actual represented endpoint loading change, positive into the bound pool.
    /// This is a net interval exchange, not a resolved adsorption/desorption history.
    pub exchange_into_sorbed_mol_m2: Vec<[f64; N]>,
    pub report: SorbingDiluteSoluteReport<N>,
}

impl LinearSoluteSorption {
    fn capacity(self) -> f64 {
        match self {
            Self::Equilibrium { capacity_m } | Self::Kinetic { capacity_m, .. } => capacity_m,
        }
    }

    fn validate(self) -> Result<()> {
        if !nonnegative(self.capacity())
            || matches!(self, Self::Kinetic { rate_s_inverse, .. } if !nonnegative(rate_s_inverse))
        {
            return Err(DiluteSoluteError::Invalid(
                "finite nonnegative sorption capacities and rates required",
            ));
        }
        Ok(())
    }

    fn weight(self, dt: f64) -> f64 {
        match self {
            Self::Equilibrium { .. } => 1.0,
            Self::Kinetic {
                rate_s_inverse: 0.0,
                ..
            } => 0.0,
            Self::Kinetic {
                rate_s_inverse: k, ..
            } => {
                let x = k * dt;
                if x <= 1.0 {
                    x / (1.0 + x)
                } else {
                    1.0 / (1.0 + (1.0 / k) / dt)
                }
            }
        }
    }

    fn retained(self, old: f64, dt: f64) -> Result<f64> {
        match self {
            Self::Equilibrium { .. } => Ok(0.0),
            Self::Kinetic {
                rate_s_inverse: k, ..
            } => {
                let x = k * dt;
                // If k*dt overflows then k>1 and dt>1; sequential division
                // cannot lose a representable final retained amount prematurely.
                let value = if x.is_finite() {
                    old / (1.0 + x)
                } else {
                    (old / k) / dt / (1.0 + (1.0 / k) / dt)
                };
                if old > 0.0 && value == 0.0 {
                    return Err(DiluteSoluteError::Invalid(
                        "positive retained sorbed amount is unrepresentable",
                    ));
                }
                finite(value)
            }
        }
    }
}

pub(super) struct SorptionStorage<'a, const N: usize> {
    original: &'a [[f64; N]],
    laws: &'a [[LinearSoluteSorption; N]],
    duration_s: f64,
    bound: Vec<[f64; N]>,
    exchange: Vec<[f64; N]>,
    maximum_defect: [f64; N],
}

impl<'a, const N: usize> SorptionStorage<'a, N> {
    pub(super) fn new(
        n: usize,
        duration_s: f64,
        original: &'a [[f64; N]],
        laws: &'a [[LinearSoluteSorption; N]],
        controls: DiluteSoluteControls<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let (_, base) = transport::admit_extent(n, controls)?;
        let bytes = base + Self::extra_bytes(n);
        if bytes > controls.maximum_storage_bytes {
            return Err(DiluteSoluteError::Workspace {
                required_bytes: bytes,
                maximum_bytes: controls.maximum_storage_bytes,
            });
        }
        if original.len() != n || laws.len() != n || !transport::positive(duration_s) {
            return Err(DiluteSoluteError::Invalid(
                "matching sorption extents and positive duration required",
            ));
        }
        for i in 0..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            for k in 0..N {
                laws[i][k].validate()?;
                if !nonnegative(original[i][k]) {
                    return Err(DiluteSoluteError::Invalid(
                        "finite nonnegative original sorbed amounts required",
                    ));
                }
                if let LinearSoluteSorption::Kinetic { rate_s_inverse, .. } = laws[i][k]
                    && rate_s_inverse > 0.0
                    && laws[i][k].weight(duration_s) == 0.0
                {
                    return Err(DiluteSoluteError::Invalid(
                        "positive sorption relaxation is unrepresentable",
                    ));
                }
            }
        }
        Ok(Self {
            original,
            laws,
            duration_s,
            bound: vec![[0.0; N]; n],
            exchange: vec![[0.0; N]; n],
            maximum_defect: [0.0; N],
        })
    }

    fn extra_bytes(n: usize) -> usize {
        // n,N already admitted to 4096,32. Two output fields, aggregation
        // arrays, retained defects and the simultaneously assembled report.
        (2 * n * N + 16 * N) * core::mem::size_of::<f64>() + 512
    }

    pub(super) fn finish(
        self,
        initial: &[[f64; N]],
        step: DiluteSoluteStep<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SorbingDiluteSoluteStep<N>> {
        let mut aqueous_before = [0.0; N];
        let mut aqueous_after = [0.0; N];
        let mut aqueous_change = [0.0; N];
        let mut bound_before = [0.0; N];
        let mut bound_after = [0.0; N];
        let mut bound_change = [0.0; N];
        for k in 0..N {
            poll(cancelled)?;
            aqueous_before[k] = sum(initial.iter().map(|a| a[k]))?;
            aqueous_after[k] = sum(step.amount_mol_m2.iter().map(|a| a[k]))?;
            aqueous_change[k] = sum(initial
                .iter()
                .zip(&step.amount_mol_m2)
                .map(|(a, b)| b[k] - a[k]))?;
            bound_before[k] = sum(self.original.iter().map(|a| a[k]))?;
            bound_after[k] = sum(self.bound.iter().map(|a| a[k]))?;
            bound_change[k] = sum(self.exchange.iter().map(|a| a[k]))?;
        }
        poll(cancelled)?;
        Ok(SorbingDiluteSoluteStep {
            aqueous_amount_mol_m2: step.amount_mol_m2,
            sorbed_amount_mol_m2: self.bound,
            concentration_mol_m3: step.concentration_mol_m3,
            exchange_into_sorbed_mol_m2: self.exchange,
            report: SorbingDiluteSoluteReport {
                total: step.report,
                aqueous_amount_before_mol_m2: aqueous_before,
                aqueous_amount_after_mol_m2: aqueous_after,
                aqueous_amount_change_mol_m2: aqueous_change,
                sorbed_amount_before_mol_m2: bound_before,
                sorbed_amount_after_mol_m2: bound_after,
                sorbed_amount_change_mol_m2: bound_change,
                maximum_sorption_equation_defect_mol_m2: self.maximum_defect,
            },
        })
    }
}

impl<const N: usize> ImmobileStorage<N> for SorptionStorage<'_, N> {
    fn additional_bytes(&self) -> usize {
        Self::extra_bytes(self.original.len())
    }

    fn admit(
        &self,
        mobile: &[[f64; N]],
        water: &[f64],
        controls: DiluteSoluteControls<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        for i in 0..mobile.len() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            for k in 0..N {
                if let LinearSoluteSorption::Equilibrium { capacity_m } = self.laws[i][k] {
                    let expected = integrated_flux(capacity_m, finite(mobile[i][k] / water[i])?)?;
                    let residual = self.original[i][k] - expected;
                    let allowance = controls.absolute_equation_tolerance_mol_m2[k];
                    if residual.abs() > allowance {
                        return Err(DiluteSoluteError::Sorption {
                            cell: i,
                            species: k,
                            residual_mol_m2: residual,
                            allowance_mol_m2: allowance,
                        });
                    }
                }
            }
        }
        Ok(())
    }

    fn terms(&self, i: usize, k: usize) -> Result<(f64, f64)> {
        let alpha = self.laws[i][k].weight(self.duration_s);
        Ok((
            integrated_flux(alpha, self.laws[i][k].capacity())?,
            integrated_flux(alpha, self.original[i][k])?,
        ))
    }

    fn finish_cell(
        &mut self,
        i: usize,
        k: usize,
        c: f64,
        controls: DiluteSoluteControls<N>,
    ) -> Result<()> {
        let (capacity, release) = self.terms(i, k)?;
        let uptake = integrated_flux(capacity, c)?;
        let retained = self.laws[i][k].retained(self.original[i][k], self.duration_s)?;
        let endpoint = finite(retained + uptake)?;
        let actual = finite(endpoint - self.original[i][k])?;
        let predicted = finite(uptake - release)?;
        let residual = finite(actual - predicted)?;
        let allowance = finite(
            controls.absolute_equation_tolerance_mol_m2[k]
                + controls.relative_equation_tolerance * actual.abs().max(predicted.abs()),
        )?;
        if !nonnegative(endpoint) || residual.abs() > allowance {
            return Err(DiluteSoluteError::Sorption {
                cell: i,
                species: k,
                residual_mol_m2: residual,
                allowance_mol_m2: allowance,
            });
        }
        self.bound[i][k] = endpoint;
        self.exchange[i][k] = actual;
        self.maximum_defect[k] = self.maximum_defect[k].max(residual.abs());
        Ok(())
    }

    fn amounts(&self, i: usize, k: usize) -> [f64; 3] {
        [self.original[i][k], self.bound[i][k], self.exchange[i][k]]
    }
}

#[cfg(test)]
mod tests;
