use super::super::gram::{GramWorkspace, clear};
use super::*;

/// Sparse original bank/frequency block. Fields already carry the shared
/// joint aperture/shutter center rotation; no source covariance is inferred.
pub(super) struct SourceGroup {
    pub control: Option<usize>,
    pub frequency: usize,
    pub fields: Vec<FieldParts>,
}

pub(super) struct EnergyAction {
    pub frequencies: Vec<f64>,
    pub gradient: Vec<C64>,
    modes: usize,
    spatial: Option<Spatial<RealKernel>>,
    temporal: Vec<RealKernel>,
    groups: Vec<SourceGroup>,
    weights: [Dd; 6],
    gram: GramWorkspace,
}

impl EnergyAction {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        controls: usize,
        modes: usize,
        frequencies: Vec<f64>,
        spatial: Option<Spatial<RealKernel>>,
        temporal: Vec<RealKernel>,
        groups: Vec<SourceGroup>,
        epsilon: f64,
        mu: f64,
    ) -> Result<Self, AngularSpectrumError> {
        let gram = GramWorkspace::new(modes, frequencies.len(), frequencies.len(), &mut || false)?;
        Ok(Self {
            frequencies,
            gradient: filled(controls, C64::ZERO)?,
            modes,
            spatial,
            temporal,
            groups,
            weights: [epsilon, epsilon, epsilon, mu, mu, mu].map(Dd::from_f64),
            gram,
        })
    }

    /// Covector W*K*r is the derivative of J=r^T*W*K*r/2. Pullback uses
    /// the complex multiplication's real Jacobian: (Fr,Fi), (-Fi,Fr).
    pub fn evaluate(
        &mut self,
        coefficients: &[C64],
        include_fixed: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        poll(cancelled)?;
        if coefficients.len() != self.gradient.len()
            || coefficients
                .iter()
                .any(|a| !a.re.is_finite() || !a.im.is_finite())
        {
            return Err(bad(
                "optical control coefficients must match the original banks and be finite",
            ));
        }
        self.gradient.fill(C64::ZERO);
        self.combine(coefficients, include_fixed, cancelled)?;
        let f = self.frequencies.len();
        if f != 0 {
            let spatial = self
                .spatial
                .as_ref()
                .ok_or_else(|| bad("optical control spatial kernel is unavailable"))?;
            self.gram.real(spatial, &self.temporal, f, cancelled)?;
            for output in &mut self.gram.covectors {
                poll(cancelled)?;
                for (component, weight) in output.iter_mut().zip(self.weights) {
                    component[0] = component[0] * weight;
                    component[1] = component[1] * weight;
                }
                if output.iter().flatten().any(|&v| !finite_dd(v)) {
                    return Err(bad("optical control field covector is unrepresentable"));
                }
            }
        }
        let mut energy = Dd::ZERO;
        if include_fixed {
            // Recombine before the scalar contraction. Recovering energy
            // from large individual source-gradient terms would lose nulls.
            self.combine(coefficients, true, cancelled)?;
            for (fields, covectors) in self.gram.residual.iter().zip(&self.gram.covectors) {
                poll(cancelled)?;
                for (value, covector) in fields.iter().zip(covectors) {
                    energy = energy
                        + (value[0] * covector[0] + value[1] * covector[1]) * Dd::from_f64(0.5);
                }
            }
            if !finite_dd(energy)
                || !energy.to_f64().is_finite()
                || energy.hi < 0.0
                || (energy.hi == 0.0 && energy.lo < 0.0)
            {
                return Err(bad("optical control energy density is unrepresentable"));
            }
        }
        // Original groups are in bank order. Accumulate each source's full
        // frequency/modal derivative in DD before publishing its C64 entry.
        let mut cursor = 0;
        for control in 0..self.gradient.len() {
            poll(cancelled)?;
            let mut derivative = [Dd::ZERO; 2];
            while cursor < self.groups.len()
                && self.groups[cursor]
                    .control
                    .is_none_or(|index| index <= control)
            {
                let group = &self.groups[cursor];
                if group.control == Some(control) {
                    for (mode, fields) in group.fields.iter().enumerate() {
                        if mode % 32 == 0 {
                            poll(cancelled)?;
                        }
                        let covectors = &self.gram.covectors[group.frequency * self.modes + mode];
                        for ([fr, fi], [gr, gi]) in fields.iter().zip(covectors) {
                            derivative[0] = derivative[0] + *fr * *gr + *fi * *gi;
                            derivative[1] = derivative[1] - *fi * *gr + *fr * *gi;
                        }
                    }
                }
                cursor += 1;
            }
            if derivative
                .iter()
                .any(|&v| !finite_dd(v) || !v.to_f64().is_finite())
            {
                return Err(bad("optical source gradient is unrepresentable"));
            }
            self.gradient[control] = C64::new(derivative[0].to_f64(), derivative[1].to_f64());
        }
        poll(cancelled)?;
        Ok(energy.to_f64())
    }

    fn combine(
        &mut self,
        coefficients: &[C64],
        include_fixed: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        clear(&mut self.gram.residual, cancelled)?;
        for group in &self.groups {
            poll(cancelled)?;
            let coefficient = match group.control {
                Some(control) => coefficients[control],
                None if include_fixed => C64::ONE,
                None => continue,
            };
            let (ar, ai) = (Dd::from_f64(coefficient.re), Dd::from_f64(coefficient.im));
            for (mode, fields) in group.fields.iter().enumerate() {
                if mode % 32 == 0 {
                    poll(cancelled)?;
                }
                let result = &mut self.gram.residual[group.frequency * self.modes + mode];
                for (value, [fr, fi]) in result.iter_mut().zip(fields) {
                    value[0] = value[0] + ar * *fr - ai * *fi;
                    value[1] = value[1] + ar * *fi + ai * *fr;
                }
            }
        }
        for fields in &self.gram.residual {
            poll(cancelled)?;
            if fields.iter().flatten().any(|&v| !finite_dd(v)) {
                return Err(bad("optical control field is unrepresentable"));
            }
        }
        Ok(())
    }
}
