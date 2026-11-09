//! Source participation in the complete cavity/wall/Maxwell realization.
//!
//! A single shared bordered factor retains nonsymmetric left/right participation.
//! Sources are evaluated at the Laplace pole: acoustic velocity injection is
//! s*b, whereas applied mechanical traction is W*t. No local wall impedance,
//! damping fit, modal mass approximation or real-frequency source is substituted.

use super::{
    CoupledPencil, MechanicalPoleError, MechanicalPoleLimits, MechanicalResonancePole,
    NonlinearPole, Query, finite, poll, wall_fields,
};
use super::super::WallTraction;
use fs_la::pole_response::{PoleResponseError, PreparedPole};
pub use fs_la::pole_response::PoleResponseOptions as MechanicalPoleResponseOptions;
use fs_math::{c64::C64, dd::Dd};

mod sensitivity;
pub use sensitivity::{
    MechanicalRelaxationParameterDirection, MechanicalResonanceParameterDerivative,
    MechanicalWallParameterDerivative, MechanicalWallParameterDirection,
};

/// Laplace residues of the original wall coordinates; each residue has one
/// additional inverse second compared with the corresponding modal field.
#[derive(Clone, Debug)]
pub struct MechanicalWallResidue {
    pub displacement_m_per_s: Vec<C64>,
    pub normal_velocity_m_s2: Vec<C64>,
    pub relaxation_memory_per_s: Vec<Vec<C64>>,
    pub relaxation_extension_per_s: Vec<Vec<C64>>,
}

/// R = lim_(s -> s0) (s-s0) x(s), with the complete analytic source at s0.
/// Near the isolated pole, x(s) = R/(s-s0) + an analytic remainder. Under this
/// owner's exp(+i*omega*t), the angular-frequency residue is -i*R.
/// This is one modal contribution, not a full forced field or a power.
#[derive(Clone, Debug)]
pub struct MechanicalResonanceResidue {
    pub rate_per_s: C64,
    /// Complete free pressure/wall/extension order, without renormalization.
    pub coordinates: Vec<C64>,
    /// Exact pressure-release zeros remain present.
    pub pressure_pa_per_s: Vec<C64>,
    pub walls: Vec<MechanicalWallResidue>,
    pub bordered_relative_residual: f64,
    pub pole_physical_residual: f64,
}

/// Immutable full bordered factor, bound to the original physical owner.
pub(crate) struct PreparedCoupledResonance<'a> {
    source: &'a CoupledPencil,
    numerical: PreparedPole,
    mode: MechanicalResonancePole,
}

fn response_error(error: PoleResponseError) -> MechanicalPoleError {
    match error {
        PoleResponseError::Cancelled => MechanicalPoleError::Cancelled,
        error => MechanicalPoleError::Response(error),
    }
}

impl CoupledPencil {
    pub(crate) fn prepare_resonance_response(
        &self,
        pole: &MechanicalResonancePole,
        options: MechanicalPoleResponseOptions,
        limits: MechanicalPoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedCoupledResonance<'_>, MechanicalPoleError> {
        poll(cancelled)?;
        if pole.rate_per_s == C64::ZERO
            || !finite(pole.rate_per_s)
            || pole.coordinates.len() != self.free_dofs.len()
            || pole.coordinates.iter().any(|&v| !finite(v))
        {
            return Err(MechanicalPoleError::Invalid(
                "complete finite nonzero physical pole and coordinates required",
            ));
        }
        let mut query = Query::new(self, limits)?;
        let seed: Vec<_> = pole.coordinates.iter().zip(&self.scales)
            .map(|(&x, &d)| x / C64::from_re(d)).collect();
        let evaluation = query.evaluate(pole.rate_per_s, cancelled)?;
        let numerical = PreparedPole::new(
            pole.rate_per_s, &seed, &evaluation, options, &mut *cancelled,
        ).map_err(response_error)?;
        // Re-evaluate every original physical row. A supplied pole report,
        // even one from a different owner, never bypasses this admission.
        let mode = query.observe(NonlinearPole {
            value: numerical.value(),
            vector: numerical.right_vector().to_vec(),
            residual: numerical.pole_residual(),
            relative_step: numerical.relative_correction(),
            residual_scale: evaluation.residual_scale,
            iterations: 0,
        }, cancelled)?;
        poll(cancelled)?;
        Ok(PreparedCoupledResonance { source: self, numerical, mode })
    }
}

impl PreparedCoupledResonance<'_> {
    pub(crate) fn mode(&self) -> &MechanicalResonancePole {
        &self.mode
    }

    /// Reuses the complete factor for each coherent source query, O(n^2).
    /// The source columns have already been bound/projected by the public FE
    /// owner. There is no assumption of positive or independently added powers.
    pub(crate) fn forcing_residue(
        &self,
        load: &[C64],
        tractions: &[(&WallTraction, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalResonanceResidue, MechanicalPoleError> {
        poll(cancelled)?;
        let source = self.source;
        let s = self.numerical.value();
        if load.len() != source.fluid_dofs || load.iter().any(|&v| !finite(v)) {
            return Err(MechanicalPoleError::Invalid(
                "original finite acoustic source at the pole required",
            ));
        }
        let mut forcing = vec![[Dd::ZERO; 2]; source.blocks[0].nrows()];
        for (row, &b) in load.iter().enumerate() {
            poll(cancelled)?;
            forcing[row] = [
                Dd::from_f64(s.re) * Dd::from_f64(b.re)
                    - Dd::from_f64(s.im) * Dd::from_f64(b.im),
                Dd::from_f64(s.re) * Dd::from_f64(b.im)
                    + Dd::from_f64(s.im) * Dd::from_f64(b.re),
            ];
        }
        for &(traction, amplitude) in tractions {
            poll(cancelled)?;
            if !finite(amplitude)
                || traction.wall >= source.walls.len()
                || traction.force.len() != source.walls[traction.wall].dofs.len()
            {
                return Err(MechanicalPoleError::Invalid(
                    "original complete wall traction and finite amplitude required",
                ));
            }
            let offset = source.offsets[traction.wall] + traction.force.len();
            for (i, &force) in traction.force.iter().enumerate() {
                poll(cancelled)?;
                forcing[offset + i][0] = forcing[offset + i][0]
                    + force * Dd::from_f64(amplitude.re);
                forcing[offset + i][1] = forcing[offset + i][1]
                    + force * Dd::from_f64(amplitude.im);
            }
        }
        // For T_hat = D*T*D and x=D*y, the physical source is D*f.
        let mut rhs = Vec::with_capacity(source.free_dofs.len());
        for (&row, &scale) in source.free_dofs.iter().zip(&source.scales) {
            poll(cancelled)?;
            let value = C64::new(
                (forcing[row][0] * Dd::from_f64(scale)).to_f64(),
                (forcing[row][1] * Dd::from_f64(scale)).to_f64(),
            );
            if !finite(value) {
                return Err(MechanicalPoleError::Invalid(
                    "complete pole source is unrepresentable",
                ));
            }
            rhs.push(value);
        }
        let mut residues = self.numerical.forcing_residues(
            &[rhs.as_slice()], &mut *cancelled,
        ).map_err(response_error)?;
        let residue = residues.pop().ok_or(MechanicalPoleError::Invalid(
            "complete source residue was not returned",
        ))?;
        let coordinates: Vec<_> = residue.vector.iter().zip(&source.scales)
            .map(|(&x, &d)| x.scale(d)).collect();
        if coordinates.iter().any(|&z| !finite(z)) {
            return Err(MechanicalPoleError::Invalid(
                "original physical pole residue is unrepresentable",
            ));
        }
        let mut full = vec![C64::ZERO; source.blocks[0].nrows()];
        for (&row, &x) in source.free_dofs.iter().zip(&coordinates) {
            full[row] = x;
        }
        let walls = wall_fields(source, s, &full, cancelled)?.into_iter()
            .map(|w| MechanicalWallResidue {
                displacement_m_per_s: w.displacement_m,
                normal_velocity_m_s2: w.normal_velocity_m_s,
                relaxation_memory_per_s: w.relaxation_memory,
                relaxation_extension_per_s: w.relaxation_extension,
            }).collect();
        poll(cancelled)?;
        Ok(MechanicalResonanceResidue {
            rate_per_s: s,
            coordinates,
            pressure_pa_per_s: full[..source.fluid_dofs].to_vec(),
            walls,
            bordered_relative_residual: residue.linear_residual,
            pole_physical_residual: self.mode.physical_residual,
        })
    }
}

#[cfg(test)]
mod tests;
