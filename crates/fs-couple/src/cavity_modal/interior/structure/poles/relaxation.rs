//! Complete coupled free modes with original generalized-Maxwell memory.
//!
//! The same original arrays enter the shared rational mode owner:
//! `M q_ddot + C q_dot + K_relaxed q + sum(b_j sigma_j) = 0`,
//! `sigma_dot_j = k_j b_j^T q_dot - sigma_j/tau_j`.
//! Here q=[u,psi], p=psi_dot, and the two original B blocks remain reciprocal.
//! No coefficient is evaluated at a guessed real resonance frequency. The
//! caller supplies the physical relaxation channels; no material law, mesh
//! derivative, frequency fit or shell constitutive image is inferred.

use super::{
    C64, Cx, Dd, StructuralInterior, StructuralInteriorError, StructuralInteriorPoleEnergy,
};
use super::{poll, response};
use fs_modal::relaxation_modes::{MaxwellModeError, PreparedMaxwellModes};

/// Original scalar channels and shared numerical accuracy/resource controls.
pub use fs_modal::relaxation_modes::{MaxwellDampingBranch, PoleOptions};

/// Original physical or native rational-mode refusal; no partial mode is returned.
#[derive(Debug)]
pub enum StructuralInteriorRelaxationError {
    Physical(StructuralInteriorError),
    Mode(MaxwellModeError),
}
impl core::fmt::Display for StructuralInteriorRelaxationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Physical(error) => error.fmt(f),
            Self::Mode(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for StructuralInteriorRelaxationError {}
impl From<StructuralInteriorError> for StructuralInteriorRelaxationError {
    fn from(error: StructuralInteriorError) -> Self {
        Self::Physical(error)
    }
}
impl From<MaxwellModeError> for StructuralInteriorRelaxationError {
    fn from(error: MaxwellModeError) -> Self {
        Self::Mode(error)
    }
}
impl From<response::ResponseError> for StructuralInteriorRelaxationError {
    fn from(error: response::ResponseError) -> Self {
        Self::Physical(error.into())
    }
}

/// One original branch in supplied source order; these are modal envelopes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StructuralInteriorRelaxationMemory {
    /// Stress conjugate to the supplied scalar strain b^T [u,psi].
    pub stress: C64,
    /// Stored internal-variable energy |sigma|²/(2 k) [J].
    pub energy_j: f64,
    /// Independent positive relaxation loss |sigma|²/(k tau) [W].
    pub dissipation_w: f64,
}

/// One local rational free mode, normalized by phi^H M phi=1 in original units.
/// Internal memory storage is separate from the relaxed spatial storage.
#[derive(Clone, Debug)]
pub struct StructuralInteriorRelaxationPole {
    pub laplace_pole_per_s: C64,
    /// |Im(s)|; a real relaxation mode has zero oscillation frequency.
    pub oscillation_frequency_rad_s: f64,
    /// -Re(s), preserving growth when the original mechanical source is active.
    pub amplitude_decay_rate_per_s: f64,
    pub structural_displacement: Vec<C64>,
    /// Original s*u, including the real decay contribution.
    pub structural_velocity: Vec<C64>,
    pub pressure_potential_pa_s: Vec<C64>,
    /// Original s*psi, including both decay and oscillation.
    pub pressure_pa: Vec<C64>,
    /// Original fluid interface force B^T*p, excluding explicit branch forces.
    pub fluid_force: Vec<C64>,
    /// Original structural interface flow B*v, excluding explicit branch flows.
    pub displaced_flow_m3_s: Vec<C64>,
    /// Hermitian envelope storage with factor 1/2 and RELAXED original K.
    pub energy: StructuralInteriorPoleEnergy,
    pub memory: Vec<StructuralInteriorRelaxationMemory>,
    /// Sum of the branch energies above, not an additional energy account.
    pub relaxation_energy_j: f64,
    /// Sum of the branch losses above, not an additional heat debit.
    pub relaxation_dissipation_w: f64,
    /// Real interface work, positive from fluid into structure.
    pub fluid_to_structure_power_w: f64,
    /// Original symmetric structural velocity work; skew work creates no heat.
    pub structural_dissipation_w: f64,
    /// Original memoryless relative-wall loss, separate from relaxation channels.
    pub wall_dissipation_w: f64,
    /// 2 Re(s)*(H_spatial+H_memory)+D_structure+D_wall+D_memory [W].
    pub decay_power_residual_w: f64,
    /// Native backward residual of all original mechanical/acoustic rows.
    pub mechanical_relative_residual: f64,
    /// Native backward residual of all original branch history equations.
    pub relaxation_relative_residual: f64,
    pub relative_correction: f64,
    pub iterations: usize,
}

/// Immutable full coupled model and original branches, with one native preparation.
pub struct PreparedStructuralInteriorRelaxationModes<'m, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    branches: &'m [MaxwellDampingBranch],
    modes: PreparedMaxwellModes<'m>,
}

impl<'a, 'f, 'c> StructuralInterior<'a, 'f, 'c> {
    /// Prepare free poles of the complete structure/cavity with explicit memory.
    ///
    /// Every branch strain has ns+na entries in original [u,psi] order. A
    /// structural-only channel is a mechanical Maxwell arm. A potential-only
    /// channel contributes the passive wall admittance k*tau/(1+s*tau) b*b^T.
    /// Mixed channels are accepted only as the caller's explicit physical law.
    /// Existing K must be RELAXED: branch stiffness must not already be in it.
    /// The source's material/domain authority remains with its caller.
    ///
    /// Native PoleOptions bound preparation and each refinement separately.
    /// This adapter reserves 160*n+64*b+512 payload bytes and work score
    /// 160*n*n+256*n+256*b for original seeds, fields and independent storage/
    /// work recovery; n=ns+na and b=branches.len(). Borrowed source arrays,
    /// previous results and allocator overhead are excluded. Native branch
    /// assembly, binary equilibration and mass factorization keep their bounds.
    ///
    /// This displacement-based rational owner excludes s=-1/tau and pure
    /// material dark modes with zero displacement. It computes seeded nearby
    /// simple modes, not a complete spectrum or a continuum error enclosure.
    pub fn prepare_relaxation_modes<'m>(
        &'m self,
        branches: &'m [MaxwellDampingBranch],
        mut options: PoleOptions,
        cx: &Cx<'_>,
    ) -> Result<
        PreparedStructuralInteriorRelaxationModes<'m, 'a, 'f, 'c>,
        StructuralInteriorRelaxationError,
    > {
        poll(cx)?;
        let n = self.structural_dofs + self.pressure_dofs();
        let (size, count) = (n as u128, branches.len() as u128);
        let bytes = 160 * size + 64 * count + 512;
        let work = 160 * size * size + 256 * size + 256 * count;
        for (quantity, requested, cap) in [
            (
                "coupled relaxation observation bytes",
                bytes,
                options.max_dense_bytes as u128,
            ),
            (
                "coupled relaxation observation work",
                work,
                options.max_work as u128,
            ),
        ] {
            if requested > cap {
                return Err(MaxwellModeError::Budget {
                    quantity,
                    requested,
                    cap,
                }
                .into());
            }
        }
        options.max_dense_bytes -= bytes as usize;
        options.max_work -= work as u64;
        for i in 0..n {
            poll(cx)?;
            if (0..i).any(|j| self.stiffness[i * n + j] != self.stiffness[j * n + i]) {
                return Err(invalid(
                    "coupled relaxation storage requires symmetric relaxed stiffness",
                ));
            }
        }
        let modes = PreparedMaxwellModes::new(
            &self.mass,
            &self.damping,
            &self.stiffness,
            n,
            branches,
            options,
            || cx.checkpoint().is_err(),
        )?;
        poll(cx)?;
        Ok(PreparedStructuralInteriorRelaxationModes {
            model: self,
            branches,
            modes,
        })
    }
}

impl PreparedStructuralInteriorRelaxationModes<'_, '_, '_, '_> {
    /// Refine one nearby original-coordinate mode without changing any source.
    /// The complete material pencil and its analytic s derivative are evaluated
    /// by the existing native owner. Pressure is s*psi, never -i*Im(s)*psi.
    /// Physical branch storage and heat are evaluated from original stresses,
    /// separately from the source-row residual and measured decay rate.
    pub fn refine(
        &self,
        laplace_seed_per_s: C64,
        structural_displacement_seed: &[C64],
        pressure_potential_seed_pa_s: &[C64],
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorRelaxationPole, StructuralInteriorRelaxationError> {
        poll(cx)?;
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        if structural_displacement_seed.len() != ns || pressure_potential_seed_pa_s.len() != na {
            return Err(invalid(
                "relaxation seeds must cover every original displacement and potential",
            ));
        }
        let mut seed = Vec::with_capacity(ns + na);
        seed.extend_from_slice(structural_displacement_seed);
        seed.extend_from_slice(pressure_potential_seed_pa_s);
        let mut mode = self
            .modes
            .refine(laplace_seed_per_s, &seed, || cx.checkpoint().is_err())?;
        let potential = mode.displacement.split_off(ns);
        let pressure = mode.velocity.split_off(ns);
        let mut cancelled = || cx.checkpoint().is_err();
        let (fluid_force, flow) = response::interface(
            &self.model.coupling,
            &mode.velocity,
            &pressure,
            &mut cancelled,
        )?;
        let half = response::energy(
            &self.model.mass,
            &self.model.stiffness,
            &mode.displacement,
            &mode.velocity,
            &potential,
            &pressure,
            &mut cancelled,
        )?;
        let energy = StructuralInteriorPoleEnergy {
            structural_kinetic_j: 2.0 * half.structural_kinetic_j,
            structural_strain_j: 2.0 * half.structural_strain_j,
            fluid_kinetic_j: 2.0 * half.fluid_kinetic_j,
            fluid_compression_j: 2.0 * half.fluid_compression_j,
        };
        let structural_dissipation_w = 2.0
            * response::dissipation(
                &self.model.damping,
                ns + na,
                0,
                &mode.velocity,
                &mut cancelled,
            )?;
        let wall_dissipation_w = 2.0
            * response::dissipation(&self.model.damping, ns + na, ns, &pressure, &mut cancelled)?;
        let exchange = 2.0 * response::work(&fluid_force, &mode.velocity);
        let mut memory = Vec::with_capacity(self.branches.len());
        let (mut stored, mut dissipated) = (Dd::from_f64(0.0), Dd::from_f64(0.0));
        for (branch, stress) in self.branches.iter().zip(mode.branch_stress) {
            poll(cx)?;
            let energy_j = stress_square_over(stress, branch.stiffness, 2.0).map_err(invalid)?;
            let dissipation_w =
                stress_square_over(stress, branch.stiffness, branch.relaxation_time_s)
                    .map_err(invalid)?;
            memory.push(StructuralInteriorRelaxationMemory {
                stress,
                energy_j,
                dissipation_w,
            });
            stored = stored + Dd::from_f64(energy_j);
            dissipated = dissipated + Dd::from_f64(dissipation_w);
        }
        let total_energy = stored
            + Dd::from_f64(energy.structural_kinetic_j)
            + Dd::from_f64(energy.structural_strain_j)
            + Dd::from_f64(energy.fluid_kinetic_j)
            + Dd::from_f64(energy.fluid_compression_j);
        let s = mode.laplace_frequency_per_s;
        let residual = (Dd::from_f64(2.0) * Dd::from_f64(s.re) * total_energy
            + Dd::from_f64(structural_dissipation_w)
            + Dd::from_f64(wall_dissipation_w)
            + dissipated)
            .to_f64();
        if [
            total_energy.to_f64(),
            stored.to_f64(),
            dissipated.to_f64(),
            structural_dissipation_w,
            wall_dissipation_w,
            exchange,
            residual,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(invalid(
                "coupled relaxation envelope storage or work is unrepresentable",
            ));
        }
        poll(cx)?;
        Ok(StructuralInteriorRelaxationPole {
            laplace_pole_per_s: s,
            oscillation_frequency_rad_s: s.im.abs(),
            amplitude_decay_rate_per_s: -s.re,
            structural_displacement: mode.displacement,
            structural_velocity: mode.velocity,
            pressure_potential_pa_s: potential,
            pressure_pa: pressure,
            fluid_force,
            displaced_flow_m3_s: flow,
            energy,
            memory,
            relaxation_energy_j: stored.to_f64(),
            relaxation_dissipation_w: dissipated.to_f64(),
            fluid_to_structure_power_w: exchange,
            structural_dissipation_w,
            wall_dissipation_w,
            decay_power_residual_w: residual,
            mechanical_relative_residual: mode.mechanical_relative_residual,
            relaxation_relative_residual: mode.relaxation_relative_residual,
            relative_correction: mode.relative_correction,
            iterations: mode.iterations,
        })
    }
}

// Compute |stress|²/(a*b) without ever forming a stress square or a*b in
// physical units. Positive normalized mantissas keep every DD operation
// bounded; only the final published value is rescaled. Storage and heat are
// calculated independently, so a very fast branch cannot lose its heat through
// an underflowed stored-energy intermediate.
fn stress_square_over(stress: C64, a: f64, b: f64) -> Result<f64, &'static str> {
    if !stress.re.is_finite()
        || !stress.im.is_finite()
        || !a.is_finite()
        || a <= 0.0
        || !b.is_finite()
        || b <= 0.0
    {
        return Err("finite original branch stress and positive coefficients required");
    }
    let magnitude = stress.re.abs().max(stress.im.abs());
    if magnitude == 0.0 {
        return Ok(0.0);
    }
    let (m, me) = fraction(magnitude);
    let (a, ae) = fraction(a);
    let (b, be) = fraction(b);
    let re = Dd::from_f64(stress.re / magnitude);
    let im = Dd::from_f64(stress.im / magnitude);
    let mut value = ((re * re + im * im) * Dd::from_f64(m) * Dd::from_f64(m)
        / Dd::from_f64(a)
        / Dd::from_f64(b))
    .to_f64();
    let mut exponent = 2 * me - ae - be;
    while exponent != 0 && value.is_finite() && value != 0.0 {
        let step = exponent.clamp(-1022, 1023);
        value *= f64::from_bits(((step + 1023) as u64) << 52);
        exponent -= step;
    }
    if !value.is_finite() || value == 0.0 {
        return Err("original branch storage or loss is unrepresentable");
    }
    Ok(value)
}

fn fraction(value: f64) -> (f64, i32) {
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    if exponent == 0 {
        let (fraction, exponent) = fraction(value * f64::from_bits(1075_u64 << 52));
        (fraction, exponent - 52)
    } else {
        (
            f64::from_bits((bits & ((1_u64 << 52) - 1)) | (1023_u64 << 52)),
            exponent - 1023,
        )
    }
}

fn invalid(what: &'static str) -> StructuralInteriorRelaxationError {
    StructuralInteriorError::Invalid(what).into()
}

#[cfg(test)]
mod arithmetic_tests {
    use super::{C64, stress_square_over};

    #[test]
    fn original_memory_energy_and_loss_survive_extreme_but_representable_units() {
        for scale in [1e-200, 1.0, 1e200] {
            let stress = C64::new(0.6 * scale, -0.8 * scale);
            let energy = stress_square_over(stress, scale, 2.0).unwrap();
            let loss = stress_square_over(stress, scale, scale).unwrap();
            assert!((energy / (0.5 * scale) - 1.0).abs() < 5e-15);
            assert!((loss - 1.0).abs() < 5e-15);
        }
        let subnormal = f64::from_bits(1);
        assert_eq!(
            stress_square_over(C64::from_re(subnormal), subnormal, subnormal).unwrap(),
            1.0
        );
        assert_eq!(stress_square_over(C64::ZERO, 2.0, 3.0).unwrap(), 0.0);
        assert!(stress_square_over(C64::ONE, 0.0, 1.0).is_err());
        assert!(stress_square_over(C64::from_re(1e200), 1e-200, 1.0).is_err());
        assert!(stress_square_over(C64::from_re(1e-200), 1e200, 1.0).is_err());
    }
}
