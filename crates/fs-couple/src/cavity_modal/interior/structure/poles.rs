//! Free resonances of the complete original structure and interior pressure field.
//!
//! The shared quadratic-mode owner solves `(s^2 M + s C + K) phi = 0` using
//! the same assembled arrays as the forced response and exact time owners.
//! The interface blocks in C stay opposite transposes; neither the bare fluid
//! nor the bare structure is inverted or reduced. For `exp(s t)`, pressure is
//! `p = s psi`, including the REAL decay contribution, not `-i Im(s) psi`.
//!
//! These are seeded local free modes, not driven amplitude peaks, a complete
//! spectrum, continuum eigenvalue enclosures or a bandwidth-derived Q.
//! The original spatial, frozen-material and small-motion assumptions remain.

use super::{StructuralInterior, StructuralInteriorError, poll, response};
use fs_exec::Cx;
use fs_math::{c64::C64, dd::Dd};
use fs_modal::quadratic_modes::{PreparedQuadraticModes, QuadraticModeError};

/// Shared accuracy, original-coordinate, payload and arithmetic limits.
pub use fs_modal::quadratic_modes::PoleOptions;

/// Hermitian modal-envelope storage [J], with factor 1/2 in each quadratic
/// form. This is NOT a finite-cycle average of the decaying real trajectory.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StructuralInteriorPoleEnergy {
    /// `v^H M_structure v / 2`.
    pub structural_kinetic_j: f64,
    /// `u^H K_structure u / 2`.
    pub structural_strain_j: f64,
    /// `psi^H K_acoustic psi / 2`.
    pub fluid_kinetic_j: f64,
    /// `p^H M_acoustic p / 2`.
    pub fluid_compression_j: f64,
}

impl StructuralInteriorPoleEnergy {
    /// Original structural kinetic plus strain envelope energy.
    #[must_use]
    pub fn structural_j(self) -> f64 {
        self.structural_kinetic_j + self.structural_strain_j
    }

    /// Original fluid kinetic plus compression envelope energy.
    #[must_use]
    pub fn acoustic_j(self) -> f64 {
        self.fluid_kinetic_j + self.fluid_compression_j
    }

    /// Complete structural and acoustic envelope storage.
    #[must_use]
    pub fn total_j(self) -> f64 {
        self.structural_j() + self.acoustic_j()
    }
}

/// One free mode in every original structural and pressure coordinate.
/// The concatenated displacement/potential shape has `phi^H M phi = 1`.
#[derive(Clone, Debug)]
pub struct StructuralInteriorPole {
    /// Laplace exponent for exp(s t) [1/s]. Growth is retained, never clipped.
    pub laplace_pole_per_s: C64,
    /// Nonnegative oscillation rate |Im(s)|; real relaxation has zero frequency.
    pub oscillation_frequency_rad_s: f64,
    /// Signed amplitude decay -Re(s) [1/s], negative for an unstable mode.
    pub amplitude_decay_rate_per_s: f64,
    /// Original free structural-coordinate displacement shape.
    pub structural_displacement: Vec<C64>,
    /// Original structural velocity `s*u`.
    pub structural_velocity: Vec<C64>,
    /// Every original pressure-potential coefficient, including constant potential.
    pub pressure_potential_pa_s: Vec<C64>,
    /// Physical pressure `s * psi`, retaining decay and oscillation together.
    pub pressure_pa: Vec<C64>,
    /// Original generalized fluid force `B^T p`.
    pub fluid_force: Vec<C64>,
    /// Original pressure-trace integrals `B v` of structural normal velocity.
    pub displaced_flow_m3_s: Vec<C64>,
    /// Independently evaluated original structural and acoustic storage forms.
    pub energy: StructuralInteriorPoleEnergy,
    /// Real Hermitian envelope work, positive from fluid into structure [W].
    pub fluid_to_structure_power_w: f64,
    /// Work of the original symmetric structural velocity operator [W].
    /// The gyroscopic part is not counted as dissipation.
    pub structural_dissipation_w: f64,
    /// Independently evaluated original relative-wall admittance work [W].
    pub wall_dissipation_w: f64,
    /// Independently evaluated `2 Re(s) H + D_structure + D_wall` [W].
    pub decay_power_residual_w: f64,
    /// Original-row backward residual of the published full quadratic mode.
    /// This is not a forward-error or eigenvalue-enclosure certificate.
    pub relative_residual: f64,
    /// Last Newton correction divided by max(reference rate, pole magnitude).
    pub relative_correction: f64,
    /// Number of original-pencil Newton evaluations.
    pub iterations: usize,
}

/// Original physical and shared numerical refusals; no partial mode is returned.
#[derive(Debug)]
pub enum StructuralInteriorPoleError {
    /// Original physical-coordinate, storage, work or context refusal.
    Physical(StructuralInteriorError),
    /// Shared quadratic mode admission, resource, convergence or cancellation.
    Mode(QuadraticModeError),
}

impl core::fmt::Display for StructuralInteriorPoleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Physical(error) => error.fmt(f),
            Self::Mode(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for StructuralInteriorPoleError {}
impl From<StructuralInteriorError> for StructuralInteriorPoleError {
    fn from(error: StructuralInteriorError) -> Self {
        Self::Physical(error)
    }
}
impl From<QuadraticModeError> for StructuralInteriorPoleError {
    fn from(error: QuadraticModeError) -> Self {
        Self::Mode(error)
    }
}
impl From<response::ResponseError> for StructuralInteriorPoleError {
    fn from(error: response::ResponseError) -> Self {
        Self::Physical(error.into())
    }
}

/// Cached shared numerical preparation borrowing the immutable full model.
/// Subsequent modes reuse equilibration and the original complete mass factor.
pub struct PreparedStructuralInteriorModes<'m, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    modes: PreparedQuadraticModes<'m>,
}

impl<'a, 'f, 'c> StructuralInterior<'a, 'f, 'c> {
    /// Prepare free modes without eliminating the constant pressure potential,
    /// off-diagonal inertia, structural gyroscopic coupling or wall losses.
    ///
    /// `options` bounds this preparation and EACH complete refinement. Native
    /// quadratic-mode resources are augmented by `128*n+256` payload bytes and
    /// work score `128*n*n+256*n` for original seeds, fields, work and storage;
    /// n counts all structural and pressure coordinates. Original model arrays,
    /// caller inputs, earlier returned modes and allocator metadata are excluded.
    ///
    /// Mass is admitted by the shared SPD owner; original stiffness must be
    /// exactly symmetric for the physical storage identity. All original C
    /// entries survive. Physical storage checks remain in force, and no decay
    /// sign, rigid mode or repeated root is repaired with an artificial floor.
    pub fn prepare_resonance_modes(
        &self,
        mut options: PoleOptions,
        cx: &Cx<'_>,
    ) -> Result<PreparedStructuralInteriorModes<'_, 'a, 'f, 'c>, StructuralInteriorPoleError> {
        poll(cx)?;
        let n = self.structural_dofs + self.pressure_dofs();
        let size = n as u128;
        let bytes = 128 * size + 256;
        let work = 128 * size * size + 256 * size;
        for (quantity, requested, cap) in [
            (
                "coupled mode observation bytes",
                bytes,
                options.max_dense_bytes as u128,
            ),
            (
                "coupled mode observation work",
                work,
                options.max_work as u128,
            ),
        ] {
            if requested > cap {
                return Err(QuadraticModeError::Budget {
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
                return Err(StructuralInteriorError::Invalid(
                    "coupled free-mode storage requires exactly symmetric stiffness",
                )
                .into());
            }
        }
        let modes = PreparedQuadraticModes::new(
            &self.mass,
            &self.damping,
            &self.stiffness,
            n,
            options,
            || cx.checkpoint().is_err(),
        )?;
        poll(cx)?;
        Ok(PreparedStructuralInteriorModes { model: self, modes })
    }
}

impl PreparedStructuralInteriorModes<'_, '_, '_, '_> {
    /// Refine one simple free pole from original structural displacement and
    /// pressure-potential seeds. Seeds need not be normalized. A nearby physical
    /// harmonic displacement/potential pair is a useful starting direction.
    ///
    /// Pressure is recovered as `s*psi`, including the decay term. Envelope
    /// storage and heat are evaluated independently in the original blocks;
    /// neither is inferred from the pole or an energy deficit. Repeated and
    /// defective roots may refuse. No accepted state, clock or material changes.
    pub fn refine(
        &self,
        laplace_seed_per_s: C64,
        structural_displacement_seed: &[C64],
        pressure_potential_seed_pa_s: &[C64],
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorPole, StructuralInteriorPoleError> {
        poll(cx)?;
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        if structural_displacement_seed.len() != ns || pressure_potential_seed_pa_s.len() != na {
            return Err(StructuralInteriorError::Invalid(
                "free-mode seeds must cover every original structural and pressure coordinate",
            )
            .into());
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
        // The existing owner uses peak-phasor cycle averages. A free complex
        // modal envelope has exactly twice that quadratic storage and work.
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
        let s = mode.laplace_frequency_per_s;
        let residual = (Dd::from_f64(2.0) * Dd::from_f64(s.re) * Dd::from_f64(energy.total_j())
            + Dd::from_f64(structural_dissipation_w)
            + Dd::from_f64(wall_dissipation_w))
        .to_f64();
        if [
            energy.total_j(),
            structural_dissipation_w,
            wall_dissipation_w,
            exchange,
            residual,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(StructuralInteriorError::Invalid(
                "coupled free-mode envelope storage or work is unrepresentable",
            )
            .into());
        }
        poll(cx)?;
        Ok(StructuralInteriorPole {
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
            fluid_to_structure_power_w: exchange,
            structural_dissipation_w,
            wall_dissipation_w,
            decay_power_residual_w: residual,
            relative_residual: mode.relative_residual,
            relative_correction: mode.relative_correction,
            iterations: mode.iterations,
        })
    }
}
