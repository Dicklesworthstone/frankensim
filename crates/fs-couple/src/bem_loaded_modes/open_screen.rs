//! Complete real-frequency reaction of a two-sided open shell.
//!
//! Source mass-normalized coordinates are unchanged. For peak phasors in
//! exp(-i omega t), the mechanical block is K - omega^2 I - i omega (C +
//! Z_screen), with every original Maxwell elastic-memory coordinate appended.
//! The same facet-local P1 work map produces the complete self/mutual Z and
//! coherent pressure-jump field. No baffled-piston law, fitted damping or
//! independent face duplication is introduced.
use super::{
    BemHarmonicPower, BemModesError, PhysicalPencil, bad, divide, finite,
    memory::RetainedMemory, poll,
};
use crate::plate_viscosity::ModalPlateViscosity;
use fs_bem::helmholtz::Medium;
use fs_bem::near_field;
use fs_bem::reciprocal_radiation::open_screen::{
    OpenOptions, OpenPanels, PreparedScreen, ScreenPortLimits, ScreenVelocityPortSolution,
};
use fs_la::eigen_complex::lu_complex_with_cancel;
use fs_la::nonlinear_eig::PoleEvaluation;
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;
use fs_plate::shell::reduction::ShellReduction;
use fs_plate::shell::reduction::radiation::open_screen::OpenScreenRadiationBudget;

/// Loaded modes of the same source/real-frequency boundary discretization.
pub mod poles;
pub use poles::{
    OpenScreenModePole, OpenScreenModeSpectrum, OpenScreenResonanceOptions,
    OpenScreenResonancePencil,
};

/// Original gas compression memory in the same two-sided shell reaction.
pub mod molecular;

/// Source geometry/tangent limits, separate from BEM and driven-solve limits.
#[derive(Clone, Copy, Debug)]
pub struct OpenScreenModesOptions {
    /// Every source mode is retained; at most 256 are admitted by this owner.
    pub maximum_modes: usize,
    pub source_surface: OpenScreenRadiationBudget,
    /// Uniform conforming subdivisions of the source's piecewise-flat surface.
    /// These are not edge grading or additional structural modes.
    pub acoustic_refinements: u32,
    /// Conservative geometry/projection/tangent payload, including refinement
    /// scratch but excluding the source ShellReduction and BEM caches.
    pub maximum_source_bytes: usize,
    /// Admission units: 128*source_facets*modes^2 +
    /// 64*final_facets*(modes+1)*(refinements+1). The existing source Hessian
    /// is evaluated once with its original cancellable facet/row loop.
    pub maximum_setup_work: u64,
}

/// Separate admission for one cold attachment of original material branches.
/// This payload is additional to the source geometry/K/C and acoustic caches.
#[derive(Clone, Copy, Debug)]
pub struct OpenScreenMemoryLimits {
    /// At most 256 branches; equal-time or dependent rows still count separately.
    pub maximum_branches: usize,
    /// 8*branches*modes + 32*branches + 256 logical retained bytes.
    /// Includes row headers/times, excludes allocator metadata and process RSS.
    pub maximum_payload_bytes: usize,
    /// 16*branches*(modes+1) units for original projection/scale validation.
    pub maximum_work: u64,
}

/// One immutable source-derived unit-modal-mass shell and exterior medium.
/// The full zero-state Hessian, viscosity and material arms retain their basis.
pub struct OpenScreenLoadedModes {
    surface: OpenPanels,
    weights: Vec<Vec<[f64; 3]>>,
    stiffness: Vec<f64>,
    viscosity: Vec<f64>,
    memory: RetainedMemory,
    medium: Medium,
}

/// Reusable original pressure-jump geometry/quadrature preparation. It borrows
/// the immutable source, so receiver fields cannot be rebound to another mesh.
pub struct PreparedOpenScreenModes<'s> {
    source: &'s OpenScreenLoadedModes,
    screen: PreparedScreen<'s>,
}

/// Explicit limits for one frequency and a complete batch of independent forces.
#[derive(Clone, Copy, Debug)]
pub struct OpenScreenHarmonicOptions {
    /// Original boundary assembly/LU, all source columns and retained jumps.
    pub boundary: ScreenPortLimits,
    /// 1..=256 complete force cases. Coherent forces must be summed FIRST.
    pub maximum_load_cases: usize,
    /// With d=modes+memories, logical live payload is
    /// cases*(16*d+32*modes+640) + 96*d^2 + 320*d + 4096 bytes.
    /// Includes full original/scaled matrices, LU, memory states and work/audit
    /// scratch. Source K/C/material rows/geometry, boundary columns/caches,
    /// allocator metadata and process RSS have separate admission.
    pub maximum_response_bytes: usize,
    /// 32*d^3 + 256*(cases+1)*d^2, d=modes+memories; boundary work is separate.
    pub maximum_work: u64,
    /// Original mechanical/material componentwise backward error, in (0,1).
    pub physical_residual_tolerance: f64,
    /// Input = viscosity + material relaxation + exact P1 screen work, relative
    /// to the sum of absolute work magnitudes; this is a numerical balance.
    pub power_relative_tolerance: f64,
}

/// Actual force-scaled source coordinates and signed reaction for one load case.
#[derive(Clone, Debug)]
pub struct OpenScreenHarmonicResponse {
    pub angular_frequency_rad_s: f64,
    /// Original source order, [m sqrt(kg)].
    pub modal_displacement: Vec<C64>,
    /// Original source order, [m sqrt(kg)/s].
    pub modal_velocity: Vec<C64>,
    /// Original Maxwell elastic extensions e=c.q-z, in branch order [sqrt(J)].
    pub elastic_memory_sqrt_j: Vec<C64>,
    /// Fluid force ON the structure, -Z*v, [N/sqrt(kg)].
    pub fluid_reaction: Vec<C64>,
    /// Source kinetic + equilibrium elastic + Maxwell elastic cycle-mean
    /// energy. No exterior storage is inferred from a single-frequency load.
    pub mean_structural_energy_j: f64,
    pub mechanical_residual: f64,
    /// Original retained material equations, independently audited.
    pub memory_residual: f64,
    /// Largest original unit-source solve residual. This is not a bound for a
    /// coherent field near cancellation, nor an estimate of spatial error.
    pub maximum_boundary_residual: f64,
    /// `radiated_w` is signed pressure-jump input work. The independently
    /// integrated coherent far-field power remains available on the batch.
    pub power: BemHarmonicPower,
}

/// One shared boundary solve and one shared mechanical LU for all load cases.
/// Fields combine the retained unit jumps with each actual loaded velocity;
/// neither receiver sampling nor angular power repeats a source solve.
pub struct OpenScreenHarmonicBatch<'p, 's> {
    radiation: ScreenVelocityPortSolution<'p, 's>,
    responses: Vec<OpenScreenHarmonicResponse>,
}

impl OpenScreenLoadedModes {
    /// Bind the supplied ShellReduction's own reference midsurface, complete
    /// mass-normalized translations and exact potential Hessian at zero.
    /// The original full ModalPlateViscosity is copied unchanged. Its producer
    /// must use this SAME basis; dimensions alone cannot prove provenance.
    ///
    /// Both faces contact one homogeneous exterior medium. This deliberately
    /// selects the zero-thickness pressure-jump chart: rotations have zero
    /// acoustic offset, while physical thickness remains in source mass/K/C.
    /// A cavity-backed face or a finite-thickness acoustic skin needs its own
    /// owner. This small-signal image contains no nonlinear/prestressed motion,
    /// hidden damping or inferred aperture. Original material relaxation is
    /// attached explicitly with with_relaxation_branches.
    pub fn from_shell(
        shell: &ShellReduction,
        medium: Medium,
        viscosity: Option<&ModalPlateViscosity>,
        options: OpenScreenModesOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BemModesError> {
        poll(&mut cancelled)?;
        let m = shell.mode_count();
        let facets = shell.facet_count();
        if m == 0
            || m > 256
            || m > options.maximum_modes
            || viscosity.is_some_and(|c| c.mode_count() != m)
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
        {
            return Err(bad(
                "open-screen source needs complete admitted shell mechanics and a positive medium",
            ));
        }
        let mut final_facets = facets;
        for _ in 0..options.acoustic_refinements {
            poll(&mut cancelled)?;
            final_facets = final_facets
                .checked_mul(4)
                .ok_or_else(|| bad("open-screen refined source size overflows"))?;
            if final_facets > options.source_surface.maximum_facets {
                return Err(bad("open-screen refined source exceeds its facet limit"));
            }
        }
        if final_facets > fs_bem::helmholtz::MAX_DENSE_PANELS {
            return Err(bad(
                "open-screen source exceeds the original dense screen panel limit",
            ));
        }
        let setup_work = 128 * (facets as u128) * (m as u128).pow(2)
            + 64 * (final_facets as u128)
                * (m as u128 + 1)
                * (options.acoustic_refinements as u128 + 1);
        // The declared vertex ceiling bounds indexed geometry and midpoint
        // scratch even before the source chart reveals its actual vertex count.
        let bytes = 48 * (options.source_surface.maximum_vertices as u128)
            + 48 * (final_facets as u128) * (m as u128 + 1)
            + 16 * (m as u128).pow(2)
            + 128 * (m as u128)
            + 1024;
        if setup_work > options.maximum_setup_work as u128
            || bytes > options.maximum_source_bytes as u128
        {
            return Err(bad(
                "open-screen complete source exceeds its setup work/payload allowance",
            ));
        }
        let chart = shell
            .open_screen_radiation(
                options.acoustic_refinements,
                options.source_surface,
                &mut cancelled,
            )
            .map_err(BemModesError::ScreenSource)?;
        let (vertices, triangles, weights) = chart.into_parts();
        poll(&mut cancelled)?;
        let surface =
            OpenPanels::new(vertices, triangles).map_err(BemModesError::ScreenBoundary)?;
        poll(&mut cancelled)?;
        let mut stiffness = vec![0.0; m * m];
        shell
            .hessian_matrix_with_cancel(
                &vec![0.0; m],
                &mut stiffness,
                &mut vec![[0.0; 3]; m],
                &mut cancelled,
            )
            .map_err(BemModesError::ScreenSource)?;
        poll(&mut cancelled)?;
        let viscosity = viscosity.map_or_else(|| vec![0.0; m * m], |c| c.matrix_per_s().to_vec());
        if stiffness.iter().chain(&viscosity).any(|v| !v.is_finite()) {
            return Err(bad(
                "open-screen original structural coefficients are not finite",
            ));
        }
        Ok(Self {
            surface,
            weights,
            stiffness,
            viscosity,
            memory: RetainedMemory::default(),
            medium,
        })
    }

    /// Attach the SAME displacement-basis Maxwell arms used by fs-phs.
    /// The source Hessian is the equilibrium stiffness: it must not already
    /// include these excess-stiffness arms. Kelvin viscosity can coexist when
    /// it represents a distinct material loss, not the same heat counted twice.
    ///
    /// Each signed c=sqrt(stiffness)*projection remains in the original basis.
    /// Equal-time and mechanically cancelling branches retain separate states.
    /// One complete cold attachment is allowed; empty input is a bitwise no-op.
    /// The caller owns matching source basis and constitutive applicability.
    pub fn with_relaxation_branches(
        mut self,
        branches: &[RelaxationBranch],
        limits: OpenScreenMemoryLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BemModesError> {
        poll(&mut cancelled)?;
        if limits.maximum_branches > 256 || branches.len() > limits.maximum_branches {
            return Err(bad("open-screen material memory exceeds its branch ceiling"));
        }
        if branches.is_empty() {
            return Ok(self);
        }
        if self.memory_count() != 0 {
            return Err(bad("open-screen memory requires one complete cold attachment"));
        }
        let m = self.mode_count() as u128;
        let r = branches.len() as u128;
        let bytes = 8 * r * m + 32 * r + 256;
        let work = 16 * r * (m + 1);
        if bytes > limits.maximum_payload_bytes as u128 || work > limits.maximum_work as u128 {
            return Err(bad("open-screen memory exceeds its source payload/work allowance"));
        }
        let memory = RetainedMemory::from_branches(
            self.mode_count(),
            branches,
            limits.maximum_branches,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        self.memory = memory;
        Ok(self)
    }

    pub fn mode_count(&self) -> usize {
        self.weights.len()
    }
    pub fn memory_count(&self) -> usize {
        self.memory.times.len()
    }
    /// Original signed c rows [1/s], in mass-normalized displacement order.
    pub fn relaxation_rows_per_s(&self) -> &[Vec<f64>] {
        &self.memory.rows
    }
    pub fn relaxation_times_s(&self) -> &[f64] {
        &self.memory.times
    }
    pub fn surface(&self) -> &OpenPanels {
        &self.surface
    }
    pub fn normal_velocity_weights(&self) -> &[Vec<[f64; 3]>] {
        &self.weights
    }
    pub fn stiffness_per_s2(&self) -> &[f64] {
        &self.stiffness
    }
    pub fn viscosity_per_s(&self) -> &[f64] {
        &self.viscosity
    }
    pub fn medium(&self) -> Medium {
        self.medium
    }

    /// Prepare the original open-screen owner for a declared positive band.
    /// Mesh, edge/triangle/angular resolution and all original work caps apply.
    /// Uniform source refinements do not establish the caller's edge-level claim.
    pub fn prepare(
        &self,
        highest_omega_rad_s: f64,
        options: OpenOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedOpenScreenModes<'_>, BemModesError> {
        poll(&mut cancelled)?;
        let screen = PreparedScreen::new(
            &self.surface,
            highest_omega_rad_s / self.medium.sound_speed,
            options,
            &mut cancelled,
        )
        .map_err(BemModesError::ScreenBoundary)?;
        Ok(PreparedOpenScreenModes {
            source: self,
            screen,
        })
    }
}

impl<'s> PreparedOpenScreenModes<'s> {
    pub fn source(&self) -> &OpenScreenLoadedModes {
        self.source
    }
    pub fn boundary(&self) -> &PreparedScreen<'s> {
        &self.screen
    }

    /// Solve every authored force against the COMPLETE self/mutual reaction at
    /// this positive REAL frequency. All unit-source BEM columns share their
    /// original assembly/factor; all force cases share one mechanical LU.
    /// Singular undamped systems refuse without artificial broadening. Raw
    /// reciprocity/passivity errors are neither symmetrized nor clipped.
    ///
    /// Fields are actual loaded structural radiation, not an observer of an
    /// unloaded motion. The result is a forced finite-dimensional response;
    /// source/mesh/triangle/angular convergence remain independent. No runtime,
    /// left-half-plane continuation or loaded-resonance claim follows from it.
    pub fn harmonic_response_batch<'p>(
        &'p self,
        omega: f64,
        forces: &[&[C64]],
        options: OpenScreenHarmonicOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<OpenScreenHarmonicBatch<'p, 's>, BemModesError> {
        poll(&mut cancelled)?;
        self.source.admit_harmonic_batch(omega, forces, options)?;
        let rows: Vec<_> = self.source.weights.iter().map(Vec::as_slice).collect();
        let radiation = self
            .screen
            .solve_velocity_ports(
                omega / self.source.medium.sound_speed,
                self.source.medium,
                &rows,
                options.boundary,
                &mut cancelled,
            )
            .map_err(BemModesError::ScreenBoundary)?;
        let responses = self.source.solve_harmonic_batch(
            omega,
            forces,
            options,
            radiation.impedance(),
            radiation.maximum_backward_residual(),
            |velocity| radiation.resisting_efforts(velocity),
            &mut cancelled,
        )?;
        Ok(OpenScreenHarmonicBatch {
            radiation,
            responses,
        })
    }
}

impl OpenScreenLoadedModes {
    /// Complete retained mechanical/material payload, also used by sound drive.
    fn harmonic_resources(&self, cases: usize) -> (u128, u128) {
        let m = self.mode_count() as u128;
        let d = m + self.memory_count() as u128;
        let c = cases as u128;
        (
            c * (16 * d + 32 * m + 640) + 96 * d * d + 320 * d + 4096,
            32 * d * d * d + 256 * (c + 1) * d * d,
        )
    }

    fn admit_harmonic_batch(
        &self,
        omega: f64,
        forces: &[&[C64]],
        options: OpenScreenHarmonicOptions,
    ) -> Result<(), BemModesError> {
        let m = self.mode_count();
        let cases = forces.len();
        let omega_sq = omega * omega;
        if !omega.is_finite()
            || omega <= 0.0
            || !omega_sq.is_finite()
            || omega_sq == 0.0
            || cases == 0
            || options.maximum_load_cases == 0
            || options.maximum_load_cases > 256
            || cases > options.maximum_load_cases
            || forces
                .iter()
                .any(|f| f.len() != m || f.iter().any(|&z| !finite(z)))
            || [
                options.physical_residual_tolerance,
                options.power_relative_tolerance,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0 || *v >= 1.0)
        {
            return Err(bad(
                "open-screen harmonic response requires a positive frequency, complete finite forces and valid limits",
            ));
        }
        let (bytes, work) = self.harmonic_resources(cases);
        if bytes > options.maximum_response_bytes as u128 || work > options.maximum_work as u128 {
            return Err(bad(
                "open-screen harmonic batch exceeds its response/work allowance",
            ));
        }
        Ok(())
    }

    /// The original full mechanical solve is shared by nondispersive and
    /// molecular exteriors. Only their source-owned boundary operator differs.
    #[allow(clippy::too_many_arguments)]
    fn solve_harmonic_batch(
        &self,
        omega: f64,
        forces: &[&[C64]],
        options: OpenScreenHarmonicOptions,
        impedance: &[C64],
        boundary_residual: f64,
        mut resisting_efforts: impl FnMut(&[C64]) -> Result<Vec<C64>, fs_bem::helmholtz::HelmholtzError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<OpenScreenHarmonicResponse>, BemModesError> {
        let m = self.mode_count();
        let d = m + self.memory_count();
        let omega_sq = omega * omega;
        let s = C64::new(0.0, -omega);
        let mut matrix = vec![C64::ZERO; d * d];
        let mut envelope = vec![0.0; d * d];
        for i in 0..m {
            poll(&mut cancelled)?;
            for j in 0..m {
                let at = i * m + j;
                let dst = i * d + j;
                matrix[dst] =
                    C64::new(self.stiffness[at], -omega * self.viscosity[at]) + s * impedance[at];
                envelope[dst] = self.stiffness[at].abs()
                    + omega * (self.viscosity[at].abs() + impedance[at].abs());
                if i == j {
                    matrix[dst] = matrix[dst] - C64::from_re(omega_sq);
                    envelope[dst] += omega_sq;
                }
            }
        }
        if self.memory_count() != 0 {
            // Keep the original elastic-memory equations, including dark and
            // equal-time histories. No dense C_memory^T*C_memory is formed.
            let mut derivative_scratch = vec![C64::ZERO; d * d];
            self.memory.add_to_pencil(
                C64::from_re(omega),
                m,
                0,
                &mut matrix,
                &mut derivative_scratch,
                &mut envelope,
                &mut cancelled,
            )?;
        }
        if matrix.iter().any(|&z| !finite(z)) || envelope.iter().any(|x| !x.is_finite()) {
            return Err(bad(
                "open-screen coupled mechanical/material operator is not representable",
            ));
        }
        // The same original-row balancing/audit used by closed-body loading
        // handles the different displacement and sqrt(J) memory coordinates.
        // No frequency derivative or loaded-pole authority is exposed here.
        let raw = PhysicalPencil {
            value: PoleEvaluation {
                matrix,
                derivative: Vec::new(),
                residual_scale: 1.0,
            },
            envelope,
        };
        let (rows, columns) = raw.fixed_scales(d, &mut cancelled)?;
        let mut balanced = raw.value.matrix.clone();
        for i in 0..d {
            poll(&mut cancelled)?;
            for j in 0..d {
                let at = i * d + j;
                balanced[at] = balanced[at].scale(columns[j] / rows[i]);
                if !finite(balanced[at]) {
                    return Err(bad("open-screen balanced material operator is not representable"));
                }
            }
        }
        let lu = lu_complex_with_cancel(&balanced, d, &mut cancelled)
            .map_err(BemModesError::ResponseSolve)?;
        let mut responses = Vec::with_capacity(forces.len());
        for &force in forces {
            poll(&mut cancelled)?;
            let mut state = vec![C64::ZERO; d];
            for i in 0..m {
                state[i] = divide(force[i], rows[i]);
            }
            if state.iter().any(|&z| !finite(z)) {
                return Err(bad("open-screen balanced force is not representable"));
            }
            lu.solve_with_cancel(&mut state, &mut cancelled)
                .map_err(BemModesError::ResponseSolve)?;
            for (value, &column) in state.iter_mut().zip(&columns) {
                *value = value.scale(column);
            }
            let velocity: Vec<_> = state[..m].iter().map(|&q| s * q).collect();
            if state.iter().chain(&velocity).any(|&z| !finite(z)) {
                return Err(bad("open-screen physical response is not finite"));
            }
            let residuals = raw.physical_residuals(
                m,
                0,
                &state,
                Some(force),
                options.physical_residual_tolerance,
                &mut cancelled,
            )?;
            let efforts = resisting_efforts(&velocity).map_err(BemModesError::ScreenBoundary)?;
            let (power, energy) = self.response_power(
                force,
                &state[..m],
                &velocity,
                &state[m..],
                &efforts,
                impedance,
                options.power_relative_tolerance,
                &mut cancelled,
            )?;
            responses.push(OpenScreenHarmonicResponse {
                angular_frequency_rad_s: omega,
                modal_displacement: state[..m].to_vec(),
                modal_velocity: velocity,
                elastic_memory_sqrt_j: state[m..].to_vec(),
                fluid_reaction: efforts.into_iter().map(|z| -z).collect(),
                mean_structural_energy_j: energy,
                mechanical_residual: residuals[0],
                memory_residual: residuals[2],
                maximum_boundary_residual: boundary_residual,
                power,
            });
        }
        poll(&mut cancelled)?;
        Ok(responses)
    }
}

impl OpenScreenLoadedModes {
    #[allow(clippy::too_many_arguments)]
    fn response_power(
        &self,
        force: &[C64],
        displacement: &[C64],
        velocity: &[C64],
        memory: &[C64],
        efforts: &[C64],
        impedance: &[C64],
        tolerance: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(BemHarmonicPower, f64), BemModesError> {
        // The original driven-BEM owner's common amplitude normalization avoids
        // squaring extreme authored forces before auditing physical work.
        let amplitude = force.iter().map(|z| z.abs()).fold(0.0, f64::max);
        let amplitude = if amplitude == 0.0 { 1.0 } else { amplitude };
        let normalize = |values: &[C64]| -> Result<Vec<C64>, BemModesError> {
            values
                .iter()
                .map(|&z| {
                    let v = divide(z, amplitude);
                    if !finite(v) || (z.re != 0.0 && v.re == 0.0) || (z.im != 0.0 && v.im == 0.0) {
                        Err(bad("open-screen work normalization is not representable"))
                    } else {
                        Ok(v)
                    }
                })
                .collect()
        };
        let (f, q, v, e, h) = (
            normalize(force)?,
            normalize(displacement)?,
            normalize(velocity)?,
            normalize(efforts)?,
            normalize(memory)?,
        );
        // Input, viscosity, screen work, absolute work, storage, relaxation.
        let mut sums = [0.0; 6];
        let mut corrections = [0.0; 6];
        let mut add = |channel: usize, value: f64| {
            let adjusted = value - corrections[channel];
            let next = sums[channel] + adjusted;
            corrections[channel] = (next - sums[channel]) - adjusted;
            sums[channel] = next;
        };
        let m = self.mode_count();
        for i in 0..m {
            poll(cancelled)?;
            add(0, 0.5 * (f[i] * v[i].conj()).re);
            add(2, 0.5 * (e[i] * v[i].conj()).re);
            add(3, 0.5 * f[i].abs() * v[i].abs());
            add(4, 0.25 * v[i].norm_sq());
            for j in 0..m {
                let at = i * m + j;
                add(1, 0.5 * self.viscosity[at] * (v[i].conj() * v[j]).re);
                add(
                    3,
                    0.5 * (self.viscosity[at].abs() + impedance[at].abs())
                        * v[i].abs()
                        * v[j].abs(),
                );
                add(4, 0.25 * self.stiffness[at] * (q[i].conj() * q[j]).re);
            }
        }
        for (&value, &tau) in h.iter().zip(&self.memory.times) {
            poll(cancelled)?;
            let norm_sq = value.norm_sq();
            let heat = 0.5 * norm_sq / tau;
            add(5, heat);
            add(3, heat);
            add(4, 0.25 * norm_sq);
        }
        let balance = sums[0] - sums[1] - sums[2] - sums[5];
        let relative = if sums[3] == 0.0 {
            balance.abs()
        } else {
            balance.abs() / sums[3]
        };
        if sums.iter().any(|x| !x.is_finite()) || !relative.is_finite() {
            return Err(bad(
                "open-screen work or structural storage audit overflows",
            ));
        }
        if relative > tolerance {
            return Err(BemModesError::PowerBalance {
                relative,
                tolerance,
            });
        }
        let restore = |x: f64| -> Result<f64, BemModesError> {
            let value = (x * amplitude) * amplitude;
            if !value.is_finite() || (x != 0.0 && value == 0.0) {
                Err(bad(
                    "open-screen physical work or structural storage is not representable",
                ))
            } else {
                Ok(value)
            }
        };
        Ok((
            BemHarmonicPower {
                input_w: restore(sums[0])?,
                kelvin_w: restore(sums[1])?,
                relaxation_w: restore(sums[5])?,
                radiated_w: restore(sums[2])?,
                balance_w: restore(balance)?,
                apparent_power_scale_w: restore(sums[3])?,
                relative_balance: relative,
            },
            restore(sums[4])?,
        ))
    }
}

impl<'p, 's> OpenScreenHarmonicBatch<'p, 's> {
    pub fn responses(&self) -> &[OpenScreenHarmonicResponse] {
        &self.responses
    }
    /// Original complete matrix/unit-source fields, without a second solve.
    pub fn radiation(&self) -> &ScreenVelocityPortSolution<'p, 's> {
        &self.radiation
    }
    fn rates(&self, case: usize) -> Result<&[C64], BemModesError> {
        self.responses
            .get(case)
            .map(|r| r.modal_velocity.as_slice())
            .ok_or_else(|| bad("open-screen load-case index is out of range"))
    }
    /// Actual coherent pressure jump [Pa] in original acoustic vertex order.
    pub fn pressure_jump(
        &self,
        case: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, BemModesError> {
        poll(&mut cancelled)?;
        self.radiation
            .pressure_jump(self.rates(case)?, cancelled)
            .map_err(BemModesError::ScreenBoundary)
    }
    /// Actual pressure [Pa] on either side, using the original adaptive P1
    /// double-layer integration and its receiver/clearance/total-work limits.
    pub fn pressure_at(
        &self,
        case: usize,
        points: &[[f64; 3]],
        separation_m: f64,
        options: near_field::Options,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, BemModesError> {
        poll(&mut cancelled)?;
        self.radiation
            .pressure_at(self.rates(case)?, points, separation_m, options, cancelled)
            .map_err(BemModesError::ScreenBoundary)
    }
    /// Coherent outgoing F [Pa m], p ~ F*exp(ikr)/r, on either side.
    pub fn far_field(
        &self,
        case: usize,
        direction: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<C64, BemModesError> {
        poll(&mut cancelled)?;
        self.radiation
            .far_field(self.rates(case)?, direction, cancelled)
            .map_err(BemModesError::ScreenBoundary)
    }
    /// Independent coherent angular power with the original selected angular
    /// rule, retaining all cross terms. Its spatial/angular convergence is not
    /// certified by the mechanical equation or pressure-work balance.
    pub fn radiated_power_w(
        &self,
        case: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, BemModesError> {
        poll(&mut cancelled)?;
        self.radiation
            .radiated_power_w(self.rates(case)?, cancelled)
            .map_err(BemModesError::ScreenBoundary)
    }
}

#[cfg(test)]
pub(super) mod tests;
