//! Thermoviscous propagation in the actual concentric annulus.
//!
//! Stinson's arbitrary-section low-reduced-frequency model (JASA 89 (1991),
//! 550–558, equations 59 and 65–69, <https://doi.org/10.1121/1.400379>) solves
//! `(s - D Laplacian) h = s`, h=0 on BOTH walls, with s=i*omega in telegraph
//! bookkeeping. The radial Laplacian is `(r h_r)_r/r`; the cross-section mean
//! is `integral(r*h dr)/integral(r dr)`. Pressure is uniform transversely, but
//! neither axial velocity nor excess temperature is replaced by its mean
//! before this radial problem is solved.
//!
//! Exact weighted P1 mass/stiffness integrals on a declared logarithmic mesh
//! retain inner curvature and both no-slip/isothermal walls. The original
//! complex tridiagonal owner solves every radial coordinate in O(N) work and
//! storage. The fluid supplies nu and thermal diffusivity separately; no
//! hydraulic radius, loss fit, gas surrogate or geometric series truncation
//! replaces the annulus. The native TwoPort supplies stable axial propagation.
//! [`molecular`] uses the original gas's complex heat capacity in that SAME
//! radial thermal solve. Its wall diffusion and molecular storage share one
//! energy equation; no separately added attenuation replaces their coupling.
//!
//! This is fixed-geometry, linear, fully developed, Newtonian continuum
//! acoustics. End fields, eccentricity, slip, mean flow, moving walls, bulk
//! viscosity and axial heat diffusion are outside this owner. The mesh and
//! diffusion-phase ceilings are explicit resolution choices. Algebraic/work
//! residuals are NOT bounds on radial discretization or 3-D model error.

use super::{FluidDuctError, bad, positive};
use crate::{DuctError, SegmentWave, compose::TwoPort};
use fs_la::eigen_complex::solve_tridiagonal_complex;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_math::{c64::C64, dd::Dd, det};

/// Molecular heat storage solved together with the original annular wall diffusion.
pub mod molecular;

/// Bounds the geometry, each native linear-solve tile and all frequency work.
pub const MAX_ANNULAR_ELEMENTS: usize = 4096;

/// Explicit numerical and plane-pressure limits; no hidden accuracy default.
#[derive(Clone, Copy, Debug)]
pub struct AnnularAdmission {
    /// max(omega/c, |k|)*outer_radius must not exceed this value in (0,1).
    /// The first nonaxisymmetric Neumann eigenvalue is at least 1/outer_radius
    /// by its Rayleigh quotient; this is a conservative transverse envelope.
    pub maximum_transverse_phase: f64,
    /// max radial cell width times sqrt(|z|), for each active diffusion argument
    /// z [1/m2]. For a constant heat capacity, |z| = omega/D.
    /// In (0,1]; a resolution ceiling, not a discretization-error certificate.
    pub maximum_cell_diffusion_phase: f64,
    /// Componentwise solve backward error and relative diffusion-work defect,
    /// in (0,1). Refinement of the radial mesh is a separate caller-owned axis.
    pub maximum_relative_residual: f64,
    /// Logical peak payload including retained geometry and native solve
    /// scratch: conservatively 192*(elements+1)+4096 bytes. Excludes allocator
    /// metadata, stacks and any other owners retained by the caller.
    pub maximum_payload_bytes: usize,
}

/// Mean radial diffusion in +i omega telegraph bookkeeping.
#[derive(Clone, Copy, Debug)]
pub struct AnnularDiffusion {
    /// Mean h, with nonnegative real/imaginary parts for passive diffusion.
    pub mean: C64,
    /// Mean 1-h, calculated through the well-scaled member of the two fields.
    pub wall: C64,
    /// Mean |grad h|² * Im(z)/|z|², independently integrated from the P1 field.
    /// Equals Im(mean) for the discrete diffusion equations. For constant
    /// diffusivity z=i*omega/D this is mean |grad h|² * D/omega. With molecular
    /// heat storage it remains a radial diffusion identity, not total absorbed heat.
    pub gradient_work: f64,
    pub maximum_cell_diffusion_phase: f64,
    pub relative_backward_residual: f64,
    pub relative_work_defect: f64,
}

/// Actual-area wave in the native exp(-i omega t) convention.
#[derive(Clone, Copy, Debug)]
pub struct AnnularWave {
    /// Pressure/mean-velocity and pressure/volume-flow impedance. The shear
    /// diagnostic uses HALF the radial gap, not an equivalent cylinder radius.
    pub wave: SegmentWave,
    pub viscous: AnnularDiffusion,
    /// Absent exactly when Cp=Cv; uncoupled thermal data do no radial work.
    pub thermal: Option<AnnularDiffusion>,
    pub transverse_phase: f64,
    area_m2: f64,
}
impl AnnularWave {
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Original stable finite-length pressure/volume-flow block. This carries
    /// no abrupt area-change, termination or three-dimensional junction law.
    pub fn two_port(&self, length_m: f64, reference: f64) -> Result<TwoPort, DuctError> {
        TwoPort::from_uniform_wave(&self.wave, self.area_m2, length_m, reference)
    }
}

/// Prepared, area-normalized radial P1 geometry with both material walls.
/// A zero inner radius is a cylinder axis, and belongs to the existing circular
/// owner rather than a Dirichlet wall. All original radial cells are retained.
#[derive(Clone, Debug)]
pub struct AnnularSection {
    inner_m: f64,
    outer_m: f64,
    gap_m: f64,
    area_m2: f64,
    fractions: Vec<f64>,
    stiffness_edges: Vec<f64>,
    mass_edges: Vec<f64>,
    mass_diagonal: Vec<f64>,
    load: Vec<f64>,
    boundary_weights: [f64; 2],
    maximum_cell_fraction: f64,
}

impl AnnularSection {
    pub fn new(inner_m: f64, outer_m: f64, elements: usize) -> Result<Self, DuctError> {
        Self::new_controlled(inner_m, outer_m, elements, || false)
    }

    /// Logarithmic radial spacing resolves the inner cylinder. Nearly equal
    /// radii use expm1 and a compensated log ratio to retain the actual gap.
    /// Cancellation is polled before allocation and every 64 radial cells.
    pub fn new_controlled(
        inner_m: f64,
        outer_m: f64,
        elements: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        poll(&mut cancelled, 0)?;
        if !positive(inner_m)
            || !positive(outer_m)
            || outer_m <= inner_m
            || !(2..=MAX_ANNULAR_ELEMENTS).contains(&elements)
        {
            return Err(bad(
                "annular section needs positive ordered radii and 2..=4096 cells",
            ));
        }
        let gap_m = outer_m - inner_m;
        let mean_radius = inner_m + 0.5 * gap_m;
        let area_m2 = core::f64::consts::TAU * mean_radius * gap_m;
        let ratio = gap_m / inner_m;
        if !positive(area_m2) || !positive(ratio) {
            return Err(bad("annular area or radius ratio is unrepresentable"));
        }
        let log_ratio = if ratio < 0.125 {
            // log(1+x) = 2 atanh(x/(2+x)); fixed terms on this short interval.
            let z = ratio / (2.0 + ratio);
            let mut term = z;
            let mut value = z;
            for n in 1..12 {
                term *= z * z;
                value += term / f64::from(2 * n + 1);
            }
            2.0 * value
        } else {
            det::ln(1.0 + ratio)
        };
        let mut fractions = filled(elements + 1, 0.0)?;
        for (index, value) in fractions.iter_mut().enumerate().skip(1).take(elements - 1) {
            if index % 64 == 0 {
                poll(&mut cancelled, index)?;
            }
            *value = det::expm1(log_ratio * index as f64 / elements as f64) / ratio;
        }
        fractions[elements] = 1.0;
        let mut stiffness_edges = filled(elements, 0.0)?;
        let mut mass_edges = filled(elements, 0.0)?;
        let mut mass_diagonal = filled(elements - 1, 0.0)?;
        let mut load = filled(elements - 1, 0.0)?;
        let mut boundary_weights = [0.0; 2];
        let mut maximum_cell_fraction = 0.0_f64;
        for i in 0..elements {
            if i % 64 == 0 {
                poll(&mut cancelled, i)?;
            }
            let dx = fractions[i + 1] - fractions[i];
            let left = inner_m / mean_radius + (gap_m / mean_radius) * fractions[i];
            let right = inner_m / mean_radius + (gap_m / mean_radius) * fractions[i + 1];
            let k = 0.5 * (left + right) / dx;
            let m = dx * ((left + right) / 12.0);
            let diagonal = [
                dx * ((3.0 * left + right) / 12.0),
                dx * ((left + 3.0 * right) / 12.0),
            ];
            if ![dx, left, right, k, m, diagonal[0], diagonal[1]]
                .into_iter()
                .all(positive)
            {
                return Err(bad("annular radial mesh is unresolved or unrepresentable"));
            }
            maximum_cell_fraction = maximum_cell_fraction.max(dx);
            stiffness_edges[i] = k;
            mass_edges[i] = m;
            for (side, diagonal) in diagonal.into_iter().enumerate() {
                let node = i + side;
                if node == 0 || node == elements {
                    boundary_weights[usize::from(node == elements)] = diagonal + m;
                } else {
                    mass_diagonal[node - 1] += diagonal;
                    load[node - 1] += diagonal + m;
                }
            }
        }
        poll(&mut cancelled, elements)?;
        Ok(Self {
            inner_m,
            outer_m,
            gap_m,
            area_m2,
            fractions,
            stiffness_edges,
            mass_edges,
            mass_diagonal,
            load,
            boundary_weights,
            maximum_cell_fraction,
        })
    }

    pub const fn inner_radius_m(&self) -> f64 {
        self.inner_m
    }
    pub const fn outer_radius_m(&self) -> f64 {
        self.outer_m
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    pub fn elements(&self) -> usize {
        self.stiffness_edges.len()
    }
    pub fn radius_m(&self, node: usize) -> Option<f64> {
        self.fractions
            .get(node)
            .map(|x| self.inner_m + self.gap_m * x)
    }

    fn admit(&self, omega: f64, admission: AnnularAdmission) -> Result<(), DuctError> {
        if !positive(omega)
            || !positive(admission.maximum_transverse_phase)
            || admission.maximum_transverse_phase >= 1.0
            || !positive(admission.maximum_cell_diffusion_phase)
            || admission.maximum_cell_diffusion_phase > 1.0
            || !positive(admission.maximum_relative_residual)
            || admission.maximum_relative_residual >= 1.0
            || 192 * (self.elements() + 1) + 4096 > admission.maximum_payload_bytes
        {
            return Err(bad(
                "annular frequency, numerical limits or payload budget is invalid",
            ));
        }
        Ok(())
    }

    pub fn evaluate<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        admission: AnnularAdmission,
    ) -> Result<AnnularWave, FluidDuctError> {
        self.evaluate_controlled(material, omega, admission, || false)
    }

    /// Query one coherent original fluid state and solve each active diffusion
    /// field once. The native O(N) solve is one bounded tile of at most 4095
    /// unknowns, with cancellation polls immediately before and after it.
    pub fn evaluate_controlled<M: ThermoviscousFluid + ?Sized>(
        &self,
        material: &M,
        omega: f64,
        admission: AnnularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AnnularWave, FluidDuctError> {
        poll(&mut cancelled, 0)?;
        self.admit(omega, admission)?;
        let state = material.thermoviscous_state()?;
        let acoustic = state.acoustic();
        let transport = state.transport();
        let k0 = omega / acoustic.sound_speed_m_s();
        if !positive(k0) || k0 * self.outer_m > admission.maximum_transverse_phase {
            return Err(bad("annulus is outside its plane-pressure frequency envelope").into());
        }
        let viscous = self.diffusion(
            omega,
            transport.kinematic_viscosity_m2_s(),
            admission,
            &mut cancelled,
        )?;
        let gamma = state.heat_capacity_ratio();
        let thermal = if gamma == 1.0 {
            None
        } else {
            Some(self.diffusion(
                omega,
                transport.thermal_diffusivity_m2_s(),
                admission,
                &mut cancelled,
            )?)
        };
        let s = C64::new(0.0, omega);
        let series = s.scale(acoustic.density_kg_m3() / self.area_m2) / viscous.mean;
        let factor = C64::ONE + thermal.map_or(C64::ZERO, |t| t.wall.scale(gamma - 1.0));
        let shunt = s.scale(self.area_m2 / acoustic.adiabatic_bulk_modulus_pa()) * factor;
        let (wavenumber, characteristic_impedance) =
            crate::helmholtz_wave_from_telegraph(series, shunt);
        let specific_impedance = characteristic_impedance.scale(self.area_m2);
        let transverse_phase = k0.max(wavenumber.abs()) * self.outer_m;
        if ![wavenumber, characteristic_impedance, specific_impedance]
            .into_iter()
            .all(finite)
            || wavenumber.re <= 0.0
            || wavenumber.im < 0.0
            || characteristic_impedance.re <= 0.0
            || !transverse_phase.is_finite()
            || transverse_phase > admission.maximum_transverse_phase
        {
            return Err(bad(
                "annular wave is unrepresentable or outside its transverse-phase envelope",
            )
            .into());
        }
        poll(&mut cancelled, self.elements())?;
        Ok(AnnularWave {
            wave: SegmentWave {
                wavenumber,
                characteristic_impedance,
                specific_impedance,
                shear_number: 0.5
                    * self.gap_m
                    * det::sqrt(omega / transport.kinematic_viscosity_m2_s()),
            },
            viscous,
            thermal,
            transverse_phase,
            area_m2: self.area_m2,
        })
    }

    fn diffusion(
        &self,
        omega: f64,
        diffusivity: f64,
        admission: AnnularAdmission,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AnnularDiffusion, DuctError> {
        let phase = self.gap_m / det::sqrt(diffusivity / omega);
        let alpha = phase * phase;
        self.diffusion_scaled(
            C64::new(0.0, alpha),
            phase,
            1.0 / phase,
            admission,
            cancelled,
        )
    }

    // One original radial solve for real transport and complex caloric storage.
    // z is dimensionless (physical argument * gap²); phase=sqrt(|z|).
    // Im(<h>) = Im(z)/|z|² * <|grad_x h|²>, x=(r-inner)/gap,
    // follows by taking the imaginary part of h* (K/z + M) h = h* b.
    fn diffusion_scaled(
        &self,
        argument: C64,
        phase: f64,
        gradient_scale: f64,
        admission: AnnularAdmission,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AnnularDiffusion, DuctError> {
        poll(cancelled, 0)?;
        let magnitude = argument.abs();
        let maximum_cell_diffusion_phase = phase * self.maximum_cell_fraction;
        if !finite(argument)
            || argument.re < 0.0
            || !positive(argument.im)
            || !positive(magnitude)
            || !positive(phase)
            || !positive(gradient_scale)
            || !maximum_cell_diffusion_phase.is_finite()
            || maximum_cell_diffusion_phase > admission.maximum_cell_diffusion_phase
        {
            return Err(bad(
                "annular radial mesh does not resolve the active diffusion layer",
            ));
        }
        // Solve h in a narrow/low-frequency gap, and w=1-h in a wide/high-
        // frequency gap. The latter retains small wall fractions without
        // subtracting a nearly uniform interior h from one after averaging.
        let solve_wall = magnitude > 1.0;
        let (k_scale, m_scale) = if solve_wall {
            (
                1.0 / magnitude,
                C64::new(argument.re / magnitude, argument.im / magnitude),
            )
        } else {
            (1.0, argument)
        };
        let n = self.load.len();
        let mut diagonal = filled(n, C64::ZERO)?;
        let mut off = filled(n - 1, C64::ZERO)?;
        let mut rhs = filled(n, C64::ZERO)?;
        let boundary = [0, self.elements() - 1].map(|i| {
            C64::new(
                -k_scale * self.stiffness_edges[i] + m_scale.re * self.mass_edges[i],
                m_scale.im * self.mass_edges[i],
            )
        });
        for i in 0..n {
            if i % 64 == 0 {
                poll(cancelled, i)?;
            }
            diagonal[i] = C64::new(
                k_scale * (self.stiffness_edges[i] + self.stiffness_edges[i + 1])
                    + m_scale.re * self.mass_diagonal[i],
                m_scale.im * self.mass_diagonal[i],
            );
            if i + 1 < n {
                off[i] = C64::new(
                    -k_scale * self.stiffness_edges[i + 1] + m_scale.re * self.mass_edges[i + 1],
                    m_scale.im * self.mass_edges[i + 1],
                );
            }
            if !solve_wall {
                rhs[i] = m_scale.scale(self.load[i]);
            }
        }
        if solve_wall {
            rhs[0] = rhs[0] - boundary[0];
            rhs[n - 1] = rhs[n - 1] - boundary[1];
        }
        if diagonal.iter().chain(&off).chain(&rhs).any(|&z| !finite(z)) {
            return Err(bad("annular diffusion matrix is unrepresentable"));
        }
        poll(cancelled, n)?;
        solve_tridiagonal_complex(&off, &diagonal, &off, &mut rhs, 1)
            .map_err(|_| DuctError::Singular)?;
        poll(cancelled, n)?;
        let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
        let mut relative_backward_residual = 0.0_f64;
        for i in 0..n {
            if i % 64 == 0 {
                poll(cancelled, i)?;
            }
            if !finite(rhs[i]) {
                return Err(bad("annular diffusion solution is unrepresentable"));
            }
            re = re + Dd::from_f64(self.load[i]) * Dd::from_f64(rhs[i].re);
            im = im + Dd::from_f64(self.load[i]) * Dd::from_f64(rhs[i].im);
            let mut drive = if solve_wall {
                C64::ZERO
            } else {
                m_scale.scale(self.load[i])
            };
            if solve_wall && i == 0 {
                drive = drive - boundary[0];
            }
            if solve_wall && i == n - 1 {
                drive = drive - boundary[1];
            }
            let mut residual = diagonal[i] * rhs[i] - drive;
            let mut scale = diagonal[i].abs() * rhs[i].abs() + drive.abs();
            if i > 0 {
                residual = residual + off[i - 1] * rhs[i - 1];
                scale += off[i - 1].abs() * rhs[i - 1].abs();
            }
            if i + 1 < n {
                residual = residual + off[i] * rhs[i + 1];
                scale += off[i].abs() * rhs[i + 1].abs();
            }
            relative_backward_residual =
                relative_backward_residual.max(relative(residual.abs(), scale));
        }
        if solve_wall {
            re = re
                + Dd::from_f64(self.boundary_weights[0])
                + Dd::from_f64(self.boundary_weights[1]);
        }
        let average = C64::new(re.to_f64(), im.to_f64());
        let complement = C64::new((Dd::from_f64(1.0) - re).to_f64(), -im.to_f64());
        let (mean, wall) = if solve_wall {
            (complement, average)
        } else {
            (average, complement)
        };
        let boundary_value = if solve_wall { C64::ONE } else { C64::ZERO };
        let mut previous = boundary_value;
        let mut work = Dd::ZERO;
        for i in 0..self.elements() {
            if i % 64 == 0 {
                poll(cancelled, i)?;
            }
            let next = if i < n { rhs[i] } else { boundary_value };
            // Scale before squaring: weak diffusion work must not disappear
            // solely because the unscaled gradient square underflows.
            let scaled = (next - previous).scale(gradient_scale);
            work = work + Dd::from_f64(self.stiffness_edges[i]) * Dd::from_f64(scaled.norm_sq());
            previous = next;
        }
        let gradient_work = work.to_f64();
        let relative_work_defect = relative(
            (gradient_work - mean.im).abs(),
            gradient_work.abs() + mean.im.abs(),
        );
        if !finite(mean)
            || !finite(wall)
            || mean.re < 0.0
            || mean.im < 0.0
            || wall.re <= 0.0
            || !gradient_work.is_finite()
            || gradient_work < 0.0
            || relative_backward_residual > admission.maximum_relative_residual
            || relative_work_defect > admission.maximum_relative_residual
        {
            return Err(bad(
                "annular diffusion solve or physical work residual exceeds its admission",
            ));
        }
        poll(cancelled, n)?;
        Ok(AnnularDiffusion {
            mean,
            wall,
            gradient_work,
            maximum_cell_diffusion_phase,
            relative_backward_residual,
            relative_work_defect,
        })
    }
}

fn filled<T: Clone>(n: usize, value: T) -> Result<Vec<T>, DuctError> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| bad("annular radial allocation exceeds available memory"))?;
    v.resize(n, value);
    Ok(v)
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn relative(error: f64, scale: f64) -> f64 {
    if !error.is_finite() || !scale.is_finite() {
        f64::INFINITY
    } else if scale > 0.0 {
        error / scale
    } else if error == 0.0 {
        0.0
    } else {
        f64::INFINITY
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool, completed_steps: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps })
    } else {
        Ok(())
    }
}
