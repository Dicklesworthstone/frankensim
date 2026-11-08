//! Incident-source diffraction by the original rigid or admittance P1 screen.
//!
//! Incident pressure is continuous across this zero-thickness sheet. Its
//! normal velocity supplies -integral(phi * v_inc.n) to the original weak
//! solve, so the total field satisfies v_total.n = -Y*[p]. In particular,
//! there is no closed-wall -Y*p_inc term in this pressure-jump equation.
use super::*;
use crate::{ComplexSum, helmholtz::bc::IncidentField};

mod receivers;
pub use receivers::PreparedScreenScatteringReceivers;

/// Complete incident bank limits in addition to the original screen limits.
#[derive(Debug, Clone, Copy)]
pub struct ScreenScatteringLimits {
    /// Complete independent incident cases, in 1..=256.
    pub maximum_cases: usize,
    /// Every retained case/vertex pressure jump, including zero edge values.
    /// One additional active-unknown weak drive is temporary solve scratch.
    pub maximum_jump_entries: usize,
    /// All analytic source terms times (quadrature points + facets), including
    /// conservative admission of finite-emitter face/edge clearance work.
    /// The same ceiling bounds terms*receivers during receiver preparation.
    pub maximum_source_point_products: u64,
    /// Conservative scalar-work units for one original assembly/LU and every
    /// incident sample, source solve and residual. Not a wall-clock bound.
    /// Caller-owned nonlocal admittance applications remain indivisible.
    pub maximum_work: u64,
    /// Positive exact distance [m] from each finite source to every retained
    /// triangle face/edge. Sources on either side of an open screen are valid.
    pub minimum_source_clearance_m: f64,
}

/// Matched incident fields and their full-wave scattered pressure jumps.
/// The source declarations and original geometry stay borrowed. One original
/// operator/LU serves all cases; no case or coherent source term is truncated.
pub struct ScreenScattering<'p, 's, 'i> {
    solutions: Vec<JumpSolution<'p, 's>>,
    incidents: Vec<IncidentField<'i>>,
    maximum_source_point_products: u64,
    frequency_work: u64,
}

impl<'s> PreparedScreen<'s> {
    /// Scatter one coherent incident field from a stationary screen with the
    /// original sheet admittance law. The returned bank contains one case.
    pub fn solve_scattering<'p, 'i>(
        &'p self,
        k: f64,
        medium: Medium,
        admittance: Admittance<'_>,
        incident: IncidentField<'i>,
        limits: ScreenScatteringLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ScreenScattering<'p, 's, 'i>, HelmholtzError> {
        self.solve_scattering_many(
            k,
            medium,
            admittance,
            std::slice::from_ref(&incident),
            limits,
            cancelled,
        )
    }

    /// Scatter every incident case using ONE original H and LU. Incident
    /// normal velocity is integrated at the retained triangle quadrature
    /// points, not centroids or free-edge-clamped velocity nodes. Close or
    /// rapidly varying sources require a separate source/mesh/quadrature
    /// refinement; the backward residual alone does not establish that error.
    ///
    /// The unknown remains the total pressure jump, since the incident field
    /// has no jump. Rigid and local/nonlocal sheet laws use the same admittance
    /// operator and passive-work checks as the original motion-driven solve.
    /// Two independently coated faces still require a second trace.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_scattering_many<'p, 'i>(
        &'p self,
        k: f64,
        medium: Medium,
        admittance: Admittance<'_>,
        incidents: &[IncidentField<'i>],
        limits: ScreenScatteringLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScreenScattering<'p, 's, 'i>, HelmholtzError> {
        poll(&mut cancelled)?;
        let rho_c = medium.density * medium.sound_speed;
        let rho_omega = k * medium.sound_speed * medium.density;
        if !k.is_finite()
            || k <= 0.0
            || k > self.highest_k
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
            || !rho_c.is_finite()
            || rho_c <= 0.0
            || !rho_omega.is_finite()
            || rho_omega <= 0.0
            || incidents.is_empty()
            || incidents.len() > limits.maximum_cases
            || limits.maximum_cases == 0
            || limits.maximum_cases > crate::helmholtz::MAX_RADIATION_FIELDS_PER_BATCH
            || !limits.minimum_source_clearance_m.is_finite()
            || limits.minimum_source_clearance_m <= 0.0
            || limits.maximum_source_point_products == 0
        {
            return Err(bad(
                "open screen scattering requires a resolved band, medium, cases and clearance",
            ));
        }
        admittance.validate(self.surface.vertices.len())?;
        for incident in incidents {
            poll(&mut cancelled)?;
            incident.validate()?;
        }
        let cases = incidents.len() as u128;
        let terms: u128 = incidents
            .iter()
            .map(|field| field.source_count() as u128)
            .sum();
        let samples: u128 = self.samples.iter().map(|row| row.len() as u128).sum();
        let n = self.n as u128;
        let vertices = self.surface.vertices.len() as u128;
        let source_work = terms * (samples + self.facets.len() as u128);
        let frequency_work =
            self.frequency_work() + 256 * source_work + 128 * cases * n * n + 64 * cases * vertices;
        if cases * vertices > limits.maximum_jump_entries as u128
            || source_work > limits.maximum_source_point_products as u128
            || frequency_work > limits.maximum_work as u128
            || frequency_work > self.options.maximum_work as u128
        {
            return Err(bad(
                "complete open screen scattering bank exceeds its work or jump limit",
            ));
        }
        let frequency_work = u64::try_from(frequency_work)
            .map_err(|_| bad("open screen scattering work overflows"))?;

        // Winding or an enclosing sphere would incorrectly exclude valid
        // emitters on an open screen. Retain its actual faces and free edges.
        let mut scale = 1.0_f64;
        for (vertex, point) in self.surface.vertices.iter().enumerate() {
            if vertex % 256 == 0 {
                poll(&mut cancelled)?;
            }
            for value in point {
                scale = scale.max(value.abs());
            }
        }
        for incident in incidents {
            for point in incident.source_positions() {
                for value in point {
                    scale = scale.max(value.abs());
                }
                if limits.minimum_source_clearance_m <= 64.0 * f64::EPSILON * scale {
                    return Err(bad(
                        "open screen source clearance is below coordinate resolution",
                    ));
                }
                for (index, facet) in self.facets.iter().enumerate() {
                    if index % 32 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let distance = crate::near_field::distance(point, facet.triangle);
                    if !distance.is_finite() || distance <= limits.minimum_source_clearance_m {
                        return Err(bad("open screen finite source touches a face or edge"));
                    }
                }
            }
        }
        let (h, lu) = self.factor_wave(k, true, admittance, &mut cancelled)?;
        let mut solutions = Vec::with_capacity(incidents.len());
        for incident in incidents {
            poll(&mut cancelled)?;
            let mut weak = vec![ComplexSum::default(); self.n];
            for (facet, samples) in self.facets.iter().zip(&self.samples) {
                poll(&mut cancelled)?;
                for sample in samples {
                    let field = incident.sample_wave(
                        k,
                        rho_omega,
                        Some(rho_c),
                        sample.point,
                        &mut cancelled,
                    )?;
                    let mut normal = ComplexSum::default();
                    for axis in 0..3 {
                        normal.add(field.velocity_m_s[axis].scale(facet.normal[axis]));
                    }
                    let drive = -normal.finish();
                    for a in 0..3 {
                        if let Some(u) = self.dof[facet.ids[a]] {
                            weak[u].add(drive.scale(sample.weight * sample.shape[a]));
                        }
                    }
                }
            }
            let drive: Vec<_> = weak.into_iter().map(ComplexSum::finish).collect();
            if drive.iter().any(|&value| !finite(value)) {
                return Err(bad("open screen incident weak drive overflows"));
            }
            solutions.push(self.finish_weak_drive(
                k,
                true,
                medium,
                admittance,
                &drive,
                &h,
                &lu,
                &mut cancelled,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(ScreenScattering {
            solutions,
            incidents: incidents.to_vec(),
            maximum_source_point_products: limits.maximum_source_point_products,
            frequency_work,
        })
    }
}

impl<'p, 's> ScreenScattering<'p, 's, '_> {
    pub fn cases(&self) -> usize {
        self.solutions.len()
    }

    /// The original scattered field for one case: pressure jump, scattered
    /// far field, absorption and quadrature remain available on JumpSolution.
    /// Its driving power is work of the blocked incident forcing, not an
    /// applied mechanical drive or the infinite incident plane-wave power.
    pub fn scattered(&self, case: usize) -> Result<&JumpSolution<'p, 's>, HelmholtzError> {
        self.solutions
            .get(case)
            .ok_or_else(|| bad("open screen scattering case is out of range"))
    }

    pub fn frequency_work(&self) -> u64 {
        self.frequency_work
    }

    /// Algebraic residual only; no source/geometry/quadrature error bound.
    pub fn maximum_backward_residual(&self) -> f64 {
        self.solutions
            .iter()
            .map(|s| s.backward_residual)
            .fold(0.0, f64::max)
    }

    /// Net Cartesian force of the fluid on the stationary screen [N],
    /// -integral([p] n dS). The analytic incident field has no pressure jump;
    /// the original solved P1 jump therefore supplies the complete load.
    pub fn fluid_force_n(
        &self,
        case: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[C64; 3], HelmholtzError> {
        poll(&mut cancelled)?;
        let source = self.scattered(case)?;
        let mut force = [ComplexSum::default(); 3];
        for facet in &source.prepared.facets {
            poll(&mut cancelled)?;
            let integral = facet
                .ids
                .into_iter()
                .fold(C64::ZERO, |sum, v| sum + source.jump[v])
                .scale(-facet.area / 3.0);
            for axis in 0..3 {
                force[axis].add(integral.scale(facet.normal[axis]));
            }
        }
        let force = force.map(ComplexSum::finish);
        if force.iter().any(|&value| !finite(value)) {
            return Err(bad("open screen fluid force overflows"));
        }
        poll(&mut cancelled)?;
        Ok(force)
    }
}
