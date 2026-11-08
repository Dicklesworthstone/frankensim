//! Prescribed finite sources radiating from the bounded fluid through its
//! original transmitting interface. All phasors use exp(-i omega t).
//!
//! The native direct transmission matrix M=I+C and its factor are unchanged.
//! An interior source adds p_s and q_s to the interior Green identities, so
//! M [p, q_o/k_o] = [p_s, q_s/k_o], q_s=i*omega*rho_i*v_s.
//! For the reference x_s=[p_s,i*Z_o*v_s], solve the correction directly:
//! M delta = -C x_s + [0, i*Z_o*(rho_i-rho_o)/rho_o*v_s].
//! This preserves weak reflected traces without total-minus-source subtraction.
//! The interior reflected field is S_i delta_q_i-D_i delta_p; the exterior
//! outgoing field is D_o p-S_o q_o. Analytic source fields are added only in
//! their original interior medium. Exactly matched media have zero reflection
//! and retain the analytic source field in both observation regions.
//!
//! These are prescribed compact volume/force/stress sources in two stationary,
//! inviscid, nondispersive fluids. Source feedback, motion, viscosity, finite
//! source size and nonlinear interface deformation are outside this owner.
//! The original P0 BEM mesh and source/receiver quadrature need independent
//! refinement. The algebraic residual is not a continuum-error certificate.

use super::*;
use crate::ComplexSum;
use crate::helmholtz::bc::IncidentSource;

/// Source geometry and original corrected-equation admission.
#[derive(Clone, Copy, Debug)]
pub struct InternalSourceControls {
    /// Exact minimum point-to-triangle distance [m], including face edges.
    pub minimum_source_clearance_m: f64,
    /// Maximum componentwise backward error of the native corrected system.
    pub maximum_relative_residual: f64,
}

/// Caller-declared physical fluid containing the receivers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InternalSourceRegion {
    Interior,
    Exterior,
}

impl<'s> FluidInterface<'s> {
    /// Radiate a coherent finite source array lying strictly INSIDE this
    /// original closed surface. Plane waves are not bounded interior sources.
    /// Each source uses the original interior density and sound speed; each
    /// exterior Green kernel uses the original exterior state at the same omega.
    ///
    /// Reuses the complete native factor and its dense/vector allowance. The
    /// native source/receiver cardinality and exact geometric work caps apply.
    /// All panel sampling, matrix rows, LU and geometry tiles are cancellable.
    /// No source material is re-queried and no bare interior inverse is formed.
    pub fn radiate_internal_field<'p>(
        &'p self,
        sources: IncidentField<'p>,
        controls: InternalSourceControls,
        gate: &CancelGate,
    ) -> Result<InternalSourceRadiation<'p, 's>, HelmholtzError> {
        poll(gate)?;
        sources.validate()?;
        if !controls.maximum_relative_residual.is_finite()
            || controls.maximum_relative_residual <= 0.0
            || controls.maximum_relative_residual >= 1.0
        {
            return Err(bad("internal source residual ceiling must lie in (0,1)"));
        }
        if matches!(sources, IncidentField::PlaneWave(_))
            || matches!(sources, IncidentField::Sources(terms)
                if terms.iter().any(|s| matches!(s, IncidentSource::PlaneWave(_))))
        {
            return Err(bad(
                "internal radiation requires only bounded finite sources",
            ));
        }
        near_field::Geometry::new_interior_cancellable(
            self.surface,
            &sources.source_positions(),
            controls.minimum_source_clearance_m,
            gate,
        )?;
        let n = self.surface.areas().len();
        let dimension = 2 * n;
        let z_o = self.exterior.characteristic_impedance_pa_s_m();
        let interior = medium(self.interior);
        let mut reference = vec![C64::ZERO; dimension];
        for i in 0..n {
            poll(gate)?;
            let sample = sources.sample_wave(
                self.k_i,
                self.omega * interior.density,
                Some(interior.density * interior.sound_speed),
                self.surface.centroids()[i],
                &mut || gate.is_requested(),
            )?;
            let mut normal = ComplexSum::default();
            for axis in 0..3 {
                normal.add(sample.velocity_m_s[axis].scale(self.surface.normals()[i][axis]));
            }
            reference[i] = sample.pressure_pa;
            reference[n + i] = C64::new(0.0, z_o) * normal.finish();
        }
        if reference.iter().any(|&z| !finite(z)) {
            return Err(bad("internal source trace scaling is unrepresentable"));
        }
        let mut correction = vec![C64::ZERO; dimension];
        let mut relative_backward_residual = 0.0_f64;
        if let Some(factor) = &self.factor {
            let density_jump = (self.interior.density_kg_m3() - self.exterior.density_kg_m3())
                / self.exterior.density_kg_m3();
            for i in 0..dimension {
                poll(gate)?;
                let mut drive = ComplexSum::default();
                if i >= n {
                    drive.add(reference[i].scale(density_jump));
                }
                for (j, &value) in reference.iter().enumerate() {
                    drive.add(-(self.contrast[i * dimension + j] * value));
                }
                correction[i] = drive.finish();
            }
            if correction.iter().any(|&z| !finite(z)) {
                return Err(bad("internal transmission forcing overflow"));
            }
            let drive = correction.clone();
            factor
                .solve_with_cancel(&mut correction, || gate.is_requested())
                .map_err(lu_error)?;
            for i in 0..dimension {
                poll(gate)?;
                let mut residual = ComplexSum::default();
                residual.add(correction[i]);
                residual.add(-drive[i]);
                let mut scale = correction[i].abs() + drive[i].abs();
                for (j, &value) in correction.iter().enumerate() {
                    let coefficient = self.contrast[i * dimension + j];
                    residual.add(coefficient * value);
                    scale += coefficient.abs() * value.abs();
                }
                let defect = residual.finish().abs();
                let relative = if !defect.is_finite() || !scale.is_finite() {
                    f64::INFINITY
                } else if scale > 0.0 {
                    defect / scale
                } else if defect == 0.0 {
                    0.0
                } else {
                    f64::INFINITY
                };
                relative_backward_residual = relative_backward_residual.max(relative);
            }
        }
        if correction.iter().any(|&z| !finite(z))
            || relative_backward_residual > controls.maximum_relative_residual
        {
            return Err(bad(
                "internal transmission solve exceeds its residual admission",
            ));
        }
        let reflected_velocity: Vec<_> = correction[n..]
            .iter()
            .map(|&w| C64::new(w.im / z_o, -w.re / z_o))
            .collect();
        let total_pressure: Vec<_> = reference[..n]
            .iter()
            .zip(&correction[..n])
            .map(|(&source, &reflected)| source + reflected)
            .collect();
        let total_velocity: Vec<_> = reference[n..]
            .iter()
            .zip(&correction[n..])
            .map(|(&source, &reflected)| {
                let w = source + reflected;
                C64::new(w.im / z_o, -w.re / z_o)
            })
            .collect();
        let transmitted = self.radiation(
            self.k_o,
            medium(self.exterior),
            total_pressure,
            total_velocity,
        )?;
        let reflected = self.radiation(
            self.k_i,
            interior,
            correction[..n].to_vec(),
            reflected_velocity,
        )?;
        poll(gate)?;
        Ok(InternalSourceRadiation {
            interface: self,
            sources,
            transmitted,
            reflected,
            relative_backward_residual,
        })
    }
}

/// Original finite sources, continuous total boundary traces and their directly
/// solved interior reflection. No total-minus-source recovery is performed.
pub struct InternalSourceRadiation<'p, 's> {
    interface: &'p FluidInterface<'s>,
    sources: IncidentField<'p>,
    transmitted: RadiationSolution,
    reflected: RadiationSolution,
    /// Native corrected-system backward error, excluding spatial/model error.
    pub relative_backward_residual: f64,
}

impl InternalSourceRadiation<'_, '_> {
    /// Original outgoing exterior traces for native far-field/directivity
    /// owners. Pressure and normal velocity are the same physical traces on
    /// both sides. The signed surface work is actual outgoing work, not a
    /// source-free defect or an inferred absorption channel. Like the native
    /// BEM far-field integral, this panel representation needs mesh refinement
    /// even when the two materials happen to match exactly.
    pub fn transmitted(&self) -> &RadiationSolution {
        &self.transmitted
    }
    pub fn reflected_pressure_pa(&self) -> &[C64] {
        &self.reflected.pressure
    }
    pub fn reflected_normal_velocity_m_s(&self) -> &[C64] {
        &self.reflected.velocity
    }

    /// Full pressure AND particle velocity in an explicitly declared region.
    /// Quadrature discrepancy estimates cover only the integrated boundary
    /// field. They exclude BEM/source error and cancellation conditioning.
    /// Matched media return the analytic field with zero quadrature work.
    pub fn pressure_velocity_at(
        &self,
        region: InternalSourceRegion,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: near_field::Options,
        gate: &CancelGate,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        poll(gate)?;
        options.validate()?;
        let model = self.interface;
        let geometry = match region {
            InternalSourceRegion::Interior => near_field::Geometry::new_interior_cancellable(
                model.surface,
                points_m,
                minimum_clearance_m,
                gate,
            )?,
            InternalSourceRegion::Exterior => near_field::Geometry::new_iterative(
                model.surface,
                points_m,
                minimum_clearance_m,
                gate,
            )?,
        };
        let mut observation = if model.factor.is_none() {
            zero_observation(points_m.len())
        } else if region == InternalSourceRegion::Interior {
            self.observe_reflection(&geometry, options, gate)?
        } else {
            geometry
                .prepare_velocity_cancellable(model.k_o, medium(model.exterior), options, gate)?
                .evaluate_velocity(&self.transmitted)?
        };
        if region == InternalSourceRegion::Interior || model.factor.is_none() {
            let state = medium(model.interior);
            for (i, &point) in points_m.iter().enumerate() {
                poll(gate)?;
                let source = self.sources.sample_wave(
                    model.k_i,
                    model.omega * state.density,
                    Some(state.density * state.sound_speed),
                    point,
                    &mut || gate.is_requested(),
                )?;
                observation.scalar.pressure[i] =
                    observation.scalar.pressure[i] + source.pressure_pa;
                for axis in 0..3 {
                    observation.particle_velocity_m_s[i][axis] =
                        observation.particle_velocity_m_s[i][axis] + source.velocity_m_s[axis];
                }
                if !finite(observation.scalar.pressure[i])
                    || observation.particle_velocity_m_s[i]
                        .iter()
                        .any(|&v| !finite(v))
                {
                    return Err(bad("internal-source total observation overflow"));
                }
            }
        }
        poll(gate)?;
        Ok(observation)
    }

    /// Weak interior reflected fields remain observable separately from the
    /// potentially much larger direct source field.
    pub fn reflected_pressure_velocity_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: near_field::Options,
        gate: &CancelGate,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        poll(gate)?;
        options.validate()?;
        let geometry = near_field::Geometry::new_interior_cancellable(
            self.interface.surface,
            points_m,
            minimum_clearance_m,
            gate,
        )?;
        self.observe_reflection(&geometry, options, gate)
    }

    fn observe_reflection(
        &self,
        geometry: &near_field::Geometry<'_>,
        options: near_field::Options,
        gate: &CancelGate,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        poll(gate)?;
        if self.interface.factor.is_none() {
            return Ok(zero_observation(geometry.points_m().len()));
        }
        let mut field = geometry
            .prepare_velocity_cancellable(
                self.interface.k_i,
                medium(self.interface.interior),
                options,
                gate,
            )?
            .evaluate_velocity(&self.reflected)?;
        // The existing receiver evaluates Dp-Sq; the regular interior field
        // has the opposite sign. Velocity changes with that same field once.
        for (pressure, velocity) in field
            .scalar
            .pressure
            .iter_mut()
            .zip(&mut field.particle_velocity_m_s)
        {
            *pressure = -*pressure;
            *velocity = velocity.map(|v| -v);
        }
        poll(gate)?;
        Ok(field)
    }
}

fn zero_observation(count: usize) -> near_field::VectorObservation {
    near_field::VectorObservation {
        scalar: near_field::Observation {
            pressure: vec![C64::ZERO; count],
            quadrature_error_estimate_pa: vec![0.0; count],
        },
        particle_velocity_m_s: vec![[C64::ZERO; 3]; count],
        quadrature_error_estimate_m_s: vec![[0.0; 3]; count],
    }
}
