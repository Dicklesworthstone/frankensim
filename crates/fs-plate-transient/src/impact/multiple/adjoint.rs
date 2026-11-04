//! Pullbacks for one shared plate and all of its simultaneous contact bodies.
//!
//! The plate/load prefix appears once. Each opted-in authored contact appends
//! its original mass/K/chi/gap coordinates in body order. Optional body-force
//! parameters follow; only the forcing pullback contributes to that suffix.
//! Every state action
//! uses the complete shared plate state and the corresponding body coordinate.

use super::PlateImpactSet;
use crate::PlateLoad;
use crate::impact::{ContactPart, finite};
use fs_time::galpha::SecondOrderProblem;
use fs_time::galpha::second_order_adjoint::SecondOrderVjp;
use fs_time::galpha::second_order_adjoint::trajectory::StructuralTrajectoryModel;

impl<L: PlateLoad + ?Sized> SecondOrderVjp for PlateImpactSet<'_, L> {
    fn parameter_count(&self) -> usize {
        self.parameters
    }

    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        if !self.io(&[seed], output) {
            return Err("invalid simultaneous impact mass transpose".into());
        }
        let n = self.plate_dofs;
        self.plate()
            .mass_transpose_apply(&seed[..n], &mut output[..n])?;
        for (i, body) in self.bodies.iter().enumerate() {
            output[n + i] = body.indenter_mass_kg * seed[n + i];
        }
        if !finite(output) {
            return Err("nonfinite simultaneous impact mass transpose".into());
        }
        Ok(())
    }

    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        if !self.io(&[seed], output) {
            return Err("invalid simultaneous impact damping transpose".into());
        }
        let n = self.plate_dofs;
        self.plate()
            .damping_transpose_apply(&seed[..n], &mut output[..n])?;
        output[n..].fill(0.0);
        Ok(())
    }

    fn tangent_transpose_apply(
        &self,
        q: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        if !self.io(&[q, seed], output) {
            return Err("invalid simultaneous impact elastic transpose".into());
        }
        let n = self.plate_dofs;
        self.plate()
            .tangent_transpose_apply(&q[..n], &seed[..n], &mut output[..n])?;
        output[n..].fill(0.0);
        self.add_contacts(q, None, Some(seed), ContactPart::ElasticTangent, output)
    }

    fn resistance_displacement_transpose_apply(
        &self,
        q: &[f64],
        v: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        if !self.io(&[q, v, seed], output) {
            return Err("invalid simultaneous impact displacement transpose".into());
        }
        let n = self.plate_dofs;
        self.plate().resistance_displacement_transpose_apply(
            &q[..n],
            &v[..n],
            &seed[..n],
            &mut output[..n],
        )?;
        output[n..].fill(0.0);
        // Each original local contact tangent is a scalar times B_i^T B_i.
        // Only that contact contribution is symmetric; the plate transpose
        // above remains explicitly owned by the original sparse operators.
        self.add_contacts(
            q,
            Some(v),
            Some(seed),
            ContactPart::DisplacementTangent,
            output,
        )
    }

    fn resistance_velocity_transpose_apply(
        &self,
        q: &[f64],
        v: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        if !self.io(&[q, v, seed], output) {
            return Err("invalid simultaneous impact velocity transpose".into());
        }
        let n = self.plate_dofs;
        self.plate().resistance_velocity_transpose_apply(
            &q[..n],
            &v[..n],
            &seed[..n],
            &mut output[..n],
        )?;
        output[n..].fill(0.0);
        self.add_contacts(q, Some(v), Some(seed), ContactPart::VelocityTangent, output)
    }

    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        a: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        self.check(&[q, v, a, seed])
            .map_err(|error| error.to_string())?;
        if output.len() != self.parameters {
            return Err("invalid simultaneous impact parameter derivative dimension".into());
        }
        let body_forcing = self.checked_body_forcing()?;
        // Admit every original source before evaluating or publishing any
        // parameter contribution, including bodies whose coefficients are fixed.
        self.validate_state(q, v)?;
        let n = self.plate_dofs;
        let prefix = self.plate().parameter_count();
        self.plate().residual_parameter_vjp(
            &q[..n],
            &v[..n],
            &a[..n],
            &seed[..n],
            &mut output[..prefix],
        )?;
        for (i, body) in self.bodies.iter().enumerate() {
            let range = self.parameter_ranges[i].clone();
            if range.is_empty() {
                continue;
            }
            let q = self.local(q, i).map_err(|error| error.to_string())?;
            let v = self.local(v, i).map_err(|error| error.to_string())?;
            let a = self.local(a, i).map_err(|error| error.to_string())?;
            let seed = self.local(seed, i).map_err(|error| error.to_string())?;
            body.contact_parameter_vjp(&q, &v, &a, &seed, &mut output[range])
                .map_err(|error| format!("impact body {i}: {error}"))?;
        }
        if let Some((_, range)) = body_forcing {
            // Applied force controls enter through forcing_vjp, never twice
            // through the mass/resistance residual at fixed state.
            output[range].fill(0.0);
        }
        if !finite(output) {
            return Err("nonfinite simultaneous impact parameter derivative".into());
        }
        Ok(())
    }
}

impl<L: PlateLoad + ?Sized> StructuralTrajectoryModel for PlateImpactSet<'_, L> {
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
        if output.len() != self.dimension || !time.is_finite() {
            return Err("invalid simultaneous impact forcing time or dimension".into());
        }
        let body_forcing = self.checked_body_forcing()?;
        self.plate().forcing(time, &mut output[..self.plate_dofs])?;
        let bodies = &mut output[self.plate_dofs..];
        if let Some((load, _)) = body_forcing {
            bodies.fill(f64::NAN);
            load.forcing(time, bodies)
                .map_err(|error| format!("body-force history: {error}"))?;
            if !finite(bodies) {
                return Err("nonfinite or unwritten body-force history output".into());
            }
        } else {
            bodies.fill(0.0);
        }
        Ok(())
    }

    fn forcing_vjp(&self, time: f64, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.check(&[seed]).map_err(|error| error.to_string())?;
        if output.len() != self.parameters || !time.is_finite() {
            return Err("invalid simultaneous impact forcing derivative time or dimension".into());
        }
        let body_forcing = self.checked_body_forcing()?;
        let prefix = self.plate().parameter_count();
        self.plate()
            .forcing_vjp(time, &seed[..self.plate_dofs], &mut output[..prefix])?;
        output[prefix..].fill(0.0);
        if let Some((load, range)) = body_forcing {
            output[range.clone()].fill(f64::NAN);
            load.forcing_vjp(time, &seed[self.plate_dofs..], &mut output[range])
                .map_err(|error| format!("body-force history derivative: {error}"))?;
            if !finite(output) {
                return Err("nonfinite or unwritten body-force derivative output".into());
            }
        }
        Ok(())
    }
}
