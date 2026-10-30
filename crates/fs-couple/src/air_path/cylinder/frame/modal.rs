//! Exact constant-inertia coordinates for the shared modal assembly.
//!
//! For the original loaded M=L L^T, q=L^-T r gives unit kinetic storage.
//! Every physical force/observation column uses L^-1 b. This is a congruence,
//! not a wet-mode fit: stiffness, solid losses and fluid memory remain coupled.

use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use crate::render::schedule::force::coupled::{ModalAttachment, ModalConnection, ModalCouplingConfig};
use fs_la::factor::Cholesky;
use fs_math::c64::C64;

mod component;
pub use component::{
    FrameFluidComponent, FrameFluidConnectionLayout, FrameFluidPlacement, FrameFluidPortMap,
};

/// Original source and complete constant-inertia coordinate change. The source
/// basis remains dry; use `original_states` before reconstructing physical motion.
/// No state, material selection, reference geometry or mode is discarded.
pub struct FrameFluidCoordinates<'m, 'f, F> {
    source: &'m MaterialFluidFrame<'f, F>,
    factor: Cholesky,
}

impl<'m, 'f, F> FrameFluidCoordinates<'m, 'f, F> {
    pub const fn source(&self) -> &'m MaterialFluidFrame<'f, F> { self.source }

    /// b_native=L^-1 b_original. The SAME transform applies to an observation
    /// column and its conjugate generalized force. It does not alter units.
    pub fn project_column(&self, original: &[f64]) -> Result<Vec<f64>> {
        let mut result = vec![0.0; self.source.basis.modes().len()];
        self.project_column_into(original, &mut result)?;
        Ok(result)
    }

    /// Allocation-free force/column projection into caller scratch. Refusal may
    /// overwrite scratch but never changes a physical state or source object.
    pub fn project_column_into(&self, original: &[f64], output: &mut [f64]) -> Result<()> {
        lower_solve_into(&self.factor, original, output, self.source.basis.modes().len())
    }

    /// r=L^T q and rdot=L^T qdot, including a nonzero initial history.
    pub fn normalized_states(&self, original: &[ModalAcousticState]) -> Result<Vec<ModalAcousticState>> {
        map_states(&self.factor, original, self.source.basis.modes().len(), false)
    }

    /// q=L^-T r and qdot=L^-T rdot in the original structural basis. Reference
    /// eigenfrequencies are not the evolving stiffness; this only changes charts.
    pub fn original_states(&self, normalized: &[ModalAcousticState]) -> Result<Vec<ModalAcousticState>> {
        map_states(&self.factor, normalized, self.source.basis.modes().len(), true)
    }

    /// Physical nodal force/moment -> dry generalized force -> native force.
    /// Fixed supports and stationary exterior only; no base-acceleration term.
    pub fn project_nodal_force(&self, nodal: &[f64], gate: &CancelGate) -> Result<Vec<f64>> {
        poll(gate)?;
        let mut force = vec![0.0; self.source.basis.modes().len()];
        self.source.basis.project_force_into(nodal, [0.0; 2], &mut force, || gate.is_requested())?;
        let result = self.project_column(&force)?;
        poll(gate)?;
        Ok(result)
    }
}

/// Free structural coordinates and ALL fluid memory on the original native
/// modal owner. This intentionally supplies NO solid stiffness or viscosity:
/// the source constitutive owner must add those ONCE using projected columns.
/// Connections refer to structural component 0 and fluid component 1.
/// Retain the coordinates/source together with every returned model and link.
pub struct FrameFluidModalParts<'m, 'f, F> {
    pub coordinates: FrameFluidCoordinates<'m, 'f, F>,
    pub structural: ModalAcousticTimeModel,
    pub memory: Option<ModalAcousticTimeModel>,
    pub connections: Vec<ModalConnection>,
}

impl<'f, F> MaterialFluidFrame<'f, F> {
    /// Lift the SAME complete physical state into the shared free-coordinate
    /// assembly, without recomputing solid spectra or fluid Bessel quadrature.
    ///
    /// The modal owner gives each fluid momentum a free position gauge; it has
    /// no spring, observation or input port. These additional coordinates count
    /// toward its complete mode/setup budgets. Gauge positions start at zero;
    /// velocities retain every original fluid momentum. No rest-fluid reset.
    ///
    /// `reserved_connections` includes ALL separately supplied solid source
    /// ports. Final assembly admission still owns the complete native budgets.
    /// The fluid is fixed at its prepared state; live density or moving geometry
    /// cannot be represented by changing this constant coordinate transform.
    #[allow(clippy::too_many_arguments)]
    pub fn modal_parts<'m>(
        &'m self,
        state: &[f64],
        sample_rate_hz: u32,
        modal_budget: ModalAcousticTimeBudget,
        reserved_connections: usize,
        config: ModalCouplingConfig,
        gate: &CancelGate,
    ) -> Result<FrameFluidModalParts<'m, 'f, F>> {
        poll(gate)?;
        self.admit_state(state)?;
        let m = self.basis.modes().len();
        let k = self.memory_masses.len();
        let total = m.checked_add(k).ok_or_else(|| invalid("complete modal fluid count overflow"))?;
        let links = self.stations.len().checked_add(k)
            .and_then(|n| n.checked_add(reserved_connections))
            .ok_or_else(|| invalid("complete modal fluid connection count overflow"))?;
        if total > config.max_modes || links > config.max_connections
            || links.checked_add(1).and_then(|n| n.checked_mul(n))
                .and_then(|n| n.checked_mul(total))
                .is_none_or(|n| n > config.max_setup_terms)
            || total.checked_mul(total).and_then(|n| n.checked_mul(links))
                .is_none_or(|n| n > config.max_setup_terms) {
            return Err(invalid("all fluid gauges, material ports and modal setup must fit the native budget"));
        }
        let factor = cholesky(&self.modal_mass, m).map_err(FrameFluidError::Factor)?;
        let coordinates = FrameFluidCoordinates { source: self, factor };
        let modes = |count| vec![ModalAcousticMode { angular_frequency_rad_s: 0.0,
            damping_ratio: 0.0, pressure_per_modal_velocity: C64::ZERO }; count];
        let modal = |error| FrameFluidError::Frame(FrameVibrationError::Modal(error));
        let mut structural = ModalAcousticTimeModel::try_new_with_free_coordinates(
            sample_rate_hz, modes(m), modal_budget).map_err(modal)?;
        structural.restore_states(&coordinates.normalized_states(&self.modal_states(state)?)?).map_err(modal)?;
        let memory = if k == 0 { None } else {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                sample_rate_hz, modes(k), modal_budget).map_err(modal)?;
            let mut states = Vec::with_capacity(k);
            for (i, mass) in self.memory_masses.iter().enumerate() {
                poll(gate)?;
                states.push(ModalAcousticState { displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: finite(state[2*m+i] / positive(mass.sqrt())?)? });
            }
            model.restore_states(&states).map_err(modal)?;
            Some(model)
        };
        let mut connections = Vec::with_capacity(self.stations.len()+k);
        for station in &self.stations {
            poll(gate)?;
            let shape = coordinates.project_column(&station.projection)?;
            let flow = &self.members[station.member].flow;
            connections.push(ModalConnection {
                left: ModalAttachment { component: 0, shapes: shape.clone() },
                right: ModalAttachment { component: 0, shapes: vec![0.0; m] },
                stiffness_n_m: 0.0,
                damping_n_s_m: positive(station.weight_m*flow.low_tail_resistance_per_length_n_s_m2())?,
                rest_extension_m: 0.0,
            });
            for (i, branch) in flow.branches().iter().enumerate() {
                poll(gate)?;
                let index = station.first_memory+i;
                let mut fluid = vec![0.0; k];
                fluid[index] = positive(1.0 / self.memory_masses[index].sqrt())?;
                connections.push(ModalConnection {
                    left: ModalAttachment { component: 0, shapes: shape.clone() },
                    right: ModalAttachment { component: 1, shapes: fluid },
                    stiffness_n_m: 0.0,
                    damping_n_s_m: positive(station.weight_m*branch.resistance_per_length_n_s_m2)?,
                    rest_extension_m: 0.0,
                });
            }
        }
        poll(gate)?;
        Ok(FrameFluidModalParts { coordinates, structural, memory, connections })
    }
}

#[cfg(test)]
fn lower_solve(l: &Cholesky, input: &[f64], n: usize) -> Result<Vec<f64>> {
    let mut out = vec![0.0; n];
    lower_solve_into(l, input, &mut out, n)?;
    Ok(out)
}
fn lower_solve_into(l: &Cholesky, input: &[f64], out: &mut [f64], n: usize) -> Result<()> {
    if input.len() != n || out.len() != n || input.iter().any(|x| !x.is_finite()) {
        return Err(invalid("complete finite inertia-coordinate column required"));
    }
    out.copy_from_slice(input);
    for i in 0..n {
        let mut value = out[i];
        for j in 0..i { value = (-l.l(i,j)).mul_add(out[j], value); }
        out[i] = finite(value/l.l(i,i))?;
    }
    Ok(())
}
fn map_states(l: &Cholesky, input: &[ModalAcousticState], n: usize, inverse: bool)
    -> Result<Vec<ModalAcousticState>> {
    if input.len() != n || input.iter().any(|s| !s.displacement_m_sqrt_kg.is_finite()
        || !s.velocity_m_sqrt_kg_per_s.is_finite()) {
        return Err(invalid("complete finite inertia-coordinate state required"));
    }
    let mut out = input.to_vec();
    for i in (0..n).rev() {
        let mut q = if inverse { input[i].displacement_m_sqrt_kg }
            else { l.l(i,i)*input[i].displacement_m_sqrt_kg };
        let mut v = if inverse { input[i].velocity_m_sqrt_kg_per_s }
            else { l.l(i,i)*input[i].velocity_m_sqrt_kg_per_s };
        for j in i+1..n {
            if inverse {
                q = (-l.l(j,i)).mul_add(out[j].displacement_m_sqrt_kg, q);
                v = (-l.l(j,i)).mul_add(out[j].velocity_m_sqrt_kg_per_s, v);
            } else {
                q = l.l(j,i).mul_add(input[j].displacement_m_sqrt_kg, q);
                v = l.l(j,i).mul_add(input[j].velocity_m_sqrt_kg_per_s, v);
            }
        }
        if inverse { q /= l.l(i,i); v /= l.l(i,i); }
        out[i] = ModalAcousticState { displacement_m_sqrt_kg: finite(q)?,
            velocity_m_sqrt_kg_per_s: finite(v)? };
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn near(a: f64, b: f64) { assert!((a-b).abs() <= 2e-12*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}"); }
    #[test]
    fn full_inertia_coordinate_roundtrip_preserves_energy_and_virtual_work() {
        let mass = [4.0,1.0,-0.5, 1.0,3.0,0.2, -0.5,0.2,2.0];
        let l = cholesky(&mass,3).unwrap();
        let x: Vec<_> = (0..3).map(|i| ModalAcousticState {
            displacement_m_sqrt_kg: [0.1,-0.3,0.2][i],
            velocity_m_sqrt_kg_per_s: [-0.2,0.5,0.7][i] }).collect();
        let r = map_states(&l,&x,3,false).unwrap();
        let restored = map_states(&l,&r,3,true).unwrap();
        let b = [0.7,-0.8,0.4];
        let native = lower_solve(&l,&b,3).unwrap();
        let mut physical_energy = 0.0;
        for i in 0..3 {
            near(x[i].displacement_m_sqrt_kg,restored[i].displacement_m_sqrt_kg);
            near(x[i].velocity_m_sqrt_kg_per_s,restored[i].velocity_m_sqrt_kg_per_s);
            for j in 0..3 { physical_energy += 0.5*mass[i*3+j]*x[i].velocity_m_sqrt_kg_per_s*x[j].velocity_m_sqrt_kg_per_s; }
        }
        near(physical_energy,r.iter().map(|s| 0.5*s.velocity_m_sqrt_kg_per_s.powi(2)).sum());
        near(b.iter().zip(&x).map(|(b,s)| b*s.velocity_m_sqrt_kg_per_s).sum(),
            native.iter().zip(&r).map(|(b,s)| b*s.velocity_m_sqrt_kg_per_s).sum());
        near(b.iter().zip(&x).map(|(b,s)| b*s.displacement_m_sqrt_kg).sum(),
            native.iter().zip(&r).map(|(b,s)| b*s.displacement_m_sqrt_kg).sum());
    }
    #[test]
    fn inertia_coordinate_conversion_refuses_incomplete_and_nonfinite_data() {
        let l = cholesky(&[2.0,0.3,0.3,1.0],2).unwrap();
        assert!(lower_solve(&l,&[1.0],2).is_err());
        assert!(lower_solve(&l,&[1.0,f64::NAN],2).is_err());
        assert!(map_states(&l,&[ModalAcousticState::default()],2,true).is_err());
        let mut x = [ModalAcousticState::default();2];
        x[1].velocity_m_sqrt_kg_per_s = f64::INFINITY;
        assert!(map_states(&l,&x,2,false).is_err());
    }
}
