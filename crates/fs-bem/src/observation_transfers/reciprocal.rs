//! Receiver transfers of the original symmetric P1/P0 radiation owner.
//! Reuse common source/declaration/FIR admission, but observe its actual mixed
//! layer densities. Panel pressure means belong to work, not to this field.
use super::*;
use crate::reciprocal_radiation::{Budget,Prepared};
use crate::near_field::Options;
use std::cell::RefCell;

/// Complete finite-band receiver of the same variational source used for load.
pub struct ReciprocalBemObservation<'a> {
    declaration: BemObservation<'a>,
    boundary: Prepared<'a>,
    highest_omega: f64,
    angular_refinement: u32,
    quadrature: Options,
}
impl<'a> ReciprocalBemObservation<'a> {
    /// Bind the complete source basis and actual receiver points to one
    /// immutable variational geometry and explicit finite band. No collocation
    /// fallback, source truncation, P1 pressure flattening or auxiliary-density
    /// relabeling occurs. Numerical budgets are admitted before source solves.
    #[allow(clippy::too_many_arguments)]
    pub fn new(surface:&'a SpherePanels,medium:Medium,sources:&[SourceId],
        traces:&[Vec<f64>],receivers:&[ReceiverSpec],minimum_clearance_m:f64,
        highest_omega:f64,budget:Budget,quadrature:Options,
        mut cancelled:impl FnMut()->bool) -> Result<Self,TransferError>
    {
        if cancelled() { return Err(TransferError::Cancelled); }
        quadrature.validate()?;
        let declaration = BemObservation::new_admitted(surface,medium,sources,traces,
            receivers,minimum_clearance_m,true)?;
        let boundary = Prepared::new(surface,highest_omega/medium.sound_speed,
            budget,&mut cancelled)?;
        if boundary.frequency_work(sources.len())? > budget.maximum_work {
            return Err(TransferError::Invalid("complete observation source exceeds reciprocal work allowance"));
        }
        Ok(Self {declaration,boundary,highest_omega,
            angular_refinement:budget.angular_refinement,quadrature})
    }
    /// Original source solve work per frequency. Adaptive receiver quadrature
    /// has its own separate maximum_kernel_evaluations allowance.
    pub fn frequency_work(&self)->Result<u64,HelmholtzError> {
        self.boundary.frequency_work(self.declaration.fields.len())
    }
    pub fn geometric_delays_s(&self)->Vec<f64> { self.declaration.geometric_delays_s() }
    pub fn operator_hash(&self)->ContentHash {
        let mut hash=DomainHasher::new("fs-bem:observation-boundary:reciprocal-p1p0:v1");
        hash.update(self.declaration.operator_hash().as_bytes());
        for value in [self.highest_omega.to_bits(),u64::from(self.angular_refinement),
            self.quadrature.relative_tolerance.to_bits(),self.quadrature.maximum_depth as u64,
            self.quadrature.maximum_kernel_evaluations as u64] {
            hash.update(&value.to_le_bytes());
        }
        hash.finalize()
    }
    /// Original exp(-iwt) pressure and Cartesian pressure-gradient samples.
    /// Sources, receiver preparation and evaluation share this SAME P1/P0 solve.
    pub fn frequency_sample(&self,omega:f64,mut cancelled:impl FnMut()->bool)
        ->Result<Vec<C64>,HelmholtzError>
    {
        if cancelled() { return Err(HelmholtzError::Cancelled); }
        if !omega.is_finite() || omega<0.0 || omega>self.highest_omega {
            return Err(HelmholtzError::BadParameter {what:"reciprocal observation frequency outside its band"});
        }
        let data=&self.declaration;
        if omega==0.0 { return Ok(vec![C64::ZERO;4*data.points.len()*data.fields.len()]); }
        let k=omega/data.medium.sound_speed;
        let geometry=self.boundary.receivers(&data.points,data.minimum_clearance_m)?;
        let rows=geometry.prepare(k,data.medium,self.quadrature,true,&mut cancelled)?;
        let fields:Vec<_>=data.fields.iter().map(Vec::as_slice).collect();
        let solutions=self.boundary.solve_batch(k,data.medium,&fields,&mut cancelled)?;
        let mut output=Vec::with_capacity(4*data.points.len()*data.fields.len());
        let momentum=C64::new(0.0,omega*data.medium.density);
        for solution in &solutions {
            if cancelled() { return Err(HelmholtzError::Cancelled); }
            let field=rows.evaluate_velocity(solution)?;
            for (pressure,velocity) in field.scalar.pressure.into_iter().zip(field.particle_velocity_m_s) {
                output.push(pressure);
                // Euler: grad(p)=i*omega*rho*u under the same exp(-iwt) convention.
                output.extend(velocity.map(|v|momentum*v));
            }
        }
        if output.iter().any(|z|!z.re.is_finite() || !z.im.is_finite()) {
            return Err(HelmholtzError::BadParameter {what:"reciprocal receiver gradient overflow"});
        }
        Ok(output)
    }
    /// Prepare p/u statistics from this SAME reciprocal radiation solution.
    pub fn prepare_harmonic(&self,omega:f64,cancelled:impl FnMut()->bool)
        ->Result<PreparedHarmonicObservation,TransferError>
    {
        let data=&self.declaration;
        harmonic::admit_frequency(omega,data.medium)?;
        let samples=self.frequency_sample(omega,cancelled)?;
        PreparedHarmonicObservation::from_samples(omega,data.medium,&data.sources,
            &data.receivers,self.operator_hash(),samples)
    }
    /// Common causal FIR realization, with cancellation propagated into the
    /// original boundary assembly/LU and adaptive receiver integration.
    pub fn bake(&self,decl:&TransferDeclaration,cancelled:impl FnMut()->bool)
        ->Result<ObservationTransferBake,TransferError>
    {
        let source=&self.declaration;
        if decl.sources!=source.sources || decl.receivers!=source.receivers
            || decl.geometric_delays_s!=self.geometric_delays_s()
            || decl.producer_hash!=self.operator_hash()
            || core::f64::consts::TAU*decl.taper_stop_hz>self.highest_omega
            || matches!(decl.producer_uncertainty,ProducerUncertainty::ExactAnalytic)
        { return Err(TransferError::Invalid("cross-wired reciprocal observation declaration")); }
        let cancel=RefCell::new(cancelled);
        bake_observation_transfers(
            |omega|self.frequency_sample(omega,||(*cancel.borrow_mut())()),
            decl,||(*cancel.borrow_mut())())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn receiver_uses_original_mixed_field_and_converts_velocity_to_pressure_gradient() {
        let surface=SpherePanels::icosphere(0.1,0).unwrap();
        let mut traces=vec![vec![0.0;20];2];traces[0][0]=1.0;
        traces[1][7]=1.0/surface.areas()[7];
        let sources=vec![SourceId::StructuralMode("body".into()),
            SourceId::ApertureVolumeFlow("mouth".into())];
        let points=[ReceiverSpec {id:"near".into(),position_m:[0.17,0.12,0.21]}];
        let budget=Budget {maximum_unknowns:64,maximum_cache_bytes:1_000_000,
            maximum_work:100_000_000,angular_refinement:1};
        let quadrature=Options::default();let medium=Medium::air();let omega=500.0;
        let observed=ReciprocalBemObservation::new(&surface,medium,&sources,&traces,
            &points,0.01,omega,budget,quadrature,||false).unwrap();
        let actual=observed.frequency_sample(omega,||false).unwrap();
        let original=Prepared::new(&surface,omega/medium.sound_speed,budget,||false).unwrap();
        let receivers=original.receivers(&[points[0].position_m],0.01).unwrap();
        let row=receivers.prepare(omega/medium.sound_speed,medium,quadrature,true,||false).unwrap();
        let refs:Vec<_>=traces.iter().map(Vec::as_slice).collect();
        let solved=original.solve_velocity_ports(omega/medium.sound_speed,medium,&refs,||false).unwrap();
        for (index,solution) in solved.solutions().iter().enumerate() {
            let field=row.evaluate_velocity(solution).unwrap();
            assert_eq!(actual[4*index],field.scalar.pressure[0]);
            for axis in 0..3 {
                assert_eq!(actual[4*index+1+axis],
                    C64::new(0.0,omega*medium.density)*field.particle_velocity_m_s[0][axis]);
            }
        }
        assert!(actual.iter().any(|z|z.abs()>0.0));
        assert!(observed.frequency_sample(omega,||true).is_err());
        assert!(observed.frequency_sample(2.0*omega,||false).is_err());
    }
}
