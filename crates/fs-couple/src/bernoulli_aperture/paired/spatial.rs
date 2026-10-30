//! Spatial interface selection and its two-coordinate nonlinear island.
//! The original material midpoint response, pressure/flow law, energy ledger,
//! sample publication and tube/network compositions remain the only owners.

use super::{PairedAperture, Trial, ApertureDrive, AcousticRealizeError, invalid};
use crate::bernoulli_aperture::plate::closure::{PlateClosureSpec, paired::PairedPlateClosure};

#[derive(Clone, Copy)]
struct SpatialTrial {
    velocity: [f64;2],
    dp: f64,
    residual: [f64;2],
    scale: [f64;2],
}
impl SpatialTrial {
    fn norm(&self) -> f64 { self.residual[0].hypot(self.residual[1]) }
}

impl PairedAperture {
    /// Select actual nodal clearance and distributed contact before stepping.
    /// This explicitly REPLACES the uniform-gap contact supplied at construction;
    /// it is not an additional lay or a live contact/energy reset. Initial blade
    /// vibration and each already-attached material memory remain unchanged.
    ///
    /// Both original plate meshes/slit edges must match in their planar interface
    /// coordinates. Their material sections, supports and reduced shapes may
    /// differ. Only the common slit/contact geometry changes: the original full
    /// uniform-pressure faces and swept-volume areas remain work conjugate.
    /// No pressure-area masking by contact, nonmatching contact or curved reed
    /// reconstruction is claimed by this small-slope matching-interface image.
    ///
    /// # Errors
    /// Repeated/late selection, missing or nonmatching source plates, invalid
    /// profile/contact law, initial local penetration/slope or energy violation.
    pub fn with_plate_closure(mut self, spec: PlateClosureSpec) -> Result<Self,AcousticRealizeError> {
        if self.samples!=0 || self.closure.is_some() {
            return Err(invalid("spatial paired closure requires one explicit selection before stepping"));
        }
        let plates=self.plates.as_ref().ok_or_else(||invalid("spatial paired closure requires both original source plates"))?;
        let closure=plates[0].compile_paired_closure(&plates[1],spec)?;
        closure.validate_displacements(self.blades.map(|s|s.displacement_m))?;
        if let Some(law) = self.viscous_slit {
            closure.viscous_flow(law, self.blades.map(|s|s.displacement_m), self.spec.density_kg_m3, 0.0)?;
        }
        self.closure=Some(closure.into());
        self.validate_state(self.blades)?;
        self.validate_energy(self.stored_energy_j())?;
        Ok(self)
    }
    /// Actual source-bound spatial interface; None keeps the legacy uniform slit.
    #[must_use]
    pub fn plate_closure(&self) -> Option<&PairedPlateClosure> { self.closure.as_deref() }

    /// Solve BOTH blade momentum equations. A single total-gap root cannot
    /// represent different spatial force projections or opposing nodal motion.
    /// Numerical Jacobian trials and every backtrack consume the existing
    /// explicit evaluation budget. There is no fallback to a uniform slit,
    /// explicit force, weaker material law or partially advanced state.
    pub(super) fn solve_profiled(&self,drive:ApertureDrive,free:[f64;2],mobility:[f64;2],
        closure:&PairedPlateClosure) -> Result<Trial,AcousticRealizeError>
    {
        let before=self.blades.map(|s|s.displacement_m);
        let dt=self.spec.time_step_s;
        let z=self.spec.impedance_pa_s_m3;
        let areas=self.spec.blades.map(|b|b.pressure_area_m2);
        let denominator=1.0+z*(mobility[0]*areas[0]*areas[0]+mobility[1]*areas[1]*areas[1]);
        if !denominator.is_finite() { return Err(invalid("spatial paired pressure mobility overflowed")); }
        let mut evaluations=0;
        let mut evaluate=|velocity:[f64;2]| -> Result<SpatialTrial,AcousticRealizeError> {
            evaluations+=1;
            if evaluations>self.spec.limits.max_evaluations {
                return Err(invalid("spatial paired aperture nonlinear work budget exhausted"));
            }
            let after=[before[0]+dt*velocity[0],before[1]+dt*velocity[1]];
            let middle=[f64::midpoint(before[0],after[0]),f64::midpoint(before[1],after[1])];
            let contact=closure.response(before,after,velocity)?;
            let base=[free[0]+mobility[0]*contact.forces_n[0],free[1]+mobility[1]*contact.forces_n[1]];
            let moving=drive.body_flow_m3_s-areas[0]*base[0]-areas[1]*base[1];
            let source=(drive.upstream_pressure_pa-2.0*drive.incoming_pressure_pa-z*moving)/denominator;
            // Eliminate the same mechanical swept response for either channel
            // law. Viscous resistance uses local gaps, not the area-equivalent
            // opening that suffices only for the Bernoulli contribution.
            let dp = match self.viscous_slit {
                Some(law) => crate::bernoulli_aperture::viscous::profiled_pressure_drop(z/denominator, source,
                    |dp| closure.viscous_flow(law, middle, self.spec.density_kg_m3, dp))?,
                None => {
                    let opening = closure.open_area_m2(middle)?/self.spec.width_m;
                    crate::bernoulli_aperture::moving::fixed_opening_pressure_drop(
                        self.spec.width_m, opening, self.spec.density_kg_m3, z/denominator, source)?
                }
            };
            let predicted=[base[0]-mobility[0]*areas[0]*dp,base[1]-mobility[1]*areas[1]*dp];
            let residual=[velocity[0]-predicted[0],velocity[1]-predicted[1]];
            let scale=[velocity[0].abs()+predicted[0].abs()+free[0].abs(),
                velocity[1].abs()+predicted[1].abs()+free[1].abs()];
            if residual.iter().chain(&scale).any(|x| !x.is_finite()) {
                return Err(invalid("spatial paired momentum trial is not finite"));
            }
            Ok(SpatialTrial {velocity,dp,residual,scale})
        };
        let resolved=|t:&SpatialTrial| (0..2).all(|i| t.residual[i].abs()
            <=self.spec.limits.velocity_tolerance_m_s+self.spec.limits.relative_tolerance*t.scale[i]);
        let mut trial=evaluate(free)?;
        loop {
            if resolved(&trial) {
                return Ok(Trial {velocity:trial.velocity,dp:trial.dp,residual:trial.norm(),
                    scale:trial.scale[0]+trial.scale[1]});
            }
            let mut jacobian=[[0.0;2];2];
            for column in 0..2 {
                let delta=1e-6*(trial.velocity[column].abs()+0.01);
                let mut left=trial.velocity; let mut right=trial.velocity;
                left[column]-=delta; right[column]+=delta;
                let width=right[column]-left[column];
                if !width.is_finite() || width<=0.0 { return Err(invalid("spatial paired Jacobian increment is unresolved")); }
                let a=evaluate(left)?; let b=evaluate(right)?;
                for row in 0..2 { jacobian[row][column]=(b.residual[row]-a.residual[row])/width; }
            }
            // Scaled two-equation algebra; no iterative structural solver or
            // new time integrator. Refuse a singular/overflowing Newton system.
            let scale=jacobian.iter().flatten().fold(0.0_f64,|m,x|m.max(x.abs()));
            if !scale.is_finite() || scale==0.0 || jacobian.iter().flatten().any(|x|!x.is_finite()) {
                return Err(invalid("spatial paired Jacobian is singular or nonfinite"));
            }
            let a=jacobian[0][0]/scale; let b=jacobian[0][1]/scale;
            let c=jacobian[1][0]/scale; let d=jacobian[1][1]/scale;
            let determinant=a*d-b*c;
            let rhs=[-trial.residual[0]/scale,-trial.residual[1]/scale];
            let direction=[(rhs[0]*d-b*rhs[1])/determinant,(a*rhs[1]-rhs[0]*c)/determinant];
            if !determinant.is_finite() || determinant==0.0 || direction.iter().any(|x|!x.is_finite()) {
                return Err(invalid("spatial paired Newton direction is unresolved"));
            }
            let mut fraction=1.0;
            let mut accepted=None;
            for _ in 0..16 {
                let velocity=[trial.velocity[0]+fraction*direction[0],trial.velocity[1]+fraction*direction[1]];
                let next=evaluate(velocity)?;
                if resolved(&next) || next.norm()<trial.norm() {
                    accepted=Some(next); break;
                }
                fraction*=0.5;
            }
            trial=accepted.ok_or_else(||invalid("spatial paired Newton solve failed to reduce the momentum residuals"))?;
        }
    }
}
