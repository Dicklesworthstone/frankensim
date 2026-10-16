//! Compact inertial openings into a zero-gauge-pressure reservoir.
//!
//! Let v be outward displaced volume and Q=v'. A neck has acoustic inertance
//! L=rho*ell_eff/S, kinetic energy L*Q^2/2, and pressure drop R*Q. In the existing
//! unit-mass coordinates q=sqrt(L)*v, p=sqrt(L)*Q, add psi_bar_j*q/sqrt(L) to
//! each cavity compression spring. Thus p'=sum_j psi_bar_j*a_j/sqrt(L)-R*p/L:
//! the SAME aperture average drives pressure feedback and displaced volume.
//! Neck dissipation is R*Q^2, with no added damping of head coordinates.
//!
//! The uniform fixed-wall limit is omega_H=c*sqrt(S/(V*ell_eff)). Effective
//! length includes caller-supplied end corrections; none are guessed here.
//! Reference: UNSW Musical Acoustics, https://phys.unsw.edu.au/jw/Helmholtz.html.
//! This is a linear compact-neck chart, NOT a distributed duct, turbulent jet,
//! frequency-dependent thermoviscous impedance or exterior radiation solver.
use super::{CavityCoupling,ImpactError,ModalAcousticState,invalid};
use fs_exec::CancelGate;

/// Separate interior neck impedance from a coupled exterior pressure field.
pub mod exterior;

/// Physical opening and its averaged pressure-basis values, in cavity order.
#[derive(Debug,Clone)]
pub struct CavityNeck {
    /// Throat cross section [m^2].
    pub area_m2:f64,
    /// Acoustic length [m], including explicitly chosen end corrections.
    pub effective_length_m:f64,
    /// Passive, frequency-independent pressure/volume-flow resistance [Pa s/m^3].
    pub resistance_pa_s_m3:f64,
    /// Integral(psi_j dA)/area, using the SAME normalization as the cavity.
    /// Finite-aperture averaging belongs to the geometry provider, not this host.
    pub pressure_shape_averages:Vec<f64>,
    /// Outward displaced volume at initialization [m^3]; zero is an unshifted slug.
    pub initial_volume_m3:f64,
    /// Outward volume flow at initialization [m^3/s].
    pub initial_flow_m3_s:f64,
}

#[derive(Debug,Clone)]
pub(super) struct CompiledNeck {
    // True only when the supplied inertia/resistance explicitly omit exterior
    // radiation. Ordinary effective-length necks keep their original meaning.
    pub exterior_load:bool,
    pub area_m2:f64,
    pub effective_length_m:f64,
    pub coordinate:usize,
    pub volume_weight:f64,
    pub drag_per_s:f64,
    pub resistance:f64,
    pub initial:ModalAcousticState,
    pub fixed_wall_omega:f64,
    pub averages:Vec<f64>,
}

/// Read-only observation port of an admitted compact neck. This does not load
/// the neck with exterior pressure or add another end correction/inertia.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct NeckRadiationPort {
    /// Original coordinate index, not an output channel.
    pub coordinate:usize,
    /// Nominal throat area [m^2].
    pub area_m2:f64,
    /// Declared acoustic length [m], including the caller's end corrections.
    pub effective_length_m:f64,
    /// Outward volume/flow per q/p: 1/sqrt(rho*ell_eff/S) [m^2/sqrt(kg)].
    pub volume_weight_m2_per_sqrt_kg:f64,
}

/// Physical quantities reconstructed from one accepted state, without stepping.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct NeckObservation {
    /// Index in the concatenated mechanical coordinate basis (not interleaved).
    pub coordinate:usize,
    /// Outward displaced volume [m^3].
    pub displaced_volume_m3:f64,
    /// Outward volume flow [m^3/s].
    pub volume_flow_m3_s:f64,
    /// Aperture-averaged interior pressure minus reservoir pressure [Pa].
    pub driving_pressure_pa:f64,
    /// Resistive pressure drop, signed in the direction of outward flow [Pa].
    pub resistive_pressure_drop_pa:f64,
    /// Stored neck kinetic energy [J]. Cavity compression is counted elsewhere.
    pub kinetic_energy_j:f64,
    /// Instantaneous nonnegative R*Q^2 [W], not a discrete step loss estimate.
    pub dissipated_power_w:f64,
}

impl CavityCoupling {
    /// Append at most eight compact necks before `build`, preserving all prior
    /// solid and acoustic addresses. Contacts, pads and external solid forcing
    /// do not act on the appended gas coordinates. Calls may append more necks;
    /// the combined count and original `new_with_mode_budget` ceiling still
    /// apply, including on repeated additions. A cancelled or invalid
    /// construction returns no partially compiled model.
    ///
    /// The retained pressure modes and each fixed-wall neck frequency must obey
    /// k*max(ell_eff,sqrt(S/pi)) <= 0.3. This is an explicit chart guard, not an
    /// error certificate or a check on a nonlinear strike's entire bandwidth.
    /// Reservoir pressure is zero gauge; no silent pressure-source work is added.
    pub fn with_necks(mut self,necks:Vec<CavityNeck>,gate:&CancelGate)->Result<Self,ImpactError> {
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        if necks.len()>8-self.necks.len()
            || self.total.checked_add(necks.len()).is_none_or(|n| n>self.maximum_modes) {
            return Err(invalid("cavity openings exceed neck or total state budget"));
        }
        let total=self.total+necks.len();
        let mut additions=Vec::with_capacity(necks.len());
        let mut columns=Vec::with_capacity(necks.len());
        for (index,neck) in necks.into_iter().enumerate() {
            if gate.is_requested() {return Err(ImpactError::Cancelled);}
            if ![neck.area_m2,neck.effective_length_m].iter().all(|v|v.is_finite() && *v>0.0)
                || !neck.resistance_pa_s_m3.is_finite() || neck.resistance_pa_s_m3<0.0
                || !neck.initial_volume_m3.is_finite() || !neck.initial_flow_m3_s.is_finite()
                || neck.pressure_shape_averages.len()!=self.cavity_modes()
                || neck.pressure_shape_averages.iter().any(|v|!v.is_finite()) {
                return Err(invalid("neck needs positive SI geometry, passive resistance and the complete finite pressure basis"));
            }
            let inertance=self.medium.rho0*neck.effective_length_m/neck.area_m2;
            let weight=1.0/inertance.sqrt();
            let drag=neck.resistance_pa_s_m3/inertance;
            let initial=ModalAcousticState {
                displacement_m_sqrt_kg:neck.initial_volume_m3/weight,
                velocity_m_sqrt_kg_per_s:neck.initial_flow_m3_s/weight,
            };
            let column:Vec<_>=neck.pressure_shape_averages.iter().map(|v|v*weight).collect();
            let stiffness=self.springs.iter().zip(&column)
                .map(|(s,b)|(s.bulk_modulus_pa/s.volume_m3)*b*b).sum::<f64>();
            if !inertance.is_finite() || inertance<=0.0 || !weight.is_finite() || weight<=0.0
                || !drag.is_finite() || (neck.resistance_pa_s_m3>0.0 && drag==0.0)
                || !stiffness.is_finite() || stiffness<=0.0
                || !initial.displacement_m_sqrt_kg.is_finite() || !initial.velocity_m_sqrt_kg_per_s.is_finite()
                || (neck.initial_volume_m3!=0.0 && initial.displacement_m_sqrt_kg==0.0)
                || (neck.initial_flow_m3_s!=0.0 && initial.velocity_m_sqrt_kg_per_s==0.0)
                || column.iter().zip(&neck.pressure_shape_averages)
                    .any(|(c,a)|!c.is_finite() || (*a!=0.0 && *c==0.0)) {
                return Err(invalid("neck inertia, coupling, drag or initial state is unrepresentable or disconnected"));
            }
            let largest=self.omegas.iter().copied().fold(stiffness.sqrt(),f64::max);
            let extent=neck.effective_length_m.max((neck.area_m2/core::f64::consts::PI).sqrt());
            let compactness=(largest/self.medium.c0)*extent;
            if !compactness.is_finite() || compactness>0.3 {
                return Err(invalid("neck is outside the declared compact acoustic chart; use a distributed duct"));
            }
            additions.push(CompiledNeck {exterior_load:false,area_m2:neck.area_m2,effective_length_m:neck.effective_length_m,
                coordinate:self.total+index,volume_weight:weight,
                drag_per_s:drag,resistance:neck.resistance_pa_s_m3,initial,
                fixed_wall_omega:stiffness.sqrt(),averages:neck.pressure_shape_averages});
            columns.push(column);
        }
        for (j,spring) in self.springs.iter_mut().enumerate() {
            spring.areas.resize(total,0.0);
            for (i,column) in columns.iter().enumerate() {spring.areas[self.total+i]=column[j];}
        }
        self.total=total;self.necks.extend(additions);
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        Ok(self)
    }

    /// Recover the exact physical port used by this mechanical neck. A pressure
    /// observer must project through this weight rather than interpreting the
    /// mass-normalized coordinate as a physical volume or surface velocity.
    pub fn neck_radiation_port(&self,index:usize)->Result<NeckRadiationPort,ImpactError> {
        let neck=self.necks.get(index).ok_or_else(||invalid("unknown cavity neck"))?;
        Ok(NeckRadiationPort {coordinate:neck.coordinate,area_m2:neck.area_m2,
            effective_length_m:neck.effective_length_m,
            volume_weight_m2_per_sqrt_kg:neck.volume_weight})
    }

    /// Number of explicit compact openings; zero retains the sealed model.
    #[must_use]
    pub fn neck_count(&self)->usize {self.necks.len()}

    /// Observe accepted neck flow and the reciprocal pressure drive. Invalid
    /// states, indices or arithmetic are refused rather than publishing NaNs.
    pub fn neck_observation(&self,state:&[f64],index:usize)->Result<NeckObservation,ImpactError> {
        let neck=self.necks.get(index).ok_or_else(||invalid("unknown cavity neck"))?;
        let driving_pressure_pa=self.pressure_at(state,&neck.averages)?;
        let q=state[2*neck.coordinate];let p=state[2*neck.coordinate+1];
        let volume=q*neck.volume_weight;let flow=p*neck.volume_weight;
        let drop=neck.resistance*flow;let energy=0.5*p*p;let power=drop*flow;
        if [volume,flow,drop,energy,power].iter().any(|v|!v.is_finite()) {
            return Err(invalid("neck observation overflow"));
        }
        Ok(NeckObservation {coordinate:neck.coordinate,displaced_volume_m3:volume,
            volume_flow_m3_s:flow,driving_pressure_pa:driving_pressure_pa,
            resistive_pressure_drop_pa:drop,kinetic_energy_j:energy,dissipated_power_w:power})
    }
}
