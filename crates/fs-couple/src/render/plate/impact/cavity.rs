//! Reciprocal time-domain realization of the existing rigid-wall cavity basis.
//!
//! With outward structural displacement q, C_rj = integral(phi_r psi_j dA),
//! and A_j = rho*c^2/Lambda_j, retain
//! H_air = sum_j A_j (C_j.q + omega_j*z_j/sqrt(A_j))^2/2 + p_j^2/2.
//! Nonzero modes get a free canonical (z,p) pair. The uniform zero mode gets
//! only its compression spring: no spurious free acoustic coordinate. Pressure
//! a_j = -A_j (C_j.q + omega_j*z_j/sqrt(A_j)) acts back through +C_j*a_j.
//! Eliminating z gives a_j'' + omega_j^2*a_j = -A_j*C_j.q'' when lossless,
//! exactly the conservative frequency-domain model in `vibroacoustic`.
//!
//! This compiles into existing ImpactBody/VolumeSpring/ContactStorage owners.
//! No second integrator, one-way forcing, synthetic resonances, or extra bulk
//! spring is introduced. Any geometrically derived CavityModes can be supplied.
use super::{BodyPotential, ImpactBody, ImpactConfig, ImpactError, ImpactSystem,
    MAX_IMPACT_MODES, VolumeSpring, felt::FeltPad, invalid};
use crate::modal_acoustic_time::ModalAcousticState;
use crate::vibroacoustic::CavityModes;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

/// Geometry-derived cylindrical pressure basis, using the existing eigensolver.
pub mod cylinder;
/// Passive inertial openings coupled to the same distributed pressure field.
pub mod neck;
/// Prepared linear bodies, simultaneous contact and distributed air.
pub mod prepared;


/// Cold-compiled acoustic storage. Interface signs are OUTWARD from the gas,
/// unlike a compression-positive drum pickup. Cavity norms and couplings must
/// come from the same basis; rescaling one without the other changes physics.
#[derive(Debug, Clone)]
pub struct CavityCoupling {
    structural: usize,
    omegas: Vec<f64>,
    damping: Vec<f64>,
    springs: Vec<VolumeSpring>,
    dynamic: Vec<Option<usize>>,
    total: usize,
    // Retain the caller's admission across later neck additions. A prepared
    // wire bank may exceed 64 modes; the reference builder still refuses it.
    maximum_modes: usize,
    medium: crate::vibroacoustic::AcousticMedium,
    necks: Vec<neck::CompiledNeck>,
}
impl CavityCoupling {
    /// Compile the existing row-major structural-by-acoustic overlap matrix.
    /// `damping_per_s` is explicit drag on each acoustic momentum, NOT a
    /// conversion of frequency-domain hysteretic loss. A nonzero loss_factor
    /// refuses rather than inventing a causal constitutive law. Zero modes
    /// cannot receive momentum drag because they have no momentum coordinate.
    pub fn new(cavity: &CavityModes, structural: usize, coupling: &[f64],
        damping_per_s: &[f64]) -> Result<Self, ImpactError>
    {
        Self::new_with_mode_budget(cavity, structural, coupling, damping_per_s, MAX_IMPACT_MODES)
    }

    /// The same cavity storage with an explicit total-coordinate budget.
    /// Prepared linear bodies may exceed the nonlinear reference's 64 modes;
    /// this does not enlarge that reference solver's admission. The hard limit
    /// is the existing modal owner's 4096 coordinates. No mode is truncated.
    pub fn new_with_mode_budget(cavity: &CavityModes, structural: usize, coupling: &[f64],
        damping_per_s: &[f64], maximum_modes: usize) -> Result<Self, ImpactError>
    {
        let count = cavity.omegas.len();
        if maximum_modes == 0
            || maximum_modes > crate::modal_acoustic_time::MAX_TIME_DOMAIN_ACOUSTIC_MODES
            || structural == 0 || structural > maximum_modes || count == 0 || count > 8
            || cavity.lambdas.len() != count || cavity.interface.len() != count
            || coupling.len() != structural*count || damping_per_s.len() != count
            || cavity.loss_factor != 0.0
            || ![cavity.rho0,cavity.c0].iter().all(|v| v.is_finite() && *v > 0.0)
            || cavity.omegas.iter().any(|w| !w.is_finite() || *w < 0.0)
            || cavity.lambdas.iter().any(|l| !l.is_finite() || *l <= 0.0)
            || coupling.iter().chain(damping_per_s).any(|v| !v.is_finite())
            || cavity.interface.iter().flatten().any(|v| !v.is_finite())
            || damping_per_s.iter().any(|d| *d < 0.0)
        { return Err(invalid("cavity requires bounded consistent basis, finite overlaps and explicit causal loss")); }
        let bulk = cavity.rho0*cavity.c0*cavity.c0;
        let total = structural + cavity.omegas.iter().filter(|&&w| w > 0.0).count();
        if !bulk.is_finite() || bulk <= 0.0 || total > maximum_modes {
            return Err(invalid("cavity exceeds mechanical state or finite bulk-modulus budget"));
        }
        let mut next = structural;
        let mut springs = Vec::with_capacity(count);
        let mut dynamic = Vec::with_capacity(count);
        for j in 0..count {
            let scale = bulk/cavity.lambdas[j];
            if !scale.is_finite() || scale <= 0.0 {
                return Err(invalid("cavity normalization is unrepresentable"));
            }
            let mut areas = vec![0.0;total];
            for r in 0..structural { areas[r] = coupling[r*count+j]; }
            if cavity.omegas[j] > 0.0 {
                let weight = cavity.omegas[j]/scale.sqrt();
                if !weight.is_finite() || weight <= 0.0 {
                    return Err(invalid("cavity inertial coordinate is unrepresentable"));
                }
                areas[next] = weight; dynamic.push(Some(next)); next += 1;
            } else {
                if damping_per_s[j] != 0.0 { return Err(invalid("uniform cavity mode has no momentum drag")); }
                dynamic.push(None);
            }
            springs.push(VolumeSpring { bulk_modulus_pa:bulk, volume_m3:cavity.lambdas[j], areas });
        }
        Ok(Self { structural, omegas:cavity.omegas.clone(), damping:damping_per_s.to_vec(),
            springs, dynamic, total, maximum_modes,
            medium: crate::vibroacoustic::AcousticMedium {rho0:cavity.rho0,c0:cavity.c0},
            necks:Vec::new() })
    }

    /// Original solid/striker coordinates; these keep their original addresses.
    #[must_use]
    pub fn structural_modes(&self) -> usize { self.structural }
    /// All mechanical coordinates, including standing-wave and neck inertia.
    #[must_use]
    pub fn total_modes(&self) -> usize { self.total }
    /// Retained cavity basis size, including uniform compression when supplied.
    #[must_use]
    pub fn cavity_modes(&self) -> usize { self.omegas.len() }

    /// Compose before execution, preserving every original body's initial motion
    /// and every pad's declared conditioning. Standing-wave coordinates start
    /// at zero; necks use their explicit initial volumes/flows. Initial pressure
    /// includes those volumes and is not silently relaxed.
    /// Cavity springs REPLACE the compact VolumeSpring; there is no argument for
    /// adding a second copy of that same gas compliance. Contact/pad rows are
    /// extended by exact zeros so solid contact never acts directly on gas modes.
    /// Caller forces and radiation projections must likewise append exact zeros.
    ///
    /// Acoustic momentum drag is passive. Its pressure equation is
    /// a''+d*a'+omega^2*a = -A*C.(q''+d*q') without necks; with necks the same
    /// equation includes their outward volume derivatives. Not hysteretic loss.
    /// No radiation load, mean flow, thermoviscous spectrum or RT claim is inferred.
    pub fn build(self, bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>,
        pads: Vec<FeltPad>, config: ImpactConfig, gate: &CancelGate)
        -> Result<(ImpactSystem, Self), ImpactError>
    {
        self.build_with_dampers(bodies, contacts, pads, Vec::new(), config, gate)
    }

    /// Attach spatial viscous loss only to the original structural coordinates.
    /// The complete resistance participates in the same implicit cavity solve.
    pub fn build_with_dampers(self, bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>,
        pads: Vec<FeltPad>, dampers: Vec<super::damping::ViscousDamper>,
        config: ImpactConfig, gate: &CancelGate) -> Result<(ImpactSystem, Self), ImpactError>
    {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        if self.total > MAX_IMPACT_MODES {
            return Err(invalid("distributed cavity exceeds the nonlinear reference mode budget"));
        }
        let (bodies, contacts, pads) = self.extend_parts(bodies, contacts, pads, config.dt_s, gate)?;
        let dampers = super::damping::extend(dampers, self.structural, self.total)?;
        let system = ImpactSystem::new_with_dampers(bodies, contacts, pads, self.springs.clone(), dampers, config)?;
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        Ok((system, self))
    }

    // Both numerical images use exactly this address/initial-state lowering.
    // Only the downstream solver and its own work admission differ.
    fn extend_parts(&self, mut bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>,
        mut pads: Vec<FeltPad>, dt_s: f64, gate: &CancelGate)
        -> Result<(Vec<ImpactBody>, Vec<Obstacle>, Vec<FeltPad>), ImpactError>
    {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        let count = bodies.iter().try_fold(0usize, |n,b| n.checked_add(b.potential.count()));
        if count != Some(self.structural) || contacts.len() > 32 || pads.len() > 16
            || !dt_s.is_finite() || dt_s <= 0.0
            || self.omegas.iter().any(|w| w*dt_s >= 0.9*core::f64::consts::PI)
            || self.necks.iter().any(|n| n.fixed_wall_omega*dt_s >= 0.9*core::f64::consts::PI)
        { return Err(invalid("cavity/body basis mismatch or acoustic Nyquist limit")); }
        let mut extended = Vec::with_capacity(contacts.len());
        for contact in contacts {
            if gate.is_requested() { return Err(ImpactError::Cancelled); }
            if contact.n_points() > 4096 || contact.collocation().len() != contact.n_points()*self.structural {
                return Err(invalid("cavity contact must use original structural coordinates"));
            }
            let mut rows = vec![0.0;contact.n_points()*self.total];
            for (source,target) in contact.collocation().chunks_exact(self.structural)
                .zip(rows.chunks_exact_mut(self.total)) {
                target[..self.structural].copy_from_slice(source);
            }
            let ob = Obstacle::new(rows,contact.n_points(),self.total,contact.gaps().to_vec(),
                contact.weights().to_vec(),contact.stiffness(),contact.alpha(),contact.provenance().to_string())
                .and_then(|ob| ob.with_internal_loss(contact.internal_loss()))
                .map_err(|e| ImpactError::Owner(e.to_string()))?;
            extended.push(ob);
        }
        for pad in &mut pads {
            pad.validate(self.structural)?;
            pad.weights.resize(self.total,0.0);
        }
        if self.total > self.structural {
            let mut initial=vec![ModalAcousticState::default();self.total-self.structural];
            let mut damping=vec![0.0;initial.len()];
            for (coordinate,&drag) in self.dynamic.iter().zip(&self.damping) {
                if let Some(i)=coordinate { damping[*i-self.structural]=drag; }
            }
            for neck in &self.necks {
                initial[neck.coordinate-self.structural]=neck.initial;
                damping[neck.coordinate-self.structural]=neck.drag_per_s;
            }
            bodies.push(ImpactBody { potential:BodyPotential::Linear(vec![0.0;self.total-self.structural]),
                initial,damping_per_s:damping });
        }
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        Ok((bodies, extended, pads))
    }

    /// Pressure coefficients [Pa] of the accepted state in the supplied basis.
    /// Reuses caller storage; bad state/dimensions/overflow leave it unchanged.
    /// Kelvin memory can follow the interleaved prefix and is not acoustic data.
    pub fn pressures_into(&self, state: &[f64], output: &mut [f64]) -> Result<(), ImpactError> {
        if state.len() < 2*self.total || output.len() != self.springs.len()
            || state[..2*self.total].iter().any(|v| !v.is_finite()) {
            return Err(invalid("cavity pressure needs the complete finite accepted state"));
        }
        let mut values = [0.0;8];
        for (j,spring) in self.springs.iter().enumerate() {
            let volume = spring.areas.iter().enumerate().map(|(r,c)| c*state[2*r]).sum::<f64>();
            values[j] = -(spring.bulk_modulus_pa/spring.volume_m3)*volume;
            if !values[j].is_finite() { return Err(invalid("cavity pressure overflow")); }
        }
        output.copy_from_slice(&values[..self.springs.len()]);
        Ok(())
    }

    /// Interior physical pressure at a point with basis values `psi_j(point)`.
    /// Not a microphone outside the enclosure; exterior observation stays BEM.
    pub fn pressure_at(&self, state: &[f64], basis_values: &[f64]) -> Result<f64, ImpactError> {
        if basis_values.len() != self.cavity_modes() || basis_values.iter().any(|v| !v.is_finite()) {
            return Err(invalid("cavity receiver must use the same finite pressure basis"));
        }
        let mut a = [0.0;8];
        self.pressures_into(state,&mut a[..self.cavity_modes()])?;
        let pressure = a.iter().zip(basis_values).map(|(a,b)|a*b).sum::<f64>();
        if !pressure.is_finite() { return Err(invalid("cavity receiver pressure overflow")); }
        Ok(pressure)
    }
}
