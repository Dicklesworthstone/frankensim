//! Homogeneous translating specimens sharing the exact bulk/interface sources.
use super::{ImpactBody, ImpactConfig, ImpactError, ImpactSystem, MaterialContact,
    SourcePointContact, invalid};
use fs_material::state_point::DENSITY_PROPERTY;

/// Explicit geometry and initial motion of one uniform, fixed-axis specimen.
/// Surface A closes the gap by moving positively; surface B opens it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TranslatingSpecimen {
    /// Entire moving specimen volume [m3], not its contact patch volume.
    pub volume_m3: f64,
    /// Initial position along the shared contact axis [m].
    pub initial_position_m: f64,
    /// Initial velocity along the shared contact axis [m/s].
    pub initial_velocity_m_per_s: f64,
}

/// Immutable cold assembly; all motion stays in the existing ImpactSystem.
///
/// The two surface bulk states also describe their corresponding WHOLE uniform
/// specimens. For layered tips or heterogeneous bodies, use resolve_body and
/// explicit MaterialContact maps instead. This compiler does not silently
/// assign surface density to an unrelated backing structure.
#[derive(Clone, Debug, PartialEq)]
pub struct TranslatingMaterialPair {
    key: String,
    source: SourcePointContact,
    specimens: [TranslatingSpecimen; 2],
    masses_kg: [f64; 2],
    gap_m: f64,
}
impl TranslatingMaterialPair {
    /// Bind volumes and motion to the SAME ordered material states that already
    /// supply contact elasticity/yield. No separate density, mass, stiffness,
    /// exponent, restitution or damping input is accepted. No re-query can pick
    /// a different claim or temperature for the inertial half of this model.
    pub fn new(key: String, source: SourcePointContact,
        specimens: [TranslatingSpecimen; 2], gap_m: f64) -> Result<Self, ImpactError> {
        if key.trim().is_empty() || key.len()>128 || !gap_m.is_finite() {
            return Err(invalid("material pair needs a bounded contact key and finite reference gap"));
        }
        let elastic=source.model().normal_state().elastic();
        let states=[elastic.surface_a_state(),elastic.surface_b_state()];
        let mut masses_kg=[0.0;2];
        for i in 0..2 {
            let specimen=specimens[i];
            if !specimen.volume_m3.is_finite() || specimen.volume_m3<=0.0 {
                return Err(invalid("material pair volume must be finite and positive"));
            }
            let density=states[i].resolved().property(DENSITY_PROPERTY)
                .ok_or_else(||invalid("material pair bulk source omitted density"))?.value_si();
            masses_kg[i]=density*specimen.volume_m3;
            // Same admission as the downstream inertial owner, including the
            // authored motion. No clipping or normalization of launch energy.
            ImpactBody::free_mass(masses_kg[i],specimen.initial_position_m,
                specimen.initial_velocity_m_per_s)?;
            let root=masses_kg[i].sqrt();
            if !(specimen.initial_position_m*root).is_finite()
                || !(specimen.initial_velocity_m_per_s*root).is_finite()
                || !(1.0/root).is_finite() {
                return Err(invalid("material pair mass-normalized motion is unrepresentable"));
            }
        }
        Ok(Self { key,source,specimens,masses_kg,gap_m })
    }
    /// Full ordered material/interface model used for both inertia and contact.
    pub const fn source(&self) -> &SourcePointContact { &self.source }
    /// Original authored geometric volumes and physical initial motion.
    pub const fn specimens(&self) -> &[TranslatingSpecimen;2] { &self.specimens }
    /// Source-derived masses in ordered A/B roles [kg].
    pub const fn masses_kg(&self) -> [f64;2] { self.masses_kg }

    /// Construct the existing joint nonlinear runtime from this frozen assembly.
    /// Rebuilding creates a new experiment; it never changes a live history.
    /// Contact limits and the original total-energy admission still apply.
    pub fn build(&self, config: ImpactConfig) -> Result<ImpactSystem, ImpactError> {
        let mut bodies=Vec::with_capacity(2);let mut weights=Vec::with_capacity(2);
        for (i,specimen) in self.specimens.iter().enumerate() {
            let (body,weight)=ImpactBody::free_mass(self.masses_kg[i],
                specimen.initial_position_m,specimen.initial_velocity_m_per_s)?;
            bodies.push(body);weights.push(if i==0 {weight}else{-weight});
        }
        let contact=MaterialContact::new(self.key.clone(),self.source.clone(),weights,self.gap_m)?;
        ImpactSystem::new_with_material_contacts(bodies,vec![contact],vec![],vec![],config)
    }
}
