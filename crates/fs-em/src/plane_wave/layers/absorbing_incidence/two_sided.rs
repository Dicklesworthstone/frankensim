//! Phase-sensitive illumination from both original absorbing half-spaces.
//!
//! Both sources share the original frequency, signed real kx and phasor clock.
//! Their peak [Ex,Ey] are referenced to their respective physical faces on
//! the SAME fixed Cartesian electric axes. Normal wavevectors have opposite
//! signs. Source amplitudes are local; no distant launch attenuation is inferred.
//!
//! Two-source superposition and absorbing-port Poynting accounting follow
//! S. J. Byrnes, "Multilayer optical calculations", §§3.1, 4.3–4.5 and
//! Appendix B, <https://arxiv.org/pdf/1603.02720>, with the parent's conjugate
//! time convention and explicit original source permeabilities.
use super::{
    AbsorbingLayerPoint, AbsorbingTwoSidedReflection, PlanarAbsorbingDielectricLayers,
    ReflectionError,
    kernel::two_sided::{PreparedTwoSided, TwoSidedIllumination},
};
use fs_math::c64::C64;

/// One additional O(films) bounded sweep of the original resolved stack.
/// Both exterior modes must carry positive normal power at unit amplitude.
/// Repeated source-phase/amplitude and depth queries allocate nothing.
pub struct PreparedAbsorbingTwoSidedLayers<'s, 'a> {
    source: &'s PlanarAbsorbingDielectricLayers<'a>,
    prepared: PreparedTwoSided<'s>,
}

/// Coherent actual fields and film heat for two prescribed face amplitudes.
pub struct AbsorbingTwoSidedIllumination<'p, 's, 'a> {
    prepared: &'p PreparedAbsorbingTwoSidedLayers<'s, 'a>,
    electric: [[C64; 2]; 2],
    solved: TwoSidedIllumination,
}

impl<'a> PlanarAbsorbingDielectricLayers<'a> {
    /// Prepare simultaneous illumination of both physical faces. The original
    /// one-sided owner still accepts evanescent/critical exits; they refuse here
    /// because such an exterior cannot supply an incoming positive-power wave.
    /// Material laws are not reevaluated, and all original sources remain bound.
    pub fn prepare_two_sided(
        &self,
    ) -> Result<PreparedAbsorbingTwoSidedLayers<'_, 'a>, ReflectionError> {
        Ok(PreparedAbsorbingTwoSidedLayers {
            source: self,
            prepared: PreparedTwoSided::new(&self.prepared)?,
        })
    }
}

impl<'s, 'a> PreparedAbsorbingTwoSidedLayers<'s, 'a> {
    pub fn source(&self) -> &'s PlanarAbsorbingDielectricLayers<'a> {
        self.source
    }
    /// Face (first,last), then channel (TM,TE), in A/V. These are the original
    /// complex physical admittances, not real power-normalization factors.
    pub fn port_admittances_s(&self) -> [[C64; 2]; 2] {
        self.prepared.admittances()
    }
    /// Channel (TM,TE), outgoing face (first,last), incoming face (first,last).
    /// Each coefficient maps actual peak fixed-axis tangential E to E.
    pub fn tangential_e_scattering(&self) -> [[[C64; 2]; 2]; 2] {
        self.prepared.scattering()
    }
    /// Actual complex peak [Ex,Ey] [V/m] at the first and last faces.
    /// Relative source phase survives in outgoing fields and positive film heat.
    /// Zero amplitude at either face is admitted; zero at both gives dark fields.
    pub fn illuminate(
        &self,
        first_face_e_v_m: [C64; 2],
        last_face_e_v_m: [C64; 2],
    ) -> Result<AbsorbingTwoSidedIllumination<'_, 's, 'a>, ReflectionError> {
        Ok(AbsorbingTwoSidedIllumination {
            prepared: self,
            electric: [first_face_e_v_m, last_face_e_v_m],
            solved: self
                .prepared
                .illuminate(first_face_e_v_m, last_face_e_v_m)?,
        })
    }
}

impl<'p, 's, 'a> AbsorbingTwoSidedIllumination<'p, 's, 'a> {
    pub fn source(&self) -> &'s PlanarAbsorbingDielectricLayers<'a> {
        self.prepared.source
    }
    pub fn response(&self) -> &AbsorbingTwoSidedReflection {
        &self.solved.response
    }
    /// Integrated dielectric heat in each original film [W/m²], in original
    /// source order. Thermal consumers own area/time and deposit each term once.
    pub fn absorbed_power_per_layer_w_m2(&self) -> &[f64] {
        &self.solved.absorbed[..self.prepared.source.layer_count()]
    }
    /// Original film depth in metres from its FIRST face. E/H and signed
    /// Poynting flux use the global +z direction, including the right source.
    pub fn layer_field(
        &self,
        layer_index: usize,
        depth_m: f64,
    ) -> Result<AbsorbingLayerPoint, ReflectionError> {
        self.prepared
            .prepared
            .field(self.electric[0], self.electric[1], layer_index, depth_m)
    }
}
