//! Piano surface projection into the shared stationary Rayleigh receiver.
//! The bank still owns the loaded coordinate map; the shared owner alone owns
//! acoustic differentiation, filtering, geometry kernels and propagation memory.
use super::board_geometry::SurfaceSample;
use super::linear::Bank;
use fs_bem::helmholtz::Medium;
use fs_couple::pcm_wav::baffled::{BaffledPressure, RayleighMedium, SurfaceSample as AcousticSample};

pub struct Microphone { receiver: BaffledPressure }
impl std::ops::Deref for Microphone {
    type Target = BaffledPressure;
    fn deref(&self) -> &Self::Target { &self.receiver }
}
impl std::ops::DerefMut for Microphone {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.receiver }
}
impl Microphone {
    fn surface(surface: &[SurfaceSample], bank: &Bank) -> Result<Vec<AcousticSample>, String> {
        surface.iter().map(|p| Ok(AcousticSample { position_m:p.position_m, area_m2:p.area_m2,
            mode_shape:bank.project_board_shape(&p.mode_shape)? })).collect()
    }
    pub fn new(surface: &[SurfaceSample], bank: &Bank, rate:u32,
        position_m:[f64;3], medium:Medium) -> Result<Self,String> {
        let loaded=Self::surface(surface,bank)?;
        Ok(Self { receiver:BaffledPressure::new(&loaded,bank.board_count,rate,position_m,
            RayleighMedium {density:medium.density,sound_speed:medium.sound_speed})? })
    }
    pub fn new_multirate(surface:&[SurfaceSample], bank:&Bank, rate:u32,
        position_m:[f64;3], medium:Medium) -> Result<Self,String> {
        let loaded=Self::surface(surface,bank)?;
        Ok(Self { receiver:BaffledPressure::new_multirate(&loaded,bank.board_count,bank.rate,rate,position_m,
            RayleighMedium {density:medium.density,sound_speed:medium.sound_speed})? })
    }
}
