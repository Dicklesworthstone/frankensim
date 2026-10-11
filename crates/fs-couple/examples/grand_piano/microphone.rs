//! Piano surface projection into the shared stationary Rayleigh receiver.
//! The bank still owns the loaded coordinate map; the shared owner alone owns
//! acoustic differentiation, filtering, geometry kernels and propagation memory.
//! Bandlimited acceleration and travel-time reconstruction retain upper partial
//! amplitude and phase without two-sample differencing or linear-delay droop.
//! At 4x oversampling the declared numerical delay is 84 output samples;
//! geometry, signed spatial cancellation and physical pressure gain are unchanged.
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
        Ok(Self { receiver:BaffledPressure::new_differentiating(&loaded,bank.board_count,rate,position_m,
            RayleighMedium {density:medium.density,sound_speed:medium.sound_speed})? })
    }
    pub fn new_multirate(surface:&[SurfaceSample], bank:&Bank, rate:u32,
        position_m:[f64;3], medium:Medium) -> Result<Self,String> {
        let loaded=Self::surface(surface,bank)?;
        Ok(Self { receiver:BaffledPressure::new_multirate_differentiating(&loaded,bank.board_count,bank.rate,rate,position_m,
            RayleighMedium {density:medium.density,sound_speed:medium.sound_speed})? })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn played_piano_selects_bandlimited_acceleration_and_travel_without_changing_the_loaded_projection() {
        let course=super::super::geometry::demonstration_scale().unwrap()[48];
        let board=super::super::board::demonstration();
        let surface=[SurfaceSample {position_m:[0.0;3],area_m2:0.2,
            mode_shape:board.iter().map(|b|b.volume/0.2).collect()}];
        let piano=super::super::engine::Instrument::new(vec![course],&board,48_000,4,8,true).unwrap();
        let mic=Microphone::new_multirate(&surface,&piano.bank,48_000,[0.0,0.0,1.0],Medium::air()).unwrap();
        assert_eq!(mic.interpolation_delay_samples(),32);
        assert_eq!(mic.filter_delay_samples(),84.0);
        let legacy=Microphone::surface(&surface,&piano.bank).unwrap();
        assert_eq!(legacy[0].mode_shape,piano.bank.project_board_shape(&surface[0].mode_shape).unwrap());
    }
}
