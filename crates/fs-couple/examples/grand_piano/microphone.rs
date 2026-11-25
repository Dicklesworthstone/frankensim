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
pub use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::baffled::{BaffledPressure, RayleighMedium, SurfaceSample as AcousticSample};

#[derive(Clone, Copy)]
pub struct Placement {
    pub position_m: [f64; 3],
    pub pattern: FirstOrder,
}

impl Placement {
    pub fn omni(position_m: [f64; 3]) -> Self { Self { position_m, pattern: FirstOrder::default() } }
}

pub fn parse_pattern(value: &str) -> Result<FirstOrder, String> {
    let invalid = "expected microphone pressure fraction and unit front axis: alpha,x,y,z";
    let mut values = value.split(',');
    let mut number = || values.next().and_then(|v| v.parse::<f64>().ok()).ok_or(invalid);
    let pattern = FirstOrder::new(number()?, [number()?, number()?, number()?]).map_err(|e| e.to_string())?;
    if values.next().is_some() { return Err(invalid.into()); }
    Ok(pattern)
}

pub fn pressure_unit(pattern: FirstOrder) -> &'static str {
    if pattern.pressure_fraction() == 1.0 { "pa" } else { "pa_equivalent" }
}

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
        Self::placed(surface, bank, rate, Placement::omni(position_m), medium)
    }
    pub fn placed(surface: &[SurfaceSample], bank: &Bank, rate: u32, placement: Placement,
        medium: Medium) -> Result<Self,String> {
        let loaded = Self::surface(surface,bank)?;
        Ok(Self { receiver: BaffledPressure::new_multirate_first_order(&loaded, bank.board_count,
            bank.rate, rate, placement.position_m,
            RayleighMedium {density:medium.density, sound_speed:medium.sound_speed}, placement.pattern)? })
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
