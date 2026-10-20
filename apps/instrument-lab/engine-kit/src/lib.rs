//! Browser binding for the live drum kit: one `KitPiece` per kit piece, each the fs-vkplate
//! solver on a coarser grid (see `fs_vkplate::live`). A Worker owns each piece and steps it ahead
//! of playback; the page never post-processes its samples.

pub use fs_vkplate::live::{Live, PIECES};
pub use fs_vkplate::materials::KitMaterials;

#[cfg(target_arch = "wasm32")]
mod wasm {
    use super::*;
    use wasm_bindgen::prelude::*;

    #[wasm_bindgen]
    pub struct KitPiece {
        live: Live,
        scale: f64,
    }

    #[wasm_bindgen]
    impl KitPiece {
        /// `piece` is one of `pieces()`; `grid_scale` coarsens the offline grid (5 runs in real time).
        #[wasm_bindgen(constructor)]
        pub fn new(piece: &str, grid_scale: f64) -> Result<KitPiece, JsError> {
            Live::new(piece, grid_scale).map(|live| KitPiece { live, scale: grid_scale }).map_err(|e| JsError::new(&e))
        }
        /// Use the shared material gas grammar: `T P RH`, or
        /// `mixture T P mole|mass COUNT SPECIES FRACTION ...` (SI units).
        #[wasm_bindgen(js_name = withGas)]
        pub fn with_gas(piece: &str, grid_scale: f64, gas_record: &str) -> Result<KitPiece, JsError> {
            Live::with_gas_record(piece, grid_scale, gas_record)
                .map(|live| KitPiece { live, scale: grid_scale })
                .map_err(|e| JsError::new(&e))
        }
        /// Canonical material-card packs for both surfaces and a spherical
        /// impactor. Mass, elasticity and Hertz stiffness follow the resolved
        /// data; preset T60/support/rim/wire laws remain independent estimates.
        #[wasm_bindgen(js_name = withMaterials)]
        pub fn with_materials(piece: &str, grid_scale: f64, gas_record: &str,
            surface_pack: &[u8], secondary_pack: &[u8], striker_pack: &[u8],
            temperature_axis: &str, temperature_k: f64, striker_volume_m3: f64,
            striker_tip_radius_m: f64, impact_loss_s_m: f64) -> Result<KitPiece, JsError> {
            let materials = KitMaterials::from_packs_at_temperature(surface_pack, secondary_pack, striker_pack,
                temperature_axis, temperature_k, striker_volume_m3, striker_tip_radius_m, impact_loss_s_m)
                .map_err(|e| JsError::new(&e))?;
            Live::with_materials_record(piece, grid_scale, gas_record, materials)
                .map(|live| KitPiece { live, scale: grid_scale })
                .map_err(|e| JsError::new(&e))
        }
        /// Explicit source-curvature-loss route. Its distinct export prevents
        /// an older WASM binary from silently ignoring the loss selection.
        #[wasm_bindgen(js_name = withMaterialBendingLoss)]
        pub fn with_material_bending_loss(piece: &str, grid_scale: f64, gas_record: &str,
            surface_pack: &[u8], secondary_pack: &[u8], striker_pack: &[u8],
            temperature_axis: &str, temperature_k: f64, striker_volume_m3: f64,
            striker_tip_radius_m: f64, impact_loss_s_m: f64, frequency_axis: &str,
            frequency_hz: f64, primary_angle_rad: f64, secondary_angle_rad: f64) -> Result<KitPiece, JsError> {
            let materials = KitMaterials::from_packs_with_bending_loss(surface_pack, secondary_pack, striker_pack,
                temperature_axis, temperature_k, striker_volume_m3, striker_tip_radius_m, impact_loss_s_m,
                frequency_axis, frequency_hz, [primary_angle_rad, secondary_angle_rad])
                .map_err(|e| JsError::new(&e))?;
            Live::with_materials_record(piece, grid_scale, gas_record, materials)
                .map(|live| KitPiece { live, scale: grid_scale })
                .map_err(|e| JsError::new(&e))
        }
        /// A second piece in exactly this one's state: it renders the same samples as this one would.
        /// Building a snare or a hi-hat is almost all its settle, so a copy is far cheaper than `new`.
        pub fn copy(&self) -> KitPiece {
            KitPiece { live: self.live.clone(), scale: self.scale }
        }
        /// Strike at (x, y) metres from the centre of the playing surface at `speed` m/s.
        pub fn strike(&mut self, x: f64, y: f64, speed: f64) {
            self.live.strike(x, y, speed);
        }
        /// The next `frames` samples at 48 kHz: microphone pressure [Pa] times `gain`.
        pub fn render(&mut self, frames: usize, gain: f64) -> Vec<f32> {
            let mut v = vec![0.0f32; frames];
            self.live.render(&mut v, gain);
            v
        }
        pub fn asleep(&self) -> bool {
            self.live.asleep()
        }
        /// Material bending loss only; contact/support/gas terms are separate.
        #[wasm_bindgen(js_name = bendingLossJ)]
        pub fn bending_loss_j(&self) -> f64 {
            self.live.ins.bodies.iter().map(|b| b.bending_loss_dissipated_j()).sum()
        }
        pub fn radius(&self) -> f64 {
            self.live.radius
        }
        #[wasm_bindgen(js_name = surfaceXY)]
        pub fn surface_xy(&self) -> Vec<f32> {
            self.live.surface_xy()
        }
        #[wasm_bindgen(js_name = surfaceH)]
        pub fn surface_h(&self) -> f64 {
            self.live.surface_h()
        }
        /// Transverse velocity of each surface node [m/s] (same order as `surfaceXY`).
        pub fn velocity(&self) -> Vec<f32> {
            self.live.surface_velocity()
        }
        pub fn displacement(&self) -> Vec<f32> {
            self.live.surface_displacement()
        }
        /// Kinetic energy of the piece after the last block [J].
        pub fn kinetic(&self) -> f64 {
            self.live.kinetic
        }
        /// Largest stick (or beater) force during the last block [N].
        pub fn force(&self) -> f64 {
            self.live.force
        }
        /// JSON: name, grid scale, surface radius and node count, every body's name and nodes.
        pub fn info(&self) -> String {
            let b: Vec<String> = self.live.ins.bodies.iter().map(|b| format!("{{\"name\":\"{}\",\"nodes\":{}}}", b.name, b.n)).collect();
            format!(
                "{{\"piece\":\"{}\",\"gridScale\":{},\"radius\":{},\"surfaceNodes\":{},\"h\":{},\"rate\":{},\"gas\":{{\"temperatureK\":{},\"pressurePa\":{},\"densityKgM3\":{},\"soundSpeedMps\":{}}},\"materials\":{},\"bodies\":[{}]}}",
                self.live.name, self.scale, self.live.radius, self.live.surface_xy().len() / 2, self.live.surface_h(), self.live.ins.rate,
                self.live.ins.gas().temperature, self.live.ins.gas().pressure,
                self.live.ins.gas().density, self.live.ins.gas().sound_speed,
                self.live.ins.materials().map_or_else(|| "null".into(), KitMaterials::info_json), b.join(",")
            )
        }
    }

    /// The kit's piece names, JSON array.
    #[wasm_bindgen]
    pub fn pieces() -> String {
        format!("[{}]", PIECES.iter().map(|p| format!("\"{p}\"")).collect::<Vec<_>>().join(","))
    }
}
