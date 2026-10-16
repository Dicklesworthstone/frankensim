//! Browser binding for the live drum kit: one `KitPiece` per kit piece, each the fs-vkplate
//! solver on a coarser grid (see `fs_vkplate::live`). A Worker owns each piece and steps it ahead
//! of playback; the page never post-processes its samples.

pub use fs_vkplate::live::{Live, PIECES};

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
                "{{\"piece\":\"{}\",\"gridScale\":{},\"radius\":{},\"surfaceNodes\":{},\"h\":{},\"rate\":{},\"bodies\":[{}]}}",
                self.live.name, self.scale, self.live.radius, self.live.surface_xy().len() / 2, self.live.surface_h(), self.live.ins.rate, b.join(",")
            )
        }
    }

    /// The kit's piece names, JSON array.
    #[wasm_bindgen]
    pub fn pieces() -> String {
        format!("[{}]", PIECES.iter().map(|p| format!("\"{p}\"")).collect::<Vec<_>>().join(","))
    }
}
