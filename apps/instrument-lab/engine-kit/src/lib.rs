//! Browser binding for the live drum kit: one `KitPiece` per kit piece, each the fs-vkplate
//! solver on a coarser grid (see `fs_vkplate::live`). A Worker owns each piece and steps it ahead
//! of playback; the page never post-processes its samples.

pub use fs_vkplate::live::{Live, PIECES};
pub use fs_vkplate::materials::KitMaterials;

/// Preset drum heads keep their recording-fitted loss in the live kit until a
/// sourced film or edge loss exists, and say so (ruling 6533).
pub const FITTED_HEAD_LOSS_RECEIPT: &str = "head loss: FITTED (declared estimate; film and edge loss not yet sourced)";

/// The head-loss receipt of one live piece, for drums only (cymbal and hi-hat
/// loss is pexfp.372): FITTED on the preset heads, the stated-unknown receipt
/// on a card route that declared it, and none when the cards source the loss.
pub fn head_loss_receipt(piece: &str, materials: Option<&KitMaterials>) -> Option<&'static str> {
    if !matches!(piece, "kick" | "snare" | "snare-off" | "tom-12" | "tom-14") {
        return None;
    }
    match materials {
        None => Some(FITTED_HEAD_LOSS_RECEIPT),
        Some(m) if m.film_loss_is_stated_unknown() => Some(fs_vkplate::materials::FILM_LOSS_STATED_UNKNOWN_RECEIPT),
        Some(_) => None,
    }
}

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
        /// One shared pair of hi-hat cymbals, settled under the declared downward
        /// initial pedal force [N]. Open, closed and foot events can share its motion.
        #[wasm_bindgen(js_name = withHiHatPedal)]
        pub fn with_hihat_pedal(grid_scale: f64, initial_pedal_n: f64) -> Result<KitPiece, JsError> {
            Live::with_hihat_pedal(grid_scale, initial_pedal_n)
                .map(|live| KitPiece { live, scale: grid_scale })
                .map_err(|e| JsError::new(&e))
        }
        /// Use the shared material gas grammar: `T P RH`, or
        /// `mixture T P mole|mass COUNT SPECIES FRACTION ...` (SI units).
        #[wasm_bindgen(js_name = withGas)]
        pub fn with_gas(piece: &str, grid_scale: f64, gas_record: &str) -> Result<KitPiece, JsError> {
            Live::with_gas_record(piece, grid_scale, gas_record)
                .map(|live| KitPiece { live, scale: grid_scale })
                .map_err(|e| JsError::new(&e))
        }
        /// Canonical species, NASA-9 and transport material packs. The packet
        /// is shared by all native physics consumers; no browser gas law.
        #[wasm_bindgen(js_name = withGasSources)]
        pub fn with_gas_sources(piece: &str, grid_scale: f64, gas_record: &str,
            bytes: &[u8], lengths: &[u32], water_species: &str, phase: &str) -> Result<KitPiece, JsError> {
            Live::with_gas_source_packet(piece, grid_scale, gas_record, bytes, lengths,
                water_species, phase, None)
                .map(|live| KitPiece { live, scale: grid_scale })
                .map_err(|e| JsError::new(&e))
        }
        /// Canonical material-card packs for both surfaces and a spherical
        /// impactor. Mass, elasticity and Hertz stiffness follow the resolved
        /// data; preset support/rim/wire laws remain independent estimates, and
        /// cymbals keep their preset T60. Drum heads refuse here: their film
        /// loss needs a sourced or stated-unknown route below.
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
        /// Material cards whose head films make no loss claim, with that loss
        /// declared a stated unknown entered as 0 (`materials::FILM_LOSS_STATED_UNKNOWN_RECEIPT`
        /// rides in `info`). Without this or a sourced bending loss, drums refuse.
        #[wasm_bindgen(js_name = withMaterialsStatedUnknownFilmLoss)]
        pub fn with_materials_stated_unknown_film_loss(piece: &str, grid_scale: f64, gas_record: &str,
            surface_pack: &[u8], secondary_pack: &[u8], striker_pack: &[u8],
            temperature_axis: &str, temperature_k: f64, striker_volume_m3: f64,
            striker_tip_radius_m: f64, impact_loss_s_m: f64) -> Result<KitPiece, JsError> {
            let materials = KitMaterials::from_packs_at_temperature(surface_pack, secondary_pack, striker_pack,
                temperature_axis, temperature_k, striker_volume_m3, striker_tip_radius_m, impact_loss_s_m)
                .and_then(KitMaterials::with_stated_unknown_film_loss)
                .map_err(|e| JsError::new(&e))?;
            Live::with_materials_record(piece, grid_scale, gas_record, materials)
                .map(|live| KitPiece { live, scale: grid_scale })
                .map_err(|e| JsError::new(&e))
        }
        /// The same declaration with an independently sourced atmosphere.
        #[wasm_bindgen(js_name = withMaterialsStatedUnknownFilmLossAndGasSources)]
        pub fn with_materials_stated_unknown_film_loss_and_gas_sources(piece: &str, grid_scale: f64, gas_record: &str,
            surface_pack: &[u8], secondary_pack: &[u8], striker_pack: &[u8],
            temperature_axis: &str, temperature_k: f64, striker_volume_m3: f64,
            striker_tip_radius_m: f64, impact_loss_s_m: f64,
            gas_bytes: &[u8], gas_lengths: &[u32], water_species: &str, phase: &str) -> Result<KitPiece, JsError> {
            let materials = KitMaterials::from_packs_at_temperature(surface_pack, secondary_pack, striker_pack,
                temperature_axis, temperature_k, striker_volume_m3, striker_tip_radius_m, impact_loss_s_m)
                .and_then(KitMaterials::with_stated_unknown_film_loss)
                .map_err(|e| JsError::new(&e))?;
            Live::with_gas_source_packet(piece, grid_scale, gas_record, gas_bytes, gas_lengths,
                water_species, phase, Some(materials))
                .map(|live| KitPiece { live, scale: grid_scale })
                .map_err(|e| JsError::new(&e))
        }
        /// Independently sourced solids and atmosphere in the same trajectory.
        #[wasm_bindgen(js_name = withMaterialsAndGasSources)]
        pub fn with_materials_and_gas_sources(piece: &str, grid_scale: f64, gas_record: &str,
            surface_pack: &[u8], secondary_pack: &[u8], striker_pack: &[u8],
            temperature_axis: &str, temperature_k: f64, striker_volume_m3: f64,
            striker_tip_radius_m: f64, impact_loss_s_m: f64,
            gas_bytes: &[u8], gas_lengths: &[u32], water_species: &str, phase: &str) -> Result<KitPiece, JsError> {
            let materials = KitMaterials::from_packs_at_temperature(surface_pack, secondary_pack, striker_pack,
                temperature_axis, temperature_k, striker_volume_m3, striker_tip_radius_m, impact_loss_s_m)
                .map_err(|e| JsError::new(&e))?;
            Live::with_gas_source_packet(piece, grid_scale, gas_record, gas_bytes, gas_lengths,
                water_species, phase, Some(materials))
                .map(|live| KitPiece { live, scale: grid_scale })
                .map_err(|e| JsError::new(&e))
        }
        /// Source viscosity and source atmosphere compose without restoring
        /// preset structural damping or replacing either material selection.
        #[wasm_bindgen(js_name = withMaterialBendingLossAndGasSources)]
        pub fn with_material_bending_loss_and_gas_sources(piece: &str, grid_scale: f64, gas_record: &str,
            surface_pack: &[u8], secondary_pack: &[u8], striker_pack: &[u8],
            temperature_axis: &str, temperature_k: f64, striker_volume_m3: f64,
            striker_tip_radius_m: f64, impact_loss_s_m: f64, frequency_axis: &str,
            frequency_hz: f64, primary_angle_rad: f64, secondary_angle_rad: f64,
            gas_bytes: &[u8], gas_lengths: &[u32], water_species: &str, phase: &str) -> Result<KitPiece, JsError> {
            let materials = KitMaterials::from_packs_with_bending_loss(surface_pack, secondary_pack, striker_pack,
                temperature_axis, temperature_k, striker_volume_m3, striker_tip_radius_m, impact_loss_s_m,
                frequency_axis, frequency_hz, [primary_angle_rad, secondary_angle_rad])
                .map_err(|e| JsError::new(&e))?;
            Live::with_gas_source_packet(piece, grid_scale, gas_record, gas_bytes, gas_lengths,
                water_species, phase, Some(materials))
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
        /// Downward pedal force [N], with an explicit transition duration [s].
        /// This moves the shared hi-hat rod; it does not truncate emitted audio.
        #[wasm_bindgen(js_name = setPedalForce)]
        pub fn set_pedal_force(&mut self, force_n: f64, transition_seconds: f64) -> Result<(), JsError> {
            self.live.set_hihat_pedal(force_n, transition_seconds).map_err(|e| JsError::new(&e))
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
        /// JSON: name, grid scale, surface radius and node count, every body's name and nodes,
        /// and a drum's head-loss receipt (`headLoss`, null for cymbals and sourced loss).
        pub fn info(&self) -> String {
            let b: Vec<String> = self.live.ins.bodies.iter().map(|b| format!("{{\"name\":\"{}\",\"nodes\":{}}}", b.name, b.n)).collect();
            format!(
                "{{\"piece\":\"{}\",\"gridScale\":{},\"radius\":{},\"surfaceNodes\":{},\"h\":{},\"rate\":{},\"gas\":{{\"temperatureK\":{},\"pressurePa\":{},\"densityKgM3\":{},\"soundSpeedMps\":{},\"sourceSpecies\":{}}},\"materials\":{},\"headLoss\":{},\"bodies\":[{}]}}",
                self.live.name, self.scale, self.live.radius, self.live.surface_xy().len() / 2, self.live.surface_h(), self.live.ins.rate,
                self.live.ins.gas().temperature, self.live.ins.gas().pressure,
                self.live.ins.gas().density, self.live.ins.gas().sound_speed,
                self.live.gas_source_species(),
                self.live.ins.materials().map_or_else(|| "null".into(), KitMaterials::info_json),
                // Both receipts are fixed ASCII constants without quotes or backslashes.
                head_loss_receipt(&self.live.name, self.live.ins.materials()).map_or_else(|| "null".into(), |r| format!("\"{r}\"")),
                b.join(",")
            )
        }
    }

    /// The kit's piece names, JSON array.
    #[wasm_bindgen]
    pub fn pieces() -> String {
        format!("[{}]", PIECES.iter().map(|p| format!("\"{p}\"")).collect::<Vec<_>>().join(","))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEAD_FILM: &[u8] =
        include_bytes!("../../../../crates/fs-vkplate/packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
    const HICKORY: &[u8] = include_bytes!("../../../../crates/fs-vkplate/packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk");

    #[test]
    fn preset_drums_label_their_fitted_head_loss_and_a_card_route_states_its_own() {
        for piece in PIECES {
            let drum = !matches!(piece, "hihat" | "hihat-closed" | "hihat-open" | "crash" | "ride");
            assert_eq!(head_loss_receipt(piece, None), drum.then_some(FITTED_HEAD_LOSS_RECEIPT), "{piece}");
        }
        let cards = KitMaterials::from_packs_at_temperature(HEAD_FILM, &[], HICKORY, "T", 293.15, 9e-6, 0.005, 0.05).unwrap();
        assert_eq!(head_loss_receipt("tom-12", Some(&cards)), None);
        let declared = cards.with_stated_unknown_film_loss().unwrap();
        assert_eq!(
            head_loss_receipt("tom-12", Some(&declared)),
            Some(fs_vkplate::materials::FILM_LOSS_STATED_UNKNOWN_RECEIPT)
        );
        assert_eq!(head_loss_receipt("crash", Some(&declared)), None);
    }
}
