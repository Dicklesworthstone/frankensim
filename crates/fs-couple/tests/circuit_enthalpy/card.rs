//! The fixtures' resistance law for the declared-graph comparisons of the
//! circuit's and the moving field's couplings: a resistor follows its card
//! as the module's `coefficients` resolves it.

use fs_couple::cavity_phs::piston::thermal::power::joule::MaterialOhmicBody;
use fs_couple::cavity_phs::piston::thermal::source::ThermalMaterialPoint;
use fs_matdb::MaterialCard;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::MaterialPropertySelection;
use fs_system::circuit::ResistanceLaw;

use crate::common::point;
use crate::kernel::enthalpy_at;

/// The fixture's resistor, resolved as the module's `coefficients` resolves
/// it: the card's conductor at the phase state of the cell, here the state
/// at the temperature the resistor's heat port receives, its geometry
/// checked against the cell's volume.
#[derive(Debug)]
pub struct CardResistance {
    pub card: &'static MaterialCard,
    pub curve: &'static ResolvedHeatCapacityCurve,
    /// The destination cell's volume \[m^3\].
    pub volume_m3: f64,
    /// Names the card in the circuit's identity: the fixture's index.
    pub tag: u64,
}

impl ResistanceLaw for CardResistance {
    fn ohms(&self, kelvin: f64) -> Result<f64, String> {
        let h = enthalpy_at(self.curve, kelvin)?;
        let phase = self
            .curve
            .curve()
            .state_at_specific_enthalpy(h)
            .map_err(|error| format!("{error:?}"))?;
        let context = point(300.0);
        let body = MaterialOhmicBody {
            material: ThermalMaterialPoint {
                material: self.card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            selection: MaterialPropertySelection::SingleClaimOnly,
            length_m: 1.0,
            area_m2: 1.0,
        };
        body.resolve_conductor(phase, self.volume_m3)
            .map(|conductor| conductor.resistance_ohm())
            .map_err(|error| format!("{error:?}"))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.extend([self.tag, self.volume_m3.to_bits()]);
    }
}
