//! Authored scalar dissipation using the shared Damper and material shift owner.

use std::sync::Arc;

use fs_material::visco::temperature::TemperatureShift;
use fs_system::elements::DissipationLaw;

use super::{J, Refusal, SceneShift, fields, positive, scene_shift};

/// Read the optional law for a damper, rotary damper or scalar resistor.
/// The reference coefficient remains in the element's own effort/flow units.
pub(super) fn law(definition: &J, what: &str) -> Result<Option<Arc<dyn DissipationLaw>>, Refusal> {
    let Some(body) = definition.get("dissipation_law") else {
        return Ok(None);
    };
    let what = format!("{what}.dissipation_law");
    fields(body, &["quadratic_flow_scale", "shift"], &what)?;
    let flow_scale = if body.get("quadratic_flow_scale").is_some() {
        Some(positive(body, "quadratic_flow_scale", &what)?)
    } else {
        None
    };
    let shift = scene_shift(body, &what)?;
    if flow_scale.is_none() && shift.is_none() {
        return Err(format!("{what}: declare quadratic_flow_scale, shift, or both").into());
    }
    Ok(Some(Arc::new(SceneDissipation { flow_scale, shift })))
}

#[derive(Debug)]
struct SceneDissipation {
    flow_scale: Option<f64>,
    shift: Option<TemperatureShift>,
}

impl DissipationLaw for SceneDissipation {
    fn factor(&self, kelvin: f64, relative_flow: f64) -> Result<(f64, [f64; 2]), String> {
        let (thermal, thermal_slope) = if let Some(shift) = self.shift {
            let state = shift.evaluate(kelvin).map_err(|error| error.to_string())?;
            (
                state.factor(),
                state.factor() * state.log_derivative_per_k(),
            )
        } else {
            (1.0, 0.0)
        };
        let (flow, flow_slope) = if let Some(scale) = self.flow_scale {
            // |u| is not differentiable at zero, but u |u| and its heat
            // u^2 |u| have continuous derivatives there. Zero is their
            // correct product-rule choice at the origin.
            let slope = if relative_flow == 0.0 {
                0.0
            } else {
                relative_flow.signum() / scale
            };
            (1.0 + relative_flow.abs() / scale, slope)
        } else {
            (1.0, 0.0)
        };
        Ok((thermal * flow, [thermal_slope * flow, thermal * flow_slope]))
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        let mut identity =
            fs_blake3::DomainHasher::new("org.frankensim.fs-cli.system-scene.dissipation-law.v1");
        match self.flow_scale {
            Some(scale) => {
                identity.update(&[1]);
                identity.update(&scale.to_bits().to_le_bytes());
            }
            None => {
                identity.update(&[0]);
            }
        }
        match self.shift {
            Some(shift) => {
                identity.update(&[1]);
                identity.update(&SceneShift(shift).identity());
            }
            None => {
                identity.update(&[0]);
            }
        }
        for chunk in identity.finalize().as_bytes().chunks_exact(8) {
            words.push(u64::from_le_bytes(
                chunk.try_into().expect("eight hash bytes"),
            ));
        }
    }
}
