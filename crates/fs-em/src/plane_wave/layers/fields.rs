//! Fields inside a selected original film, including critical and opaque limits.
//!
//! The same bounded slab/state recursion supplies both boundary traces. The
//! interior uses waves decaying from opposite faces, or an entire E/H
//! propagator near kz=0. No division by a vanishing normal wave number and no
//! growing evanescent transfer matrix is needed at a critical film.
//! Depth-cell heat reuses the same positive slab loss modes between its two
//! field traces. It needs neither depth quadrature nor a flux subtraction.

use super::*;

pub(super) mod kernel;
use kernel::LayerFieldKernel;

mod preparation;
use preparation::LayerFieldPreparation;

mod plasma;
pub use plasma::{
    IonizedGasLayerExcitation, IonizedGasLayerHeating, IonizedGasLayerPoint,
    PreparedIonizedGasLayerIllumination,
};

pub(in crate::plane_wave) const PREPARATION_WORKSPACE_BYTES: usize =
    core::mem::size_of::<LayerFieldPreparation<'static, 'static>>();

/// Total peak tangential fields for ONE V/m of incident tangential electric
/// field. Multiply both by the actual complex incident amplitude. The H axis
/// is power-conjugate: -H_x for TE (E_y), H_y for TM (E_x), with +z into the
/// stack. Normal flux is Re(E conj(H))/2, including standing-wave terms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerTangentialField {
    pub electric_v_per_m: C64,
    pub conjugate_magnetic_a_per_m: C64,
}

/// One source-bound selected film and polarization. Preparation is O(layers)
/// with fixed stack storage; repeated depth queries allocate nothing and do
/// not reevaluate a material law or solve the stack again.
pub struct PreparedLayerField<'s, 'm> {
    source: &'s PlanarDielectricLayers<'m>,
    layer_index: usize,
    polarization: Polarization,
    incidence_cosine: f64,
    admittance_s: f64,
    kernel: LayerFieldKernel,
    absorbed_power_w_m2: f64,
}

impl<'m> PlanarDielectricLayers<'m> {
    /// Prepare the field of a unit tangential-E plane wave incident from the
    /// first face, using the original exp(+i omega t) convention. Lossy and
    /// evanescent exits remain admitted; zero-thickness films are identities.
    pub fn prepare_layer_field(
        &self,
        layer_index: usize,
        incidence_cosine: f64,
        polarization: Polarization,
    ) -> Result<PreparedLayerField<'_, 'm>, ReflectionError> {
        if self.layers.get(layer_index).is_none() {
            return Err(bad("selected field layer is outside the original stack"));
        }
        LayerFieldPreparation::new(self, incidence_cosine, polarization)?.field(layer_index)
    }
}

impl<'s, 'm> PreparedLayerField<'s, 'm> {
    pub fn source(&self) -> &'s PlanarDielectricLayers<'m> {
        self.source
    }
    pub fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub fn polarization(&self) -> Polarization {
        self.polarization
    }
    /// Physical heat integrated through this film for unit incident peak E.
    pub fn absorbed_power_w_m2(&self) -> f64 {
        self.absorbed_power_w_m2
    }
    /// Physical cycle-mean heat in a depth interval [W/m²] for unit incident
    /// peak tangential E. Both depths are measured from this film's left face;
    /// 0 <= start <= end <= thickness. An empty interval has exactly zero heat.
    ///
    /// The original slab's positive dissipative modes integrate the local
    /// dielectric loss exactly, including standing waves and longitudinal TM
    /// fields. No depth mesh, flux subtraction, new material evaluation or
    /// stack solve is needed. The query is allocation-free and O(1). Multiply
    /// by the actual incident peak amplitude squared. A thermal consumer owns
    /// the area/time conversion and deposits this energy only once.
    pub fn absorbed_power_in_depth_interval_w_m2(
        &self,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<f64, ReflectionError> {
        let layer = self.source.layers[self.layer_index];
        let width = admit_depth_interval(layer.source.thickness_m, start_depth_m, end_depth_m)?;
        if width == 0.0 || layer.wave.epsilon_f_m.im == 0.0 {
            return Ok(0.0);
        }
        if start_depth_m == 0.0 && end_depth_m == layer.source.thickness_m {
            return Ok(self.absorbed_power_w_m2);
        }
        let slab = Slab::new_with_phase_thickness(
            layer.wave,
            self.source.wave_number,
            self.incidence_cosine,
            self.polarization,
            width,
            Dd::from_f64(end_depth_m) - Dd::from_f64(start_depth_m),
        )?;
        let trace = |depth| {
            self.kernel
                .field(depth)
                .ok_or_else(|| bad("depth-interval field trace is unrepresentable"))
        };
        let modes = interval_loss_transfer(
            slab,
            self.admittance_s,
            trace(start_depth_m)?,
            trace(end_depth_m)?,
        )?;
        let heat = modes[0].norm_sq() + modes[1].norm_sq();
        if !heat.is_finite() {
            return Err(bad("depth-interval absorbed power is unrepresentable"));
        }
        Ok(heat)
    }
    /// Depth in metres from this film's incident face, including both faces.
    pub fn field(&self, depth_m: f64) -> Result<LayerTangentialField, ReflectionError> {
        let [electric_v_per_m, conjugate_magnetic_a_per_m] = self
            .kernel
            .field(depth_m)
            .ok_or_else(|| bad("field depth is outside the layer or is unrepresentable"))?;
        Ok(LayerTangentialField {
            electric_v_per_m,
            conjugate_magnetic_a_per_m,
        })
    }
}

pub(super) fn admit_depth_interval(
    thickness_m: f64,
    start_depth_m: f64,
    end_depth_m: f64,
) -> Result<f64, ReflectionError> {
    if !start_depth_m.is_finite()
        || !end_depth_m.is_finite()
        || start_depth_m < 0.0
        || end_depth_m < start_depth_m
        || end_depth_m > thickness_m
    {
        return Err(bad(
            "depth interval must be ordered and inside the original film",
        ));
    }
    Ok(end_depth_m - start_depth_m)
}

// Restrict the already solved film to these two traces. The artificial cut
// faces use the SAME real reference admittance as every original slab; they
// add no material interface. Contract incoming reference waves through the
// positive slab eigenmodes, preserving weak heat without subtracting fluxes.
pub(super) fn interval_loss_transfer(
    slab: Slab,
    admittance_s: f64,
    left: [C64; 2],
    right: [C64; 2],
) -> Result<[C64; 2], ReflectionError> {
    let incoming = (left[0] + divided(left[1], admittance_s)).scale(0.5);
    let opposite = (right[0] - divided(right[1], admittance_s)).scale(0.5);
    let root_reference = 0.5 * det::sqrt(admittance_s);
    let modes = [
        (incoming + opposite).scale(root_reference * det::sqrt(slab.even_loss)),
        (incoming - opposite).scale(root_reference * det::sqrt(slab.odd_loss)),
    ];
    if modes.iter().any(|&z| !finite(z)) {
        return Err(bad("depth-interval loss transfer is unrepresentable"));
    }
    Ok(modes)
}
