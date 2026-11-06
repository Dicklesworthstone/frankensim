//! Amplitude-admitted coherent films containing original collisional gases.
//!
//! The native Maxwell layer solution supplies both field traces and positive
//! integrated loss. An analytic whole-film electric envelope is admitted by the
//! original carrier owner before any wave is published. This retains standing
//! waves, longitudinal TM fields and all pair-relative drift limits; testing
//! only an incoming amplitude or a few depths would not establish that domain.

use super::*;

mod one_sided;
mod thermal;
use fs_material::{
    field_storage::linear_electric_permittivity,
    ionized_gas::{IonizedGasError, IonizedGasHarmonicResponse},
};

/// Simultaneous complex peak tangential E sources at the two original faces.
/// Both share the native exp(+i omega t) clock and transverse wavevector.
#[derive(Clone, Copy, Debug)]
pub struct IonizedGasLayerExcitation {
    pub incidence_cosine: f64,
    pub polarization: Polarization,
    pub incoming_tangential_e_v_m: [C64; 2],
}

/// Original carrier friction integrated through one finite film, in W/m².
/// The signed exchange correction is already included in electron/heavy heat;
/// zero-drift thermal exchange is a separate material process.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IonizedGasLayerHeating {
    pub total_w_m2: f64,
    pub electron_w_m2: f64,
    pub heavy_w_m2: f64,
    pub exchange_correction_to_electrons_w_m2: f64,
}

struct AdmittedLayer<'m> {
    source: DielectricLayer<'m>,
    wave: WaveMedium,
    kernel: Option<LayerFieldKernel>,
    longitudinal: C64,
    root_loss: f64,
    permeability_h_m: f64,
    background_permittivity_f_m: Option<f64>,
    electric_envelope_v_m: Option<f64>,
    carriers: Option<IonizedGasHarmonicResponse<'m>>,
    absorbed_w_m2: f64,
    carrier_heating: Option<IonizedGasLayerHeating>,
}

/// One admitted finite-layer wave with fixed original physical face sources.
///
/// Simultaneous sources require two lossless propagating exterior ports. The
/// one-sided constructor also admits an absorbing or evanescent exit. Finite
/// films retain their original ionized gas, Debye, Ohmic or optical owners. Every
/// positive-thickness gas film is admitted at a conservative analytic bound
/// on its full vector field, including coherent counterpropagating sources.
/// A bound may reject an actually admissible field; it never substitutes a
/// depth grid for the original weak-drift or Coulomb-relative-drift domain.
///
/// Preparation is bounded by MAX_LAYERS and each original carrier budget.
/// Repeated depth/heat/carrier queries are allocation-free and do not resolve
/// the stack or carrier modes. The unrestricted normalized stack is private
/// to preparation, so this owner cannot apply a new unadmitted source amplitude.
/// Frozen populations, local unmagnetized collisions and the original source
/// applicability remain. No sheath, ionization evolution, transient reservoir,
/// arbitrary angular beam or measured broadband accuracy is supplied.
pub struct PreparedIonizedGasLayerIllumination<'m> {
    incident: DielectricResponse<'m>,
    exit: DielectricResponse<'m>,
    omega_rad_s: f64,
    wave_number_rad_m: f64,
    reference_admittance_s: f64,
    excitation: IonizedGasLayerExcitation,
    boundary: CoherentLayerIllumination,
    layers: Vec<AdmittedLayer<'m>>,
    one_sided_exit: Option<one_sided::ExitWave>,
}

/// Native local Maxwell fields plus the original free-carrier response.
/// All depths are measured from the selected film's original left face.
pub struct IonizedGasLayerPoint<'w, 'm> {
    wave: &'w PreparedIonizedGasLayerIllumination<'m>,
    layer_index: usize,
    depth_m: f64,
    envelope_ratio: f64,
    background_energy_j_m3: Option<f64>,
    /// Peak E/H, their RMS squares, signed Poynting flux and local optical loss.
    fields: LayerPointField,
    /// Free-carrier current only; zero in a film with no ionized-gas owner.
    pub current_density_a_m2: [C64; 3],
}

impl<'m> PreparedIonizedGasLayerIllumination<'m> {
    /// Prepare the actual simultaneous face sources. A zero right input is
    /// ordinary one-sided illumination. Cancellation or any material-domain
    /// failure returns no partial wave. Zero-thickness films retain their
    /// source and identity scattering but have no volume-field queries.
    pub fn try_new(
        incident: DielectricResponse<'m>,
        exit: DielectricResponse<'m>,
        omega_rad_s: f64,
        layers: &[DielectricLayer<'m>],
        excitation: IonizedGasLayerExcitation,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReflectionError> {
        poll(&mut cancelled)?;
        if excitation
            .incoming_tangential_e_v_m
            .iter()
            .any(|&e| !finite(e))
        {
            return Err(bad("plasma layer sources must be finite peak fields"));
        }
        let stack =
            PlanarDielectricLayers::try_new_for_fields(incident, exit, omega_rad_s, layers, true)?;
        poll(&mut cancelled)?;
        let c = excitation.incidence_cosine;
        let polarization = excitation.polarization;
        let illumination = stack.prepare_illumination(c, polarization)?;
        let mut heat = [0.0; MAX_LAYERS];
        let boundary = illumination.evaluate(
            excitation.incoming_tangential_e_v_m,
            &mut heat[..layers.len()],
        )?;
        let parallel_over_omega = det::sqrt(stack.incident_epsilon_f_m)
            * (det::sqrt(stack.incident_mu_h_m) * det::sqrt((1.0 - c) * (1.0 + c)));
        if !parallel_over_omega.is_finite() {
            return Err(bad(
                "plasma layer transverse propagation is unrepresentable",
            ));
        }
        let mut admitted = Vec::new();
        admitted
            .try_reserve_exact(layers.len())
            .map_err(|_| bad("plasma layer field storage is unavailable"))?;
        for (index, selected) in stack.layers.iter().enumerate() {
            poll(&mut cancelled)?;
            let thickness = selected.source.thickness_m;
            let kernel = if thickness == 0.0 {
                None
            } else {
                let trace = |depth| -> Result<[C64; 2], ReflectionError> {
                    let point = illumination
                        .prepare_point(index, depth)?
                        .evaluate(excitation.incoming_tangential_e_v_m)?;
                    Ok(match polarization {
                        Polarization::Te => [point.electric_v_m[1], -point.magnetic_a_m[0]],
                        Polarization::Tm => [point.electric_v_m[0], point.magnetic_a_m[1]],
                    })
                };
                let mu = linear_magnetic_permeability(selected.source.response.material())?;
                let normal = selected.wave.normal(c)?;
                let kz = normal.scale(stack.wave_number);
                let omega_mu = positive(omega_rad_s * mu)?;
                let omega_epsilon = selected.wave.epsilon_f_m.scale(omega_rad_s);
                let (a, b) = match polarization {
                    Polarization::Te => (C64::from_re(omega_mu), (kz * kz).scale(1.0 / omega_mu)),
                    Polarization::Tm => ((kz * kz) / omega_epsilon, omega_epsilon),
                };
                Some(
                    LayerFieldKernel::new(
                        thickness,
                        kz,
                        Dd::from_f64(normal.re) * Dd::from_f64(stack.wave_number),
                        a,
                        b,
                        trace(0.0)?,
                        trace(thickness)?,
                    )
                    .ok_or_else(|| bad("plasma layer field is unrepresentable"))?,
                )
            };
            admitted.push(admit_layer(
                &stack,
                index,
                polarization,
                parallel_over_omega,
                kernel,
                heat[index],
                &mut cancelled,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            incident,
            exit,
            omega_rad_s,
            wave_number_rad_m: stack.wave_number,
            reference_admittance_s: illumination.port_admittance_s()[0],
            excitation,
            boundary,
            layers: admitted,
            one_sided_exit: None,
        })
    }

    pub const fn incident_response(&self) -> DielectricResponse<'m> {
        self.incident
    }
    pub const fn exit_response(&self) -> DielectricResponse<'m> {
        self.exit
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    pub const fn excitation(&self) -> IonizedGasLayerExcitation {
        self.excitation
    }
    pub const fn boundary(&self) -> &CoherentLayerIllumination {
        &self.boundary
    }
    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    pub fn layer_source(&self, index: usize) -> Option<DielectricLayer<'m>> {
        self.layers.get(index).map(|l| l.source)
    }
    /// Physical integrated film loss, including coherent source interference.
    pub fn layer_absorbed_power_w_m2(&self, index: usize) -> Option<f64> {
        self.layers.get(index).map(|l| l.absorbed_w_m2)
    }
    /// Original electron/heavy partitions; None for a non-gas or invalid film.
    pub fn layer_carrier_heating(&self, index: usize) -> Option<IonizedGasLayerHeating> {
        self.layers.get(index).and_then(|l| l.carrier_heating)
    }
    /// Conservative preparation envelope; None for absent gas volume.
    pub fn layer_electric_envelope_v_m(&self, index: usize) -> Option<f64> {
        self.layers.get(index).and_then(|l| l.electric_envelope_v_m)
    }

    /// Native vector field at x=y=0. The shared transverse phase at any other
    /// x multiplies every phasor equally and does not change drift/heat bounds.
    pub fn field(
        &self,
        layer_index: usize,
        depth_m: f64,
    ) -> Result<IonizedGasLayerPoint<'_, 'm>, ReflectionError> {
        let layer = self
            .layers
            .get(layer_index)
            .ok_or_else(|| bad("plasma film index is outside the original stack"))?;
        let [electric, magnetic] = layer
            .kernel
            .ok_or_else(|| bad("zero-thickness film has no interior volume"))?
            .field(depth_m)
            .ok_or_else(|| bad("plasma field depth is outside its film or unrepresentable"))?;
        let (e, h) = match self.excitation.polarization {
            Polarization::Te => (
                [C64::ZERO, electric, C64::ZERO],
                [-magnetic, C64::ZERO, layer.longitudinal * electric],
            ),
            Polarization::Tm => (
                [electric, C64::ZERO, layer.longitudinal * magnetic],
                [C64::ZERO, magnetic, C64::ZERO],
            ),
        };
        let fields = point_field(e, h, layer.root_loss)?;
        let norm = vector_norm(e);
        let (current, ratio) = if let Some(carriers) = &layer.carriers {
            let bound = layer.electric_envelope_v_m.unwrap();
            let ratio = if bound == 0.0 { 0.0 } else { norm / bound };
            if !ratio.is_finite() || norm > bound * (1.0 + 64.0 * f64::EPSILON) {
                return Err(bad("plasma field escaped its admitted analytic envelope"));
            }
            let heat = &carriers.mean_heating_w_m3;
            if [
                heat.total_power_density_w_m3(),
                heat.electron_power_density_w_m3(),
                heat.heavy_power_density_w_m3(),
                heat.thermal_exchange_correction_to_electrons_w_m3(),
                carriers.mean_drift_kinetic_energy_j_m3,
                carriers.reactive_power_var_m3,
            ]
            .iter()
            .any(|&value| !((value * ratio) * ratio).is_finite())
            {
                return Err(bad(
                    "plasma layer heat or carrier storage is unrepresentable",
                ));
            }
            (e.map(|value| carriers.conductivity_s_m * value), ratio)
        } else {
            ([C64::ZERO; 3], 0.0)
        };
        if current.iter().any(|&j| !finite(j)) {
            return Err(bad("plasma layer current is unrepresentable"));
        }
        let background_energy_j_m3 = layer.background_permittivity_f_m.map(|epsilon| {
            let e_scale = 0.5 * det::sqrt(epsilon);
            let h_scale = 0.5 * det::sqrt(layer.permeability_h_m);
            e.iter()
                .map(|v| v.scale(e_scale).norm_sq())
                .chain(h.iter().map(|v| v.scale(h_scale).norm_sq()))
                .fold(Dd::ZERO, |sum, x| sum + Dd::from_f64(x))
                .to_f64()
        });
        if background_energy_j_m3.is_some_and(|v| !v.is_finite() || v < 0.0) {
            return Err(bad("plasma background field storage is unrepresentable"));
        }
        Ok(IonizedGasLayerPoint {
            wave: self,
            layer_index,
            depth_m,
            envelope_ratio: ratio,
            background_energy_j_m3,
            fields,
            current_density_a_m2: current,
        })
    }
}

// Both physical source configurations pass their actual whole-film fields
// through the same original carrier domain and thermal partition admission.
fn admit_layer<'m>(
    stack: &PlanarDielectricLayers<'m>,
    index: usize,
    polarization: Polarization,
    parallel_over_omega: f64,
    kernel: Option<LayerFieldKernel>,
    absorbed_w_m2: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AdmittedLayer<'m>, ReflectionError> {
    let selected = stack.layers[index];
    let source = selected.source;
    let epsilon = selected.wave.epsilon_f_m;
    let mu = linear_magnetic_permeability(source.response.material())?;
    let longitudinal = match polarization {
        Polarization::Te => C64::from_re(parallel_over_omega / mu),
        Polarization::Tm => C64::from_re(-parallel_over_omega) / epsilon,
    };
    let root_loss =
        det::sqrt(stack.omega) * det::sqrt(-epsilon.im) * core::f64::consts::FRAC_1_SQRT_2;
    if !finite(longitudinal)
        || !root_loss.is_finite()
        || !absorbed_w_m2.is_finite()
        || absorbed_w_m2 < 0.0
    {
        return Err(bad("plasma layer field scale or power is unrepresentable"));
    }
    let mut layer = AdmittedLayer {
        source,
        wave: selected.wave,
        kernel,
        longitudinal,
        root_loss,
        permeability_h_m: mu,
        background_permittivity_f_m: None,
        electric_envelope_v_m: None,
        carriers: None,
        absorbed_w_m2,
        carrier_heating: None,
    };
    if let DielectricResponse::IonizedGas(plasma) = source.response {
        let Some(kernel) = kernel else {
            layer.carrier_heating = Some(IonizedGasLayerHeating::default());
            return Ok(layer);
        };
        let [e, h] = kernel
            .amplitude_bounds()
            .ok_or_else(|| bad("plasma layer field envelope is unrepresentable"))?;
        let envelope = match polarization {
            Polarization::Te => e,
            Polarization::Tm => det::hypot(e, longitudinal.abs() * h),
        };
        if !envelope.is_finite() {
            return Err(bad("plasma vector field envelope is unrepresentable"));
        }
        // The complex vector norm bounds every instantaneous vector norm.
        // Isotropy makes this real-axis drive conservative for every original
        // carrier and pair-relative drift, including longitudinal TM fields.
        let carriers = plasma.transport().response(
            stack.omega,
            [C64::from_re(envelope), C64::ZERO, C64::ZERO],
            cancelled,
        )?;
        if carriers.conductivity_s_m != plasma.conductivity_s_m() {
            return Err(bad("plasma film changed its original carrier response"));
        }
        layer.carrier_heating = Some(partition_heat(&carriers, absorbed_w_m2)?);
        layer.background_permittivity_f_m = Some(linear_electric_permittivity(plasma.material())?);
        layer.electric_envelope_v_m = Some(envelope);
        layer.carriers = Some(carriers);
    }
    Ok(layer)
}

impl<'w, 'm> IonizedGasLayerPoint<'w, 'm> {
    pub const fn wave(&self) -> &'w PreparedIonizedGasLayerIllumination<'m> {
        self.wave
    }
    pub const fn layer_index(&self) -> usize {
        self.layer_index
    }
    pub const fn depth_m(&self) -> f64 {
        self.depth_m
    }
    /// Original admitted peak fields and their derived local statistics.
    /// The shared reference preserves the field used by carrier rescaling.
    pub const fn fields(&self) -> &LayerPointField {
        &self.fields
    }
    /// Carrier order and physical units are those of the original gas owner.
    pub fn carrier_velocity_m_s(&self, carrier: usize) -> Option<[C64; 3]> {
        let layer = &self.wave.layers[self.layer_index];
        let velocity = layer.carriers.as_ref()?.carrier_velocity_m_s.get(carrier)?[0];
        let bound = layer.electric_envelope_v_m?;
        Some(if bound == 0.0 {
            [C64::ZERO; 3]
        } else {
            self.fields
                .electric_v_m
                .map(|e| velocity * C64::new(e.re / bound, e.im / bound))
        })
    }
    pub fn mean_carrier_heating_w_m3(&self) -> Option<f64> {
        self.carrier_power(|r| r.mean_heating_w_m3.total_power_density_w_m3())
    }
    pub fn mean_electron_heating_w_m3(&self) -> Option<f64> {
        self.carrier_power(|r| r.mean_heating_w_m3.electron_power_density_w_m3())
    }
    pub fn mean_heavy_heating_w_m3(&self) -> Option<f64> {
        self.carrier_power(|r| r.mean_heating_w_m3.heavy_power_density_w_m3())
    }
    pub fn mean_exchange_correction_to_electrons_w_m3(&self) -> Option<f64> {
        self.carrier_power(|r| {
            r.mean_heating_w_m3
                .thermal_exchange_correction_to_electrons_w_m3()
        })
    }
    pub fn mean_drift_kinetic_energy_j_m3(&self) -> Option<f64> {
        self.carrier_power(|r| r.mean_drift_kinetic_energy_j_m3)
    }
    pub fn reactive_carrier_power_var_m3(&self) -> Option<f64> {
        self.carrier_power(|r| r.reactive_power_var_m3)
    }
    /// Lossless bound-charge/background magnetic storage; carrier inertia is
    /// reported separately. Complex epsilon is never used as stored energy.
    pub fn mean_background_field_energy_j_m3(&self) -> Option<f64> {
        self.background_energy_j_m3
    }
    fn carrier_power(
        &self,
        select: impl FnOnce(&IonizedGasHarmonicResponse<'m>) -> f64,
    ) -> Option<f64> {
        let value = select(self.wave.layers[self.layer_index].carriers.as_ref()?);
        Some((value * self.envelope_ratio) * self.envelope_ratio)
    }
}

fn partition_heat(
    carriers: &IonizedGasHarmonicResponse<'_>,
    total: f64,
) -> Result<IonizedGasLayerHeating, ReflectionError> {
    let heat = &carriers.mean_heating_w_m3;
    let full = heat.total_power_density_w_m3();
    if full == 0.0 {
        if total != 0.0 {
            return Err(bad("nonzero plasma film heat with zero admitted field"));
        }
        return Ok(IonizedGasLayerHeating::default());
    }
    let result = IonizedGasLayerHeating {
        total_w_m2: total,
        electron_w_m2: total * (heat.electron_power_density_w_m3() / full),
        heavy_w_m2: total * (heat.heavy_power_density_w_m3() / full),
        exchange_correction_to_electrons_w_m2: total
            * (heat.thermal_exchange_correction_to_electrons_w_m3() / full),
    };
    if [
        result.total_w_m2,
        result.electron_w_m2,
        result.heavy_w_m2,
        result.exchange_correction_to_electrons_w_m2,
    ]
    .iter()
    .any(|v| !v.is_finite())
    {
        return Err(bad("plasma film thermal partition is unrepresentable"));
    }
    Ok(result)
}

fn vector_norm(v: [C64; 3]) -> f64 {
    det::hypot(v[0].abs(), det::hypot(v[1].abs(), v[2].abs()))
}

fn point_field(
    e: [C64; 3],
    h: [C64; 3],
    root_loss: f64,
) -> Result<LayerPointField, ReflectionError> {
    let rms = |v: [C64; 3]| v.map(|x| x.scale(core::f64::consts::FRAC_1_SQRT_2).norm_sq());
    let poynting = core::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        (e[j] * h[k].conj() - e[k] * h[j].conj()).scale(0.5)
    });
    let heat = e
        .iter()
        .fold(Dd::ZERO, |sum, x| {
            sum + Dd::from_f64(x.scale(root_loss).norm_sq())
        })
        .to_f64();
    let stats = LayerPointStatistics {
        electric_rms_squared_v2_m2: rms(e),
        magnetic_rms_squared_a2_m2: rms(h),
        complex_poynting_w_m2: poynting,
        absorbed_power_w_m3: heat,
    };
    if e.iter().chain(&h).chain(&poynting).any(|&v| !finite(v))
        || stats
            .electric_rms_squared_v2_m2
            .iter()
            .chain(&stats.magnetic_rms_squared_a2_m2)
            .any(|v| !v.is_finite())
        || !heat.is_finite()
        || heat < 0.0
    {
        return Err(bad("plasma film field or local power is unrepresentable"));
    }
    Ok(LayerPointField {
        electric_v_m: e,
        magnetic_a_m: h,
        statistics: stats,
    })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ReflectionError> {
    if cancelled() {
        Err(ReflectionError::IonizedGas(IonizedGasError::Cancelled))
    } else {
        Ok(())
    }
}
