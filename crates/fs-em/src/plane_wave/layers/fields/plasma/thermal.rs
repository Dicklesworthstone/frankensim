//! Conservative depth-cell loss from the already admitted physical film fields.
use super::*;

impl PreparedIonizedGasLayerIllumination<'_> {
    /// Physical cycle-mean dielectric loss through a film depth interval [W/m²].
    /// Depths are measured from this film's left face and must satisfy
    /// `0 <= start <= end <= thickness`. Empty intervals, including a
    /// zero-thickness identity film, have exactly zero heat.
    ///
    /// The already admitted coherent E/H traces feed the native positive slab
    /// loss modes restricted to this interval. Standing waves, longitudinal TM
    /// loss, and weak/opaque limits retain their original Maxwell solution.
    /// This O(1), allocation-free query uses no depth quadrature, flux subtraction,
    /// new material evaluation, carrier solve or stack solve. It cannot change
    /// the prepared physical source amplitudes or their whole-film drift domain.
    ///
    /// Disjoint cells sum to the original film power. A thermal consumer applies
    /// area/time and deposits each cell's energy once; these spectral cycle means
    /// do not evolve a temperature, transient reservoir or material population.
    pub fn layer_absorbed_power_in_depth_interval_w_m2(
        &self,
        layer_index: usize,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<f64, ReflectionError> {
        let layer = self
            .layers
            .get(layer_index)
            .ok_or_else(|| bad("plasma depth-interval film index is out of range"))?;
        let width = admit_depth_interval(layer.source.thickness_m, start_depth_m, end_depth_m)?;
        if width == 0.0 || layer.wave.epsilon_f_m.im == 0.0 {
            return Ok(0.0);
        }
        if start_depth_m == 0.0 && end_depth_m == layer.source.thickness_m {
            return Ok(layer.absorbed_w_m2);
        }
        let slab = Slab::new_with_phase_thickness(
            layer.wave,
            self.wave_number_rad_m,
            self.excitation.incidence_cosine,
            self.excitation.polarization,
            width,
            Dd::from_f64(end_depth_m) - Dd::from_f64(start_depth_m),
        )?;
        let kernel = layer
            .kernel
            .ok_or_else(|| bad("positive film interval requires its admitted field"))?;
        let trace = |depth| {
            kernel
                .field(depth)
                .ok_or_else(|| bad("plasma depth-interval field trace is unrepresentable"))
        };
        let modes = interval_loss_transfer(
            slab,
            self.reference_admittance_s,
            trace(start_depth_m)?,
            trace(end_depth_m)?,
        )?;
        let heat = (Dd::from_f64(modes[0].norm_sq()) + Dd::from_f64(modes[1].norm_sq())).to_f64();
        if !heat.is_finite() || heat < 0.0 {
            return Err(bad(
                "plasma depth-interval absorbed power is unrepresentable",
            ));
        }
        Ok(heat)
    }

    /// Original electron/heavy friction partitions over one depth cell [W/m²].
    /// The signed electron exchange correction is already included in these
    /// partitions; it is not an additional heat deposit. Zero-drift thermal
    /// exchange remains a separate material process.
    ///
    /// Returns `None` for an ordinary non-gas film and a zero partition for an
    /// empty gas interval. Index/depth admission is identical to the total-loss
    /// query, including for non-gas films. The original cached carrier fractions
    /// are applied to the same positive cell loss without resolving carriers.
    pub fn layer_carrier_heating_in_depth_interval(
        &self,
        layer_index: usize,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<Option<IonizedGasLayerHeating>, ReflectionError> {
        let total = self.layer_absorbed_power_in_depth_interval_w_m2(
            layer_index,
            start_depth_m,
            end_depth_m,
        )?;
        let layer = &self.layers[layer_index];
        if start_depth_m == 0.0 && end_depth_m == layer.source.thickness_m {
            return Ok(layer.carrier_heating);
        }
        match &layer.carriers {
            Some(carriers) => Ok(Some(partition_heat(carriers, total)?)),
            // A gas identity film has no volume carrier solve; its admitted
            // integrated heat is exactly zero. Ordinary films retain None.
            None => Ok(layer
                .carrier_heating
                .map(|_| IonizedGasLayerHeating::default())),
        }
    }
}
