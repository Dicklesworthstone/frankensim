//! The original one-sided Maxwell solution with admitted finite gas films.
use super::*;

pub(super) struct ExitWave {
    electric_v_m: [C64; 3],
    magnetic_a_m: [C64; 3],
    kz: C64,
    phase_wavenumber: Dd,
    root_loss: f64,
    entering_w_m2: f64,
}

impl<'m> PreparedIonizedGasLayerIllumination<'m> {
    /// Admit one actual incoming plane wave and an outgoing passive half-space.
    /// The right source must be exactly zero. The original incident medium is
    /// lossless; the original non-ionized exit may be absorbing, evanescent or
    /// exactly critical. This does not create an incoming power channel in a
    /// lossy medium or relax the two-sided constructor's port requirements.
    ///
    /// All finite gas films retain their original whole-vector drift admission
    /// and carrier heat partitions. A shared O(layers) native field preparation
    /// supplies every film; each selected field then costs O(1). The actual
    /// amplitude is fixed for this owner and cannot be replaced after admission.
    pub fn try_new_one_sided(
        incident: DielectricResponse<'m>,
        exit: DielectricResponse<'m>,
        omega_rad_s: f64,
        layers: &[DielectricLayer<'m>],
        excitation: IonizedGasLayerExcitation,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ReflectionError> {
        poll(&mut cancelled)?;
        let [amplitude, right] = excitation.incoming_tangential_e_v_m;
        if !finite(amplitude) || right != C64::ZERO {
            return Err(bad(
                "one-sided plasma films require one finite left source and zero right source",
            ));
        }
        let stack =
            PlanarDielectricLayers::try_new_for_fields(incident, exit, omega_rad_s, layers, true)?;
        let c = excitation.incidence_cosine;
        let polarization = excitation.polarization;
        let response = stack.response(c, polarization)?;
        poll(&mut cancelled)?;
        let prepared = LayerFieldPreparation::new(&stack, c, polarization)?;
        let reference = prepared.admittance_s();
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
        let mut total_heat = Dd::ZERO;
        for (index, source) in layers.iter().enumerate() {
            poll(&mut cancelled)?;
            let (kernel, heat) = if source.thickness_m == 0.0 {
                (None, 0.0)
            } else {
                let unit = prepared.field(index)?;
                let heat = amplitude
                    .scale(det::sqrt(unit.absorbed_power_w_m2))
                    .norm_sq();
                let kernel = unit
                    .kernel
                    .scaled(amplitude)
                    .ok_or_else(|| bad("physical plasma film field is unrepresentable"))?;
                (Some(kernel), heat)
            };
            admitted.push(admit_layer(
                &stack,
                index,
                polarization,
                parallel_over_omega,
                kernel,
                heat,
                &mut cancelled,
            )?);
            total_heat = total_heat + Dd::from_f64(heat);
        }
        let traces = prepared
            .boundary_traces()
            .map(|face| face.map(|v| v * amplitude));
        let incoming = amplitude
            .scale(det::sqrt(reference) * core::f64::consts::FRAC_1_SQRT_2)
            .norm_sq();
        let reflected = incoming * response.reflected_power_fraction;
        let transmitted = incoming * response.transmitted_power_fraction;
        let absorbed = total_heat.to_f64();
        let outgoing = [response.reflected_tangential_e * amplitude, traces[1][0]];
        let power_scale = incoming.max(reflected).max(transmitted).max(absorbed);
        let imbalance = (Dd::from_f64(incoming)
            - Dd::from_f64(reflected)
            - Dd::from_f64(transmitted)
            - total_heat)
            .to_f64();
        let tolerance = 512.0 * f64::EPSILON * (layers.len() + 1) as f64;
        if traces
            .iter()
            .flatten()
            .chain(&outgoing)
            .any(|&v| !finite(v))
            || [incoming, reflected, transmitted, absorbed]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
            || (power_scale > 0.0 && imbalance.abs() / power_scale > tolerance)
        {
            return Err(bad(
                "one-sided plasma field or Poynting balance is unrepresentable",
            ));
        }
        let boundary = CoherentLayerIllumination {
            outgoing_tangential_e_v_m: outgoing,
            boundary_tangential_e_v_m: traces.map(|face| face[0]),
            boundary_conjugate_h_a_m: traces.map(|face| face[1]),
            power: LayerIlluminationPower {
                incoming_w_m2: [incoming, 0.0],
                outgoing_w_m2: [reflected, transmitted],
                absorbed_w_m2: absorbed,
            },
        };
        let epsilon = stack.exit_wave.epsilon_f_m;
        let mu = linear_magnetic_permeability(exit.material())?;
        let [et, ht] = traces[1];
        // Keep independent H at a critical TM exit: E_t=k_z=0 while H_t
        // and E_z remain finite. No E_t/k_z reconstruction is used.
        let (electric_v_m, magnetic_a_m) = match polarization {
            Polarization::Te => (
                [C64::ZERO, et, C64::ZERO],
                [-ht, C64::ZERO, et.scale(parallel_over_omega / mu)],
            ),
            Polarization::Tm => (
                [
                    et,
                    C64::ZERO,
                    ht * (C64::from_re(-parallel_over_omega) / epsilon),
                ],
                [C64::ZERO, ht, C64::ZERO],
            ),
        };
        let normal = stack.exit_wave.normal(c)?;
        let exit_wave = ExitWave {
            electric_v_m,
            magnetic_a_m,
            kz: response.exit_normal_wavenumber_rad_m,
            phase_wavenumber: Dd::from_f64(normal.re) * Dd::from_f64(stack.wave_number),
            root_loss: det::sqrt(omega_rad_s)
                * det::sqrt(-epsilon.im)
                * core::f64::consts::FRAC_1_SQRT_2,
            entering_w_m2: transmitted,
        };
        point_field(electric_v_m, magnetic_a_m, exit_wave.root_loss)?;
        poll(&mut cancelled)?;
        Ok(Self {
            incident,
            exit,
            omega_rad_s,
            wave_number_rad_m: stack.wave_number,
            reference_admittance_s: reference,
            excitation,
            boundary,
            layers: admitted,
            one_sided_exit: Some(exit_wave),
        })
    }

    /// Original outgoing E/H and local dielectric heat in the exit half-space,
    /// measured from the last film's exit face at x=y=0. Requires one-sided
    /// preparation. Evaluation is allocation-free and retains the native
    /// passive complex phase, including evanescent and critical TM fields.
    pub fn exit_field(&self, depth_m: f64) -> Result<LayerPointField, ReflectionError> {
        let wave = self.admitted_exit(depth_m)?;
        let phase = super::super::kernel::decay(
            wave.kz,
            wave.phase_wavenumber,
            depth_m,
            Dd::from_f64(depth_m),
        )
        .ok_or_else(|| bad("exit half-space propagation is unrepresentable"))?;
        point_field(
            wave.electric_v_m.map(|v| v * phase),
            wave.magnetic_a_m.map(|v| v * phase),
            wave.root_loss,
        )
    }

    /// Physical loss from the original exit face to depth [W/m²]. Finite-film
    /// heat is separate. The native attenuation integrates loss without a
    /// depth mesh or subtraction of nearly equal Poynting powers. The zero
    /// loss of an evanescent/critical lossless exit remains exact.
    pub fn absorbed_exit_power_to_depth_w_m2(&self, depth_m: f64) -> Result<f64, ReflectionError> {
        let wave = self.admitted_exit(depth_m)?;
        let heat = wave.entering_w_m2 * -det::expm1((wave.kz.im * depth_m) * 2.0);
        if !heat.is_finite() || heat < 0.0 {
            return Err(bad("exit half-space absorption is unrepresentable"));
        }
        Ok(heat)
    }

    fn admitted_exit(&self, depth_m: f64) -> Result<&ExitWave, ReflectionError> {
        if !depth_m.is_finite() || depth_m < 0.0 {
            return Err(bad("exit half-space depth must be finite and nonnegative"));
        }
        self.one_sided_exit
            .as_ref()
            .ok_or_else(|| bad("exit half-space queries require one-sided preparation"))
    }
}
