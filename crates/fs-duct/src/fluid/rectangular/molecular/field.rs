//! Local fields and separately positive heat losses of the original molecular
//! rectangular/slit wave. The same dynamic heat capacity and wall diffusion
//! solve `rho cp(-i omega) (-i omega) T - kappa Laplacian(T) = -i omega p`.
//!
//! Each original reservoir obeys `(-i omega) T_j = (T - T_j)/tau_j`.
//! With peak phasors, irreversible work is `mu |grad u|²/2`,
//! `kappa |grad T|²/(2 T0)`, and
//! `rho C_j |T - T_j|²/(2 T0 tau_j)`. These are three independent positive
//! channels, not an attenuation residual or an additional bulk-viscosity law.
//! The shared rectangular modes supply both local gradients and analytic
//! whole-section integrals. No second molecular or wall attenuation is added.

use super::super::field::{DiffusionField, amplitudes};
use super::super::{
    FluidDuctError, RectangularAdmission, RectangularSection, RectangularWave, bad, finite, poll,
    positive,
};
use crate::DuctError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::relaxation::{EvaluatedProcess, RelaxationState};
use fs_math::{c64::C64, dd::Dd, det};

/// Actual local fields under peak exp(-i omega t), with coordinates/gradients
/// ordered [short side, long side], independently of constructor argument order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularRectangularField {
    /// Distance from the first short-side wall [m].
    pub short_coordinate_m: f64,
    /// Distance from the adjacent long-side wall [m].
    pub long_coordinate_m: f64,
    /// Cross-section-uniform pressure phasor [Pa].
    pub pressure_pa: C64,
    /// Actual axial particle velocity [m/s].
    pub axial_velocity_m_s: C64,
    /// Both transverse derivatives of axial velocity [1/s].
    pub transverse_velocity_gradient_per_s: [C64; 2],
    /// Translational/equilibrated-background excess temperature [K].
    pub gas_temperature_k: C64,
    /// Both transverse derivatives of that temperature [K/m].
    pub transverse_temperature_gradient_k_m: [C64; 2],
    /// Original transverse shear work per volume [W/m³].
    pub viscous_dissipation_w_m3: f64,
    /// Original transverse Fourier heat-conduction work per volume [W/m³].
    pub conduction_dissipation_w_m3: f64,
    /// Sum of original internal-reservoir exchange work per volume [W/m³].
    pub relaxation_dissipation_w_m3: f64,
}

/// One original species/reservoir's local thermal lag and positive loss.
/// The borrowed process preserves its original species, index, strength and time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularRectangularReservoirField<'material> {
    /// The source material's evaluated process, not a fitted pole.
    pub process: &'material EvaluatedProcess,
    /// Excess temperature of the gas's equilibrated background [K].
    pub gas_temperature_k: C64,
    /// Excess temperature of this internal reservoir [K].
    pub reservoir_temperature_k: C64,
    /// T - T_j [K], evaluated without cancellation of equilibrated temperatures.
    pub temperature_difference_k: C64,
    /// This reservoir's cycle-mean exchange loss per fluid volume [W/m³].
    pub relaxation_dissipation_w_m3: f64,
}

/// Cross-section-integrated irreversible work per axial metre [W/m].
/// These are analytic positive field integrals, never flux differences.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MolecularRectangularDissipationPerLength {
    /// Original transverse shear work.
    pub viscous_w_m: f64,
    /// Original transverse Fourier heat-conduction work.
    pub conduction_w_m: f64,
    /// Sum of original internal-reservoir exchange work.
    pub relaxation_w_m: f64,
    /// Sum of the three independently positive channels.
    pub total_w_m: f64,
}

#[derive(Clone, Copy, Debug)]
struct Reservoir {
    temperature_ratio: C64,
    temperature_difference_ratio: C64,
    /// sqrt(rho C_j/(2 T0 tau_j)) * |1 - 1/(1 - i omega tau_j)|.
    loss_root: f64,
}

/// Source-bound resolved fields of the SAME molecular mean-wave owner.
/// Preparation is O(N + R), a point query O(N), and whole-section work O(1).
/// Two bounded N-mode arrays and at most the original 32 reservoir records are
/// allocated fallibly; queries allocate nothing and do not evaluate material laws.
///
/// The original mean-series tail admission does not certify pointwise fields.
/// Refine the authored term count to resolve gradients/corners. Rigid no-slip
/// isothermal walls, linear continuum motion, uniform section pressure and no
/// axial thermal diffusion retain the original model's domain. These thermal
/// observations are not a background-temperature or species-evolution solver.
#[derive(Debug)]
pub struct PreparedMolecularRectangularFields<'geometry, 'material> {
    section: &'geometry RectangularSection,
    source: &'material RelaxationState,
    state: ThermoviscousFluidState,
    omega: f64,
    wave: RectangularWave,
    heat_capacity: C64,
    viscous: DiffusionField,
    thermal: DiffusionField,
    reservoirs: Vec<Reservoir>,
    velocity_scale: C64,
    temperature_scale: C64,
    loss_roots: [f64; 3],
    total_loss_roots: [f64; 3],
}

impl RectangularSection {
    /// Resolve gas/reservoir temperatures, shear and all three positive loss
    /// channels of the original joint wall/molecular wave at one frequency.
    pub fn prepare_molecular_fields<'g, 'm>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        admission: RectangularAdmission,
    ) -> Result<PreparedMolecularRectangularFields<'g, 'm>, FluidDuctError> {
        self.prepare_molecular_fields_controlled(material, omega, admission, || false)
    }

    /// Preparation with the original bounded cancellation checks. Every
    /// reservoir and each block of 32 transverse modes polls before work.
    /// No partially prepared response is returned on cancellation or refusal.
    pub fn prepare_molecular_fields_controlled<'g, 'm>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        admission: RectangularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularRectangularFields<'g, 'm>, FluidDuctError> {
        // Reuse the original mean-wave admission and its exact empty-reservoir
        // path. RelaxationState is immutable: the cached caloric/transport views
        // below refer to precisely the same source and frequency.
        let wave =
            self.evaluate_molecular_controlled(material, omega, admission, &mut cancelled)?;
        let gas = material.gas();
        let state = gas.thermoviscous_state()?;
        let transport = state.transport();
        let (heat_capacity, thermal_argument) = if material.processes().is_empty() {
            (
                C64::from_re(transport.specific_heat_j_kg_k()),
                C64::new(0.0, omega / transport.thermal_diffusivity_m2_s()),
            )
        } else {
            let cp = material.caloric_response(omega)?.cp_j_kg_k;
            (
                cp,
                C64::new(cp.im, cp.re).scale(
                    omega * state.acoustic().density_kg_m3()
                        / transport.thermal_conductivity_w_m_k(),
                ),
            )
        };
        poll(&mut cancelled, 0)?;
        let viscous = self.prepare_diffusion_fields(
            C64::new(0.0, omega / transport.kinematic_viscosity_m2_s()),
            &mut cancelled,
        )?;
        let thermal = self.prepare_diffusion_fields(thermal_argument, &mut cancelled)?;
        let velocity_scale = wave.viscous.mean.conj().recip() / C64::from_re(self.area_m2);
        let temperature_scale = heat_capacity.recip() / C64::from_re(gas.density);
        let mut reservoirs = Vec::new();
        reservoirs
            .try_reserve_exact(material.processes().len())
            .map_err(|_| {
                bad("rectangular molecular reservoir allocation exceeds available memory")
            })?;
        let mut molecular_energy = Dd::ZERO;
        for (index, process) in material.processes().iter().enumerate() {
            poll(&mut cancelled, index)?;
            let x = omega * process.relaxation_time_s;
            // Stable Debye factors, including omega*tau overflowing to infinity.
            // Retain the tiny fast-reservoir lag instead of subtracting T_j from T.
            let (ratio, difference, rate_root) = if x <= 1.0 {
                let denominator = 1.0 + x * x;
                (
                    C64::new(1.0 / denominator, x / denominator),
                    C64::new(x * x / denominator, -x / denominator),
                    det::sqrt(omega) * det::sqrt(x / denominator),
                )
            } else {
                let inverse = 1.0 / x;
                let denominator = 1.0 + inverse * inverse;
                (
                    C64::new(inverse * inverse / denominator, inverse / denominator),
                    C64::new(1.0 / denominator, -inverse / denominator),
                    (1.0 / det::sqrt(process.relaxation_time_s)) / det::sqrt(denominator),
                )
            };
            let loss_root = (det::sqrt(gas.density) * det::sqrt(process.heat_capacity_j_kg_k)
                / det::sqrt(gas.temperature))
                * core::f64::consts::FRAC_1_SQRT_2
                * rate_root;
            if !finite(ratio) || !finite(difference) || !loss_root.is_finite() || loss_root < 0.0 {
                return Err(bad("rectangular molecular reservoir field is unrepresentable").into());
            }
            molecular_energy = molecular_energy + Dd::from_f64(loss_root * loss_root);
            reservoirs.push(Reservoir {
                temperature_ratio: ratio,
                temperature_difference_ratio: difference,
                loss_root,
            });
        }
        let loss_roots = [
            det::sqrt(transport.dynamic_viscosity_pa_s()) * core::f64::consts::FRAC_1_SQRT_2,
            (det::sqrt(transport.thermal_conductivity_w_m_k()) / det::sqrt(gas.temperature))
                * core::f64::consts::FRAC_1_SQRT_2,
            det::sqrt(molecular_energy.to_f64()),
        ];
        let total_loss_roots = [
            velocity_scale.scale(loss_roots[0]).abs() * viscous.gradient_norm,
            temperature_scale.scale(loss_roots[1]).abs() * thermal.gradient_norm,
            temperature_scale.scale(loss_roots[2]).abs() * thermal.field_norm,
        ];
        if !finite(velocity_scale)
            || velocity_scale == C64::ZERO
            || !finite(temperature_scale)
            || temperature_scale == C64::ZERO
            || !loss_roots.into_iter().all(|r| r.is_finite() && r >= 0.0)
            || !total_loss_roots
                .into_iter()
                .all(|r| r.is_finite() && r >= 0.0)
            || !positive(total_loss_roots[0])
            || !positive(total_loss_roots[1])
        {
            return Err(bad(
                "rectangular molecular field normalization or work is unrepresentable",
            )
            .into());
        }
        poll(&mut cancelled, self.terms())?;
        Ok(PreparedMolecularRectangularFields {
            section: self,
            source: material,
            state,
            omega,
            wave,
            heat_capacity,
            viscous,
            thermal,
            reservoirs,
            velocity_scale,
            temperature_scale,
            loss_roots,
            total_loss_roots,
        })
    }
}

impl<'g, 'm> PreparedMolecularRectangularFields<'g, 'm> {
    /// Original immutable gas and reservoir source.
    pub const fn source(&self) -> &'m RelaxationState {
        self.source
    }
    /// Original geometry and authored transverse resolution.
    pub const fn section(&self) -> &'g RectangularSection {
        self.section
    }
    /// The original gas's coherent equilibrium acoustic/transport state.
    pub const fn state(&self) -> &ThermoviscousFluidState {
        &self.state
    }
    /// Bound angular frequency [rad/s].
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// Exactly the same mean response as RectangularSection::evaluate_molecular.
    pub const fn wave(&self) -> &RectangularWave {
        &self.wave
    }
    /// Original dynamic isobaric heat capacity, under exp(-i omega t) [J/(kg K)].
    pub const fn heat_capacity_j_kg_k(&self) -> C64 {
        self.heat_capacity
    }

    /// Local velocity, actual gas temperature, both transverse gradients and
    /// independently positive viscous/conduction/relaxation power densities.
    pub fn field_at(
        &self,
        short_coordinate_m: f64,
        long_coordinate_m: f64,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<MolecularRectangularField, DuctError> {
        self.field_at_controlled(
            short_coordinate_m,
            long_coordinate_m,
            pressure_pa,
            volume_flow_m3_s,
            || false,
        )
    }

    /// As field_at, with cancellation at every block of 32 modes.
    pub fn field_at_controlled(
        &self,
        short_coordinate_m: f64,
        long_coordinate_m: f64,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularRectangularField, DuctError> {
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let fraction = self.fraction(short_coordinate_m, long_coordinate_m)?;
        let velocity = self
            .viscous
            .at(
                fraction,
                long_coordinate_m,
                self.section.long_m,
                &mut cancelled,
            )?
            .map(|v| v.conj() * self.velocity_scale * volume_flow_m3_s);
        let temperature = self
            .thermal
            .at(
                fraction,
                long_coordinate_m,
                self.section.long_m,
                &mut cancelled,
            )?
            .map(|t| t.conj() * self.temperature_scale * pressure_pa);
        let viscous = velocity[1..]
            .iter()
            .map(|v| v.scale(self.loss_roots[0]).norm_sq())
            .sum::<f64>();
        let conduction = temperature[1..]
            .iter()
            .map(|t| t.scale(self.loss_roots[1]).norm_sq())
            .sum::<f64>();
        let relaxation = temperature[0].scale(self.loss_roots[2]).norm_sq();
        if !velocity.into_iter().all(finite)
            || !temperature.into_iter().all(finite)
            || ![viscous, conduction, relaxation]
                .into_iter()
                .all(f64::is_finite)
        {
            return Err(bad(
                "rectangular molecular local field or work is unrepresentable",
            ));
        }
        Ok(MolecularRectangularField {
            short_coordinate_m,
            long_coordinate_m,
            pressure_pa,
            axial_velocity_m_s: velocity[0],
            transverse_velocity_gradient_per_s: [velocity[1], velocity[2]],
            gas_temperature_k: temperature[0],
            transverse_temperature_gradient_k_m: [temperature[1], temperature[2]],
            viscous_dissipation_w_m3: viscous,
            conduction_dissipation_w_m3: conduction,
            relaxation_dissipation_w_m3: relaxation,
        })
    }

    /// One source reservoir's temperature, thermal lag and own loss. The index
    /// addresses source().processes(), preserving that canonical source order.
    /// Only the thermal modes are evaluated; velocity modes are not needed.
    pub fn reservoir_field_at(
        &self,
        process_index: usize,
        short_coordinate_m: f64,
        long_coordinate_m: f64,
        pressure_pa: C64,
    ) -> Result<MolecularRectangularReservoirField<'m>, DuctError> {
        self.reservoir_field_at_controlled(
            process_index,
            short_coordinate_m,
            long_coordinate_m,
            pressure_pa,
            || false,
        )
    }

    /// As reservoir_field_at, with bounded cancellation during the thermal sum.
    pub fn reservoir_field_at_controlled(
        &self,
        process_index: usize,
        short_coordinate_m: f64,
        long_coordinate_m: f64,
        pressure_pa: C64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularRectangularReservoirField<'m>, DuctError> {
        let reservoir = self.reservoir(process_index)?;
        amplitudes(pressure_pa, C64::ZERO)?;
        let fraction = self.fraction(short_coordinate_m, long_coordinate_m)?;
        let temperature = self.thermal.at(
            fraction,
            long_coordinate_m,
            self.section.long_m,
            &mut cancelled,
        )?[0]
            .conj()
            * self.temperature_scale
            * pressure_pa;
        let internal = temperature * reservoir.temperature_ratio;
        let difference = temperature * reservoir.temperature_difference_ratio;
        let work = temperature.scale(reservoir.loss_root).norm_sq();
        if ![temperature, internal, difference].into_iter().all(finite) || !work.is_finite() {
            return Err(bad(
                "rectangular internal-reservoir observation is unrepresentable",
            ));
        }
        Ok(MolecularRectangularReservoirField {
            process: &self.source.processes()[process_index],
            gas_temperature_k: temperature,
            reservoir_temperature_k: internal,
            temperature_difference_k: difference,
            relaxation_dissipation_w_m3: work,
        })
    }

    /// Analytic cross-section integrals of all three positive work channels.
    /// Pressure/flow may include coherent forward and backward waves.
    pub fn dissipation_per_axial_length(
        &self,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<MolecularRectangularDissipationPerLength, DuctError> {
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let viscous = volume_flow_m3_s.scale(self.total_loss_roots[0]).norm_sq();
        let conduction = pressure_pa.scale(self.total_loss_roots[1]).norm_sq();
        let relaxation = pressure_pa.scale(self.total_loss_roots[2]).norm_sq();
        let total = viscous + conduction + relaxation;
        if ![viscous, conduction, relaxation, total]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(bad(
                "rectangular molecular integrated work is unrepresentable",
            ));
        }
        Ok(MolecularRectangularDissipationPerLength {
            viscous_w_m: viscous,
            conduction_w_m: conduction,
            relaxation_w_m: relaxation,
            total_w_m: total,
        })
    }

    /// One original reservoir's analytic cross-section work [W/m], enabling a
    /// species/process audit without integrating a point grid or changing poles.
    pub fn reservoir_dissipation_per_axial_length(
        &self,
        process_index: usize,
        pressure_pa: C64,
    ) -> Result<f64, DuctError> {
        let reservoir = self.reservoir(process_index)?;
        amplitudes(pressure_pa, C64::ZERO)?;
        let root =
            self.temperature_scale.scale(reservoir.loss_root).abs() * self.thermal.field_norm;
        let work = pressure_pa.scale(root).norm_sq();
        if !work.is_finite() {
            return Err(bad(
                "rectangular reservoir integrated work is unrepresentable",
            ));
        }
        Ok(work)
    }

    fn reservoir(&self, process_index: usize) -> Result<&Reservoir, DuctError> {
        self.reservoirs
            .get(process_index)
            .ok_or_else(|| bad("rectangular reservoir index is outside the original material"))
    }

    fn fraction(&self, short: f64, long: f64) -> Result<f64, DuctError> {
        if !short.is_finite()
            || !(0.0..=self.section.short_m).contains(&short)
            || !long.is_finite()
            || !(0.0..=self.section.long_m).contains(&long)
        {
            return Err(bad(
                "rectangular molecular field point is outside the original section",
            ));
        }
        Ok(short / self.section.short_m)
    }
}
