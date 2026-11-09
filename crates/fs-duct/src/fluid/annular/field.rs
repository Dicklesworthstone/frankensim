//! Original radial fields and positive acoustic loss in an annular sound passage.
//!
//! The existing Stinson diffusion solutions h_v and h_t are retained, not
//! re-solved or replaced by an equal-area circle. With exp(-i omega t),
//! `u(r) = (U/A) conj(h_v(r)/<h_v>)` and the temperature response is proportional
//! to `p conj(h_t(r))`. The velocity and thermal gradients are those of the
//! original radial P1 field, including BOTH stationary walls.
//!
//! Local cycle-mean loss is `mu |du/dr|²/2` and
//! `(gamma-1) D_t |p dh_t/dr|²/(2 K_s)`. The latter follows from
//! `kappa_T-kappa_s = alpha² T/(rho Cp) = (gamma-1)/K_s` and
//! `D_t = conductivity/(rho Cp)`. It therefore applies to the original
//! admitted liquid or gas without inventing a thermal expansion coefficient
//! or temperature. A complete ring integral uses its actual `2 pi r dr`.
//! Source: Stinson (1991), equations 55--69, <https://doi.org/10.1121/1.400379>.
//! For the original molecular gas, complex Cp enters that SAME radial equation:
//! `T'=p conj(h_t)/(rho Cp)`. Fourier availability loss is
//! `conductivity |grad T'|²/(2 T0)`, and the original reservoir loss is
//! `rho omega Im(Cp) |T'|²/(2 T0)`. Only viscous work is a direct caloric
//! source. Thermal availability/entropy loss does not by itself specify local
//! mean heat transport or wall/fluid deposition.
//!
//! These are the existing fully developed, uniform-pressure model's fields,
//! not a new 3-D solution. Radial velocity, axial-gradient dissipation, bulk
//! viscosity, entrance effects and wall/fluid heat partition are absent.
//! Molecular relaxation is present only for `prepare_molecular_fields`.
//! The caller supplies coherent pressure and volume
//! flow at this prepared frequency from its original axial acoustic solve.

use super::{
    AnnularAdmission, AnnularSection, AnnularWave, FluidDuctError, RadialSolution, bad, finite,
    poll, positive, radial_mass,
};
use crate::DuctError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::gas::GasState;
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, dd::Dd, det};

/// Fields inside ONE original radial P1 cell, with derivatives taken on that
/// cell. Pressure and velocity are peak phasors under exp(-i omega t).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnnularField {
    /// Physical radial coordinate from the common cylinder axis [m].
    pub radius_m: f64,
    /// Cross-section-uniform acoustic pressure [Pa].
    pub pressure_pa: C64,
    /// Actual axial particle velocity, not the cross-section mean U/A.
    pub axial_velocity_m_s: C64,
    /// Piecewise-constant radial derivative of axial velocity [1/s].
    pub radial_velocity_gradient_per_s: C64,
    /// p*h_t [Pa]. Multiplying by the original (possibly complex) dT/dp gives T' [K].
    /// No expansion coefficient or background temperature is inferred for a
    /// general fluid. Absent exactly when Cp=Cv and thermal coupling is zero.
    pub thermal_pressure_response_pa: Option<C64>,
    /// Cycle-mean irreversible shear work per fluid volume [W/m³].
    pub viscous_dissipation_w_m3: f64,
    /// Fourier thermal availability loss [W/m³], not a local caloric source.
    pub thermal_dissipation_w_m3: f64,
    /// Original molecular-reservoir availability loss [W/m³].
    /// Zero for a nonrelaxing fluid; no added bulk viscosity is inferred.
    pub molecular_dissipation_w_m3: f64,
}

/// Irreversible acoustic loss per unit axial length [W/m].
/// This is a local cross-section/ring integral, not a finite-length integral
/// or a net boundary-power difference. It does not assign heat to a wall.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AnnularDissipationPerLength {
    /// Original radial shear work [W/m].
    pub viscous_w_m: f64,
    /// Original radial Fourier availability loss [W/m].
    pub thermal_w_m: f64,
    /// Original molecular-reservoir availability loss [W/m].
    pub molecular_w_m: f64,
    /// Sum of all independently positive acoustic loss channels [W/m].
    pub total_w_m: f64,
}

/// One source-bound, frequency-specific radial solution. Every subsequent
/// field/ring/whole-section query uses fixed work without allocation, material
/// evaluation or another radial solve. Preparation retains two native solution
/// vectors at most, within the original AnnularAdmission payload cap.
#[derive(Debug)]
pub struct PreparedAnnularFields<'geometry, 'material, M: ?Sized> {
    section: &'geometry AnnularSection,
    source: &'material M,
    state: ThermoviscousFluidState,
    omega: f64,
    wave: AnnularWave,
    viscous: RadialSolution,
    thermal: Option<RadialSolution>,
    velocity_scale: C64,
    loss_scales: [f64; 3],
    total_loss_roots: [f64; 3],
    molecular_cp: Option<C64>,
}

impl AnnularSection {
    /// Retain the native fields that produce this section's pressure/flow wave.
    /// One original material query supplies all acoustic/caloric/transport data.
    pub fn prepare_fields<'g, 'm, M: ThermoviscousFluid + ?Sized>(
        &'g self,
        material: &'m M,
        omega: f64,
        admission: AnnularAdmission,
    ) -> Result<PreparedAnnularFields<'g, 'm, M>, FluidDuctError> {
        self.prepare_fields_controlled(material, omega, admission, || false)
    }

    /// As `prepare_fields`, polling at the original solve boundaries and every
    /// 64 cells during field-weight preparation. Refusal publishes no partial
    /// field owner. Mesh resolution remains an independent accuracy axis.
    pub fn prepare_fields_controlled<'g, 'm, M: ThermoviscousFluid + ?Sized>(
        &'g self,
        material: &'m M,
        omega: f64,
        admission: AnnularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAnnularFields<'g, 'm, M>, FluidDuctError> {
        poll(&mut cancelled, 0)?;
        self.admit(omega, admission)?;
        let state = material.thermoviscous_state()?;
        let acoustic = state.acoustic();
        let transport = state.transport();
        let k0 = omega / acoustic.sound_speed_m_s();
        if !positive(k0) || k0 * self.outer_m > admission.maximum_transverse_phase {
            return Err(bad("annulus is outside its plane-pressure frequency envelope").into());
        }
        let viscous = self.diffusion_solution(
            omega,
            transport.kinematic_viscosity_m2_s(),
            admission,
            &mut cancelled,
        )?;
        let thermal = if state.heat_capacity_ratio() == 1.0 {
            None
        } else {
            Some(self.diffusion_solution(
                omega,
                transport.thermal_diffusivity_m2_s(),
                admission,
                &mut cancelled,
            )?)
        };
        let wave = self.wave_from_diffusion(
            &state,
            omega,
            viscous.response,
            thermal.as_ref().map(|t| t.response),
            admission,
        )?;
        Ok(PreparedAnnularFields::from_solutions(
            self,
            material,
            state,
            omega,
            wave,
            viscous,
            thermal,
            None,
            &mut cancelled,
        )?)
    }
}

impl<'g, 'm, M: ?Sized> PreparedAnnularFields<'g, 'm, M> {
    // Retain already-admitted native solutions. The ordinary and molecular
    // preparations share all field normalization and positive-form integration.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn from_solutions(
        section: &'g AnnularSection,
        source: &'m M,
        state: ThermoviscousFluidState,
        omega: f64,
        wave: AnnularWave,
        viscous: RadialSolution,
        thermal: Option<RadialSolution>,
        molecular: Option<(C64, f64)>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        let acoustic = state.acoustic();
        let transport = state.transport();
        let velocity_scale = viscous.response.mean.conj().recip() / C64::from_re(section.area_m2);
        let mut loss_scales = [
            det::sqrt(transport.dynamic_viscosity_pa_s()) * core::f64::consts::FRAC_1_SQRT_2,
            (det::sqrt(state.heat_capacity_ratio() - 1.0)
                / det::sqrt(acoustic.adiabatic_bulk_modulus_pa()))
                * det::sqrt(transport.thermal_diffusivity_m2_s())
                * core::f64::consts::FRAC_1_SQRT_2,
            0.0,
        ];
        if let Some((cp, temperature)) = molecular {
            // T'/p = conj(h_t)/(rho Cp). Weight before squaring and keep
            // Fourier-gradient and reservoir-temperature work independent.
            loss_scales[1] = det::sqrt(transport.thermal_conductivity_w_m_k())
                / det::sqrt(temperature)
                / acoustic.density_kg_m3()
                / cp.abs()
                * core::f64::consts::FRAC_1_SQRT_2;
            loss_scales[2] = det::sqrt(omega) * det::sqrt(cp.im)
                / det::sqrt(acoustic.density_kg_m3())
                / det::sqrt(temperature)
                / cp.abs()
                * core::f64::consts::FRAC_1_SQRT_2;
        }
        if !finite(velocity_scale)
            || velocity_scale == C64::ZERO
            || !positive(loss_scales[0])
            || !loss_scales[1].is_finite()
            || (thermal.is_some() && loss_scales[1] <= 0.0)
            || !loss_scales[2].is_finite()
            || loss_scales[2] < 0.0
        {
            return Err(bad(
                "annular field normalization or loss weights are unrepresentable",
            ));
        }
        let mut prepared = PreparedAnnularFields {
            section,
            source,
            state,
            omega,
            wave,
            viscous,
            thermal,
            velocity_scale,
            loss_scales,
            total_loss_roots: [0.0; 3],
            molecular_cp: molecular.map(|(cp, _)| cp),
        };
        let mut sums = [Dd::ZERO; 3];
        for cell in 0..section.elements() {
            if cell % 64 == 0 {
                poll(cancelled, cell)?;
            }
            for (sum, coefficient) in sums.iter_mut().zip(prepared.cell_loss_roots(cell)?) {
                *sum = *sum + Dd::from_f64(coefficient.norm_sq());
            }
        }
        prepared.total_loss_roots = sums.map(|sum| det::sqrt(sum.to_f64()));
        if !positive(prepared.total_loss_roots[0])
            || !prepared.total_loss_roots[1].is_finite()
            || (prepared.thermal.is_some() && prepared.total_loss_roots[1] <= 0.0)
            || !prepared.total_loss_roots[2].is_finite()
        {
            return Err(bad(
                "annular total physical loss weights are unrepresentable",
            ));
        }
        poll(cancelled, section.elements())?;
        Ok(prepared)
    }

    pub(super) fn with_source<'n, N: ?Sized>(
        self,
        source: &'n N,
    ) -> PreparedAnnularFields<'g, 'n, N> {
        PreparedAnnularFields {
            section: self.section,
            source,
            state: self.state,
            omega: self.omega,
            wave: self.wave,
            viscous: self.viscous,
            thermal: self.thermal,
            velocity_scale: self.velocity_scale,
            loss_scales: self.loss_scales,
            total_loss_roots: self.total_loss_roots,
            molecular_cp: self.molecular_cp,
        }
    }

    /// Retained original material owner, including its own source receipts.
    pub const fn source(&self) -> &'m M {
        self.source
    }
    /// Original geometry and radial resolution; observation never remeshes it.
    pub const fn section(&self) -> &'g AnnularSection {
        self.section
    }
    /// Original coherent equilibrium acoustic/transport snapshot. Molecular
    /// preparation additionally retains its original frequency-dependent Cp.
    pub const fn state(&self) -> &ThermoviscousFluidState {
        &self.state
    }
    /// Angular frequency binding the local pressure/flow amplitudes [rad/s].
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    /// Exactly the same average wave as the original ordinary/molecular evaluator.
    pub const fn wave(&self) -> &AnnularWave {
        &self.wave
    }

    pub(super) fn cell_loss_roots(&self, cell: usize) -> Result<[C64; 3], DuctError> {
        let stiffness = self
            .section
            .stiffness_edges
            .get(cell)
            .ok_or_else(|| bad("annular radial cell is outside the prepared section"))?;
        // Integral 2*pi*r*|dh/dr|² dr = A/gap² * K_cell * |h1-h0|².
        // Use original stored differences, including when h has rounded to one.
        let geometry = det::sqrt(self.section.area_m2) * det::sqrt(*stiffness) / self.section.gap_m;
        let coefficients = [
            (self.viscous.difference(cell).conj() * self.velocity_scale)
                .scale(geometry * self.loss_scales[0]),
            self.thermal.as_ref().map_or(C64::ZERO, |t| {
                t.difference(cell)
                    .conj()
                    .scale(geometry * self.loss_scales[1])
            }),
            self.molecular_root(cell),
        ];
        if !positive(geometry) || !coefficients.into_iter().all(finite) {
            return Err(bad("annular radial cell heat is unrepresentable"));
        }
        Ok(coefficients)
    }

    fn molecular_root(&self, cell: usize) -> C64 {
        if self.loss_scales[2] == 0.0 {
            return C64::ZERO;
        }
        let thermal = self
            .thermal
            .as_ref()
            .expect("molecular gas has a thermal field");
        let section = self.section;
        let width = section.fractions[cell + 1] - section.fractions[cell];
        let mean_radius = section.inner_m + 0.5 * section.gap_m;
        let radii = [cell, cell + 1].map(|node| {
            section.inner_m / mean_radius + (section.gap_m / mean_radius) * section.fractions[node]
        });
        let (diagonal, off) = radial_mass(width, radii[0], radii[1]);
        // Positive Cholesky form of the SAME weighted P1 mass used in the
        // diffusion solve. Unlike gradient work, reservoir loss varies inside
        // each ring; a midpoint or area-average temperature square is not exact.
        let a = det::sqrt(diagonal[0]);
        let b = off / a;
        let c = det::sqrt(diagonal[1] - b * b);
        let scale = det::sqrt(section.area_m2) * self.loss_scales[2];
        let left = thermal.node(cell);
        let right = thermal.node(cell + 1);
        let first = left.scale(a * scale) + right.scale(b * scale);
        let second = right.scale(c * scale);
        C64::from_re(det::hypot(first.abs(), second.abs()))
    }

    pub(super) fn total_loss_roots(&self) -> [C64; 3] {
        self.total_loss_roots.map(C64::from_re)
    }

    /// Actual ring cross-sectional area [m²], retaining both physical radii.
    pub fn cell_area_m2(&self, cell: usize) -> Result<f64, DuctError> {
        let stiffness = self
            .section
            .stiffness_edges
            .get(cell)
            .ok_or_else(|| bad("annular radial cell is outside the prepared section"))?;
        let width = self.section.fractions[cell + 1] - self.section.fractions[cell];
        let area = self.section.area_m2 * (*stiffness * width) * width;
        if !positive(area) {
            return Err(bad("annular radial cell area is unrepresentable"));
        }
        Ok(area)
    }

    /// Observe an original radial cell at fraction 0..=1. At shared nodes,
    /// velocity is continuous but its P1 gradient depends on the selected cell.
    /// Neither the axial solver nor the radial mesh is refined by observation.
    pub fn field_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<AnnularField, DuctError> {
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(bad(
                "annular radial interpolation fraction must be in [0,1]",
            ));
        }
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let area = self.cell_area_m2(cell)?;
        let left = self.section.fractions[cell];
        let width = self.section.fractions[cell + 1] - left;
        let interpolate = |field: &RadialSolution| {
            (field.node(cell).scale(1.0 - fraction) + field.node(cell + 1).scale(fraction)).conj()
        };
        let velocity = interpolate(&self.viscous) * self.velocity_scale * volume_flow_m3_s;
        let gradient = (self.viscous.difference(cell).conj() * self.velocity_scale)
            / C64::from_re(self.section.gap_m * width)
            * volume_flow_m3_s;
        let thermal_pressure = self.thermal.as_ref().map(|t| interpolate(t) * pressure_pa);
        let roots = self
            .cell_loss_roots(cell)?
            .map(|c| c.scale(1.0 / det::sqrt(area)));
        let heat = losses(roots, pressure_pa, volume_flow_m3_s)?;
        let molecular = thermal_pressure
            .unwrap_or(C64::ZERO)
            .scale(self.loss_scales[2])
            .norm_sq();
        if !finite(velocity)
            || !finite(gradient)
            || thermal_pressure.is_some_and(|t| !finite(t))
            || !molecular.is_finite()
        {
            return Err(bad("annular physical field is unrepresentable"));
        }
        Ok(AnnularField {
            radius_m: self.section.inner_m + self.section.gap_m * (left + fraction * width),
            pressure_pa,
            axial_velocity_m_s: velocity,
            radial_velocity_gradient_per_s: gradient,
            thermal_pressure_response_pa: thermal_pressure,
            viscous_dissipation_w_m3: heat.viscous_w_m,
            thermal_dissipation_w_m3: heat.thermal_w_m,
            molecular_dissipation_w_m3: molecular,
        })
    }

    /// Separate shear, Fourier and molecular work in one complete radial cell,
    /// per metre of duct length. Cross-section pressure and flow are coherent.
    pub fn cell_dissipation_per_axial_length(
        &self,
        cell: usize,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<AnnularDissipationPerLength, DuctError> {
        losses(self.cell_loss_roots(cell)?, pressure_pa, volume_flow_m3_s)
    }

    /// The same positive work over the whole cross section, with O(1) cost.
    /// Prepared weights sum the ORIGINAL gradient and reservoir mass forms;
    /// no boundary-flux subtraction or attenuation fit supplies these losses.
    pub fn dissipation_per_axial_length(
        &self,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<AnnularDissipationPerLength, DuctError> {
        losses(
            self.total_loss_roots.map(C64::from_re),
            pressure_pa,
            volume_flow_m3_s,
        )
    }
}

impl PreparedAnnularFields<'_, '_, GasState> {
    /// Temperature perturbation for the original ideal-gas material [K].
    /// Its adiabatic derivative is 1/(rho Cp); no liquid derivative is guessed.
    pub fn gas_temperature_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        pressure_pa: C64,
    ) -> Result<C64, DuctError> {
        let field = self.field_in_cell(cell, fraction, pressure_pa, C64::ZERO)?;
        let temperature = field.thermal_pressure_response_pa.unwrap_or(C64::ZERO)
            / C64::from_re(self.state.acoustic().density_kg_m3())
            / C64::from_re(self.state.transport().specific_heat_j_kg_k());
        if !finite(temperature) {
            return Err(bad("annular gas temperature field is unrepresentable"));
        }
        Ok(temperature)
    }
}

impl PreparedAnnularFields<'_, '_, RelaxationState> {
    /// Original dynamic Cp in exp(-i omega t), cached at this frequency.
    /// An empty reservoir list retains the original constant-Cp path exactly.
    pub fn complex_cp_j_kg_k(&self) -> C64 {
        self.molecular_cp
            .unwrap_or_else(|| C64::from_re(self.state.transport().specific_heat_j_kg_k()))
    }

    /// Translational gas temperature with the original molecular heat storage.
    /// This is a first-order phasor, not a mean temperature or deposited heat.
    pub fn gas_temperature_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        pressure_pa: C64,
    ) -> Result<C64, DuctError> {
        let field = self.field_in_cell(cell, fraction, pressure_pa, C64::ZERO)?;
        let temperature = field.thermal_pressure_response_pa.unwrap_or(C64::ZERO)
            / C64::from_re(self.state.acoustic().density_kg_m3())
            / self.complex_cp_j_kg_k();
        if !finite(temperature) {
            return Err(bad(
                "annular molecular gas temperature field is unrepresentable",
            ));
        }
        Ok(temperature)
    }
}

fn amplitudes(pressure: C64, flow: C64) -> Result<(), DuctError> {
    if !finite(pressure) || !finite(flow) {
        Err(bad(
            "annular section pressure and volume flow must be finite",
        ))
    } else {
        Ok(())
    }
}

fn losses(
    roots: [C64; 3],
    pressure: C64,
    flow: C64,
) -> Result<AnnularDissipationPerLength, DuctError> {
    amplitudes(pressure, flow)?;
    let viscous = (roots[0] * flow).norm_sq();
    let thermal = (roots[1] * pressure).norm_sq();
    let molecular = (roots[2] * pressure).norm_sq();
    let total = viscous + thermal + molecular;
    if ![viscous, thermal, molecular, total]
        .into_iter()
        .all(f64::is_finite)
    {
        return Err(bad("annular dissipated power is unrepresentable"));
    }
    Ok(AnnularDissipationPerLength {
        viscous_w_m: viscous,
        thermal_w_m: thermal,
        molecular_w_m: molecular,
        total_w_m: total,
    })
}
