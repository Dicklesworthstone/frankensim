//! Quasisteady viscous flow through a uniform rectangular moving slit.
//!
//! For width w, gap h, streamwise length L and mean speed v, the admitted
//! reduced law is dp = 12 mu L v / h^2 + rho v |v| / 2, with U = w h v.
//! The first term is fully developed, no-slip, parallel-plate resistance;
//! the second is the existing separated-jet loss, with no pressure recovery.
//! Width is the physical slit width, not a fitted discharge coefficient.
//!
//! A P1 gap varying ACROSS the width is integrated as parallel local channels,
//! including its exact positive portion through closure. No gap is replaced by
//! area / width, which would change the h^-3 resistance. The streamwise length
//! is common; this does not resolve transverse flow, unsteady boundary layers,
//! entrance flow, pressure recovery, compressibility or turbulence.
//! Material parameters are explicit; existing inviscid callers are unchanged.

use crate::acoustic_realize::AcousticRealizeError;
use fs_material::gas::GasState;
use fs_math::det;

/// Explicit fluid viscosity and uniform slit length. Both are strictly positive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousSlit {
    dynamic_viscosity_pa_s: f64,
    length_m: f64,
    gas: Option<GasState>,
}

/// Signed flow, local tangents and nonnegative loss of the reduced slit law.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ViscousSlitFlow {
    /// Signed volume flow from upstream to downstream [m^3/s].
    pub flow_m3_s: f64,
    /// dU/d(dp) at fixed physical gap and width [m^3/(s Pa)].
    pub flow_per_pressure_m3_s_pa: f64,
    /// dU/dh at fixed pressure, width, viscosity and length [m^2/s].
    pub flow_per_opening_m2_s: f64,
    /// Fully developed wall-friction loss [W].
    pub viscous_dissipation_w: f64,
    /// Separated-jet loss [W], without pressure recovery.
    pub separation_dissipation_w: f64,
}

/// Integrated flow of parallel channels across a nonuniform physical slit.
/// Its single pressure drop supplies nonnegative loss `dp * flow_m3_s`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ViscousProfileFlow {
    /// Signed physical volume flow [m³/s].
    pub flow_m3_s: f64,
    /// Nonnegative derivative of integrated flow with pressure [m³/(s Pa)].
    pub flow_per_pressure_m3_s_pa: f64,
}

impl ViscousSlitFlow {
    /// Sum of the two physical losses [W]. Equal to dp * U up to rounding.
    #[must_use]
    pub fn dissipated_power_w(self) -> f64 {
        self.viscous_dissipation_w + self.separation_dissipation_w
    }
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

fn numerical() -> AcousticRealizeError {
    AcousticRealizeError::Reed {
        what: "viscous slit coefficients or result cannot be represented as finite physical values",
    }
}

impl ViscousSlit {
    /// Admit SI viscosity and streamwise slit length, without inferred defaults.
    ///
    /// # Errors
    /// Nonpositive/nonfinite inputs or an unrepresentable resistance numerator.
    pub fn new(dynamic_viscosity_pa_s: f64, length_m: f64)
        -> Result<Self, AcousticRealizeError>
    {
        if !dynamic_viscosity_pa_s.is_finite() || dynamic_viscosity_pa_s <= 0.0
            || !length_m.is_finite() || length_m <= 0.0
        {
            return Err(invalid("viscous slit requires positive finite viscosity and streamwise length"));
        }
        let numerator = 12.0 * (dynamic_viscosity_pa_s * length_m);
        if !numerator.is_finite() || numerator <= 0.0 { return Err(numerical()); }
        Ok(Self { dynamic_viscosity_pa_s, length_m, gas: None })
    }

    /// Bind viscosity and density to one complete material gas snapshot.
    /// Composition, pressure and temperature therefore reach the channel
    /// without an air constant; its streamwise length remains explicit geometry.
    pub fn from_gas(gas: &GasState, length_m: f64) -> Result<Self, AcousticRealizeError> {
        gas.validate_coherence().map_err(|_| invalid("viscous slit requires one coherent source gas"))?;
        let mut law = Self::new(gas.dynamic_viscosity, length_m)?;
        law.gas = Some(*gas);
        Ok(law)
    }

    /// Original material state, when constructed from a complete source gas.
    #[must_use]
    pub const fn gas(&self) -> Option<&GasState> { self.gas.as_ref() }

    // A sourced channel follows the same new complete gas as its acoustic
    // owner. An explicitly supplied viscosity remains an explicit override.
    // The quasisteady law carries no fluid momentum or independent storage.
    pub(super) fn rebind_source_gas(self, gas: &GasState)
        -> Result<Self, AcousticRealizeError>
    {
        if self.gas.is_some() { Self::from_gas(gas, self.length_m) } else { Ok(self) }
    }

    fn admit_density(self, density: f64) -> Result<(), AcousticRealizeError> {
        if self.gas.is_some_and(|gas| (gas.density - density).abs()
            > 64.0 * f64::EPSILON * gas.density.abs().max(density.abs())) {
            return Err(invalid("viscous slit viscosity and density must use the same source gas"));
        }
        Ok(())
    }

    /// Supplied Newtonian dynamic viscosity [Pa s].
    #[must_use]
    pub const fn dynamic_viscosity_pa_s(self) -> f64 { self.dynamic_viscosity_pa_s }

    /// Supplied streamwise length [m].
    #[must_use]
    pub const fn length_m(self) -> f64 { self.length_m }

    fn linear_pressure(self, opening: f64) -> Result<f64, AcousticRealizeError> {
        let value = (12.0 * (self.dynamic_viscosity_pa_s * self.length_m) / opening) / opening;
        if !value.is_finite() || value <= 0.0 { return Err(numerical()); }
        Ok(value)
    }

    /// Evaluate dp = 12 mu L U / (w h^3) + rho U |U| / (2 w^2 h^2).
    /// The open-slit pressure tangent is finite even at dp = 0. A geometrically
    /// closed slit returns exactly zero flow and both zero tangents.
    ///
    /// # Errors
    /// Invalid geometry, density or pressure, or an unrepresentable result.
    pub fn evaluate(self, width_m: f64, opening_m: f64, density_kg_m3: f64, pressure_drop_pa: f64)
        -> Result<ViscousSlitFlow, AcousticRealizeError>
    {
        let area = admitted_area(width_m, opening_m, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if !pressure_drop_pa.is_finite() {
            return Err(invalid("viscous slit pressure drop must be finite"));
        }
        if opening_m == 0.0 { return Ok(ViscousSlitFlow::default()); }
        let linear = self.linear_pressure(opening_m)?;
        let velocity = solve_velocity(linear, density_kg_m3, pressure_drop_pa)?;
        let denominator = linear + density_kg_m3 * velocity.abs();
        if !denominator.is_finite() || denominator <= 0.0 { return Err(numerical()); }
        let flow = area * velocity;
        // Differentiating in mean velocity avoids divisions by U at rest and
        // avoids forming R ~ h^-3 and the cancellation-prone quadratic root.
        let result = ViscousSlitFlow {
            flow_m3_s: flow,
            flow_per_pressure_m3_s_pa: area / denominator,
            flow_per_opening_m2_s: (width_m * velocity) * (1.0 + 2.0 * (linear / denominator)),
            viscous_dissipation_w: (linear * velocity) * flow,
            separation_dissipation_w: ((0.5 * density_kg_m3 * velocity) * velocity.abs()) * flow,
        };
        if ![result.flow_m3_s, result.flow_per_pressure_m3_s_pa,
            result.flow_per_opening_m2_s, result.viscous_dissipation_w,
            result.separation_dissipation_w, result.dissipated_power_w()]
            .iter().all(|v| v.is_finite())
        {
            return Err(numerical());
        }
        Ok(result)
    }

    /// Solve the SAME-sample characteristic junction, including supplied swept
    /// and body flow, rather than filtering an inviscid result afterward.
    ///
    /// With A = w h and d = Pmouth - 2 Pminus - Z Umoving, solve
    /// (12 mu L / h^2 + Z A) v + rho v |v| / 2 = d. Return Pplus [Pa].
    /// No nonlinear iteration, state mutation or allocation occurs here.
    ///
    /// # Errors
    /// Invalid physical inputs, nonfinite intermediate pressure or overflow.
    #[allow(clippy::too_many_arguments)]
    pub fn characteristic_pressure(
        self, width_m: f64, opening_m: f64, density_kg_m3: f64,
        impedance_pa_s_m3: f64, incoming_pressure_pa: f64,
        upstream_pressure_pa: f64, moving_flow_m3_s: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let area = admitted_area(width_m, opening_m, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if ![impedance_pa_s_m3, incoming_pressure_pa, upstream_pressure_pa, moving_flow_m3_s]
            .iter().all(|v| v.is_finite()) || impedance_pa_s_m3 <= 0.0
        {
            return Err(invalid("viscous slit junction requires finite pressures/flow and positive impedance"));
        }
        let no_jet = incoming_pressure_pa + impedance_pa_s_m3 * moving_flow_m3_s;
        let no_drop = upstream_pressure_pa - incoming_pressure_pa;
        let drive = no_drop - no_jet;
        if ![no_jet, no_drop, drive].iter().all(|v| v.is_finite()) { return Err(numerical()); }
        if opening_m == 0.0 || drive == 0.0 { return Ok(no_jet); }
        let linear = self.linear_pressure(opening_m)?;
        let wave_load = impedance_pa_s_m3 * area;
        let velocity = solve_velocity(linear + wave_load, density_kg_m3, drive)?;
        let jet_pressure = wave_load * velocity;
        // Keep the constitutive drop: reconstructing it as drive - Z U would
        // erase small pressure drops when the wave load dominates resistance.
        let dp = linear * velocity + (0.5 * density_kg_m3 * velocity) * velocity.abs();
        let outgoing = if no_drop.abs().max(dp.abs()) <= no_jet.abs().max(jet_pressure.abs()) {
            no_drop - dp
        } else {
            no_jet + jet_pressure
        };
        if ![dp, jet_pressure, outgoing].iter().all(|v| v.is_finite()) { return Err(numerical()); }
        Ok(outgoing)
    }

    /// Eliminate dp + Z U(dp) = drive after condensing mechanical swept flow.
    /// Return dp directly: subtracting large wave pressures can erase it.
    pub(super) fn fixed_opening_pressure_drop(
        self, width_m: f64, opening_m: f64, density_kg_m3: f64,
        impedance_pa_s_m3: f64, drive_pa: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let area = admitted_area(width_m, opening_m, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if !impedance_pa_s_m3.is_finite() || impedance_pa_s_m3 <= 0.0 || !drive_pa.is_finite() {
            return Err(invalid("viscous slit condensation requires positive impedance and finite drive"));
        }
        if opening_m == 0.0 || drive_pa == 0.0 { return Ok(drive_pa); }
        let linear = self.linear_pressure(opening_m)?;
        let velocity = solve_velocity(linear + impedance_pa_s_m3 * area, density_kg_m3, drive_pa)?;
        let dp = linear * velocity + (0.5 * density_kg_m3 * velocity) * velocity.abs();
        if !dp.is_finite() { return Err(numerical()); }
        Ok(dp)
    }

    /// Integrate this same law over a linearly varying widthwise gap. Negative
    /// gaps contribute no flow; a partially closed edge retains its actual open
    /// fraction. The pressure derivative is integrated analytically as well.
    ///
    /// With K=12 mu L and C=2 rho |dp|, flow per width is
    /// `(sqrt(K²+C h⁴)-K)/(rho h)`. Its primitive on h>0 is
    /// `[r-K ln(K+r)]/(2 rho)`, r=sqrt(K²+C h⁴). The divided difference below
    /// rationalizes delta-r and log1p so nearly equal gaps and the laminar limit
    /// do not subtract nearly equal primitives. This is exact for the declared
    /// parallel-channel P1 model, not a spatial Navier--Stokes solution.
    pub fn linear_gap_flow(self, width_m: f64, gaps_m: [f64; 2], density_kg_m3: f64,
        pressure_drop_pa: f64) -> Result<ViscousProfileFlow, AcousticRealizeError>
    {
        if gaps_m.iter().any(|g| !g.is_finite()) || !pressure_drop_pa.is_finite() {
            return Err(invalid("viscous slit profile requires finite gaps and pressure"));
        }
        let lower = gaps_m[0].min(gaps_m[1]);
        let upper = gaps_m[0].max(gaps_m[1]).max(0.0);
        admitted_area(width_m, upper, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if upper == 0.0 { return Ok(ViscousProfileFlow::default()); }
        if lower == upper {
            let flow = self.evaluate(width_m, upper, density_kg_m3, pressure_drop_pa)?;
            return Ok(ViscousProfileFlow { flow_m3_s: flow.flow_m3_s,
                flow_per_pressure_m3_s_pa: flow.flow_per_pressure_m3_s_pa });
        }
        let open_fraction = if lower >= 0.0 { 1.0 } else {
            let scale = upper.max(-lower);
            (upper / scale) / (upper / scale + (-lower) / scale)
        };
        let ratio = lower.max(0.0) / upper;
        let linear = self.linear_pressure(upper)?;
        let quadratic = (std::f64::consts::SQRT_2 * det::sqrt(density_kg_m3))
            * det::sqrt(pressure_drop_pa.abs());
        if !quadratic.is_finite() { return Err(numerical()); }
        let scale = linear.max(quadratic);
        let (l, s) = (linear / scale, quadratic / scale);
        let r0 = det::hypot(l, s * ratio * ratio);
        let r1 = det::hypot(l, s);
        let polynomial = (1.0 + ratio) * (1.0 + ratio * ratio);
        let tangent = (width_m * open_fraction) * (upper / scale)
            * (polynomial / (2.0 * (r0 + r1)));
        let delta = s * s * (1.0 - ratio) * polynomial / (r0 + r1);
        let log_argument = delta / (l + r0);
        let log_ratio = if !log_argument.is_finite() { 0.0 } else {
            let sum = 1.0 + log_argument;
            if sum == 1.0 { 1.0 } else { det::ln(sum) / (sum - 1.0) }
        };
        let factor = if l == 0.0 { 1.0 } else { 1.0 - (l / (l + r0)) * log_ratio };
        let flow = (pressure_drop_pa * tangent) * (2.0 * factor);
        if !flow.is_finite() || !tangent.is_finite() || tangent < 0.0 {
            return Err(numerical());
        }
        Ok(ViscousProfileFlow { flow_m3_s: flow, flow_per_pressure_m3_s_pa: tangent })
    }
}

/// Shared monotone passive pressure solve for a supplied spatial slit integral.
/// The constitutive law is odd and increasing. All 96 bounded evaluations are
/// at the same mechanical trial; no accepted state or pressure history changes.
pub(super) fn profiled_pressure_drop(
    impedance: f64, drive: f64,
    flow: impl Fn(f64) -> Result<ViscousProfileFlow, AcousticRealizeError>,
) -> Result<f64, AcousticRealizeError> {
    if !impedance.is_finite() || impedance <= 0.0 || !drive.is_finite() {
        return Err(invalid("profiled slit requires finite drive and positive passive impedance"));
    }
    if drive == 0.0 { return Ok(0.0); }
    let source = drive.abs();
    let (mut lo, mut hi, mut p) = (0.0, source, source);
    for _ in 0..96 {
        let value = flow(p)?;
        let load = impedance * value.flow_m3_s;
        let derivative = 1.0 + impedance * value.flow_per_pressure_m3_s_pa;
        let residual = p + load - source;
        if ![load, derivative, residual].iter().all(|x| x.is_finite())
            || value.flow_m3_s < 0.0 || value.flow_per_pressure_m3_s_pa < 0.0 {
            return Err(numerical());
        }
        if residual.abs() <= 32.0 * f64::EPSILON * source {
            return Ok(drive.signum() * p);
        }
        if residual > 0.0 { hi = p; } else { lo = p; }
        let newton = p - residual / derivative;
        p = if newton > lo && newton < hi { newton } else { f64::midpoint(lo, hi) };
        if p == lo || p == hi { break; }
    }
    Err(invalid("profiled viscous slit pressure balance did not converge within its bounded work"))
}

/// Same spatial pressure solve before mechanical swept-flow condensation.
pub(super) fn profiled_characteristic_pressure(
    impedance: f64, incoming: f64, upstream: f64, moving: f64,
    flow: impl Fn(f64) -> Result<ViscousProfileFlow, AcousticRealizeError>,
) -> Result<f64, AcousticRealizeError> {
    let no_jet = incoming + impedance * moving;
    let no_drop = upstream - incoming;
    let dp = profiled_pressure_drop(impedance, no_drop - no_jet, &flow)?;
    let jet_pressure = impedance * flow(dp)?.flow_m3_s;
    let outgoing = if no_drop.abs().max(dp.abs()) <= no_jet.abs().max(jet_pressure.abs()) {
        no_drop - dp
    } else { no_jet + jet_pressure };
    if !outgoing.is_finite() { return Err(numerical()); }
    Ok(outgoing)
}

fn admitted_area(width: f64, opening: f64, rho: f64) -> Result<f64, AcousticRealizeError> {
    if ![width, opening, rho].iter().all(|v| v.is_finite())
        || width <= 0.0 || opening < 0.0 || rho <= 0.0
    {
        return Err(invalid("viscous slit requires positive finite width/density and a nonnegative finite gap"));
    }
    let area = width * opening;
    if !area.is_finite() || (opening > 0.0 && area == 0.0) { return Err(numerical()); }
    Ok(area)
}

fn solve_velocity(linear: f64, rho: f64, drive: f64) -> Result<f64, AcousticRealizeError> {
    if !linear.is_finite() || linear <= 0.0 { return Err(numerical()); }
    if drive == 0.0 { return Ok(0.0); }
    // v = 2 |d| / (a + sqrt(a^2 + 2 rho |d|)), rationalized and scaled.
    // sqrt(|d|) keeps sub-threshold pressures continuous; a is never squared.
    let quadratic = (std::f64::consts::SQRT_2 * det::sqrt(rho)) * det::sqrt(drive.abs());
    if !quadratic.is_finite() { return Err(numerical()); }
    let scale = linear.max(quadratic);
    let an = linear / scale;
    let bn = quadratic / scale;
    let denominator = 0.5 * (an + det::sqrt(an * an + bn * bn));
    let velocity = drive.signum() * ((drive.abs() / scale) / denominator);
    if !velocity.is_finite() { return Err(numerical()); }
    Ok(velocity)
}

#[cfg(test)]
#[allow(clippy::float_cmp)] // Exact closure, rest, symmetry and refusal pins.
mod tests {
    use super::*;

    fn law() -> ViscousSlit { ViscousSlit::new(1.81e-5, 1e-3).unwrap() }

    fn close(actual: f64, expected: f64, relative: f64) {
        assert!((actual - expected).abs() <= relative * expected.abs().max(1e-250),
            "actual={actual:e}, expected={expected:e}");
    }

    #[test]
    fn g3_material_gas_supplies_viscosity_and_cannot_be_mixed_with_another_density() {
        use fs_material::gas::GasSpec;
        let cold = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
        let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 393.15, 101325.0).unwrap();
        let a = ViscousSlit::from_gas(&cold, 0.001).unwrap();
        let b = ViscousSlit::from_gas(&hot, 0.001).unwrap();
        assert_eq!(a.gas(), Some(&cold));
        assert_eq!(a.dynamic_viscosity_pa_s(), cold.dynamic_viscosity);
        let flow = |law: ViscousSlit, gas: GasState| law.evaluate(0.01, 1e-5, gas.density, 100.0).unwrap().flow_m3_s;
        assert!(flow(b, hot) < flow(a, cold), "a thin channel must respond to material viscosity, not only density");
        assert!(a.evaluate(0.01, 1e-5, hot.density, 100.0).is_err());
        let mut corrupt = cold; corrupt.dynamic_viscosity *= 2.0;
        assert!(ViscousSlit::from_gas(&corrupt, 0.001).is_err());
    }

    #[test]
    fn g1_spatial_laminar_limit_integrates_cubic_gap_and_partial_closure() {
        let (width, h, dp) = (0.013, 4e-5, 1e-8);
        let k = 12.0 * 1.81e-5 * 1e-3;
        let triangular = law().linear_gap_flow(width, [0.0, 2.0*h], 1.2, dp).unwrap();
        let mean_gap = law().evaluate(width, h, 1.2, dp).unwrap();
        close(triangular.flow_m3_s, width*dp*(2.0*h).powi(3)/(4.0*k), 1e-10);
        close(triangular.flow_m3_s/mean_gap.flow_m3_s, 2.0, 1e-10);
        let partial = law().linear_gap_flow(width, [-2.0*h, 2.0*h], 1.2, dp).unwrap();
        close(partial.flow_m3_s, 0.5*triangular.flow_m3_s, 1e-14);
        let tangent = law().linear_gap_flow(width, [0.0, 2.0*h], 1.2, 0.0).unwrap();
        close(tangent.flow_per_pressure_m3_s_pa, width*(2.0*h).powi(3)/(4.0*k), 1e-14);
        assert_eq!(law().linear_gap_flow(width, [-h, 0.0], 1.2, dp).unwrap(), ViscousProfileFlow::default());
    }

    #[test]
    fn g1_spatial_integral_matches_independent_point_quadrature_and_pressure_tangent() {
        let (width, rho) = (0.013, 1.2);
        for gaps in [[-4e-5, 8e-5], [1e-5, 1e-4], [1e-4, 1e-4*(1.0+1e-14)]] {
            for dp in [-2000.0, 0.1, 2000.0] {
                let exact = law().linear_gap_flow(width, gaps, rho, dp).unwrap();
                let mut numerical = 0.0;
                for i in 0..8192 {
                    let t = (i as f64 + 0.5)/8192.0;
                    let h = ((1.0-t)*gaps[0]+t*gaps[1]).max(0.0);
                    numerical += law().evaluate(width/8192.0, h, rho, dp).unwrap().flow_m3_s;
                }
                close(exact.flow_m3_s, numerical, 4e-8);
                let delta = dp.abs()*1e-5;
                let left = law().linear_gap_flow(width, gaps, rho, dp-delta).unwrap();
                let right = law().linear_gap_flow(width, gaps, rho, dp+delta).unwrap();
                close(exact.flow_per_pressure_m3_s_pa, (right.flow_m3_s-left.flow_m3_s)/(2.0*delta), 2e-9);
                let mid = f64::midpoint(gaps[0], gaps[1]);
                let a = law().linear_gap_flow(width/2.0, [gaps[0],mid], rho, dp).unwrap();
                let b = law().linear_gap_flow(width/2.0, [mid,gaps[1]], rho, dp).unwrap();
                close(exact.flow_m3_s, a.flow_m3_s+b.flow_m3_s, 2e-13);
                let z = 1e6;
                let solved = profiled_pressure_drop(z, dp+z*exact.flow_m3_s,
                    |p| law().linear_gap_flow(width, gaps, rho, p)).unwrap();
                close(solved, dp, 1e-12);
            }
        }
    }

    #[test]
    fn g1_manufactured_signed_flows_recover_both_pressure_losses() {
        let (width, rho) = (0.013, 1.2);
        for h in [1e-7, 1e-5, 4e-4] {
            for velocity in [-20.0_f64, -0.25, 0.0, 0.25, 20.0] {
                let pv = 12.0 * 1.81e-5 * 1e-3 * velocity / (h * h);
                let ps = 0.5 * rho * velocity * velocity.abs();
                let actual = law().evaluate(width, h, rho, pv + ps).unwrap();
                let expected = width * h * velocity;
                close(actual.flow_m3_s, expected, 3e-14);
                close(actual.viscous_dissipation_w, pv * expected, 6e-14);
                close(actual.separation_dissipation_w, ps * expected, 6e-14);
                close(actual.dissipated_power_w(), (pv + ps) * expected, 6e-14);
                assert!(actual.viscous_dissipation_w >= 0.0);
                assert!(actual.separation_dissipation_w >= 0.0);
            }
        }
    }

    #[test]
    fn g1_laminar_limit_has_finite_rest_tangent_and_cubic_gap_scaling() {
        let (width, h, rho) = (0.013, 1e-5, 1.2);
        let at_rest = law().evaluate(width, h, rho, 0.0).unwrap();
        close(at_rest.flow_per_pressure_m3_s_pa, width * h * h * h / (12.0 * 1.81e-5 * 1e-3), 3e-14);
        assert_eq!(at_rest.flow_m3_s, 0.0);
        assert_eq!(at_rest.flow_per_opening_m2_s, 0.0);
        let first = law().evaluate(width, h, rho, 1e-6).unwrap();
        let doubled = law().evaluate(width, 2.0 * h, rho, 1e-6).unwrap();
        close(doubled.flow_m3_s / first.flow_m3_s, 8.0, 1e-8);
    }

    #[test]
    fn g1_low_viscosity_recovers_the_separated_jet_limit() {
        let (width, h, rho, dp) = (0.013, 4e-4, 1.2_f64, 250.0_f64);
        let inviscid = width * h * (2.0 * dp / rho).sqrt();
        let nearly_inviscid = ViscousSlit::new(1e-12, 1e-5).unwrap().evaluate(width, h, rho, dp).unwrap();
        close(nearly_inviscid.flow_m3_s, inviscid, 1e-8);
        let viscous = law().evaluate(width, h, rho, dp).unwrap();
        assert!(viscous.flow_m3_s > 0.0 && viscous.flow_m3_s < inviscid);
    }

    #[test]
    fn g3_pressure_and_gap_tangents_match_independent_differences() {
        let (width, h, rho) = (0.013, 4e-5, 1.2);
        for dp in [-250.0_f64, -1e-5, 1e-5, 250.0] {
            let value = law().evaluate(width, h, rho, dp).unwrap();
            let ep = dp.abs() * 1e-5;
            let eh = h * 1e-5;
            let flow = |gap, pressure| law().evaluate(width, gap, rho, pressure).unwrap().flow_m3_s;
            close(value.flow_per_pressure_m3_s_pa, (flow(h, dp + ep) - flow(h, dp - ep)) / (2.0 * ep), 2e-8);
            close(value.flow_per_opening_m2_s, (flow(h + eh, dp) - flow(h - eh, dp)) / (2.0 * eh), 2e-8);
        }
    }

    #[test]
    fn g1_characteristic_junction_recovers_manufactured_bidirectional_flow() {
        let (width, rho, impedance) = (0.013, 1.2, 1e6);
        for h in [1e-7, 4e-5, 4e-4] {
            for velocity in [-3.0_f64, -1e-5, 1e-5, 3.0] {
                let dp = 12.0 * 1.81e-5 * 1e-3 * velocity / (h * h)
                    + 0.5 * rho * velocity * velocity.abs();
                for incoming in [-200.0, 0.0, 200.0] {
                    for moving in [-1e-6, 0.0, 1e-6] {
                        let expected = incoming + impedance * (width * h * velocity + moving);
                        let upstream = dp + expected + incoming;
                        let actual = law().characteristic_pressure(width, h, rho, impedance, incoming, upstream, moving).unwrap();
                        assert!((actual - expected).abs() <= 3e-12 * (1.0 + expected.abs()));
                    }
                }
            }
        }
    }

    #[test]
    fn g0_exact_closure_and_arbitrarily_small_signed_drives() {
        let closed = law().evaluate(0.013, 0.0, 1.2, 500.0).unwrap();
        assert_eq!(closed, ViscousSlitFlow::default());
        assert_eq!(law().characteristic_pressure(0.013, 0.0, 1.2, 1e6, 75.0, 1000.0, 2e-7).unwrap(), 75.2);
        for dp in [1e-20, 1e-100, 1e-200] {
            let forward = law().evaluate(0.013, 1e-4, 1.2, dp).unwrap();
            let reverse = law().evaluate(0.013, 1e-4, 1.2, -dp).unwrap();
            assert!(forward.flow_m3_s > 0.0);
            assert_eq!(reverse.flow_m3_s, -forward.flow_m3_s);
            assert_eq!(reverse.flow_per_pressure_m3_s_pa, forward.flow_per_pressure_m3_s_pa);
        }
    }

    #[test]
    fn g0_invalid_and_unrepresentable_inputs_are_refused() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(ViscousSlit::new(bad, 1e-3).is_err());
            assert!(ViscousSlit::new(1.81e-5, bad).is_err());
            assert!(law().evaluate(bad, 1e-4, 1.2, 100.0).is_err());
            assert!(law().evaluate(0.013, 1e-4, bad, 100.0).is_err());
            assert!(law().characteristic_pressure(0.013, 1e-4, 1.2, bad, 0.0, 100.0, 0.0).is_err());
        }
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(law().evaluate(0.013, bad, 1.2, 100.0).is_err());
        }
        assert!(law().evaluate(0.013, 0.0, 1.2, f64::NAN).is_err());
        assert!(law().evaluate(0.013, 1e-200, 1.2, 100.0).is_err());
        assert!(ViscousSlit::new(f64::MAX, f64::MAX).is_err());
    }
}
