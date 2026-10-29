//! Cooper/CFAST exchange through a shallow horizontal ceiling/floor opening.
//!
//! Dense gas above lighter gas produces counterflow even at zero pressure
//! difference. Pressure-driven flow and density-driven exchange remain separate
//! in each direction; donor density, species and enthalpy must not be replaced
//! by a net-flow average. The original [`DischargeCoefficient`] lowers the
//! published pressure coefficient into the existing quadratic loss owner.
//!
//! Source: NIST TN 1889v1 (2015), section 4.2, equations 4.8--4.10,
//! <https://doi.org/10.6028/NIST.TN.1889v1>, following Cooper, NISTIR 89-4052,
//! <https://doi.org/10.6028/NIST.IR.89-4052>. Shape lengths follow the reference
//! implementation `firemodels/cfast`, `Source/CFAST/flowvertical.f90::ventcf`,
//! last changed at `180799b62d42a770ae86fb1a7e7becf49a4ac7dc`: circle diameter
//! `2 sqrt(A/pi)`, square side `sqrt(A)`. This distinction is ambiguous in the
//! older guide's common-D prose. No CFAST small-pressure damping is introduced.
//!
//! This is an estimated low-Mach, shallow-opening correlation. Only circular
//! and square unobstructed holes are admitted; ducts, shafts, louvers, thick
//! openings, choking, plug-holing and vent-generated plume entrainment are not
//! resolved. Transfers use the layer immediately adjacent to each side. The
//! declared resistance allowance is not propagated, and total model discrepancy
//! remains unquantified. A caller must enforce inventory positivity when using
//! these rates to advance a finite room.

use crate::stratified::{
    DischargeCoefficient, Layer, LayerEnthalpyRates, LayerMassRates, SpecificEnthalpy,
    TwoLayerColumn, VentError,
};
use crate::{SourceProvenance, ToleranceBasis};
use fs_math::det;
use fs_qty::{
    Acceleration, Area, Dimensionless, Length, MassFlowRate, Power, Pressure, VolumetricFlowRate,
};

/// Shapes supported by the original shallow-opening correlation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CeilingShape {
    /// Circular hole; the characteristic length is its diameter.
    Round,
    /// Square hole; the characteristic length is its side length.
    Square,
}
impl CeilingShape {
    fn geometry(self, area: f64) -> (f64, f64) {
        let root_area = det::sqrt(area);
        let (shape, length_ratio) = match self {
            Self::Round => (0.754, 2.0 / det::sqrt(core::f64::consts::PI)),
            Self::Square => (0.942, 1.0),
        };
        // Cs^2 D^5/(2 A^2) = Cs^2 (D/sqrt(A))^5 sqrt(A)/2.
        // This avoids huge fifth powers of an otherwise representable length.
        let ratio_squared = length_ratio * length_ratio;
        let flooding_length =
            0.5 * shape * shape * ratio_squared * ratio_squared * length_ratio * root_area;
        (root_area, flooding_length)
    }
}

/// Horizontal opening between the ceiling of the lower room and the floor of
/// the upper room. A zero area is an explicitly closed opening.
#[derive(Clone, Debug, PartialEq)]
pub struct CeilingOpening {
    height: Length,
    area: Area,
    shape: CeilingShape,
    root_area: f64,
    flooding_length: f64,
    discharge: DischargeCoefficient,
}
impl CeilingOpening {
    /// Bind the published coefficient `Cd=0.68` and actual geometric area.
    ///
    /// `relative_resistance_allowance` is a caller-chosen engineering allowance
    /// in `[0,1)`, not a published confidence interval. Zero means no additional
    /// allowance was declared; it does not assert a physically exact model.
    pub fn cooper(
        height: Length,
        area: Area,
        shape: CeilingShape,
        relative_resistance_allowance: f64,
    ) -> Result<Self, VentError> {
        if !height.value().is_finite() || !area.value().is_finite() || area.value() < 0.0 {
            return Err(VentError::Invalid(
                "finite ceiling height and nonnegative area required",
            ));
        }
        let (root_area, flooding_length) = shape.geometry(area.value());
        if !root_area.is_finite()
            || !flooding_length.is_finite()
            || (area.value() > 0.0 && (root_area <= 0.0 || flooding_length <= 0.0))
        {
            return Err(VentError::NonFinite("horizontal opening geometry"));
        }
        Ok(Self {
            height,
            area,
            shape,
            root_area,
            flooding_length,
            discharge: DischargeCoefficient::new(
                Dimensionless::new(0.68),
                SourceProvenance::new(
                    "Cooper ceiling/floor vent correlation; NIST TN 1889v1, section 4.2, equations 4.8-4.10",
                    "NIST-TN-1889v1:4.2:Cooper-horizontal-vent",
                ),
                relative_resistance_allowance,
                ToleranceBasis::EngineeringAllowance,
            )?,
        })
    }
    /// World elevation of the horizontal opening.
    pub const fn height(&self) -> Length {
        self.height
    }
    /// Actual free area, without an extra discharge-coefficient multiplier.
    pub const fn area(&self) -> Area {
        self.area
    }
    /// Original correlation's admitted aperture shape.
    pub const fn shape(&self) -> CeilingShape {
        self.shape
    }
    /// Original loss owner with the published coefficient and explicit allowance.
    pub const fn discharge(&self) -> &DischargeCoefficient {
        &self.discharge
    }
}

/// Complete directional flow at one horizontal opening; all rates are
/// nonnegative and retain both directions even if net transfer vanishes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CeilingFlow {
    /// Layer immediately below the opening, in the physically lower room.
    pub below_layer: Layer,
    /// Layer immediately above the opening, in the physically upper room.
    pub above_layer: Layer,
    /// Below-minus-above pressure at the aperture [Pa]. Zero when closed.
    pub pressure_difference: Pressure,
    /// Pressure magnitude suppressing density exchange [Pa]. Absent for stable
    /// density ordering, zero gravity, or a closed opening.
    pub flooding_pressure: Option<Pressure>,
    /// Equal volume exchange in each direction before donor density [m3/s].
    pub exchange_volume_rate: VolumetricFlowRate,
    /// Total upward volume flow, including forced flow and exchange [m3/s].
    pub upward_volume_rate: VolumetricFlowRate,
    /// Total downward volume flow, including forced flow and exchange [m3/s].
    pub downward_volume_rate: VolumetricFlowRate,
    /// Actual donor mass transported from below to above [kg/s].
    pub upward_mass_rate: MassFlowRate,
    /// Actual donor mass transported from above to below [kg/s].
    pub downward_mass_rate: MassFlowRate,
}
impl CeilingFlow {
    /// Paired layer source rates. `left` means the physically lower column and
    /// `right` the upper column; this is not a left/right geometric convention.
    pub fn mass_rates(self) -> LayerMassRates {
        let net = self.upward_mass_rate.value() - self.downward_mass_rate.value();
        let mut left = [MassFlowRate::new(0.0); 2];
        let mut right = left;
        left[self.below_layer.index()] = MassFlowRate::new(-net);
        right[self.above_layer.index()] = MassFlowRate::new(net);
        LayerMassRates { left, right }
    }

    /// Advect each direction's actual donor enthalpy in one declared common
    /// energy gauge. Pressure work is included by transporting enthalpy.
    pub fn advect_enthalpy(
        self,
        below: [SpecificEnthalpy; 2],
        above: [SpecificEnthalpy; 2],
    ) -> Result<LayerEnthalpyRates, VentError> {
        if !below.iter().chain(&above).all(|h| h.value().is_finite()) {
            return Err(VentError::Invalid(
                "finite ceiling donor enthalpies required",
            ));
        }
        let up = finite(
            self.upward_mass_rate.value() * below[self.below_layer.index()].value(),
            "upward ceiling enthalpy",
        )?;
        let down = finite(
            self.downward_mass_rate.value() * above[self.above_layer.index()].value(),
            "downward ceiling enthalpy",
        )?;
        let net = finite(up - down, "net ceiling enthalpy")?;
        let mut left = [Power::new(0.0); 2];
        let mut right = left;
        left[self.below_layer.index()] = Power::new(-net);
        right[self.above_layer.index()] = Power::new(net);
        Ok(LayerEnthalpyRates { left, right })
    }
}

/// Evaluate shallow horizontal flow between the supplied lower and upper spaces.
///
/// Hydrostatic pressures come from the original column owner. The opening must
/// be a ceiling for the lower finite room and a floor for the upper finite room;
/// topology admission belongs to the network, because a column carries no room
/// bounds. A homogeneous ambient may occupy either side. No allocation occurs.
pub fn ceiling_vent_flow(
    opening: &CeilingOpening,
    below: TwoLayerColumn,
    above: TwoLayerColumn,
    gravity: Acceleration,
    cancelled: impl FnMut() -> bool,
) -> Result<CeilingFlow, VentError> {
    // Traces approach the lower space from below and the upper from above.
    let below_layer = if opening.height.value() <= below.interface_height().value() {
        Layer::Lower
    } else {
        Layer::Upper
    };
    let above_layer = if opening.height.value() < above.interface_height().value() {
        Layer::Lower
    } else {
        Layer::Upper
    };
    flow_from_layers(
        opening,
        below,
        above,
        below_layer,
        above_layer,
        gravity,
        cancelled,
    )
}

/// The finite-room network owns the exact layer inventories. Its admitted
/// ceiling/floor topology can therefore select the present adjacent layer even
/// when a world-coordinate interface rounds onto the boundary. Both density
/// traces and transported inventories then use that same source-owned layer.
pub(super) fn flow_from_layers(
    opening: &CeilingOpening,
    below: TwoLayerColumn,
    above: TwoLayerColumn,
    below_layer: Layer,
    above_layer: Layer,
    gravity: Acceleration,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CeilingFlow, VentError> {
    poll(&mut cancelled)?;
    if !gravity.value().is_finite() || gravity.value() < 0.0 {
        return Err(VentError::Invalid(
            "finite nonnegative ceiling gravity required",
        ));
    }
    let mut result = CeilingFlow {
        below_layer,
        above_layer,
        pressure_difference: Pressure::new(0.0),
        flooding_pressure: None,
        exchange_volume_rate: VolumetricFlowRate::new(0.0),
        upward_volume_rate: VolumetricFlowRate::new(0.0),
        downward_volume_rate: VolumetricFlowRate::new(0.0),
        upward_mass_rate: MassFlowRate::new(0.0),
        downward_mass_rate: MassFlowRate::new(0.0),
    };
    if opening.area.value() == 0.0 {
        poll(&mut cancelled)?;
        return Ok(result);
    }
    let difference = below
        .pressure_difference_at(above, opening.height, gravity)?
        .value();
    result.pressure_difference = Pressure::new(difference);
    let rho_below = below.density(below_layer).value();
    let rho_above = above.density(above_layer).value();
    let density_difference = rho_above - rho_below;
    let mut exchange = 0.0;
    if density_difference > 0.0 && gravity.value() > 0.0 {
        let buoyant_head = positive(
            gravity.value() * density_difference,
            "ceiling buoyant density head",
        )?;
        let flood = positive(
            buoyant_head * opening.flooding_length,
            "ceiling flooding pressure",
        )?;
        result.flooding_pressure = Some(Pressure::new(flood));
        if difference.abs() < flood {
            let average_density = 0.5 * rho_below + 0.5 * rho_above;
            let velocity_squared = positive(
                buoyant_head * opening.root_area / average_density,
                "ceiling exchange velocity squared",
            )?;
            let zero_pressure_exchange = positive(
                0.1 * opening.area.value() * det::sqrt(velocity_squared),
                "ceiling zero-pressure exchange",
            )?;
            exchange = positive(
                zero_pressure_exchange * (1.0 - difference.abs() / flood),
                "ceiling pressure-limited exchange",
            )?;
        }
    }
    poll(&mut cancelled)?;
    let mut up = exchange;
    let mut down = exchange;
    if difference != 0.0 {
        let donor_density = if difference > 0.0 {
            rho_below
        } else {
            rho_above
        };
        let inverse_root_resistance =
            1.0 / det::sqrt(opening.discharge.reference_loss().resistance.value());
        let forced = positive(
            opening.area.value()
                * inverse_root_resistance
                * (det::sqrt(difference.abs()) / det::sqrt(donor_density)),
            "ceiling pressure-driven volume flow",
        )?;
        if difference > 0.0 {
            up = finite(up + forced, "upward ceiling volume flow")?;
        } else {
            down = finite(down + forced, "downward ceiling volume flow")?;
        }
    }
    let upward_mass = finite(up * rho_below, "upward ceiling mass flow")?;
    let downward_mass = finite(down * rho_above, "downward ceiling mass flow")?;
    if (up > 0.0 && upward_mass == 0.0) || (down > 0.0 && downward_mass == 0.0) {
        return Err(VentError::NonFinite(
            "nonzero ceiling mass flow underflowed",
        ));
    }
    result.exchange_volume_rate = VolumetricFlowRate::new(exchange);
    result.upward_volume_rate = VolumetricFlowRate::new(up);
    result.downward_volume_rate = VolumetricFlowRate::new(down);
    result.upward_mass_rate = MassFlowRate::new(upward_mass);
    result.downward_mass_rate = MassFlowRate::new(downward_mass);
    poll(&mut cancelled)?;
    Ok(result)
}

fn finite(value: f64, stage: &'static str) -> Result<f64, VentError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(VentError::NonFinite(stage))
    }
}
fn positive(value: f64, stage: &'static str) -> Result<f64, VentError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(VentError::NonFinite(stage))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), VentError> {
    if cancelled() {
        Err(VentError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_qty::Density;

    fn column(p: f64, height: f64, interface: f64, densities: [f64; 2]) -> TwoLayerColumn {
        TwoLayerColumn::new(
            Length::new(height),
            Pressure::new(p),
            Length::new(interface),
            Density::new(densities[0]),
            Density::new(densities[1]),
        )
        .unwrap()
    }
    fn aperture(area: f64, shape: CeilingShape) -> CeilingOpening {
        CeilingOpening::cooper(Length::new(0.0), Area::new(area), shape, 0.0).unwrap()
    }
    fn flow(area: f64, shape: CeilingShape, dp: f64, density: [f64; 2], g: f64) -> CeilingFlow {
        ceiling_vent_flow(
            &aperture(area, shape),
            column(100_000.0 + dp, 0.0, -1.0, [density[0]; 2]),
            column(100_000.0, 0.0, 1.0, [density[1]; 2]),
            Acceleration::new(g),
            || false,
        )
        .unwrap()
    }
    fn close(actual: f64, expected: f64, relative: f64) {
        assert!(
            (actual - expected).abs() <= relative * expected.abs().max(1e-30),
            "actual={actual:.17e}; expected={expected:.17e}"
        );
    }

    #[test]
    fn unstable_equal_pressure_exchanges_equal_volumes_with_actual_donor_masses() {
        let result = flow(1.0, CeilingShape::Square, 0.0, [0.6, 1.2], 10.0);
        let q = 0.258_198_889_747_161_1;
        close(result.exchange_volume_rate.value(), q, 2e-15);
        assert_eq!(result.upward_volume_rate, result.downward_volume_rate);
        close(result.upward_mass_rate.value(), 0.6 * q, 2e-15);
        close(result.downward_mass_rate.value(), 1.2 * q, 2e-15);
        close(result.flooding_pressure.unwrap().value(), 2.662_092, 2e-15);
        assert_eq!(result.below_layer, Layer::Upper);
        assert_eq!(result.above_layer, Layer::Lower);
        let m = result.mass_rates();
        assert_eq!(m.left[0], MassFlowRate::new(0.0));
        assert_eq!(m.right[1], MassFlowRate::new(0.0));
        assert_eq!(m.left[1].value(), -m.right[0].value());
    }

    #[test]
    fn pressure_flow_reuses_original_coefficient_and_actual_upstream_density() {
        let opening = aperture(0.7, CeilingShape::Round);
        assert_eq!(opening.discharge().value().value(), 0.68);
        assert_eq!(
            opening.discharge().reference_loss().tolerance_basis,
            ToleranceBasis::EngineeringAllowance
        );
        for dp in [-100.0_f64, 100.0] {
            let result = flow(0.7, CeilingShape::Round, dp, [1.2, 0.6], 10.0);
            let (density, mass, other) = if dp > 0.0 {
                (1.2, result.upward_mass_rate, result.downward_mass_rate)
            } else {
                (0.6, result.downward_mass_rate, result.upward_mass_rate)
            };
            close(
                mass.value(),
                0.68 * 0.7 * (2.0 * density * dp.abs()).sqrt(),
                3e-15,
            );
            assert_eq!(other, MassFlowRate::new(0.0));
            assert_eq!(result.flooding_pressure, None);
        }
    }

    #[test]
    fn stable_order_equal_density_zero_gravity_and_closed_hole_have_no_exchange() {
        for density in [[1.2, 0.6], [1.2, 1.2]] {
            let result = flow(1.0, CeilingShape::Round, 0.0, density, 9.81);
            assert_eq!(result.upward_mass_rate.value(), 0.0);
            assert_eq!(result.downward_mass_rate.value(), 0.0);
            assert_eq!(result.flooding_pressure, None);
        }
        for (area, g) in [(0.0, 9.81), (1.0, 0.0)] {
            let result = flow(area, CeilingShape::Square, 0.0, [0.6, 1.2], g);
            assert_eq!(result.exchange_volume_rate.value(), 0.0);
            assert_eq!(result.flooding_pressure, None);
        }
        let closed = flow(0.0, CeilingShape::Square, 100.0, [0.6, 1.2], 9.81);
        assert_eq!(closed.upward_mass_rate.value(), 0.0);
        assert_eq!(closed.downward_mass_rate.value(), 0.0);
    }

    #[test]
    fn counterflow_is_continuously_suppressed_by_pressure_in_both_directions() {
        let zero = flow(1.0, CeilingShape::Square, 0.0, [0.6, 1.2], 10.0);
        let flood = zero.flooding_pressure.unwrap().value();
        for sign in [-1.0, 1.0] {
            let half = flow(
                1.0,
                CeilingShape::Square,
                sign * 0.5 * flood,
                [0.6, 1.2],
                10.0,
            );
            close(
                half.exchange_volume_rate.value(),
                0.5 * zero.exchange_volume_rate.value(),
                1e-11,
            );
            assert!(half.upward_mass_rate.value() > 0.0 && half.downward_mass_rate.value() > 0.0);
            let beyond = flow(
                1.0,
                CeilingShape::Square,
                sign * 1.001 * flood,
                [0.6, 1.2],
                10.0,
            );
            assert_eq!(beyond.exchange_volume_rate.value(), 0.0);
            if sign > 0.0 {
                assert_eq!(beyond.downward_mass_rate.value(), 0.0);
            } else {
                assert_eq!(beyond.upward_mass_rate.value(), 0.0);
            }
        }
    }

    #[test]
    fn circle_and_square_use_the_original_shape_specific_lengths() {
        let square = flow(1.0, CeilingShape::Square, 0.0, [0.6, 1.2], 10.0);
        let circle = flow(1.0, CeilingShape::Round, 0.0, [0.6, 1.2], 10.0);
        let diameter = 2.0 / core::f64::consts::PI.sqrt();
        close(
            circle.flooding_pressure.unwrap().value(),
            0.5 * 0.754_f64.powi(2) * 10.0 * 0.6 * diameter.powi(5),
            3e-15,
        );
        assert_eq!(circle.exchange_volume_rate, square.exchange_volume_rate);
        assert_ne!(circle.flooding_pressure, square.flooding_pressure);
    }

    #[test]
    fn area_scaling_preserves_the_five_halves_power_without_large_intermediates() {
        let small = flow(1.0, CeilingShape::Square, 0.0, [0.6, 1.2], 10.0);
        let large = flow(4.0, CeilingShape::Square, 0.0, [0.6, 1.2], 10.0);
        close(
            large.exchange_volume_rate.value() / small.exchange_volume_rate.value(),
            4.0 * 2.0_f64.sqrt(),
            3e-15,
        );
        close(
            large.flooding_pressure.unwrap().value() / small.flooding_pressure.unwrap().value(),
            2.0,
            2e-15,
        );
    }

    #[test]
    fn traces_at_an_interface_select_existing_layers_and_ignore_absent_properties() {
        let result = ceiling_vent_flow(
            &aperture(1.0, CeilingShape::Square),
            column(100_000.0, 0.0, 0.0, [0.6, 50.0]),
            column(100_000.0, 0.0, 0.0, [40.0, 1.2]),
            Acceleration::new(10.0),
            || false,
        )
        .unwrap();
        assert_eq!(result.below_layer, Layer::Lower);
        assert_eq!(result.above_layer, Layer::Upper);
        close(
            result.upward_mass_rate.value(),
            0.6 * 0.258_198_889_747_161_1,
            2e-15,
        );
        close(
            result.downward_mass_rate.value(),
            1.2 * 0.258_198_889_747_161_1,
            2e-15,
        );
    }

    #[test]
    fn hydrostatic_reference_and_world_translation_leave_flow_invariant() {
        let vent = aperture(1.0, CeilingShape::Square);
        let base = ceiling_vent_flow(
            &vent,
            column(100_000.0, 0.0, -1.0, [0.6; 2]),
            column(100_000.0, 0.0, 1.0, [1.2; 2]),
            Acceleration::new(10.0),
            || false,
        )
        .unwrap();
        let moved = CeilingOpening::cooper(
            Length::new(200.0),
            Area::new(1.0),
            CeilingShape::Square,
            0.0,
        )
        .unwrap();
        let changed = ceiling_vent_flow(
            &moved,
            column(100_012.0, 198.0, 199.0, [0.6; 2]),
            column(99_964.0, 203.0, 201.0, [1.2; 2]),
            Acceleration::new(10.0),
            || false,
        )
        .unwrap();
        assert_eq!(base, changed);
    }

    #[test]
    fn enthalpy_uses_both_donors_and_is_exactly_equal_and_opposite() {
        let result = flow(1.0, CeilingShape::Square, 0.0, [0.6, 1.2], 10.0);
        let h = result
            .advect_enthalpy(
                [SpecificEnthalpy::new(700_000.0); 2],
                [SpecificEnthalpy::new(200_000.0); 2],
            )
            .unwrap();
        let net = result.upward_mass_rate.value() * 700_000.0
            - result.downward_mass_rate.value() * 200_000.0;
        assert_eq!(h.left[1].value(), -net);
        assert_eq!(h.right[0].value(), net);
        assert_eq!(h.left[0].value(), 0.0);
        assert_eq!(h.right[1].value(), 0.0);
        assert!(net > 0.0);
        assert!(result.upward_mass_rate.value() < result.downward_mass_rate.value());
    }

    #[test]
    fn invalid_geometry_and_nonfinite_derived_rates_refuse() {
        for area in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(
                CeilingOpening::cooper(Length::new(0.0), Area::new(area), CeilingShape::Round, 0.0)
                    .is_err()
            );
        }
        assert!(
            CeilingOpening::cooper(
                Length::new(f64::NAN),
                Area::new(1.0),
                CeilingShape::Round,
                0.0
            )
            .is_err()
        );
        assert!(
            CeilingOpening::cooper(Length::new(0.0), Area::new(1.0), CeilingShape::Round, 1.0)
                .is_err()
        );
        let overflow = ceiling_vent_flow(
            &aperture(f64::MAX, CeilingShape::Round),
            column(100_001.0, 0.0, -1.0, [0.6; 2]),
            column(100_000.0, 0.0, 1.0, [1.2; 2]),
            Acceleration::new(10.0),
            || false,
        );
        assert!(matches!(overflow, Err(VentError::NonFinite(_))));
    }

    #[test]
    fn cancellation_refuses_before_any_result_is_published() {
        for stop in 1..=3 {
            let mut calls = 0;
            let result = ceiling_vent_flow(
                &aperture(1.0, CeilingShape::Square),
                column(100_000.0, 0.0, -1.0, [0.6; 2]),
                column(100_000.0, 0.0, 1.0, [1.2; 2]),
                Acceleration::new(10.0),
                || {
                    calls += 1;
                    calls == stop
                },
            );
            assert_eq!(result, Err(VentError::Cancelled));
        }
    }
}
