//! The current-surface fixtures as declared graphs, in the gate lane of
//! `gate/graph.rs`.
//!
//! Each switch's channel follows its package card's complete I(Vds, Vgs)
//! surface times the card's thermal multiplier at the package's temperature
//! ([`CardSurface`]). Both are resolved and evaluated as the retired module's
//! `GateCurrentSurfaceMaterial` did, through public fs-material and fs-matdb
//! calls. The rest of the lane is the gate fixture's.

use std::sync::Arc;

use fs_matdb::PreparedSurfaceResponse;
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_phs::PhsError;
use fs_phs::circuit::converter::commutation::gate::{
    GateChannelLaw, GateChannelPoint, GateChannelResponse,
};
use fs_qty::{Current, Dims};
use fs_system::commutated_drive::{ChannelLaw, GateStorage};

use super::super::graph::{Card, Case, Row, held_to_owner, kernel_rows, leak, text};
use super::*;
use super::{CELLS, gate_case, gate_owner};

struct CurrentFixture {
    gate: GateFixture,
    current: [PropertyKey; 4],
    scale: PropertyKey,
}
fn current_fixture(miller: f64, slope: f64) -> CurrentFixture {
    let mut gate = gate_fixture(miller, false, 0.0);
    let current = core::array::from_fn(|i| {
        PropertyKey::new(format!("channel_current_surface_{i}"), Current::DIMS)
    });
    let scale = PropertyKey::new("channel_temperature_multiplier", Dims::NONE);
    let mut claims = ClaimSet::new();
    for (_, claim) in gate.base.card.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    let thermal = || {
        ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            299.0,
            350.0,
        )
    };
    let mut add = |key: &PropertyKey, value, validity| {
        claims.insert_claim(PropertyClaim{key:key.clone(),value,validity,
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,observations:vec![],
            provenance:Provenance{source:"manufactured original triangular channel-only current; not qualified transistor data".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
    };
    for (device, key) in current.iter().enumerate() {
        let drain: [f64; 5] = [-100.0, -1.0, 0.0, 1.0, 100.0];
        let gate_v: [f64; 4] = [-4.0, 0.0, 1.0, 4.0];
        let reversed = device % 2 == 1;
        let mut samples = Vec::new();
        for &vg in &gate_v {
            for &vd in &drain {
                let amps = vd.signum()
                    * vd.abs().min(1.0)
                    * vg.max(0.0)
                    * if vd < 0.0 { 0.4 } else { 1.0 };
                samples.push((if reversed { [vg, vd] } else { [vd, vg] }, amps));
            }
        }
        let mut triangles = Vec::new();
        for row in 0..gate_v.len() - 1 {
            for col in 0..drain.len() - 1 {
                let a = (row * drain.len() + col) as u32;
                let b = a + 1;
                let d = a + drain.len() as u32;
                let e = d + 1;
                for mut triangle in [[a, b, d], [b, e, d]] {
                    if reversed {
                        triangle.swap(1, 2);
                    }
                    triangles.push(triangle);
                }
            }
        }
        add(
            key,
            PropertyValue::Surface2D {
                axes: if reversed {
                    ["Vgs".into(), "Vds".into()]
                } else {
                    ["Vds".into(), "Vgs".into()]
                },
                axis_dims: [Voltage::DIMS; 2],
                samples,
                triangles,
                dims: Current::DIMS,
            },
            thermal()
                .with_quantity(
                    "Vds",
                    QuantitySpec::dimensional(Voltage::DIMS),
                    -100.0,
                    100.0,
                )
                .with_quantity("Vgs", QuantitySpec::dimensional(Voltage::DIMS), -4.0, 4.0),
        );
    }
    add(
        &scale,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(299.0, 1.0 - slope), (350.0, 1.0 + 50.0 * slope)],
            dims: Dims::NONE,
        },
        thermal(),
    );
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-device".into(),
            phase: "solid".into(),
            process: "original signed triangular current and explicit thermal multiplier".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &common::point(299.0),
            upper: &common::point(350.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    gate.base.card = card;
    gate.base.caloric = caloric;
    CurrentFixture {
        gate,
        current,
        scale,
    }
}
fn current_initial(steady: &SteadyDrive, p: &GateCommutationParameters) -> [f64; 15] {
    let mut state = initial(steady, p);
    state[3..12].copy_from_slice(
        &p.state(0.3, [12.0, 6.0], [10.0, 7.5], [2.0, 0.2, 2.0, 0.2])
            .unwrap(),
    );
    state
}
fn command(steady: &SteadyDrive) -> BatteryMotorSwitchCommand {
    BatteryMotorSwitchCommand {
        auxiliary_load: None,
        ..steady.command
    }
}
fn approx(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-10, "{a:.17e} != {b:.17e}");
}

/// A switch's channel from its package's card, as the module resolved it:
/// the complete surface of `current`, under the claim its reference point
/// (0 V drain, -4 V gate) selects, times the multiplier `scale`.
#[derive(Debug, Clone, Copy)]
struct CardSurface {
    card: Card,
    current: &'static PropertyKey,
    scale: &'static PropertyKey,
    tag: u64,
}

impl CardSurface {
    /// The surface and its multiplier at `kelvin`.
    fn at(&self, kelvin: f64) -> Result<Scaled, String> {
        let query = self.card.query(kelvin)?;
        let reference = query.with("Vds", 0.0).and_then(|q| q.with("Vgs", -4.0));
        let rule = ScalarAdmissibility::Finite;
        let requirement = ScalarPropertyRequirement::try_with_key(self.current, rule);
        let snapshot = resolve_material_state_point(
            self.card.material,
            &reference.map_err(text)?,
            &[requirement.map_err(text)?],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(text)?;
        let property = snapshot.property_by_key(self.current).ok_or("no channel")?;
        let claims = self.card.material.claims();
        let claim = claims.claim(property.answer().receipt.selected);
        let surface = claim.ok_or("no claim")?.value.prepare_surface_response();
        let surface = surface.map_err(text)?;
        let drain = surface.axes().iter().position(|axis| axis == "Vds");
        let rule = ScalarAdmissibility::StrictlyPositive;
        Ok(Scaled {
            drain: drain.ok_or("no drain axis")?,
            surface: Arc::new(surface),
            factor: self.card.values(kelvin, &[(self.scale, rule)])?[0],
        })
    }
}

impl ChannelLaw for CardSurface {
    fn channel(&self, kelvin: f64) -> Result<GateChannelLaw, String> {
        GateChannelLaw::new(self.at(kelvin)?).map_err(text)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.tag);
    }
}

/// A surface times `factor` on the module's support, -20 V to 20 V drain
/// and -4 V to 4 V gate, evaluated as the module's response evaluated it.
#[derive(Debug, Clone)]
struct Scaled {
    surface: Arc<PreparedSurfaceResponse>,
    drain: usize,
    factor: f64,
}

impl GateChannelResponse for Scaled {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        [-4.0, 4.0]
    }

    fn evaluate(&self, vds: f64, vgs: f64) -> Result<GateChannelPoint, PhsError> {
        let what = "channel voltages outside the surface's support";
        if !((-20.0..=20.0).contains(&vds) && (-4.0..=4.0).contains(&vgs)) {
            return Err(PhsError::Unadmitted { what });
        }
        let mut at = [vgs; 2];
        at[self.drain] = vds;
        let value = self.surface.evaluate(at);
        let value = value.map_err(|_| PhsError::Unadmitted { what })?;
        Ok(GateChannelPoint {
            current_a: value.value * self.factor,
            drain_derivative_s: value.gradient[self.drain] * self.factor,
            gate_derivative_s: value.gradient[1 - self.drain] * self.factor,
        })
    }
}

/// The fixture's surface in closed form, times the multiplier
/// `1 + slope (T - 300 K)`: I = sign(Vds) min(|Vds|, 1 V) max(Vgs, 0), and
/// 0.4 of that for negative Vds. Where the module's steps run (|Vds| above
/// 1 V, Vgs positive), each of the card's triangles carries this affine law.
#[derive(Debug, Clone)]
struct ClosedSurface(f64);

impl GateChannelResponse for ClosedSurface {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        [-4.0, 4.0]
    }

    fn evaluate(&self, vds: f64, vgs: f64) -> Result<GateChannelPoint, PhsError> {
        let k = if vds < 0.0 { 0.4 } else { 1.0 } * self.0;
        let (d, g) = (vds.signum() * vds.abs().min(1.0), vgs.max(0.0));
        Ok(GateChannelPoint {
            current_a: k * d * g,
            drain_derivative_s: if vds.abs() < 1.0 { k * g } else { 0.0 },
            gate_derivative_s: if vgs > 0.0 { k * d } else { 0.0 },
        })
    }
}

/// The module's authored parameters for Cgd `miller`, each channel
/// `channel(i, T)` at its package's temperature.
fn parameters_at<C: GateChannelResponse + Clone + 'static>(
    miller: f64,
    channel: impl Fn(usize, f64) -> C,
) -> impl Fn(&[f64]) -> GateCommutationParameters<C> {
    move |kelvin: &[f64]| {
        let a = authored(miller, false);
        GateCommutationParameters {
            converter: a.converter,
            intrinsic_drain_charge_curves: a.intrinsic_drain_charge_curves,
            gate_source_capacitance_f: a.gate_source_capacitance_f,
            gate_drain_capacitance_f: a.gate_drain_capacitance_f,
            channel_conductance: core::array::from_fn(|i| channel(i, kelvin[i + 1])),
            gate_drivers: a.gate_drivers,
            maximum_cells_per_chord: a.maximum_cells_per_chord,
        }
    }
}

/// A current-surface case of the module's fixture `f` (Cgd 0.05 F, the
/// multiplier rising 0.2 per K), authored as `p`: the gate case with each
/// channel following [`CardSurface`], started from the module's current
/// state, every cell at `start_k`. One step of 1/1024 s under the steady
/// command without its auxiliary load. Also returns the four laws.
fn surface_case(
    f: &'static CurrentFixture,
    p: &GateCommutationParameters,
    start_k: f64,
) -> (Case<'static, GateStorage>, [CardSurface; 4]) {
    let steady = steady_drive(false);
    let command = command(&steady);
    let steady = leak(SteadyDrive { command, ..steady });
    let mut case = gate_case(&f.gate, steady, p, start_k, 1.0 / 1024.0);
    let card = Card {
        material: &f.gate.base.card,
        curve: &f.gate.base.caloric,
        context: case.context,
    };
    let laws = core::array::from_fn(|i| CardSurface {
        card,
        current: &f.current[i],
        scale: &f.scale,
        tag: i as u64,
    });
    case.storage.channels = laws.map(|law| Arc::new(law) as Arc<dyn ChannelLaw>);
    case.start = current_initial(steady, p);
    (case, laws)
}

/// The module test
/// `original_triangles_retain_signed_saturation_cutoff_zero_mesh_and_named_partials`
/// on the graph's copy of the module's surface law. Until the module
/// retired, the copy was also held bit for bit to the module's bound
/// sources, at 63 points per channel.
#[test]
fn a_current_graph_channel_keeps_the_signed_saturation_cutoff_zero_row_and_partials() {
    let f = leak(current_fixture(0.05, 0.2));
    let p = authored(0.05, false);
    let (_, laws) = surface_case(f, &p, 300.0);
    for (i, law) in laws.iter().enumerate() {
        // The claim's own surface, drain first on the even packages, and the
        // multiplier 1 at 300 K.
        let scaled = law.at(300.0).unwrap();
        assert_eq!(scaled.drain, i % 2);
        near(scaled.factor, 1.0);
        let mut claims = law.card.material.claims().claims_ordered();
        let mut current = claims
            .by_ref()
            .filter(|(_, claim)| claim.key == f.current[i]);
        let (_, claim) = current.next().unwrap();
        assert!(current.next().is_none(), "one current claim per package");
        let PropertyValue::Surface2D {
            samples, triangles, ..
        } = &claim.value
        else {
            panic!("original surface lost")
        };
        assert_eq!(scaled.surface.samples(), samples.as_slice());
        assert_eq!(scaled.surface.triangles(), triangles.as_slice());
        // Saturation, signed reverse current, the linear cell and cutoff,
        // each with its partials and their finite differences.
        let channel = law.channel(300.0).unwrap();
        for (vd, vg, expected) in [
            (2.0, 2.0, [2.0, 0.0, 1.0]),
            (-2.0, 2.0, [-0.8, 0.0, -0.4]),
            (0.75, 0.75, [0.5, 1.0, 1.0]),
            (2.0, -1.0, [0.0, 0.0, 0.0]),
        ] {
            let point = channel.evaluate(vd, vg).unwrap();
            let values = [
                point.current_a,
                point.drain_derivative_s,
                point.gate_derivative_s,
            ];
            for (a, b) in values.into_iter().zip(expected) {
                approx(a, b);
            }
            let e = 1e-5;
            let amps = |vd, vg| channel.evaluate(vd, vg).unwrap().current_a;
            approx((amps(vd + e, vg) - amps(vd - e, vg)) / (2.0 * e), values[1]);
            approx((amps(vd, vg + e) - amps(vd, vg - e)) / (2.0 * e), values[2]);
        }
        // No current along the 0 V drain row, and none outside the support.
        for vg in [-4.0, -0.2, 0.0, 0.7, 2.0, 4.0] {
            assert_eq!(channel.evaluate(0.0, vg).unwrap().current_a, 0.0);
        }
        for (vd, vg) in [(20.01, 1.0), (-20.01, 1.0), (1.0, 4.01)] {
            assert!(channel.evaluate(vd, vg).is_err());
        }
    }
}

#[test]
fn a_current_graph_scales_its_channel_power_at_the_accepted_temperature() {
    // The module's fixture: Cgd 0.05 F, the multiplier rising 0.2 per K, two
    // steps of 1/1024 s from the module's current state.
    let f = leak(current_fixture(0.05, 0.2));
    let p = authored(0.05, false);
    let (case, _) = surface_case(f, &p, 300.0);
    let graph: Vec<Row<_, [f64; CELLS]>> = kernel_rows(&case, 2);
    // The owner by hand with the surface in closed form.
    let closed = parameters_at(0.05, |_, t| ClosedSurface(1.0 + 0.2 * (t - 300.0)));
    held_to_owner(&case, &graph, gate_owner(&case, closed), [1e-14, 2e-13]);
    // The module's closed forms for the first step: each package warms by
    // its channel's and diode's heat at the midpoint voltages, each driver
    // cell by its resistor's.
    let h = case.dt_s;
    let (d0, g0) = voltages(&case.start, 0.05);
    let (d1, g1) = voltages(&graph[0].state, 0.05);
    for (i, driver) in drivers().into_iter().enumerate() {
        let t = graph[0].kelvin[i + 1];
        let (vd, vg) = (f64::midpoint(d0[i], d1[i]), f64::midpoint(g0[i], g1[i]));
        assert!(vd.abs() > 1.0 && vg > 0.0);
        let amps = vd.signum() * if vd < 0.0 { 0.4 } else { 1.0 } * vg * (1.0 + 0.2 * (t - 300.0));
        near(
            t - 300.0,
            h * vd * amps + h * (-vd) * ((-vd - 0.5) / 0.2).max(0.0),
        );
        let (u, ohms) = if case.gates(0).on[i] {
            (driver.on_voltage_v, driver.on_resistance_ohm)
        } else {
            (driver.off_voltage_v, driver.off_resistance_ohm)
        };
        let amps = (u - vg) / ohms;
        near(
            graph[0].kelvin[IDS.len() + i] - 300.0,
            h * ohms * amps * amps,
        );
    }
}

/// The module's drive, bound at one accepted state, refused a step from
/// another (`current_surface_drive_bound_at_one_state_refuses_a_step_from_another_accepted_state`).
/// The graph binds nothing: from that other state it takes the step at the
/// laws' answers there.
#[test]
fn a_current_graph_steps_from_another_accepted_state_at_the_answers_there() {
    // The refusal's fixture from 320 K: only the channels' multiplier
    // varies with temperature.
    let f = leak(current_fixture(0.05, 0.2));
    let p = authored(0.05, false);
    let (case, laws) = surface_case(f, &p, 320.0);
    let graph: Vec<Row<_, [f64; CELLS]>> = kernel_rows(&case, 1);
    // The owner by hand at the laws' answers at the graph's end
    // temperatures, its new state bit for bit.
    let owner = gate_owner(&case, parameters_at(0.05, |i, t| laws[i].at(t).unwrap()));
    let (next, _) = owner(0, &case.start, &graph[0].kelvin);
    assert_eq!(next.map(f64::to_bits), graph[0].state.map(f64::to_bits));
    held_to_owner(&case, &graph, owner, [0.0, 2e-13]);
    // The multipliers at the step's end differ in their bits from those at
    // 300 K, the state the module's drive was bound at.
    let end: [f64; CELLS] = graph[0].kelvin;
    assert!(end.iter().all(|&kelvin| kelvin >= 320.0), "{end:?}");
    for (law, kelvin) in laws.iter().zip(&end[1..]) {
        let factor = |kelvin| law.at(kelvin).unwrap().factor.to_bits();
        assert_ne!(factor(*kelvin), factor(300.0));
    }
}
