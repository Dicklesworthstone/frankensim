//! The current-surface fixtures as declared graphs, in the gate lane of
//! `gate/graph.rs`.
//!
//! Each switch's channel follows its package card's complete I(Vds, Vgs)
//! surface times the card's thermal multiplier at the package's temperature
//! ([`CardSurface`]). Both are resolved and evaluated as the module's
//! `GateCurrentSurfaceMaterial` does, through public fs-material and
//! fs-matdb calls. The rest of the lane is the gate fixture's.

use std::sync::Arc;

use fs_matdb::PreparedSurfaceResponse;
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_phs::PhsError;
use fs_phs::circuit::converter::commutation::gate::{
    GateChannelLaw, GateChannelPoint, GateChannelResponse,
};
use fs_system::commutated_drive::{ChannelLaw, GateStorage};

use super::super::super::graph::{Card, Case, compare, held_to_owner, kernel_run, leak, text};
use super::super::graph::{CELLS, gate_case, gate_owner, module_steps};
use super::*;

/// A switch's channel from its package's card, as the module resolves it:
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
/// and -4 V to 4 V gate, evaluated as the module's response evaluates it.
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

/// `run` on the module's thermal network for the fixture, its channels the
/// module's current surfaces.
fn module_network<R>(
    f: &CurrentFixture,
    run: impl FnOnce(&GateBatteryMotorThermalNetwork<'_, GateCurrentSurfaceBinding<'_>>) -> R,
) -> R {
    let context = gate_context();
    let cells = f.gate.cells(&context);
    let (r, d) = (drive_resistors(&f.gate.base), f.gate.base.devices());
    let (q, cgs, cgd) = (f.gate.charges(), f.gate.caps(false), f.gate.caps(true));
    let (current, driver) = (f.bindings(), f.gate.driver_resistors());
    let bindings = material_bindings(&q, &cgs, &cgd, &current, &driver);
    run(&current_network(&cells, &r, &d, bindings))
}

fn bits(point: GateChannelPoint) -> [u64; 3] {
    let values = [
        point.current_a,
        point.drain_derivative_s,
        point.gate_derivative_s,
    ];
    values.map(f64::to_bits)
}

#[test]
fn a_current_graph_channel_answers_as_the_module_surface_does() {
    // The module's sources bound at 300 K against the graph's laws there:
    // a grid over the support, bit for bit, and the support's edges.
    let f = leak(current_fixture(0.05, 0.2, Defect::None));
    let p = authored(0.05, false);
    let (case, laws) = surface_case(f, &p, 300.0);
    let bound: CurrentDrive = module_network(f, |network| {
        let thermal = f.gate.thermal();
        let native = &case.steady.native;
        network
            .bind(&thermal, native, &p, case.start, || false)
            .unwrap()
    });
    let module = bound.gate_material().channel_conductance();
    for (law, module) in laws.iter().zip(module) {
        let (graph, module) = (law.channel(300.0).unwrap(), module.native_response());
        for vd in [-20.0, -2.0, -1.0, -0.75, 0.0, 0.75, 1.0, 2.0, 20.0] {
            for vg in [-4.0, -1.0, 0.0, 0.75, 1.0, 2.0, 4.0] {
                let [a, b] = [&graph, &module].map(|law| law.evaluate(vd, vg).unwrap());
                assert_eq!(bits(a), bits(b), "{vd} V, {vg} V");
            }
        }
        for (vd, vg) in [(20.01, 1.0), (-20.01, 1.0), (1.0, 4.01)] {
            assert!(graph.evaluate(vd, vg).is_err());
        }
    }
}

#[test]
fn a_current_graph_scales_its_channel_power_at_the_accepted_temperature() {
    // The module's fixture: Cgd 0.05 F, the multiplier rising 0.2 per K, two
    // steps of 1/1024 s from the module's current state.
    let f = leak(current_fixture(0.05, 0.2, Defect::None));
    let p = authored(0.05, false);
    let (case, _) = surface_case(f, &p, 300.0);
    let graph = kernel_run(&case, 2);
    let module = module_network(f, |network| module_steps(network, &case, &p, 2));
    compare(&module, &graph, [1e-14, 1e-12]);
    // Without the module: the owner by hand with the surface in closed form.
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

/// The module's drive bound at one accepted state refuses a step from
/// another (`current_surface_drive_bound_at_one_state_refuses_a_step_from_another_accepted_state`).
/// The graph binds nothing: from that other state it takes the step at the
/// laws' answers there.
#[test]
fn a_current_graph_steps_from_another_accepted_state_at_the_answers_there() {
    // The refusal's fixture from 320 K: only the channels' multiplier
    // varies with temperature.
    let f = leak(current_fixture(0.05, 0.2, Defect::None));
    let p = authored(0.05, false);
    let (case, laws) = surface_case(f, &p, 320.0);
    let graph = kernel_run(&case, 1);
    // The module bound at that state takes the same step.
    let module = module_network(f, |network| module_steps(network, &case, &p, 1));
    compare(&module, &graph, [1e-14, 1e-12]);
    // Without the module: the owner by hand at the laws' answers at the
    // graph's end temperatures, its new state bit for bit.
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
