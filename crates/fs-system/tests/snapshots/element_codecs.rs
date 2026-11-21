//! Native interface-law state, parameter binding and friction-history replay.

use super::*;
use fs_ports::component::LoadLine;
use fs_ports::{CoordinateBinding, JunctionLaw, JunctionSchema, PortOrientation, PortTimestamp};
use fs_system::elements::{
    AffineSignal, ControlledSource, DiracJunction, HeatNode, ModulatedTransformer, Resistance,
    Slider, SliderEnd, Slip,
};

fn junction(coupling: f64, frame: &str) -> DiracJunction {
    let id = |text: &str| fs_ports::StableId::new(text).unwrap();
    let ports = ["a", "b"]
        .map(|name| {
            PortKind::ElectricalVoltageCurrent
                .scalar_seed_schema(
                    id(name),
                    CoordinateBinding::new(
                        id("basis"),
                        id(frame),
                        PortOrientation::OutwardFromOwner,
                    ),
                    PortTimestamp::new(id("clock"), 0),
                )
                .unwrap()
        })
        .to_vec();
    DiracJunction::new(
        JunctionSchema::try_new(
            id("junction"),
            ports,
            JunctionLaw::Dirac {
                rows: 2,
                incidence: vec![0.0, -coupling, coupling, 0.0],
            },
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn native_interface_codecs_bind_laws_domains_and_declared_coordinates() {
    let mechanical = PortKind::MechanicalForceVelocity;
    let electrical = PortKind::ElectricalVoltageCurrent;
    parameters_refuse(
        HeatNode::new(&[2.0, 3.0]).unwrap(),
        HeatNode::new(&[2.0, 3.0_f64.next_up()]).unwrap(),
        (),
    );
    parameters_refuse(
        Resistance::new(mechanical, vec![2.0, 0.5, 0.5, 3.0]).unwrap(),
        Resistance::new(mechanical, vec![2.0, 0.5, 0.5, 3.0_f64.next_up()]).unwrap(),
        (),
    );
    parameters_refuse(
        junction(2.0, "frame"),
        junction(2.0_f64.next_up(), "frame"),
        (),
    );
    parameters_refuse(junction(2.0, "frame"), junction(2.0, "another-frame"), ());
    parameters_refuse(
        ModulatedTransformer::new(mechanical, electrical),
        ModulatedTransformer::new(electrical, electrical),
        (),
    );
    let mut account = SourceAccount::starting_at(1.0e16);
    account.clock_s.add(1.0);
    account.energy_in_j = -0.0;
    parameters_refuse(
        ControlledSource::effort(mechanical),
        ControlledSource::flow(mechanical),
        account,
    );
    parameters_refuse(
        Slider::new(10.0, 0.5, 0.3).unwrap(),
        Slider::new(10.0, 0.5, 0.3_f64.next_up()).unwrap(),
        Slip::sliding(-0.0),
    );
    parameters_refuse(
        Slider::new(10.0, 0.5, 0.3).unwrap(),
        Slider::new(10.0, 0.5_f64.next_up(), 0.3).unwrap(),
        Slip::stuck(),
    );
    parameters_refuse(
        Slider::new(10.0, 0.5, 0.3).unwrap(),
        Slider::new(10.0, 0.5, 0.3).unwrap().told(SliderEnd::B),
        Slip::stuck(),
    );
}

#[test]
fn signal_feedback_codecs_bind_every_gain_bias_and_published_observation() {
    let make = |matrix, bias| AffineSignal::new(2, 1, matrix, bias).unwrap();
    parameters_refuse(
        make(vec![3.0, -2.0], vec![0.25]),
        make(vec![3.0, (-2.0_f64).next_up()], vec![0.25]),
        (),
    );
    parameters_refuse(
        make(vec![3.0, -2.0], vec![0.25]),
        make(vec![3.0, -2.0], vec![0.25_f64.next_up()]),
        (),
    );
    parameters_refuse(
        Mass::new(2.0).unwrap(),
        Mass::new(2.0).unwrap().publishing(),
        Motion {
            position_m: -0.0,
            velocity_m_s: -0.25,
        },
    );
    let map = make(vec![3.0, -2.0], vec![0.25]);
    let mut words = vec![11];
    assert!(map.encode_state(&(), &mut words));
    assert_eq!(words, [11]);
    assert_eq!(map.decode_state(&[]), Some(()));
    assert_eq!(map.decode_state(&[0]), None);
}

#[test]
fn slider_codec_retains_prescribed_end_memory_and_refusal_history() {
    let slider = Slider::new(10.0, 0.5, 0.3).unwrap().told(SliderEnd::B);
    let mut state = Slip::stuck();
    let mut output = [0.0; 3];
    with_cx(|cx| {
        // A wrongly signed prescribed span is admitted for two observations
        // before refusal. Resuming must retain that history, not clear it.
        for _ in 0..2 {
            let mut next = Slip::stuck();
            slider
                .propose_on_line(
                    &state,
                    0.1,
                    &[0.0, 1.0, 300.0, 1.0, 1.0],
                    &[1.0, 0.0, 0.0, 0.0],
                    &mut next,
                    &mut output,
                    cx,
                )
                .unwrap();
            state = next;
        }
        let mut words = Vec::new();
        assert!(slider.encode_state(&state, &mut words));
        assert_eq!(words[7], 2);
        let resumed = slider.decode_state(&words).unwrap();
        let mut again = Vec::new();
        assert!(slider.encode_state(&resumed, &mut again));
        assert_eq!(words, again);
        for old in [&state, &resumed] {
            let mut next = Slip::stuck();
            assert!(
                slider
                    .propose_on_line(
                        old,
                        0.1,
                        &[0.0, 1.0, 300.0, 1.0, 1.0],
                        &[1.0, 0.0, 0.0, 0.0],
                        &mut next,
                        &mut output,
                        cx,
                    )
                    .is_err()
            );
        }
        for (index, invalid) in [
            (0, f64::NAN.to_bits()),
            (1, 2),
            (5, 2),
            (6, f64::INFINITY.to_bits()),
            (7, 3),
        ] {
            let mut bad = words.clone();
            bad[index] = invalid;
            assert!(slider.decode_state(&bad).is_none(), "word {index}");
        }
        assert!(
            Slider::new(10.0, 0.5, 0.3)
                .unwrap()
                .decode_state(&words)
                .is_none()
        );
    });
    let mut words = vec![123];
    assert!(!slider.encode_state(&Slip::sliding(f64::NAN), &mut words));
    assert_eq!(words, [123]);
}

struct Belt {
    system: System,
    mass: NodeId,
    slider: NodeId,
    belt: NodeId,
    drive: NodeId,
    bodies: [NodeId; 2],
    ports: Vec<PortRef>,
}

impl Belt {
    fn new() -> Self {
        let mut builder = SystemBuilder::new();
        let mass = builder
            .add_with_tangent(
                "mass",
                Mass::new(2.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 0.0,
                },
            )
            .unwrap();
        let slider = builder
            .add_on_line(
                "slider",
                Slider::new(19.6, 0.5, 0.3).unwrap().told(SliderEnd::B),
                Slip::stuck(),
            )
            .unwrap();
        let belt = builder
            .add_with_tangent(
                "belt",
                FlowSource::new(
                    Waveform::table(vec![
                        (0.0, 0.0),
                        (0.2, 1.0),
                        (0.4, 1.0),
                        (0.6, -0.5),
                        (0.8, 0.0),
                    ])
                    .unwrap(),
                )
                .publishing(),
                SourceAccount::default(),
            )
            .unwrap();
        let drive = builder
            .add_with_tangent(
                "drive",
                ControlledSource::effort(PortKind::MechanicalForceVelocity),
                SourceAccount::default(),
            )
            .unwrap();
        let plane = builder
            .add_with_tangent("plane", HeatNode::new(&[2.0, 4.0]).unwrap(), ())
            .unwrap();
        let bodies = [
            builder
                .add_with_tangent("pad", HeatCapacity::new(10.0).unwrap(), 300.0)
                .unwrap(),
            builder
                .add_with_tangent("disc", HeatCapacity::new(20.0).unwrap(), 300.0)
                .unwrap(),
        ];
        let port = |builder: &SystemBuilder, node, name| builder.port(node, name).unwrap();
        builder
            .common_flow(
                port(&builder, mass, "p"),
                &[port(&builder, slider, "a"), port(&builder, drive, "p")],
            )
            .unwrap();
        builder
            .bond(port(&builder, belt, "p"), port(&builder, slider, "b"))
            .unwrap();
        builder
            .signal(
                port(&builder, belt, "span"),
                port(&builder, slider, "span_b"),
            )
            .unwrap();
        builder
            .boundary(
                port(&builder, drive, "command"),
                Box::new(|_, _, values| values[0] = 0.4),
            )
            .unwrap();
        builder
            .bond(
                port(&builder, plane, "source"),
                port(&builder, slider, "heat"),
            )
            .unwrap();
        for (body, name) in bodies.into_iter().zip(["body0", "body1"]) {
            builder
                .bond(port(&builder, body, "heat"), port(&builder, plane, name))
                .unwrap();
        }
        let ports = [
            (mass, "p"),
            (slider, "a"),
            (slider, "b"),
            (slider, "heat"),
            (slider, "span_b"),
            (belt, "p"),
            (belt, "span"),
            (drive, "p"),
            (drive, "command"),
            (plane, "source"),
            (plane, "body0"),
            (plane, "body1"),
            (bodies[0], "heat"),
            (bodies[1], "heat"),
        ]
        .map(|(node, name)| port(&builder, node, name))
        .to_vec();
        let mut system = builder.seal().unwrap();
        assert!(system.guess(ports[3], &[300.0]));
        assert!(system.guess(ports[10], &[300.0]));
        assert!(system.guess(ports[11], &[300.0]));
        Self {
            system,
            mass,
            slider,
            belt,
            drive,
            bodies,
            ports,
        }
    }

    fn view(&self) -> String {
        format!(
            "{:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}",
            self.system.time_s().to_bits(),
            self.system.steps(),
            self.system.state::<Motion>(self.mass).unwrap(),
            self.system.state::<Slip>(self.slider).unwrap(),
            self.system.state::<SourceAccount>(self.belt).unwrap(),
            self.system.state::<SourceAccount>(self.drive).unwrap(),
            self.system.state::<f64>(self.bodies[0]).unwrap(),
            self.system.state::<f64>(self.bodies[1]).unwrap(),
            wires(&self.system, &self.ports),
            self.system.ledger()
        )
    }
}

#[test]
fn a_driven_friction_heat_network_resumes_with_every_native_history_intact() {
    with_cx(|cx| {
        let mut uninterrupted = Belt::new();
        for _ in 0..47 {
            uninterrupted.system.step(0.01, cx).unwrap();
        }
        let sealed = uninterrupted
            .system
            .snapshot(&context(&uninterrupted.system), limits(), || false)
            .unwrap();
        let mut resumed = Belt::new();
        resumed.system.step(0.03, cx).unwrap();
        resumed.system.resume(&opened_here(&sealed)).unwrap();
        assert_eq!(uninterrupted.view(), resumed.view());
        for _ in 0..80 {
            uninterrupted.system.step(0.01, cx).unwrap();
            resumed.system.step(0.01, cx).unwrap();
            assert_eq!(uninterrupted.view(), resumed.view());
        }
        assert!(
            *uninterrupted
                .system
                .state::<f64>(uninterrupted.bodies[0])
                .unwrap()
                > 300.0
        );
    });
}
