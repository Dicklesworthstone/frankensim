//! G0/G3 tests of the actual StreamPort/enthalpy-chart path.

use super::*;
use crate::{
    ChemicalEnergyInput, CoordinateBinding, PortOrientation, PortTimestamp,
    StreamStressWorkConvention,
};
use fs_qty::chemistry::SpeciesId;

fn id(value: &str) -> StableId {
    StableId::new(value).unwrap()
}
fn species() -> StreamConstituentId {
    StreamConstituentId::Species(SpeciesId::new("fixture-species").unwrap())
}
fn binding(frame: &str) -> StreamChartBinding {
    StreamChartBinding::try_new(
        id("surface/face"),
        id("source/state"),
        id("source/constituent-basis"),
        [species()],
        id("source/chemical-reference"),
        CoordinateBinding::new(
            id("basis/cartesian"),
            id(frame),
            PortOrientation::OutwardFromOwner,
        ),
        PortTimestamp::new(id("clock/surface"), 3),
        id("gravity/datum"),
        StreamStressWorkConvention::CauchyTensionPositiveOutwardPower,
    )
    .unwrap()
}
fn chemistry(binding: &StreamChartBinding) -> ChemicalEnergyAccounting {
    ChemicalEnergyAccounting::try_new(
        binding,
        ChemicalEnergyInput::IncludedInStatePotential {
            reference_state: binding.chemical_reference_state().clone(),
        },
    )
    .unwrap()
}
fn assemble_in(
    frame: &str,
    surface: StreamControlSurface,
    state: StreamSurfaceState,
) -> Result<StreamPort, CoupleError> {
    let bound = binding(frame);
    let chemical = chemistry(&bound);
    StreamPort::try_from_uniform_surface(
        bound,
        surface,
        state,
        [AdvectedConstituent {
            id: species(),
            amount_per_mass: SpecificAmount::new(0.25),
        }],
        chemical,
        id("source/stress-operator"),
        id("source/stress-evidence"),
    )
}
fn assemble(surface: StreamControlSurface, state: StreamSurfaceState) -> StreamPort {
    assemble_in("frame/lab", surface, state).unwrap()
}
fn surface() -> StreamControlSurface {
    StreamControlSurface {
        area: Area::new(2.0),
        outward_normal: [1.0, 0.0, 0.0],
        velocity: [Velocity::new(0.0); 3],
        position: [Length::new(0.0); 3],
    }
}
fn state(velocity: [f64; 3]) -> StreamSurfaceState {
    StreamSurfaceState {
        density: Density::new(2.0),
        pressure: Pressure::new(4.0),
        specific_enthalpy: SpecificEnergy::new(50.0),
        kinematics: StreamKinematics::try_new(
            velocity.map(Velocity::new),
            SpecificEnergy::new(6.0),
        )
        .unwrap(),
        specific_entropy: SpecificEntropy::new(7.0),
        deviatoric_traction: [Pressure::new(0.0); 3],
    }
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * (1.0 + expected.abs()),
        "actual={actual:e} expected={expected:e}"
    );
}
fn vector_near(actual: [f64; 3], expected: [f64; 3]) {
    for i in 0..3 {
        near(actual[i], expected[i]);
    }
}
fn ordinary_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[test]
fn moving_impermeable_wall_supplies_pressure_impulse_and_conjugate_work() {
    let wall_mass = 2.0;
    let duration = 0.25;
    let force = 8.0;
    for before in [-2.0, 1.0] {
        let after = before + force * duration / wall_mass;
        let midpoint = 0.5 * (before + after);
        let mut cut = surface();
        cut.velocity[0] = Velocity::new(midpoint);
        let port = assemble(cut, state([midpoint, 0.0, 0.0]));
        assert_eq!(port.mass_flow().value(), 0.0);
        assert_eq!(port.constituent_flows()[0].amount_flow().value(), 0.0);
        assert_eq!(port.entropy_flow().value(), 0.0);
        assert_eq!(port.momentum_flow().map(Force::value), [force, 0.0, 0.0]);
        assert_eq!(port.energy_flow().value(), force * midpoint);
        let transport = port.surface_transport().unwrap();
        assert_eq!(transport.pressure_motion_work().value(), force * midpoint);
        assert_eq!(transport.stress_power_on_surface(), port.energy_flow());
        // The same outward fluid impulse acts on the finite wall. Its midpoint
        // velocity makes work equal its independently evaluated kinetic change.
        near(
            wall_mass * (after - before),
            port.momentum_flow()[0].value() * duration,
        );
        near(
            0.5 * wall_mass * (after * after - before * before),
            port.energy_flow().value() * duration,
        );
        assert_eq!(
            port.energy_flow().value().is_sign_positive(),
            midpoint > 0.0
        );
    }
}

#[test]
fn open_moving_cut_obeys_the_inertial_energy_and_momentum_boost() {
    let mut cut = surface();
    cut.outward_normal = [0.6, 0.8, 0.0];
    cut.velocity = [1.0, -2.0, 4.0].map(Velocity::new);
    let mut material = state([5.0, -3.0, 2.0]);
    material.density = Density::new(2.5);
    material.pressure = Pressure::new(3.0);
    material.specific_enthalpy = SpecificEnergy::new(40.0);
    material.specific_entropy = SpecificEntropy::new(-2.0);
    material.deviatoric_traction = [0.5, -0.75, 1.25].map(Pressure::new);
    let original = assemble(cut, material);
    let boost = [7.0, -1.0, -4.0];
    let u = material.kinematics.velocity.map(Velocity::value);
    let w = cut.velocity.map(Velocity::value);
    let mut moved_material = material;
    moved_material.kinematics = StreamKinematics::try_new(
        core::array::from_fn(|i| Velocity::new(u[i] - boost[i])),
        SpecificEnergy::new(6.0),
    )
    .unwrap();
    let mut moved_cut = cut;
    moved_cut.velocity = core::array::from_fn(|i| Velocity::new(w[i] - boost[i]));
    let shifted = assemble_in("frame/boosted", moved_cut, moved_material).unwrap();
    let mass = original.mass_flow().value();
    let momentum = original.momentum_flow().map(Force::value);
    near(shifted.mass_flow().value(), mass);
    vector_near(
        shifted.momentum_flow().map(Force::value),
        core::array::from_fn(|i| momentum[i] - mass * boost[i]),
    );
    near(
        shifted.energy_flow().value(),
        original.energy_flow().value() - ordinary_dot(boost, momentum)
            + 0.5 * mass * ordinary_dot(boost, boost),
    );
    near(
        shifted.entropy_flow().value(),
        original.entropy_flow().value(),
    );
    near(
        original.constituent_flows()[0].amount_flow().value(),
        0.25 * mass,
    );
    near(original.entropy_flow().value(), -2.0 * mass);

    // Independent internal-energy/Cauchy form: no pressure-flow term is added
    // to enthalpy a second time, and stress work uses material velocity.
    let pressure_force = cut
        .outward_normal
        .map(|n| cut.area.value() * material.pressure.value() * n);
    let extra_force = material
        .deviatoric_traction
        .map(|t| cut.area.value() * t.value());
    let internal =
        material.specific_enthalpy.value() - material.pressure.value() / material.density.value();
    let expected = mass * (internal + 0.5 * ordinary_dot(u, u) + 6.0)
        + ordinary_dot(pressure_force, u)
        - ordinary_dot(extra_force, u);
    near(original.energy_flow().value(), expected);
}

#[test]
fn uniform_flow_on_a_translating_expanding_cell_retains_geometric_conservation() {
    let u = [3.0, -2.0, 1.0];
    let expansion_rate = 0.25;
    let mut material = state(u);
    material.pressure = Pressure::new(3.0);
    material.specific_enthalpy = SpecificEnergy::new(20.0);
    material.kinematics =
        StreamKinematics::try_new(u.map(Velocity::new), SpecificEnergy::new(0.5)).unwrap();
    let mut total_mass = 0.0;
    let mut total_momentum = [0.0; 3];
    let mut total_energy = 0.0;
    let mut volume_rate = 0.0;
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let mut cut = surface();
            cut.area = Area::new(4.0);
            cut.outward_normal = [0.0; 3];
            cut.outward_normal[axis] = sign;
            cut.position[axis] = Length::new(sign);
            cut.velocity = u.map(Velocity::new);
            cut.velocity[axis] = Velocity::new(u[axis] + sign * expansion_rate);
            let port = assemble(cut, material);
            total_mass += port.mass_flow().value();
            for (d, component) in total_momentum.iter_mut().enumerate() {
                *component += port.momentum_flow()[d].value();
            }
            total_energy += port.energy_flow().value();
            volume_rate += port
                .surface_transport()
                .unwrap()
                .swept_volume_flow()
                .value();
        }
    }
    assert_eq!(volume_rate, 6.0);
    let density = material.density.value();
    near(total_mass, -density * volume_rate);
    vector_near(total_momentum, u.map(|v| -density * v * volume_rate));
    let internal = material.specific_enthalpy.value() - material.pressure.value() / density;
    near(
        total_energy,
        -density * (internal + 0.5 * ordinary_dot(u, u) + 0.5) * volume_rate,
    );
}

#[test]
fn angular_flux_changes_origin_and_stress_wrench_pairs_with_rigid_surface_velocity() {
    let origin = [2.0, -1.0, 3.0];
    let arm = [1.0, 2.0, -1.0];
    let translation = [2.0, -1.0, 0.5];
    let omega = [0.5, -0.25, 1.0];
    let rotation_velocity = cross(omega, arm);
    let w = core::array::from_fn(|i| translation[i] + rotation_velocity[i]);
    let mut cut = surface();
    cut.position = core::array::from_fn(|i| Length::new(origin[i] + arm[i]));
    cut.velocity = w.map(Velocity::new);
    let mut material = state([w[0], w[1] + 1.0, w[2] - 2.0]);
    material.pressure = Pressure::new(5.0);
    material.deviatoric_traction = [1.0, 2.0, -3.0].map(Pressure::new);
    let port = assemble(cut, material);
    let transport = port.surface_transport().unwrap();
    let force = transport.stress_momentum_flow().map(Force::value);
    let moment = transport
        .stress_moment_about(origin.map(Length::new))
        .unwrap()
        .map(AngularMomentumFlow::value);
    vector_near(force, [8.0, -4.0, 6.0]);
    vector_near(moment, [8.0, -14.0, -20.0]);
    near(
        transport.stress_power_on_surface().value(),
        ordinary_dot(force, translation) + ordinary_dot(moment, omega),
    );
    near(transport.stress_power_on_surface().value(), 10.5);
    near(transport.stress_power_out_of_material().value(), -5.5);
    near(port.energy_flow().value(), -5.5);
    // Slip work is exposed, not relabelled as surface work or inferred heat.
    assert_ne!(transport.stress_power_on_surface(), port.energy_flow());

    // Open the same face: transported matter must carry its own moment too.
    material.kinematics = StreamKinematics::try_new(
        [w[0] + 2.0, w[1] + 1.0, w[2] - 2.0].map(Velocity::new),
        SpecificEnergy::new(6.0),
    )
    .unwrap();
    let open = assemble(cut, material);
    let total = open.momentum_flow().map(Force::value);
    let transport = open.surface_transport().unwrap();
    let about_origin = transport
        .angular_momentum_flow_about(origin.map(Length::new))
        .unwrap()
        .map(AngularMomentumFlow::value);
    vector_near(about_origin, cross(arm, total));
    let shift = [-3.0, 2.0, 4.0];
    let new_origin = core::array::from_fn(|i| Length::new(origin[i] + shift[i]));
    let about_shifted = transport
        .angular_momentum_flow_about(new_origin)
        .unwrap()
        .map(AngularMomentumFlow::value);
    let shift_moment = cross(shift, total);
    vector_near(
        about_shifted,
        core::array::from_fn(|i| about_origin[i] - shift_moment[i]),
    );
}

#[test]
fn reversing_the_same_surface_reverses_every_transport_without_losing_pressure() {
    let mut cut = surface();
    cut.velocity = [-1.0, 0.5, 0.0].map(Velocity::new);
    let mut material = state([3.0, 2.0, -1.0]);
    material.deviatoric_traction = [0.25, -1.0, 0.5].map(Pressure::new);
    let forward = assemble(cut, material);
    cut.outward_normal = cut.outward_normal.map(|x| -x);
    material.deviatoric_traction = material
        .deviatoric_traction
        .map(|x| Pressure::new(-x.value()));
    let reverse = assemble(cut, material);
    near(reverse.mass_flow().value(), -forward.mass_flow().value());
    vector_near(
        reverse.momentum_flow().map(Force::value),
        forward.momentum_flow().map(|x| -x.value()),
    );
    near(
        reverse.energy_flow().value(),
        -forward.energy_flow().value(),
    );
    near(
        reverse.entropy_flow().value(),
        -forward.entropy_flow().value(),
    );
    near(
        reverse.constituent_flows()[0].amount_flow().value(),
        -forward.constituent_flows()[0].amount_flow().value(),
    );
    near(
        reverse
            .surface_transport()
            .unwrap()
            .pressure_motion_work()
            .value(),
        -forward
            .surface_transport()
            .unwrap()
            .pressure_motion_work()
            .value(),
    );
}

#[test]
fn stationary_cut_reuses_the_existing_enthalpy_chart_without_pressure_double_counting() {
    let mut material = state([3.0, -1.0, 2.0]);
    material.deviatoric_traction = [0.5, -0.25, 1.0].map(Pressure::new);
    let port = assemble(surface(), material);
    let bound = binding("frame/lab");
    let legacy = MovingStreamEnthalpyChart::try_new(
        bound.clone(),
        material.specific_enthalpy,
        material.kinematics,
        DeviatoricStressWork::try_new(
            &bound,
            Power::new(-7.5),
            id("source/stress-operator"),
            id("source/stress-evidence"),
        )
        .unwrap(),
        chemistry(&bound),
    )
    .unwrap();
    assert_eq!(
        port.surface_transport()
            .unwrap()
            .pressure_motion_work()
            .value(),
        0.0
    );
    assert_eq!(port.energy_flow(), legacy.energy_rate(port.mass_flow()));
}

#[test]
fn invalid_geometry_and_replaced_mass_or_momentum_are_refused() {
    let mut bad = surface();
    bad.outward_normal = [2.0, 0.0, 0.0];
    assert!(matches!(
        assemble_in("frame/lab", bad, state([1.0, 0.0, 0.0])),
        Err(CoupleError::InvalidStreamSurface { .. })
    ));
    bad = surface();
    bad.area = Area::new(-1.0);
    assert!(matches!(
        assemble_in("frame/lab", bad, state([1.0, 0.0, 0.0])),
        Err(CoupleError::InvalidStreamSurface { .. })
    ));
    bad = surface();
    bad.velocity[1] = Velocity::new(f64::INFINITY);
    assert!(matches!(
        assemble_in("frame/lab", bad, state([1.0, 0.0, 0.0])),
        Err(CoupleError::NonFiniteStreamValue { .. })
    ));
    let mut material = state([1.0, 0.0, 0.0]);
    material.density = Density::new(0.0);
    assert!(matches!(
        assemble_in("frame/lab", surface(), material),
        Err(CoupleError::NonPositiveStreamDensity)
    ));
    let port = assemble(surface(), state([1.0, 0.0, 0.0]));
    let rebuild = |mass, momentum| {
        StreamPort::try_new(
            port.binding().clone(),
            mass,
            port.constituent_flows().to_vec(),
            momentum,
            port.energy_flow(),
            port.entropy_flow(),
            port.energy_chart().clone(),
        )
    };
    assert!(matches!(
        rebuild(MassFlowRate::new(2.0), port.momentum_flow()),
        Err(CoupleError::StreamSurfaceFlowMismatch {
            field: "mass_flow",
            index: None
        })
    ));
    let mut momentum = port.momentum_flow();
    momentum[2] = Force::new(1.0);
    assert!(matches!(
        rebuild(port.mass_flow(), momentum),
        Err(CoupleError::StreamSurfaceFlowMismatch {
            field: "momentum_flow",
            index: Some(2)
        })
    ));
    assert!(
        port.surface_transport()
            .unwrap()
            .angular_momentum_flow_about([Length::new(f64::NAN); 3])
            .is_err()
    );
}
