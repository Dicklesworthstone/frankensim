//! Coherent evolving turbulence in the real wing/blade/shaft/body force path.
use fs_airscrew::dynamic_inflow::first_harmonic::HarmonicInduction;
use fs_atmo::FlatSiteLogLaw;
use fs_atmo::mann::MANN_TARGET_V1;
use fs_atmo::mann_field::fixed_step::FixedStepMannField;
use fs_atmo::mann_field::{HomogeneousMannField, MannQuadrature, MannSpectralModel};
use fs_flyer::aircraft::{OpenLoopDesign, wright_openloop_v1};
use fs_flyer::atmosphere::mann_flight::{MannFlightBoundary, MannFlightWind};
use fs_flyer::spine::SixDofState;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor};
use fs_flyer::spine::rotors::propulsion::dynamic::first_harmonic::{
    HarmonicPropulsion, HarmonicState,
};
use fs_flyer::spine::rotors::propulsion::dynamic::{
    BladeIncident, DynamicPropulsion, DynamicState,
};
use fs_flyer::spine::rotors::propulsion::{
    GearedPropulsion, PropellerHandedness, PropellerMount, PropulsionInputs,
};
use fs_ga::Mat3;
use std::cell::Cell;

#[path = "mann_flight/unsteady.rs"]
mod unsteady;

const MODES: usize = 16;
const BUDGET: usize = 3 * 88 * MODES;
const WORK: usize = 8 * MODES;
const RHO: f64 = 1.225;
const ORIGIN: [f64; 3] = [-5.0, 3.0, 2.0];

fn wind() -> MannFlightWind {
    let model = MannSpectralModel::build(
        17,
        fs_atmo::mann::MannParams {
            alpha_eps23: 0.01,
            ..MANN_TARGET_V1
        },
        MannQuadrature {
            n_modes: MODES,
            min_kl: 0.1,
            max_kl: 16.0,
        },
        [7.0, 0.0, 0.0],
    )
    .unwrap();
    let mut snapshot = HomogeneousMannField::at_anchor(model, 29, vec![0.02; MODES]).unwrap();
    snapshot.advance_to(0).unwrap();
    let source = FixedStepMannField::from_snapshot(snapshot, 91, 0.002, 88 * MODES).unwrap();
    MannFlightWind::prepare(
        source,
        FlatSiteLogLaw {
            scenario_effective_z0_m: 0.001,
            displacement_height_m: 0.0,
            reference_height_m: 3.0,
            reference_speed_mps: 1.0,
        },
        MannFlightBoundary::GroundPlane,
        0.0,
        BUDGET,
    )
    .unwrap()
}

fn assembly(design: &OpenLoopDesign) -> GearedAssembly<2> {
    GearedAssembly::new(
        design.gross_mass_kg,
        Mat3 {
            m: [1800.0, 0.0, 0.0, 0.0, 420.0, 0.0, 0.0, 0.0, 1850.0],
        },
        core::array::from_fn(|k| GearedRotor {
            axis_body: [1.0, 0.0, 0.0],
            ratio: [1.0, -1.0][k],
            axial_inertia_kgm2: [1.5, 2.5][k],
        }),
    )
    .unwrap()
}

fn mounts() -> [PropellerMount; 2] {
    core::array::from_fn(|k| PropellerMount {
        rotor_index: k,
        thrust_direction_body: [1.0, 0.0, 0.0],
        handedness: if k == 0 {
            PropellerHandedness::Positive
        } else {
            PropellerHandedness::Negative
        },
    })
}

fn state(assembly: &GearedAssembly<2>) -> DynamicState {
    DynamicState {
        motion: assembly
            .state_from_speed(
                SixDofState {
                    pos_m: [1.0, -2.0, -10.0],
                    vel_mps: [12.0, 0.1, 0.8],
                    quat: [1.0, 0.0, 0.0, 0.0],
                    omega_body: [0.04, -0.02, 0.03],
                },
                0.2,
                35.0,
            )
            .unwrap(),
        induced_mps: [0.4, 0.8],
    }
}

fn inputs() -> PropulsionInputs {
    PropulsionInputs {
        density_kg_m3: RHO,
        canard_rad: 0.03,
        motor_nm: 100.0,
        gravity_world_mps2: [0.0, 0.0, 9.80665],
    }
}

#[test]
fn native_dynamic_flight_uses_evolved_midpoint_wind_at_every_blade_and_replays() {
    let design = wright_openloop_v1();
    let assembly = assembly(&design);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts()).unwrap();
    let model = DynamicPropulsion::new_oblique(&base, RHO, 16)
        .unwrap()
        .with_blade_incident(BladeIncident::WindField)
        .unwrap();
    let state = state(&assembly);
    let mut atmosphere = wind();
    let start = atmosphere.time_s();
    let dt = atmosphere.dt_s();
    let checkpoint = atmosphere.checkpoint();
    let held = atmosphere.clone();

    // Independent time ownership: advance the original producer once, then
    // feed its genuine start and midpoint states directly to the real solver.
    let mut source = FixedStepMannField::resume(checkpoint.field.clone(), 88 * MODES).unwrap();
    source
        .advance_to(source.step() + 1, WORK / 2, || false)
        .unwrap();
    let mut midpoint_checkpoint = checkpoint.clone();
    midpoint_checkpoint.field = source.checkpoint();
    let midpoint = MannFlightWind::resume(midpoint_checkpoint, BUDGET).unwrap();
    let queries = [Cell::new(0), Cell::new(0)];
    let direct = model
        .step(
            &state,
            start,
            dt,
            ORIGIN,
            |_, _| Ok(inputs()),
            |time, point| {
                let stage = if time == start {
                    0
                } else {
                    assert_eq!(time, start + 0.5 * dt);
                    1
                };
                queries[stage].set(queries[stage].get() + 1);
                [&held, &midpoint][stage].sample(point).map(|flow| flow.u)
            },
        )
        .unwrap();
    let actual = atmosphere
        .step_dynamic(&model, &state, ORIGIN, WORK, || false, |_, _| Ok(inputs()))
        .unwrap();
    assert_eq!(actual, direct);
    assert!(
        queries.iter().all(|count| count.get() > 128),
        "native blade quadrature must query the shared field"
    );
    source
        .advance_to(source.step() + 1, WORK / 2, || false)
        .unwrap();
    assert_eq!(atmosphere.checkpoint().field, source.checkpoint());

    // Holding the accepted field still samples moving material points, but it
    // loses the actual temporal evolution. Both start solves must agree.
    let frozen = model
        .step(
            &state,
            start,
            dt,
            ORIGIN,
            |_, _| Ok(inputs()),
            |_, point| held.sample(point).map(|flow| flow.u),
        )
        .unwrap();
    assert_eq!(actual.initial, frozen.initial);
    assert_ne!(
        actual.midpoint.aerodynamics.wing_lift_n.to_bits(),
        frozen.midpoint.aerodynamics.wing_lift_n.to_bits()
    );
    assert_ne!(
        actual.midpoint.aerodynamics.propellers[0]
            .torque_nm
            .to_bits(),
        frozen.midpoint.aerodynamics.propellers[0]
            .torque_nm
            .to_bits()
    );
    assert_ne!(actual.state, frozen.state);
    assert!(
        actual
            .midpoint
            .aerodynamics
            .oblique
            .unwrap()
            .iter()
            .any(|load| load.incident_variation_work_w.abs() > 1e-10)
    );

    let mut restored = MannFlightWind::resume(atmosphere.checkpoint(), BUDGET).unwrap();
    let next = atmosphere
        .step_dynamic(
            &model,
            &actual.state,
            ORIGIN,
            WORK,
            || false,
            |_, _| Ok(inputs()),
        )
        .unwrap();
    let replay = restored
        .step_dynamic(
            &model,
            &actual.state,
            ORIGIN,
            WORK,
            || false,
            |_, _| Ok(inputs()),
        )
        .unwrap();
    assert_eq!(next, replay);
    assert_eq!(atmosphere.checkpoint(), restored.checkpoint());
}

#[test]
fn harmonic_flight_and_wind_refuse_together_then_retry_the_same_turbulent_path() {
    let design = wright_openloop_v1();
    let assembly = assembly(&design);
    let base = GearedPropulsion::new_spatial(&design, &assembly, mounts()).unwrap();
    let model = HarmonicPropulsion::new(&base, RHO, 16)
        .unwrap()
        .with_blade_incident(BladeIncident::WindField)
        .unwrap();
    let state = HarmonicState {
        motion: state(&assembly).motion,
        induction: [
            HarmonicInduction {
                mean_mps: 0.4,
                first_mps: [0.2, -0.15],
            },
            HarmonicInduction {
                mean_mps: 0.8,
                first_mps: [-0.25, 0.12],
            },
        ],
    };
    let before_state = state;
    let mut atmosphere = wind();
    let before_wind = atmosphere.checkpoint();
    let start = atmosphere.time_s();
    let calls = Cell::new(0);
    // The real wake owner refuses a density substitution at the midpoint,
    // after the first complete wing/blade/mechanical prediction succeeded.
    assert!(
        atmosphere
            .step_harmonic(
                &model,
                &state,
                ORIGIN,
                WORK,
                || false,
                |time, _| {
                    calls.set(calls.get() + 1);
                    let mut input = inputs();
                    if time != start {
                        input.density_kg_m3 *= 1.1;
                    }
                    Ok(input)
                }
            )
            .is_err()
    );
    assert_eq!(calls.get(), 2);
    assert_eq!(atmosphere.checkpoint(), before_wind);
    assert_eq!(state, before_state);

    // Cancellation after the actual second force solve also withholds the
    // endpoint source and the completed mechanical result.
    let finish = Cell::new(false);
    assert_eq!(
        atmosphere
            .step_harmonic(
                &model,
                &state,
                ORIGIN,
                WORK,
                || finish.get(),
                |time, _| {
                    if time != start {
                        finish.set(true);
                    }
                    Ok(inputs())
                }
            )
            .unwrap_err()
            .code,
        "mann-flight-cancelled"
    );
    assert_eq!(atmosphere.checkpoint(), before_wind);
    assert_eq!(state, before_state);

    let mut restored = MannFlightWind::resume(before_wind, BUDGET).unwrap();
    let frozen = model
        .step(
            &state,
            start,
            atmosphere.dt_s(),
            ORIGIN,
            |_, _| Ok(inputs()),
            |_, point| atmosphere.sample(point).map(|flow| flow.u),
        )
        .unwrap();
    let actual = atmosphere
        .step_harmonic(&model, &state, ORIGIN, WORK, || false, |_, _| Ok(inputs()))
        .unwrap();
    let replay = restored
        .step_harmonic(&model, &state, ORIGIN, WORK, || false, |_, _| Ok(inputs()))
        .unwrap();
    assert_eq!(actual, replay);
    assert_eq!(atmosphere.checkpoint(), restored.checkpoint());
    assert_eq!(actual.initial, frozen.initial);
    assert_ne!(
        actual.midpoint.aerodynamics.propellers[0]
            .torque_nm
            .to_bits(),
        frozen.midpoint.aerodynamics.propellers[0]
            .torque_nm
            .to_bits()
    );
    assert_ne!(actual.state, frozen.state);
}
