//! G1/G3: one original saturating fluid drives porous drag, compliance,
//! thermal relaxation and the native duct pressure/volume-flow response.
//! Both fluid curves and pore dimensions below are numerical fixtures. They
//! make no claim that a measured dry foam remains unchanged in a liquid.

use super::*;
use fs_duct::network::NetworkLoad;
use fs_material::fluid_transport::FluidTransport;
use fs_material::porous::{PorousCard, PorousFluidLaw};
use std::cell::Cell;

const RADIUS: f64 = 0.01;
const LENGTH: f64 = 0.35;
const DEPTH: f64 = 0.1;

fn pore_card() -> PorousCard {
    PorousCard {
        name: "manufactured rigid porous skeleton",
        porosity: 0.7,
        flow_resistivity: 1e6,
        tortuosity: 1.6,
        viscous_length_m: 1e-4,
        thermal_length_m: 2e-4,
        thermal_permeability_m2: Some(6e-9),
        provenance: "manufactured fully saturated rigid-frame regression",
        valid_hz: (50.0, 5000.0),
    }
}

fn boundary(card: PorousCard) -> NetworkLoad {
    NetworkLoad::PorousRigid {
        material: PorousFluidLaw::try_new(card, 0.001).unwrap(),
        thickness_m: DEPTH,
        area_m2: core::f64::consts::PI * RADIUS * RADIUS,
    }
}

// Independent permeability/diffusivity form of JCA(L), followed by the
// analytic rigid-backed layer solution. No production porous, wave-basis,
// surface-impedance, or network conversion helper supplies this reference.
fn reference_load(card: &PorousCard, state: ThermoviscousFluidState, omega: f64) -> C64 {
    let rho = state.acoustic().density_kg_m3();
    let sound = state.acoustic().sound_speed_m_s();
    let transport = state.transport();
    let mu = transport.dynamic_viscosity_pa_s();
    let diffusivity =
        transport.thermal_conductivity_w_m_k() / (rho * transport.specific_heat_j_kg_k());
    let ratio = state.heat_capacity_ratio();
    let phi = card.porosity;
    let tortuosity = card.tortuosity;
    let permeability = 0.001 / card.flow_resistivity;
    let viscous = C64::new(
        1.0,
        4.0 * tortuosity.powi(2) * rho * omega * permeability.powi(2)
            / (mu * card.viscous_length_m.powi(2) * phi.powi(2)),
    )
    .sqrt();
    let density = C64::from_re(rho * tortuosity / phi)
        + C64::new(0.0, -mu / (permeability * omega)) * viscous;
    let thermal_permeability = card
        .thermal_permeability_m2
        .unwrap_or(phi * card.thermal_length_m.powi(2) / 8.0);
    let thermal = C64::new(
        1.0,
        4.0 * thermal_permeability.powi(2) * omega
            / (diffusivity * card.thermal_length_m.powi(2) * phi.powi(2)),
    )
    .sqrt();
    let relaxation =
        C64::ONE + C64::new(0.0, -phi * diffusivity / (omega * thermal_permeability)) * thermal;
    let bulk = C64::from_re(rho * sound * sound / phi)
        / (C64::from_re(ratio) - C64::from_re(ratio - 1.0) / relaxation);
    let characteristic = (density * bulk).sqrt();
    let phase = (density / bulk).sqrt().scale(omega * DEPTH);
    let cosine = C64::new(
        phase.re.cos() * phase.im.cosh(),
        -phase.re.sin() * phase.im.sinh(),
    );
    let sine = C64::new(
        phase.re.sin() * phase.im.cosh(),
        phase.re.cos() * phase.im.sinh(),
    );
    let surface = C64::new(0.0, -1.0) * characteristic * cosine / sine;
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    C64::new(surface.re / area, -surface.im / area)
}

#[test]
fn original_fluid_and_pore_data_change_actual_pressure_flow_and_absorbed_power() {
    let (source, cold) = source_liquid();
    let original = cold.clone();
    let hot = cold.at_query_point(&source, &point(400.0)).unwrap();
    let duct = cylinder(RADIUS, LENGTH);
    let mut changed = pore_card();
    changed.flow_resistivity *= 2.0;
    let mut signatures = Vec::new();
    for (liquid, card) in [(&cold, pore_card()), (&hot, pore_card()), (&cold, changed)] {
        let state = liquid.thermoviscous_state().unwrap();
        let acoustic = state.acoustic();
        let zc =
            acoustic.characteristic_impedance_pa_s_m() / (core::f64::consts::PI * RADIUS * RADIUS);
        let load = boundary(card.clone());
        for frequency in [100.0, 700.0, 2000.0] {
            let omega = core::f64::consts::TAU * frequency;
            let expected_load = reference_load(&card, state, omega);
            assert!(expected_load.re > 0.0);
            let theta = omega * LENGTH / acoustic.sound_speed_m_s();
            let denominator =
                C64::new(0.0, -theta.sin() / zc) * expected_load + C64::from_re(theta.cos());
            let expected_flow = C64::ONE / denominator;
            let expected_pressure = expected_load * expected_flow;
            let expected_input =
                (expected_load.scale(theta.cos()) + C64::new(0.0, -zc * theta.sin())) / denominator;
            let actual = fluid::input_impedance_boundary(
                &duct,
                liquid,
                omega,
                LossModel::Lossless,
                &load,
                WallModel::Rigid,
                ConeSlices::Auto,
                options(),
                || false,
            )
            .unwrap();
            near(actual.z_in, expected_input, zc, 3e-11);
            near(actual.u_out_over_u_in, expected_flow, 1.0, 3e-11);
            near(actual.p_out_over_u_in, expected_pressure, zc, 3e-11);
            let porous_power = 0.5 * expected_load.re * actual.u_out_over_u_in.norm_sq();
            let input_power = 0.5 * actual.z_in.re;
            assert!(porous_power > 0.0);
            assert!((input_power - porous_power).abs() < 3e-10 * input_power);

            let lossy = fluid::input_impedance_boundary(
                &duct,
                liquid,
                omega,
                LossModel::Bessel,
                &load,
                WallModel::Rigid,
                ConeSlices::Auto,
                options(),
                || false,
            )
            .unwrap();
            let absorbed = 0.5 * expected_load.re * lossy.u_out_over_u_in.norm_sq();
            assert!(
                0.5 * lossy.z_in.re > absorbed,
                "native wall diffusion must remain a separate loss from the porous load"
            );
            near(
                lossy.p_out_over_u_in,
                expected_load * lossy.u_out_over_u_in,
                zc,
                3e-11,
            );
            if frequency == 700.0 {
                signatures.push((actual.z_in, actual.u_out_over_u_in, porous_power));
            }
        }
    }
    for changed in &signatures[1..] {
        assert!((changed.0 - signatures[0].0).abs() > 1e-4 * signatures[0].0.abs());
        assert!((changed.1 - signatures[0].1).abs() > 1e-4);
        assert!((changed.2 / signatures[0].2 - 1.0).abs() > 1e-4);
    }
    assert!(cold.at_query_point(&source, &point(401.0)).is_err());
    assert_eq!(cold, original);
}

struct SingleRead {
    state: ThermoviscousFluidState,
    reads: Cell<usize>,
}
impl FluidAcoustics for SingleRead {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        panic!("consumer must request one complete fluid snapshot")
    }
}
impl FluidTransport for SingleRead {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        panic!("consumer must request one complete fluid snapshot")
    }
}
impl ThermoviscousFluid for SingleRead {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        if self.reads.get() == 1 {
            Ok(self.state)
        } else {
            Err(MaterialError::State {
                what: "original source was read twice".into(),
            })
        }
    }
}

#[test]
fn porous_boundary_uses_one_snapshot_and_refuses_before_partial_results() {
    let (_, liquid) = source_liquid();
    let state = liquid.thermoviscous_state().unwrap();
    let duct = cylinder(RADIUS, LENGTH);
    let load = boundary(pore_card());
    let omega = core::f64::consts::TAU * 700.0;
    let once = SingleRead {
        state,
        reads: Cell::new(0),
    };
    let solve = |medium: &SingleRead, load: &NetworkLoad, omega| {
        fluid::input_impedance_boundary(
            &duct,
            medium,
            omega,
            LossModel::Bessel,
            load,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || false,
        )
    };
    let complete = solve(&once, &load, omega).unwrap();
    assert_eq!(once.reads.get(), 1);
    assert!(matches!(
        solve(&once, &load, omega),
        Err(FluidDuctError::Material(_))
    ));
    let retry = SingleRead {
        state,
        reads: Cell::new(0),
    };
    assert_eq!(solve(&retry, &load, omega).unwrap().z_in, complete.z_in);

    let cancelled = SingleRead {
        state,
        reads: Cell::new(0),
    };
    assert!(matches!(
        fluid::input_impedance_boundary(
            &duct,
            &cancelled,
            omega,
            LossModel::Bessel,
            &load,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || true
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert_eq!(cancelled.reads.get(), 0);
    let mut calls = 0;
    assert!(matches!(
        fluid::input_impedance_boundary(
            &duct,
            &cancelled,
            omega,
            LossModel::Bessel,
            &load,
            WallModel::Rigid,
            ConeSlices::Auto,
            options(),
            || {
                calls += 1;
                calls >= 3
            }
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    assert_eq!(cancelled.reads.get(), 1);

    for (area, thickness, frequency) in [
        (0.0, DEPTH, 700.0),
        (1.0, -DEPTH, 700.0),
        (1.0, DEPTH, 20.0),
    ] {
        let bad = NetworkLoad::PorousRigid {
            material: PorousFluidLaw::try_new(pore_card(), 0.001).unwrap(),
            thickness_m: thickness,
            area_m2: area,
        };
        let medium = SingleRead {
            state,
            reads: Cell::new(0),
        };
        assert!(solve(&medium, &bad, core::f64::consts::TAU * frequency).is_err());
    }
    let medium = SingleRead {
        state,
        reads: Cell::new(0),
    };
    assert!(
        solve(
            &medium,
            &NetworkLoad::Radiation {
                radius: RADIUS,
                termination: fs_duct::Termination::FlangedOpen,
            },
            omega
        )
        .is_err()
    );
}
