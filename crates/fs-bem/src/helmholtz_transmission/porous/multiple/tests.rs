use super::*;
use crate::helmholtz::bc::{IncidentSource, Monopole, PlaneWave};
use crate::helmholtz::molecular::MolecularNearFieldOptions;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::PorousCard;
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};

#[path = "tests/reference.rs"]
mod reference;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn law() -> PorousFluidLaw {
    PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap()
}
fn plane(amplitude: C64) -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: amplitude,
    })
}
fn controls() -> MolecularNearFieldOptions {
    MolecularNearFieldOptions {
        quadrature: near_field::Options {
            relative_tolerance: 2e-7,
            maximum_kernel_evaluations: 2_000_000,
            ..near_field::Options::default()
        },
        ..MolecularNearFieldOptions::default()
    }
}
fn sphere(radius: f64, z: f64, level: u32) -> SpherePanels {
    let original = SpherePanels::icosphere(radius, level).unwrap();
    let triangles = original
        .triangles()
        .unwrap()
        .iter()
        .map(|t| t.map(|v| [v[0], v[1], v[2] + z]))
        .collect();
    SpherePanels::from_triangles(triangles).unwrap()
}
fn porous<'a>(surface: &'a SpherePanels, law: &'a PorousFluidLaw, gas: &'a GasState) -> Body<'a> {
    Body {
        surface,
        material: BodyMaterial::Porous {
            law,
            filling: PoreFluid::Fluid(gas),
        },
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-23),
        "{a:?} != {b:?}"
    );
}
fn relative(actual: &[C64], expected: &[C64]) -> f64 {
    let scale = expected.iter().map(|v| v.norm_sq()).sum::<f64>().sqrt();
    actual
        .iter()
        .zip(expected)
        .map(|(a, b)| (*a - *b).norm_sq())
        .sum::<f64>()
        .sqrt()
        / scale
}
fn relaxing_gas(with_reservoir: bool) -> RelaxationState {
    let text = if with_reservoir {
        "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let medium = RelaxingMedium::new(gas(), vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    law.evaluate(&medium).unwrap()
}

#[test]
fn one_porous_body_preserves_original_solution_receivers_and_heating() {
    let (gas, law, gate) = (gas(), law(), CancelGate::new());
    let (omega, surface) = (core::f64::consts::TAU * 500.0, sphere(0.02, 0.0, 1));
    let old = PorousFluidInterface::new(
        &surface,
        &gas,
        &law,
        PoreFluid::Fluid(&gas),
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    let new = PorousFluidBodies::new(
        &[porous(&surface, &law, &gas)],
        &gas,
        omega,
        1e-4,
        Options::default(),
        &gate,
    )
    .unwrap();
    assert!(std::ptr::eq(
        new.exterior_source(),
        &gas as &dyn FluidAcoustics
    ));
    assert_eq!(
        new.porous_material_response(0),
        Some(old.material_response())
    );
    let a = old
        .scatter_field(plane(C64::new(0.7, -0.3)), 1e-4, 10000, &gate)
        .unwrap();
    let b = new
        .scatter_field(plane(C64::new(0.7, -0.3)), 1e-4, 10000, &gate)
        .unwrap();
    assert!(relative(b.total_pressure_pa(), a.total_pressure_pa()) < 3e-13);
    assert!(relative(b.total_normal_velocity_m_s(), a.total_normal_velocity_m_s()) < 3e-13);
    close(
        b.interface_work()[0].0.surface_power,
        a.interface_work.surface_power,
        3e-13,
    );
    let point = [[0.002, 0.003, 0.006]];
    let ia = a
        .interior_absorption_adaptive_at(&point, 1e-4, 10000, controls(), &gate)
        .unwrap();
    let ib = b
        .interior_absorption_adaptive_at(0, &point, 1e-4, 10000, controls(), &gate)
        .unwrap();
    close(
        ib.fields.scalar.pressure[0],
        ia.fields.scalar.pressure[0],
        3e-13,
    );
    for j in 0..3 {
        close(
            ib.fields.particle_velocity_m_s[0][j],
            ia.fields.particle_velocity_m_s[0][j],
            3e-12,
        );
    }
    assert!((ib.viscous_loss_w_m3[0] / ia.viscous_loss_w_m3[0] - 1.0).abs() < 5e-13);
    assert!((ib.compression_loss_w_m3[0] / ia.compression_loss_w_m3[0] - 1.0).abs() < 5e-13);
    let point = [[0.04, 0.01, 0.06]];
    let ea = a
        .total_pressure_velocity_adaptive_at(&point, 1e-4, 10000, controls(), &gate)
        .unwrap();
    let eb = b
        .total_pressure_velocity_adaptive_at(&point, 1e-4, 10000, controls(), &gate)
        .unwrap();
    close(eb.scalar.pressure[0], ea.scalar.pressure[0], 3e-13);
    for j in 0..3 {
        close(
            eb.particle_velocity_m_s[0][j],
            ea.particle_velocity_m_s[0][j],
            3e-12,
        );
    }
}

#[test]
fn fluid_reduction_transparent_neighbors_and_selected_receiver_regions_are_preserved() {
    let (gas, law, gate) = (gas(), law(), CancelGate::new());
    let omega = core::f64::consts::TAU * 500.0;
    let surfaces = [sphere(0.02, -0.04, 0), sphere(0.02, 0.04, 0)];
    let fluid = FluidAcousticState::try_new(0.7 * gas.density, 1.3 * gas.sound_speed).unwrap();
    let all_fluid = PorousFluidBodies::new(
        &[
            Body {
                surface: &surfaces[0],
                material: BodyMaterial::Fluid(&fluid),
            },
            Body {
                surface: &surfaces[1],
                material: BodyMaterial::Fluid(&gas),
            },
        ],
        &gas,
        omega,
        1e-4,
        Options::default(),
        &gate,
    )
    .unwrap();
    let expected = assemble(
        &surfaces[0],
        omega / gas.sound_speed,
        omega / fluid.sound_speed_m_s(),
        fluid.density_kg_m3() / gas.density,
        (fluid.density_kg_m3() - gas.density) / gas.density,
        &gate,
    )
    .unwrap();
    for (a, b) in all_fluid.bodies[0].contrast.iter().zip(expected) {
        assert_eq!(a.re.to_bits(), b.re.to_bits());
        assert_eq!(a.im.to_bits(), b.im.to_bits());
    }
    assert!(all_fluid.bodies[1].contrast.is_empty());
    let owner = PorousFluidBodies::new(
        &[
            porous(&surfaces[0], &law, &gas),
            Body {
                surface: &surfaces[1],
                material: BodyMaterial::Fluid(&gas),
            },
        ],
        &gas,
        omega,
        1e-4,
        Options::default(),
        &gate,
    )
    .unwrap();
    let single = PorousFluidInterface::new(
        &surfaces[0],
        &gas,
        &law,
        PoreFluid::Fluid(&gas),
        omega,
        Options::default(),
        &gate,
    )
    .unwrap();
    let a = owner
        .scatter_field(plane(C64::ONE), 1e-4, 10000, &gate)
        .unwrap();
    let b = single
        .scatter_field(plane(C64::ONE), 1e-4, 10000, &gate)
        .unwrap();
    assert!(
        a.outgoing_pressure_pa()[20..]
            .iter()
            .all(|&p| p == C64::ZERO)
    );
    assert!(
        a.outgoing_normal_velocity_m_s()[20..]
            .iter()
            .all(|&v| v == C64::ZERO)
    );
    assert!(relative(&a.total_pressure_pa()[..20], b.total_pressure_pa()) < 3e-13);
    let p = [[0.003, 0.0, 0.045]];
    let inside = a
        .interior_absorption_adaptive_at(1, &p, 1e-4, 40, controls(), &gate)
        .unwrap();
    let expected = b
        .total_pressure_velocity_adaptive_at(&p, 1e-4, 20, controls(), &gate)
        .unwrap();
    close(
        inside.fields.scalar.pressure[0],
        expected.scalar.pressure[0],
        3e-13,
    );
    assert_eq!(inside.viscous_loss_w_m3, [0.0]);
    assert_eq!(inside.compression_loss_w_m3, [0.0]);
    assert!(
        a.interior_pressure_velocity_adaptive_at(1, &p, 1e-4, 20, controls(), &gate)
            .is_err()
    );
    assert!(
        a.interior_pressure_velocity_adaptive_at(0, &p, 1e-4, 1000, controls(), &gate)
            .is_err()
    );
    assert!(
        a.total_pressure_velocity_adaptive_at(&p, 1e-4, 1000, controls(), &gate)
            .is_err()
    );
    let transparent = PorousFluidBodies::new(
        &[
            Body {
                surface: &surfaces[0],
                material: BodyMaterial::Fluid(&gas),
            },
            Body {
                surface: &surfaces[1],
                material: BodyMaterial::Fluid(&gas),
            },
        ],
        &gas,
        omega,
        1e-4,
        Options::default(),
        &gate,
    )
    .unwrap();
    assert!(transparent.factor.is_none());
    let zero = transparent
        .scatter_field(plane(C64::ONE), 1e-4, 1000, &gate)
        .unwrap();
    let observed = zero
        .interior_pressure_velocity_adaptive_at(1, &p, 1e-4, 1000, controls(), &gate)
        .unwrap();
    let direct = plane(C64::ONE)
        .sample_wave(
            omega / gas.sound_speed,
            omega * gas.density,
            None,
            p[0],
            &mut || false,
        )
        .unwrap();
    assert_eq!(observed.scalar.pressure[0], direct.pressure_pa);
    assert_eq!(observed.particle_velocity_m_s[0], direct.velocity_m_s);
    let internal = [Monopole {
        position_m: p[0],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        owner
            .scatter_field(IncidentField::Monopoles(&internal), 1e-4, 1000, &gate)
            .is_err()
    );
}

#[test]
fn mixed_body_sources_coherent_heat_and_mutual_work_keep_original_order() {
    let (gas, law, gate, relaxing) = (gas(), law(), CancelGate::new(), relaxing_gas(true));
    let omega = core::f64::consts::TAU * 500.0;
    let surfaces = [sphere(0.02, -0.04, 0), sphere(0.02, 0.04, 0)];
    let bodies = [
        Body {
            surface: &surfaces[0],
            material: BodyMaterial::Porous {
                law: &law,
                filling: PoreFluid::RelaxingGas(&relaxing),
            },
        },
        Body {
            surface: &surfaces[1],
            material: BodyMaterial::RelaxingGas(&relaxing),
        },
    ];
    let owner =
        PorousFluidBodies::new(&bodies, &gas, omega, 1e-4, Options::default(), &gate).unwrap();
    let reverse = PorousFluidBodies::new(
        &[bodies[1], bodies[0]],
        &gas,
        omega,
        1e-4,
        Options::default(),
        &gate,
    )
    .unwrap();
    assert_eq!(
        *owner.porous_material_response(0).unwrap(),
        law.effective_relaxing(&relaxing, omega).unwrap()
    );
    assert_eq!(
        owner.interior_dynamic_bulk_modulus_pa(1),
        Some(relaxing.caloric_response(omega).unwrap().bulk_modulus_pa)
    );
    let first = IncidentSource::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    });
    let second = IncidentSource::PlaneWave(PlaneWave {
        direction: [1.0, 0.0, 0.0],
        amplitude_pa: C64::new(0.4, 0.7),
    });
    let sources = [first, second];
    let field = owner
        .scatter_field(IncidentField::Sources(&sources), 1e-4, 10000, &gate)
        .unwrap();
    let reordered = reverse
        .scatter_field(IncidentField::Sources(&sources), 1e-4, 10000, &gate)
        .unwrap();
    let first_sources = [first];
    let second_sources = [second];
    let a = owner
        .scatter_field(IncidentField::Sources(&first_sources), 1e-4, 10000, &gate)
        .unwrap();
    let b = owner
        .scatter_field(IncidentField::Sources(&second_sources), 1e-4, 10000, &gate)
        .unwrap();
    let mut incoherent_defect = 0.0_f64;
    for body in 0..2 {
        let point = [[0.003, 0.002, if body == 0 { -0.034 } else { 0.046 }]];
        let f = field
            .interior_absorption_adaptive_at(body, &point, 1e-4, 10000, controls(), &gate)
            .unwrap();
        let r = reordered
            .interior_absorption_adaptive_at(1 - body, &point, 1e-4, 10000, controls(), &gate)
            .unwrap();
        let x = a
            .interior_absorption_adaptive_at(body, &point, 1e-4, 10000, controls(), &gate)
            .unwrap();
        let y = b
            .interior_absorption_adaptive_at(body, &point, 1e-4, 10000, controls(), &gate)
            .unwrap();
        close(
            f.fields.scalar.pressure[0],
            r.fields.scalar.pressure[0],
            3e-12,
        );
        close(
            f.fields.scalar.pressure[0],
            x.fields.scalar.pressure[0] + y.fields.scalar.pressure[0],
            3e-12,
        );
        let heat = f.viscous_loss_w_m3[0] + f.compression_loss_w_m3[0];
        assert!(heat > 0.0);
        let separated = x.viscous_loss_w_m3[0]
            + x.compression_loss_w_m3[0]
            + y.viscous_loss_w_m3[0]
            + y.compression_loss_w_m3[0];
        incoherent_defect = incoherent_defect.max((heat - separated).abs() / heat);
        if body == 1 {
            assert_eq!(f.viscous_loss_w_m3[0], 0.0);
        }
        close(
            field.interface_work()[body].0.surface_power,
            reordered.interface_work()[1 - body].0.surface_power,
            3e-12,
        );
    }
    assert!(
        incoherent_defect > 0.1,
        "fixture must expose coherent heat cross terms"
    );
    let dot = |p: &[C64], v: &[C64]| {
        p.iter()
            .zip(v)
            .zip(owner.surface.areas())
            .fold(C64::ZERO, |sum, ((&p, &v), &area)| {
                sum + (p * v.conj()).scale(0.5 * area)
            })
    };
    let full = dot(
        field.scattered_pressure_pa(),
        field.scattered_normal_velocity_m_s(),
    );
    close(full, field.scattered_work().surface_power, 3e-13);
    let local = dot(
        field.outgoing_pressure_pa(),
        field.outgoing_normal_velocity_m_s(),
    );
    assert!(
        (full - local).abs() > 1e-4 * full.abs(),
        "mutual work must differ from local-source work"
    );
    eprintln!(
        "mixed-body mutual scattered-work fraction={:e}",
        (full - local).abs() / full.abs()
    );
    let (lo, hi) = field.scattered_work_roundoff_interval();
    assert!(lo <= full.re && full.re <= hi);
    let cancellation = CancelGate::new();
    cancellation.request();
    assert!(matches!(
        field.interior_absorption_adaptive_at(
            0,
            &[[0.0, 0.0, -0.04]],
            1e-4,
            1000,
            controls(),
            &cancellation
        ),
        Err(HelmholtzError::Cancelled)
    ));
}

#[test]
fn interacting_porous_spheres_converge_to_independent_multiple_partial_waves_and_heat() {
    let (gas, law, gate) = (gas(), law(), CancelGate::new());
    let (omega, radius) = (core::f64::consts::TAU * 500.0, 0.05);
    let centers = [-1.5 * radius, 1.5 * radius];
    let response = law.effective_fluid(&gas, omega).unwrap();
    let exact = reference::Pair::new(&gas, [response, response], radius, centers, 13, true);
    let refined = reference::Pair::new(&gas, [response, response], radius, centers, 15, true);
    let independent = reference::Pair::new(&gas, [response, response], radius, centers, 13, false);
    let external = [0.0, 0.0, 3.0 * radius];
    let expected = exact.exterior_axis(external[2]);
    let check = refined.exterior_axis(external[2]);
    eprintln!(
        "independent mutual-sphere retained-order pressure discrepancy={:e}",
        (expected.pressure_pa - check.pressure_pa).abs() / check.pressure_pa.abs()
    );
    close(expected.pressure_pa, check.pressure_pa, 2e-8);
    close(expected.velocity_m_s[2], check.velocity_m_s[2], 2e-8);
    assert!(
        (independent.exterior_axis(external[2]).pressure_pa - expected.pressure_pa).abs() > 1e-4,
        "fixture must resolve mutual reflection in total pressure"
    );
    for body in 0..2 {
        let loss = exact.volume_losses(body, 48);
        let more = exact.volume_losses(body, 64);
        eprintln!(
            "mutual sphere body={body}: viscous={} W, compression={} W, inward={} W",
            loss[0], loss[1], exact.inward_power[body]
        );
        assert!(loss.iter().all(|&v| v > 0.0));
        for channel in 0..2 {
            assert!((loss[channel] - more[channel]).abs() < 2e-12 * more[channel]);
        }
        assert!(
            (loss[0] + loss[1] - exact.inward_power[body]).abs() < 3e-10 * exact.inward_power[body]
        );
    }
    let (mut previous_field, mut previous_power) = (f64::INFINITY, f64::INFINITY);
    for level in 0..=2 {
        let surfaces = centers.map(|z| sphere(radius, z, level));
        let owner = PorousFluidBodies::new(
            &[
                porous(&surfaces[0], &law, &gas),
                porous(&surfaces[1], &law, &gas),
            ],
            &gas,
            omega,
            1e-4,
            Options::default(),
            &gate,
        )
        .unwrap();
        let field = owner
            .scatter_field(plane(C64::ONE), 1e-4, 100000, &gate)
            .unwrap();
        let exterior = field
            .total_pressure_velocity_adaptive_at(&[external], 1e-4, 100000, controls(), &gate)
            .unwrap();
        let mut error = ((exterior.scalar.pressure[0] - expected.pressure_pa).abs()
            / expected.pressure_pa.abs())
        .max(
            (exterior.particle_velocity_m_s[0][2] - expected.velocity_m_s[2]).abs()
                / expected.velocity_m_s[2].abs(),
        );
        let mut power_error = 0.0_f64;
        for body in 0..2 {
            let r = 0.31 * radius;
            let inside = field
                .interior_absorption_adaptive_at(
                    body,
                    &[[0.0, 0.0, centers[body] + r]],
                    1e-4,
                    100000,
                    controls(),
                    &gate,
                )
                .unwrap();
            let target = exact.interior_axis(body, r);
            error = error.max(
                (inside.fields.scalar.pressure[0] - target.pressure_pa).abs()
                    / target.pressure_pa.abs(),
            );
            error = error.max(
                (inside.fields.particle_velocity_m_s[0][2] - target.velocity_m_s[2]).abs()
                    / target.velocity_m_s[2].abs(),
            );
            let inward = -field.interface_work()[body].0.surface_power.re;
            power_error = power_error
                .max((inward - exact.inward_power[body]).abs() / exact.inward_power[body]);
            if level == 2 {
                let rho = response.dynamic_density_kg_m3().conj();
                let bulk = response.dynamic_bulk_modulus_pa().conj();
                let viscous = 0.5 * omega * rho.im * target.velocity_m_s[2].norm_sq();
                let compression = 0.5 * omega * bulk.recip().im * target.pressure_pa.norm_sq();
                assert!((inside.viscous_loss_w_m3[0] - viscous).abs() < 0.15 * viscous);
                assert!((inside.compression_loss_w_m3[0] - compression).abs() < 0.15 * compression);
            }
        }
        eprintln!(
            "mutual porous spheres level={level}: field_error={error:e}, body_power_error={power_error:e}"
        );
        assert!(error < previous_field && power_error < previous_power);
        (previous_field, previous_power) = (error, power_error);
    }
    assert!(previous_field < 0.08 && previous_power < 0.15);
}
