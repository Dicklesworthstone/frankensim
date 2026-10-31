//! Source -> DKT tensor -> original coupled modal dynamics, not a decay filter.
use super::*;
use crate::bowed_string::runtime::reciprocal::bank::{MaterialStringSource,
    ReciprocalStringBank, StringBankOptions, StringBowLocation,
    chart::SharedMaterialPlateChartBody, plate::SharedPlateBridge};
use crate::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity;
use crate::cavity_modal::{HelmholtzModalCavity, CavityNeckLoss};
use crate::modal_acoustic_time::{ModalAcousticTimeBudget, ModalAcousticTimeModel, ModalAcousticMode};
use crate::modal_acoustic_time::ModalAcousticState;
use crate::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingConfig};
use crate::render::schedule::force::coupled::contact::friction::{FrictionDrive, ModalFrictionConfig};
use crate::render::schedule::force::coupled::render::friction::FrictionGestureConfig;
use crate::string_specimen::with_uniform_circular_material_state;
use fs_exec::Cx;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, InterpolationPolicy, UncertaintyModel};
use fs_material::state_point::{ScalarPropertyRequirement, ScalarAdmissibility,
    MaterialPropertySelection, resolve_material_state_point};
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES as KEYS;
use fs_material::gas::{GasSpec,GasState};
use fs_scenario::{PrestressedString,HelmholtzCavity,
    gesture::{GestureSchedule,GestureTarget,GestureTrack,GestureValue}};
use fs_tribo::FrictionLaw;

fn state(scale:f64, omit:Option<usize>)->ResolvedMaterialStatePoint {
    let mut values = vec![("density",Density::DIMS,600.0),("young_modulus",Pressure::DIMS,8e9),
        ("young_modulus_1",Pressure::DIMS,8e9),("young_modulus_2",Pressure::DIMS,2e9),
        ("poisson_ratio",Dims::NONE,0.2),("poisson_ratio_12",Dims::NONE,0.2),
        ("shear_modulus_12",Pressure::DIMS,0.6e9),
        (KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,DynViscosity::DIMS,10.0)];
    for (i,(key,value)) in KEYS.into_iter().zip([3e7,6e6,4e6,1e6,4e5,-3e5]).enumerate() {
        if omit != Some(i) { values.push((key,DynViscosity::DIMS,value*scale)); }
    }
    let mut claims=ClaimSet::new();let mut requirements=Vec::new();
    for (key,dims,value) in values {
        claims.insert_claim(PropertyClaim { key:PropertyKey::new(key,dims),
            value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,310.0).with("omega",0.0,100000.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"synthetic anisotropic viscous body; not measured wood".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"tensor test".into(),revision:0},claims,vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap()
        .with("omega",1000.0).unwrap(),&requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn body(scale:f64)->SharedMaterialPlateChartBody {
    let mesh=PlateMesh::rectangle(0.4,0.3,5,4);
    let regions=["left","right"].into_iter().enumerate().map(|(r,name)|PlateRegionMaterial {
        region:PlateRegion {name:name.into(),triangle_indices:mesh.tris.iter().enumerate().filter_map(|(e,t)| {
            let right=t.iter().map(|&i|mesh.nodes[i].0).sum::<f64>()/3.0>0.2;
            (right==(r==1)).then_some(e)
        }).collect()},material:state(if r==0 {scale} else {7.0*scale},None),model:PlateMaterialModel::Orthotropic12,
        thickness_constraint:PlateThicknessConstraint::FixedThickness(0.002),
        material_angle_rad:if r==0 {0.17} else {-0.23},
    }).collect();
    let mut force=vec![0.0;mesh.nodes.len()];force[14]=0.75;force[15]=0.25;
    SharedMaterialPlateChartBody {
        specimen:compile_plate_material_chart(mesh,PlateMesh::rectangle_boundary(5,4),regions).unwrap(),
        reduction:PlateChartRadiation {assembly:AssemblyOptions {pretension:0.0,support:EdgeSupport::SimplySupported},
            eigenvalue_window:(1.0,1e8),n_modes:3,damping_ratio:0.0,unit_force_weights:force},
        bridges:vec![SharedPlateBridge {station_fraction:0.8,stiffness_n_m:3000.0,damping_n_s_m:0.1,rest_extension_m:0.0};2],
        max_plate_dofs:90,max_plate_elements:40,modal_budget:ModalAcousticTimeBudget::audible_reference(),listener_m:1.0,
    }
}
fn limits()->PlateViscosityLimits { PlateViscosityLimits {max_modes:8,max_projection_terms:1_000_000,reference_resistance_n_s_m:1.0} }
fn tensor(b:&SharedMaterialPlateChartBody)->PlateViscosity {
    PlateViscosity::from_chart(&b.specimen,KEYS,&CancelGate::new_clock_free()).unwrap()
}
fn config()->ModalCouplingConfig {
    ModalCouplingConfig {max_modes:32,max_connections:8,max_setup_terms:16384,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:100.0,
        solve_relative_tolerance:1e-10,energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8}
}
fn scoped<R>(f:impl FnOnce(&Cx<'_>,&CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|a|f(&Cx::new(&gate,a,fs_exec::StreamKey {seed:7,kernel_id:83,tile:0,iteration:0},
        fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic),&gate))
}
fn build(cx:&Cx<'_>,gate:&CancelGate,b:SharedMaterialPlateChartBody,selected:bool)->ReciprocalStringBank {
    let sources=[4,3].map(|n_modes|MaterialStringSource::Homogeneous(with_uniform_circular_material_state(
        PrestressedString {length_m:0.5,tension_n:20.0,lin_density_kg_m:1.0,axial_stiffness_n:1.0,width_m:1.0,
            n_modes,damping_ratio:0.0,rayleigh:None,bending_stiffness_n_m2:0.0,kelvin_voigt_bending:None,
            relaxation_bending:None,polarization_detune:0.0,moving_end:false,
            tension_loss_factor:None,bending_loss_factor:None},0.0003,&state(1.0,None)).unwrap()
        .with_kelvin_voigt_bending_loss().unwrap())).to_vec();
    let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap();
    let cavity=Some(MaterialBodyCavity {cavity:HelmholtzModalCavity::new(HelmholtzCavity {
        volume_m3:0.004,neck_radius_m:0.01,neck_length_m:0.02,inner_end:fs_scenario::NeckEnd::Flanged,
        outer_end:fs_scenario::NeckEnd::Flanged},gas,CavityNeckLoss::Resistance(100.0)).unwrap(),
        modal_budget:ModalAcousticTimeBudget::audible_reference(),listener_m:1.5});
    let options=StringBankOptions {sample_rate_hz:48000,bow:StringBowLocation {string:0,station_fraction:0.2},
        string_budget:ModalAcousticTimeBudget::audible_reference(),network:config(),
        friction:ModalFrictionConfig {max_iterations:96,maximum_force_n:100.0,maximum_slip_speed_m_s:10.0,
            force_absolute_tolerance_n:1e-10,force_relative_tolerance:1e-9,sticking_tolerance_m:1e-12}};
    let law=FrictionLaw::Coulomb {static_mu:0.8,kinetic_mu:0.4};
    if selected { let v=tensor(&b);ReciprocalStringBank::from_viscous_material_plate_chart(cx,sources,gas,b,v,
        limits(),cavity,law,options,gate).unwrap() }
    else { ReciprocalStringBank::from_material_plate_chart(cx,sources,gas,b,cavity,law,options,gate).unwrap() }
}
fn drive()->FrictionDrive {FrictionDrive {speed_m_s:0.01,normal_force_n:0.01}}
fn close(a:f64,b:f64) {assert!((a-b).abs()<=2e-9*(a.abs()+b.abs()).max(1e-24),"{a} != {b}");}

// Independent continuous spectral solution in energy coordinates y=[W*q,v]:
// y'=A*y+G*f, A=[0,W;-W,-C]. For distinct poles lambda_k, its projectors are
// P_k=product_{j!=k}(A-lambda_j*I)/(lambda_k-lambda_j). The held-force
// equilibrium is [W^-1*f,0], so delta_y=sum expm1(lambda_k*h)*P_k*(y-y_eq).
// Complex QR supplies only the poles; no production exponential/ZOH is used.
fn tensor_exact_increment(
    omega: &[f64],
    c: &[f64],
    dt: f64,
) -> impl Fn(&[ModalAcousticState], &[f64]) -> Vec<f64> {
    let omega = omega.to_vec();
    let n = omega.len();
    let d = 2 * n;
    let mut a = vec![C64::ZERO; d * d];
    for i in 0..n {
        a[(2 * i) * d + 2 * i + 1] = C64::from_re(omega[i]);
        a[(2 * i + 1) * d + 2 * i] = C64::from_re(-omega[i]);
        for j in 0..n {
            a[(2 * i + 1) * d + 2 * j + 1] = C64::from_re(-c[n * i + j]);
        }
    }
    let poles = fs_la::eigen_complex::eig(&a, d).unwrap();
    let mut projectors = Vec::with_capacity(d);
    let mut increments = Vec::with_capacity(d);
    for (k, &lambda) in poles.iter().enumerate() {
        let mut p = vec![C64::ZERO; d * d];
        for i in 0..d {
            p[d * i + i] = C64::ONE;
        }
        for (j, &mu) in poles.iter().enumerate() {
            if j == k {
                continue;
            }
            assert!(
                (lambda - mu).abs() > 1e-8 * lambda.abs().max(mu.abs()),
                "distinct physical poles required"
            );
            let mut next = vec![C64::ZERO; d * d];
            for row in 0..d {
                for col in 0..d {
                    let mut value = C64::ZERO;
                    for inner in 0..d {
                        let entry = a[d * row + inner] - if row == inner { mu } else { C64::ZERO };
                        value = value + entry * p[d * inner + col];
                    }
                    next[d * row + col] = value / (lambda - mu);
                }
            }
            p = next;
        }
        // Stable complex expm1 avoids losing the small forced displacement
        // by subtracting an equilibrium from two nearly identical endpoints.
        let z = lambda.scale(dt);
        let cosine = fs_math::det::cos(z.im);
        let real = fs_math::det::expm1(z.re) * cosine - 2.0 * fs_math::det::sin(0.5 * z.im).powi(2);
        let imaginary = fs_math::det::exp(z.re) * fs_math::det::sin(z.im);
        increments.push(C64::new(real, imaginary));
        projectors.push(p);
    }
    move |old: &[ModalAcousticState], force: &[f64]| {
        let mut initial = vec![0.0; d];
        for i in 0..n {
            initial[2 * i] = omega[i] * old[i].displacement_m_sqrt_kg - force[i] / omega[i];
            initial[2 * i + 1] = old[i].velocity_m_sqrt_kg_per_s;
        }
        let mut delta = vec![0.0; d];
        for row in 0..d {
            let mut value = C64::ZERO;
            let mut scale = 0.0;
            for (p, &increment) in projectors.iter().zip(&increments) {
                let projected =
                    (0..d).fold(C64::ZERO, |sum, j| sum + p[d * row + j].scale(initial[j]));
                let term = increment * projected;
                value = value + term;
                scale += term.abs();
            }
            assert!(
                value.im.abs() <= 2e-9 * scale.max(1e-24),
                "conjugate pole reconstruction"
            );
            delta[row] = value.re / if row % 2 == 0 { omega[row / 2] } else { 1.0 };
        }
        delta
    }
}

#[test]
fn full_dkt_tensor_drives_cross_mode_reactions_and_physical_dissipation() {
    let gate = CancelGate::new_clock_free();
    let b = body(1.0);
    let v = tensor(&b);
    let (image, ports) = b
        .specimen
        .certified_viscous_radiators_with_force_ports(
            &b.reduction,
            &[b.reduction.unit_force_weights.clone()],
            1000,
            &v,
            limits(),
            &gate,
        )
        .unwrap();
    for (i, r) in image.radiators().iter().enumerate() {
        close(ports[0][i], r.drive_participation / r.mass_kg.sqrt());
    }
    let n = image.damping().mode_count();
    assert_eq!(n, 3);
    let c = image.damping().matrix_per_s();
    assert!(
        (0..n * n).any(|i| i / n != i % n && c[i].abs() > 1e-8),
        "cross terms must not vanish"
    );
    let links = image.damping().connections(0);
    assert_eq!(links.len(), n);
    for i in 0..n {
        for j in 0..n {
            close(
                c[n * i + j],
                links
                    .iter()
                    .map(|l| l.damping_n_s_m * l.left.shapes[i] * l.left.shapes[j])
                    .sum(),
            );
        }
    }
    // The selected full tensor owns all loss; no intrinsic modal loss is added.
    assert!(image.radiators().iter().all(|r| r.zeta == 0.0));
    let omega: Vec<_> = image.radiators().iter().map(|r| r.omega).collect();
    for factor in [1_u32, 2, 4] {
        let rate = 48_000 * factor;
        let dt = f64::from(rate).recip();
        let exact = tensor_exact_increment(&omega, c, dt);
        let model = ModalAcousticTimeModel::try_new(
            rate,
            image
                .radiators()
                .iter()
                .map(|r| ModalAcousticMode {
                    angular_frequency_rad_s: r.omega,
                    damping_ratio: r.zeta,
                    pressure_per_modal_velocity: C64::ZERO,
                })
                .collect(),
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap();
        let mut system =
            CoupledModalSystem::new(vec![model], links.clone(), config(), &gate).unwrap();
        let mut total_work = 0.0;
        let mut total_loss = 0.0;
        for k in 0..512 * factor {
            let old = system.components()[0].states().to_vec();
            let force = [if k < 32 * factor { 0.001 } else { 0.0 }, 0.0, 0.0];
            let delta = exact(&old, &force);
            let expected_force: Vec<_> = (0..n)
                .map(|i| {
                    -(0..n)
                        .map(|j| {
                            c[n * i + j] * (old[j].velocity_m_sqrt_kg_per_s + delta[2 * j + 1])
                        })
                        .sum::<f64>()
                })
                .collect();
            if k == 0 {
                // Plant a true trapezoidal step of this same tensor pencil:
                // (I+h*C/2+h^2*W^2/4)*v1=h*f from q0=v0=0.
                let mut midpoint: Vec<_> = c.iter().map(|v| 0.5 * dt * v).collect();
                for i in 0..n {
                    midpoint[n * i + i] += 1.0 + 0.25 * dt * dt * omega[i] * omega[i];
                }
                let mut velocity: Vec<_> = force.iter().map(|f| dt * f).collect();
                fs_la::factor::cholesky(&midpoint, n)
                    .unwrap()
                    .solve(&mut velocity);
                assert!(
                    (0..n).any(|i| {
                        let planted = -(0..n).map(|j| c[n * i + j] * velocity[j]).sum::<f64>();
                        (planted - expected_force[i]).abs()
                            > 2e-9 * (planted.abs() + expected_force[i].abs()).max(1e-24)
                    }),
                    "trapezoid must fail the unchanged exact reaction tolerance at {rate} Hz"
                );
            }
            let f = system.step(&force).unwrap().clone();
            for i in 0..n {
                let actual = links
                    .iter()
                    .zip(&f.connection_forces_n)
                    .map(|(l, f)| l.left.shapes[i] * f)
                    .sum();
                close(expected_force[i], actual);
            }
            // Exact integral(v^T*C*v)=f.delta(q)-delta(H). Both increments
            // here come from the independent spectral trajectory above, not the
            // production endpoint or ledger. Use factored energy differences.
            let heat = (0..n)
                .map(|i| {
                    force[i] * delta[2 * i]
                        - omega[i].powi(2)
                            * (old[i].displacement_m_sqrt_kg + 0.5 * delta[2 * i])
                            * delta[2 * i]
                        - (old[i].velocity_m_sqrt_kg_per_s + 0.5 * delta[2 * i + 1])
                            * delta[2 * i + 1]
                })
                .sum();
            close(heat, f.connection_dissipation_j);
            assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
            total_work += f.external_work_j;
            total_loss += f.component_dissipation_j + f.connection_dissipation_j;
        }
        assert!((system.total_energy_j().unwrap() + total_loss - total_work).abs() < 1e-14);
        assert!(
            system.components()[0].states()[1]
                .velocity_m_sqrt_kg_per_s
                .abs()
                > 1e-15,
            "a force in mode one must produce the cross-mode viscous reaction"
        );
    }
}

#[test]
fn zero_source_is_exact_legacy_and_only_viscosity_changes_bow_body_and_sound() {
    scoped(|cx,gate| {
        let mut old=build(cx,gate,body(0.0),false);let mut zero=build(cx,gate,body(0.0),true);
        assert_eq!(zero.source().material_chart_viscosity().unwrap().damping().channel_count(),0);
        let external=vec![0.0;old.system().mode_count()];
        for _ in 0..128 {assert_eq!(old.step(&external,drive()).unwrap(),zero.step(&external,drive()).unwrap());}
        let mut a=build(cx,gate,body(1.0),true);let mut b=build(cx,gate,body(4.0),true);
        assert_eq!(a.cards()[0].fundamental_hz(),b.cards()[0].fundamental_hz());
        assert_eq!(a.system().components()[2].modes(),b.system().components()[2].modes(),"no retuning or diagonal decay replacement");
        let mut motion=false;let mut sound=false;
        for _ in 0..1024 {
            let fa=a.step(&external,drive()).unwrap().clone();let fb=b.step(&external,drive()).unwrap().clone();
            assert!(fa.energy_residual_j.abs()<=fa.energy_tolerance_j);assert!(fb.energy_residual_j.abs()<=fb.energy_tolerance_j);
            sound |= (fa.observer_pressure_pa-fb.observer_pressure_pa).abs()>1e-12;
            motion |= (a.system().components()[0].states()[0].displacement_m_sqrt_kg
                -b.system().components()[0].states()[0].displacement_m_sqrt_kg).abs()>1e-16;
        }
        assert!(motion && sound,"material damping must react on the source, not only output");
        let source=b.source().material_chart_viscosity().unwrap().source();
        assert_eq!(source.regions()[1].law().source(),&b.source().material_chart().unwrap().specimen.regions()[1].input.material);
        let states=b.system().components().iter().map(|m|m.states().to_vec()).collect::<Vec<_>>();
        let stopped=CancelGate::new_clock_free();stopped.request();
        assert!(b.step_under_gate(&external,drive(),&stopped).is_err());
        assert_eq!(states,b.system().components().iter().map(|m|m.states().to_vec()).collect::<Vec<_>>());
    });
}

#[test]
fn scheduled_pressure_retains_full_tensor_sources_and_callback_equivalence() {
    scoped(|cx,gate| {
        let schedule=GestureSchedule::try_new(700,vec![GestureTrack {id:"bow".into(),
            target:GestureTarget::BowStroke {string:0},initial:GestureValue::Bow {
                velocity_m_per_s:0.01,normal_force_n:0.01,station:0.2},events:vec![]}]).unwrap();
        let config=FrictionGestureConfig {sample_rate_hz:48000,samples:257,max_block:257,max_work:1000,max_events:16,
            station_fraction:0.2};
        let mut outputs=Vec::new();
        for size in [1,37,257] {
            let mut p=build(cx,gate,body(1.0),true).into_scheduled(&schedule,"bow",config).unwrap();
            let c=p.source().material_chart_viscosity().unwrap().damping().matrix_per_s().to_vec();
            let mut out=vec![0.0;257];
            for block in out.chunks_mut(size) {p.renderer_mut().block(block).unwrap();}
            assert_eq!(p.source().material_chart_viscosity().unwrap().damping().matrix_per_s(),c);
            assert!(out.iter().any(|v|v.abs()>1e-12));outputs.push(out);
        }
        assert_eq!(outputs[0],outputs[1]);assert_eq!(outputs[0],outputs[2]);
    });
}

#[test]
fn cold_rebind_requires_new_source_tensors_and_never_discards_an_existing_loss() {
    let gate=CancelGate::new_clock_free();let mut b=body(1.0);let v=tensor(&b);
    b.specimen=b.specimen.with_region_material("right",&state(4.0,None)).unwrap();
    assert!(b.specimen.certified_viscous_radiators(&b.reduction,&v,limits(),&gate).is_err());
    let replacement=v.rebind(&b.specimen,&gate).unwrap();
    assert_eq!(replacement.regions()[0].law().source(),v.regions()[0].law().source());
    assert_ne!(replacement.regions()[1].law().source(),v.regions()[1].law().source());
    let missing=b.specimen.with_region_material("right",&state(4.0,Some(5))).unwrap();
    assert!(replacement.rebind(&missing,&gate).is_err());
    let mut low=limits();low.max_projection_terms=1;
    assert!(b.specimen.certified_viscous_radiators(&b.reduction,&replacement,low,&gate).is_err());
    b.reduction.damping_ratio=0.001;
    assert!(b.specimen.certified_viscous_radiators(&b.reduction,&replacement,limits(),&gate).is_err());
    b.reduction.damping_ratio=0.0;gate.request();
    assert!(b.specimen.certified_viscous_radiators(&b.reduction,&replacement,limits(),&gate).is_err());
}

// IRP M0.2-W5: certified modal damping on the plate's own projected matrix.
// No new projection: `matrix_per_s()` is the reference, and fs-modal decides
// which modes keep one decay rate and which need a full block.

/// Isotropic state whose six reduced viscosities are `beta` [s] times its own
/// plane-stress stiffness, so the projected damping is exactly `beta*K`.
fn isotropic_state(beta: f64) -> ResolvedMaterialStatePoint {
    isotropic_state_valid_to(beta, 100000.0)
}

/// The same state, admitted up to `omega_max` [rad/s].
fn isotropic_state_valid_to(beta: f64, omega_max: f64) -> ResolvedMaterialStatePoint {
    let (e, nu) = (8e9, 0.2);
    let q11 = e / (1.0 - nu * nu);
    let g = e / (2.0 * (1.0 + nu));
    let mut values = vec![("density",Density::DIMS,600.0),("young_modulus",Pressure::DIMS,e),
        ("young_modulus_1",Pressure::DIMS,e),("young_modulus_2",Pressure::DIMS,e),
        ("poisson_ratio",Dims::NONE,nu),("poisson_ratio_12",Dims::NONE,nu),
        ("shear_modulus_12",Pressure::DIMS,g),
        (KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,DynViscosity::DIMS,10.0)];
    for (key,value) in KEYS.into_iter().zip([q11,q11,g,nu*q11,0.0,0.0]) {
        values.push((key,DynViscosity::DIMS,beta*value));
    }
    let mut claims=ClaimSet::new();let mut requirements=Vec::new();
    for (key,dims,value) in values {
        claims.insert_claim(PropertyClaim { key:PropertyKey::new(key,dims),
            value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,310.0).with("omega",0.0,omega_max),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            observations:vec![],provenance:Provenance {source:"synthetic isotropic viscous body; not a measured material".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"certified damping test".into(),revision:0},claims,vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap()
        .with("omega",1000.0).unwrap(),&requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}

/// A simply supported 0.3 m x 0.3*aspect m plate on a 6x6 right-triangle
/// mesh, which is mirror-symmetric about y = x when square. The lower-left
/// quadrant and the rest carry their own viscosity times `beta`. A quadrant,
/// unlike a left/right split, couples (2,1) with (1,2): their y-shapes are
/// orthogonal over any region spanning all of y.
fn two_region_plate(beta_corner: f64, beta_rest: f64, aspect: f64) -> (Vec<f64>, Vec<f64>) {
    let (chart, mut reduction, viscosity) = two_region_chart(beta_corner, beta_rest, aspect, 6, 100000.0);
    reduction.n_modes = 4;
    reduction.eigenvalue_window = (1.0, 1e8);
    let gate = CancelGate::new_clock_free();
    let image = chart.certified_viscous_radiators(&reduction, &viscosity, limits(), &gate).unwrap();
    let omega = image.radiators().iter().map(|r| r.omega).collect();
    (omega, image.damping().matrix_per_s().to_vec())
}

/// The same quadrant/rest chart on a `cells` x `cells` mesh, its materials
/// admitted up to `omega_max` [rad/s], asking for 64 modes in (1, omega_max^2].
fn two_region_chart(beta_corner: f64, beta_rest: f64, aspect: f64, cells: usize, omega_max: f64)
    -> (ResolvedPlateChart, PlateChartRadiation, PlateViscosity)
{
    let (nx, ny) = (cells, cells);
    let (a, b) = (0.3, 0.3 * aspect);
    let mesh = PlateMesh::rectangle(a, b, nx, ny);
    let regions = [("corner", beta_corner), ("rest", beta_rest)].into_iter().enumerate()
        .map(|(r, (name, beta))| PlateRegionMaterial {
            region: PlateRegion { name: name.into(), triangle_indices: mesh.tris.iter().enumerate()
                .filter_map(|(e, t)| {
                    let x = t.iter().map(|&i| mesh.nodes[i].0).sum::<f64>() / 3.0;
                    let y = t.iter().map(|&i| mesh.nodes[i].1).sum::<f64>() / 3.0;
                    let corner = x < 0.5 * a && y < 0.5 * b;
                    (corner == (r == 0)).then_some(e)
                }).collect() },
            material: isotropic_state_valid_to(beta, omega_max), model: PlateMaterialModel::Orthotropic12,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.002),
            material_angle_rad: 0.0,
        }).collect();
    // Node (2,2) on the 6x6 mesh: off every nodal line of the four lowest modes.
    let mut force = vec![0.0; mesh.nodes.len()];
    let at = cells / 3;
    force[at * (nx + 1) + at] = 1.0;
    let chart = compile_plate_material_chart(mesh, PlateMesh::rectangle_boundary(nx, ny), regions).unwrap();
    let reduction = PlateChartRadiation {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::SimplySupported },
        eigenvalue_window: (1.0, omega_max * omega_max), n_modes: 64, damping_ratio: 0.0, unit_force_weights: force,
    };
    let viscosity = PlateViscosity::from_chart(&chart, KEYS, &CancelGate::new_clock_free()).unwrap();
    (chart, reduction, viscosity)
}

/// One pole of `lambda^2 I + lambda C + Omega^2`, refined from a dense seed by
/// inverse iteration and the root of `x^T Q(lambda) x` (Q is complex
/// symmetric). Raw `quadratic_eigenvalues` can sit ~1e-8 relative off.
fn refined_pole(omega: &[f64], c: &[f64], seed: C64) -> C64 {
    let n = omega.len();
    let (mut lambda, mut x) = (seed, vec![C64::ONE; n]);
    for _ in 0..8 {
        let mut q = vec![C64::ZERO; n * n];
        for i in 0..n {
            for j in 0..n { q[n * i + j] = lambda * C64::from_re(c[n * i + j]); }
            q[n * i + i] = q[n * i + i] + lambda * lambda + C64::from_re(omega[i] * omega[i]);
        }
        let Ok(lu) = fs_la::eigen_complex::lu_complex(&q, n) else { break };
        lu.solve(&mut x);
        let scale = x.iter().map(|v| v.abs()).fold(0.0, f64::max);
        for v in &mut x { *v = v.scale(1.0 / scale); }
        let (mut a2, mut a1, mut a0) = (C64::ZERO, C64::ZERO, C64::ZERO);
        for i in 0..n {
            a2 = a2 + x[i] * x[i];
            a0 = a0 + x[i] * x[i] * C64::from_re(omega[i] * omega[i]);
            for j in 0..n { a1 = a1 + x[i] * C64::from_re(c[n * i + j]) * x[j]; }
        }
        let root = (a1 * a1 - C64::from_re(4.0) * a2 * a0).sqrt();
        let two = C64::from_re(2.0) * a2;
        let (r1, r2) = ((C64::ZERO - a1 + root) / two, (C64::ZERO - a1 - root) / two);
        let next = if (r1 - lambda).abs() <= (r2 - lambda).abs() { r1 } else { r2 };
        let step = (next - lambda).abs();
        lambda = next;
        if !(step > 4.0 * f64::EPSILON * lambda.abs()) { break; }
    }
    lambda
}

fn damping_budget() -> fs_modal::ModalDampingBudget {
    fs_modal::ModalDampingBudget { max_rel_decay_error: 1e-3, max_freq_error_cents: 0.01,
        max_cluster_size: fs_modal::MAX_DAMPING_CLUSTER }
}

#[test]
fn g1_proportional_plate_is_all_diagonal_and_matches_energy_quotient() {
    // zeta ~ 5e-3..2e-2: the certified radius cannot fall below the closed
    // form's own rounding (2*eps*nu), which must sit well under 1e-12*|Re|.
    let beta = 2e-5;
    let (omega, c) = two_region_plate(beta, beta, 1.07);
    let n = omega.len();
    let result = fs_modal::certify_modal_damping(&omega, &c, &damping_budget()).unwrap();
    assert_eq!(result.blocks.len(), n, "every mode a singleton");
    let (mut worst_radius, mut worst_quotient) = (0.0_f64, 0.0_f64);
    for (i, block) in result.blocks.iter().enumerate() {
        let fs_modal::ModalDampingBlock::Diagonal { mode, zeta, eigenvalue, radius } = block else {
            panic!("proportional damping must stay diagonal: {block:?}");
        };
        assert_eq!(*mode, i);
        assert!(*radius <= 1e-12 * eigenvalue.re.abs(), "mode {i}: radius {radius:e}");
        worst_radius = worst_radius.max(radius / eigenvalue.re.abs());
        // The energy quotient of the projected matrix, bit for bit.
        assert_eq!(zeta.to_bits(), (c[n * i + i] / (2.0 * omega[i])).to_bits());
        // Independent oracle: C_FE = beta*K gives zeta = beta*omega/2.
        let quotient = (zeta / (0.5 * beta * omega[i]) - 1.0).abs();
        worst_quotient = worst_quotient.max(quotient);
        assert!(quotient <= 1e-6, "mode {i}: zeta off beta*omega/2 by {quotient:e}");
    }
    println!("{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g1_proportional_plate_is_all_diagonal_and_matches_energy_quotient\",\"tier\":\"G1\",\"inputs\":{{\"fixture\":\"two-region isotropic plate, C_FE = beta*K\",\"n_modes\":{n},\"beta_s\":{beta:e}}},\"provenance\":\"named fixture\",\"observed\":{{\"max_radius_over_decay\":{worst_radius:e},\"max_zeta_vs_beta_omega_half\":{worst_quotient:e}}},\"tolerance\":{{\"radius_over_decay\":1e-12,\"zeta\":1e-6}},\"verdict\":\"pass\"}}");
}

#[test]
fn g1_clustered_pair_triggers_full_block_fallback() {
    // Viscosity 50:1 between the quadrant and the rest; a 0.2 % aspect
    // detuning splits the (2,1)/(1,2) pair by under 0.5 % (1.004 measured
    // 0.61 % on this mesh).
    let beta_rest = 1.3e-7;
    let (omega, c) = two_region_plate(50.0 * beta_rest, beta_rest, 1.002);
    let n = omega.len();
    let split = (omega[2] - omega[1]) / omega[1];
    assert!(split > 0.0 && split < 5e-3, "pair split {split:e}");
    let result = fs_modal::certify_modal_damping(&omega, &c, &damping_budget()).unwrap();
    let mut omega_sq = vec![0.0; n * n];
    let mut identity = vec![0.0; n * n];
    for i in 0..n {
        omega_sq[n * i + i] = omega[i] * omega[i];
        identity[n * i + i] = 1.0;
    }
    let raw: Vec<C64> = fs_modal::quadratic_eigenvalues(&identity, &c, &omega_sq, n).unwrap()
        .into_iter().filter(|p| p.im > 0.0).collect();
    let truth: Vec<C64> = raw.iter().map(|&p| refined_pole(&omega, &c, p)).collect();
    let raw_correction = raw.iter().zip(&truth).map(|(r, t)| (*r - *t).abs() / t.abs()).fold(0.0, f64::max);
    assert!(raw_correction <= 1e-6, "refinement wandered: {raw_correction:e}");
    let nearest = |z: C64| truth.iter().copied()
        .min_by(|a, b| (*a - z).abs().total_cmp(&(*b - z).abs())).unwrap();
    let nearest_raw = |z: C64| raw.iter().copied()
        .min_by(|a, b| (*a - z).abs().total_cmp(&(*b - z).abs())).unwrap();
    let mut pair_block = None;
    let (mut worst_decay, mut worst_cents) = (0.0_f64, 0.0_f64);
    let (mut raw_decay, mut raw_cents) = (0.0_f64, 0.0_f64);
    for block in &result.blocks {
        if let fs_modal::ModalDampingBlock::Full { modes, eigenvalues, radii, .. } = block {
            if modes.contains(&1) && modes.contains(&2) { pair_block = Some(modes.clone()); }
            for (&mu, &r) in eigenvalues.iter().zip(radii) {
                let t = nearest(mu);
                assert!((t - mu).abs() <= r + 8.0 * f64::EPSILON * mu.abs(),
                    "block pole off by {:e}, radius {r:e}", (t - mu).abs());
                worst_decay = worst_decay.max((t.re - mu.re).abs() / t.re.abs());
                worst_cents = worst_cents.max(1200.0 / core::f64::consts::LN_2 * ((t.im - mu.im) / t.im).abs());
                // The unrefined `quadratic_eigenvalues` poles as truth, too.
                let t = nearest_raw(mu);
                raw_decay = raw_decay.max((t.re - mu.re).abs() / t.re.abs());
                raw_cents = raw_cents.max(1200.0 / core::f64::consts::LN_2 * ((t.im - mu.im) / t.im).abs());
            }
        }
    }
    let pair_block = pair_block.expect("the near-degenerate pair must be one full block");
    assert!(worst_decay <= 1e-3 && worst_cents <= 0.01, "block poles {worst_decay:e} / {worst_cents:e} c");
    assert!(raw_decay <= 1e-3 && raw_cents <= 0.01, "block poles vs raw dense truth {raw_decay:e} / {raw_cents:e} c");
    // Forcing the diagonal on the same matrix.
    let loose = fs_modal::ModalDampingBudget { max_rel_decay_error: 1e9, max_freq_error_cents: 1e9,
        max_cluster_size: fs_modal::MAX_DAMPING_CLUSTER };
    let forced = fs_modal::certify_modal_damping(&omega, &c, &loose).unwrap();
    let (mut forced_decay, mut forced_cents) = (0.0_f64, 0.0_f64);
    for block in &forced.blocks {
        if let fs_modal::ModalDampingBlock::Diagonal { mode, eigenvalue, .. } = block
            && (*mode == 1 || *mode == 2)
        {
            let t = nearest(*eigenvalue);
            forced_decay = forced_decay.max((t.re - eigenvalue.re).abs() / t.re.abs());
            forced_cents = forced_cents.max(1200.0 / core::f64::consts::LN_2 * ((t.im - eigenvalue.im) / t.im).abs());
        }
    }
    let shape: Vec<String> = result.blocks.iter().map(|b| match b {
        fs_modal::ModalDampingBlock::Diagonal { mode, .. } => format!("{mode}"),
        fs_modal::ModalDampingBlock::Full { modes, .. } => format!("{modes:?}"),
    }).collect();
    println!("{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g1_clustered_pair_triggers_full_block_fallback\",\"tier\":\"G1\",\"inputs\":{{\"fixture\":\"quadrant/rest isotropic plate, viscosity 50:1, aspect 1.002\",\"n_modes\":{n},\"pair_split\":{split:e}}},\"provenance\":\"named fixture\",\"observed\":{{\"blocks\":\"{}\",\"pair_block\":{pair_block:?},\"block_rel_decay_error\":{worst_decay:e},\"block_freq_error_cents\":{worst_cents:e},\"block_rel_decay_error_vs_raw_dense\":{raw_decay:e},\"block_freq_error_cents_vs_raw_dense\":{raw_cents:e},\"forced_diagonal_rel_decay_error\":{forced_decay:e},\"forced_diagonal_freq_error_cents\":{forced_cents:e},\"raw_dense_truth_max_rel_correction\":{raw_correction:e}}},\"tolerance\":{{\"decay\":1e-3,\"cents\":0.01}},\"verdict\":\"pass\"}}", shape.join(" "));
    assert!(forced_decay > 1e-3 || forced_cents > 0.01,
        "forcing the diagonal must break the budget: {forced_decay:e} / {forced_cents:e} c");
}

#[test]
fn g0_bessel_far_coupling_bound_dominates_true_norm() {
    // The windowed plate projection forms entries only inside the window and
    // bounds the rest of each row by Bessel's inequality on the mass-norm
    // residual. Truth is the dense `matrix_per_s()` of the same modes.
    // - "partial": the 6x6 clustered chart at its 64-mode ceiling (123 DOFs),
    //   so the bound also covers modes that were never retained.
    // - "complete": the same quadrant/rest chart on a 4x4 mesh, all 59 DOFs
    //   retained. Bessel is then Parseval, and with the window covering
    //   everything the bound must vanish to 1e-10 of each column norm.
    // Viscosities are the clustered test's scaled down 130x, so the stiffest
    // rotary modes stay underdamped.
    let gate = CancelGate::new_clock_free();
    let limits = PlateViscosityLimits { max_modes: 64, max_projection_terms: 100_000_000,
        reference_resistance_n_s_m: 1.0 };
    let loose = fs_modal::ModalDampingBudget { max_rel_decay_error: 1e9, max_freq_error_cents: 1e9,
        max_cluster_size: fs_modal::MAX_DAMPING_CLUSTER };
    let mut rows = Vec::new();
    let mut complete_floor = None;
    for (label, cells) in [("partial", 6), ("complete", 4)] {
        let (chart, reduction, viscosity) = two_region_chart(5e-8, 1e-9, 1.002, cells, 1e8);
        let free = chart.chart().assemble(&[], &reduction.assembly).unwrap().free;
        let dense = chart.certified_viscous_radiators(&reduction, &viscosity, limits, &gate).unwrap();
        let omega: Vec<f64> = dense.radiators().iter().map(|r| r.omega).collect();
        let n = omega.len();
        assert_eq!(n, free.min(64), "{label}: retained modes");
        let c = dense.damping().matrix_per_s();
        let column: Vec<f64> = (0..n).map(|j| (0..n).map(|k| c[n * j + k] * c[n * j + k]).sum::<f64>().sqrt()).collect();
        let span = omega[n - 1] - omega[0];
        for window in [0.02 * omega[0], 0.5 * omega[0], 0.1 * span, 2.0 * span] {
            let (radiators, windowed) = chart.certified_viscous_radiators_windowed(&reduction, &viscosity,
                limits, window, &gate).unwrap();
            assert!(radiators.iter().zip(&omega).all(|(r, w)| r.omega.to_bits() == w.to_bits()),
                "the same single eigensolve and modes");
            let neighbours = fs_modal::damping_window_neighbours(&omega, window).unwrap();
            let everything = neighbours.iter().all(|row| row.len() == n - 1);
            let (mut margin, mut far_rel, mut column_excess, mut entry_error) = (f64::INFINITY, 0.0_f64, 0.0_f64, 0.0_f64);
            for j in 0..n {
                for &k in neighbours[j].iter().chain(std::iter::once(&j)) {
                    let formed = windowed.entry(j, k).expect("formed inside the window");
                    entry_error = entry_error.max((formed - c[n * j + k]).abs() / column[j]);
                }
                assert!(neighbours[j].iter().all(|&k| windowed.entry(k, j) == windowed.entry(j, k)), "symmetric");
                let true_far = (0..n).filter(|&k| k != j && neighbours[j].binary_search(&k).is_err())
                    .map(|k| c[n * j + k] * c[n * j + k]).sum::<f64>().sqrt();
                let bound = windowed.far_row_norm_bound()[j];
                assert!(bound >= true_far, "{label}, window {window:e}: mode {j} far bound {bound:e} < true {true_far:e}");
                let column_bound = windowed.column_norm_bound()[j];
                assert!(column_bound >= column[j], "{label}: mode {j} column bound {column_bound:e} < true {:e}", column[j]);
                margin = margin.min((bound - true_far) / column[j]);
                far_rel = far_rel.max(bound / column[j]);
                column_excess = column_excess.max(column_bound / column[j] - 1.0);
            }
            assert!(entry_error <= 1e-10, "{label}: formed entries off the dense projection by {entry_error:e}");
            // The certifier consumes exactly this bound.
            let certified = match windowed.certify(&loose) {
                Ok(result) => {
                    assert_eq!(result.far_coupling.far_row_norm, windowed.far_row_norm_bound());
                    true
                }
                Err(fs_modal::ModalError::DampingBlockSpansWindow { .. }) if !everything => false,
                Err(e) => panic!("{label}, window {window:e}: {e}"),
            };
            // The same column bound through a difference of squares, for the record.
            let near = |i: usize, k: usize| windowed.entry(i, k).unwrap_or(f64::NAN);
            let squares = fs_modal::certify_modal_damping_windowed(&omega, &near, windowed.column_norm_bound(),
                window, &loose).map(|r| (0..n).map(|j| r.far_coupling.far_row_norm[j] / column[j]).fold(0.0, f64::max));
            let squares = squares.map_or_else(|e| format!("\"{e}\""), |v| format!("{v:e}"));
            if everything && label == "complete" { complete_floor = Some((far_rel, column_excess)); }
            rows.push(format!("{{\"chart\":\"{label}\",\"cells\":{cells},\"free_dofs\":{free},\"n_modes\":{n},\"window_rad_s\":{window:e},\"full_window\":{everything},\"certified\":{certified},\"basis_defect\":{:e},\"min_rel_margin\":{margin:e},\"max_rel_far_bound\":{far_rel:e},\"max_rel_column_excess\":{column_excess:e},\"max_rel_entry_error\":{entry_error:e},\"squares_path_max_rel_far_bound\":{squares}}}", windowed.basis_defect()));
        }
    }
    let (floor, excess) = complete_floor.expect("the complete chart's full window ran");
    println!("{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g0_bessel_far_coupling_bound_dominates_true_norm\",\"tier\":\"G0\",\"inputs\":{{\"fixture\":\"quadrant/rest isotropic plate, viscosity 50:1, aspect 1.002\",\"beta_s\":[5e-8,1e-9]}},\"provenance\":\"named fixture\",\"observed\":[{}],\"tolerance\":{{\"complete_full_window_rel_far_bound\":1e-10,\"complete_rel_column_excess\":1e-10}},\"verdict\":\"{}\"}}",
        rows.join(","), if floor <= 1e-10 && excess <= 1e-10 { "pass" } else { "fail" });
    assert!(floor <= 1e-10, "complete space, full window: far bound {floor:e} of the column norm");
    assert!(excess <= 1e-10, "complete space: column bound exceeds the true norm by {excess:e}");
}

#[test]
fn g1_certified_plate_zeta_is_the_w3_evaluator() {
    // One owner: the certified Diagonal blocks of the full tensor image report
    // DiagonalPlateViscosity's zeta, the same bits a consumer gets by selecting
    // scalar regional loss (M0.2-W3) on the same chart and single eigensolve.
    // The routed disc grows by its centre shift, and it must still enclose the
    // refined dense truth of the full modal pencil.
    let mut rows = Vec::new();
    for (label, corner, rest, aspect, budget) in [
        ("proportional", 2e-5, 2e-5, 1.07, damping_budget()),
        ("clustered", 6.5e-6, 1.3e-7, 1.002, fs_modal::ModalDampingBudget { max_rel_decay_error: 1e9,
            max_freq_error_cents: 1e9, max_cluster_size: fs_modal::MAX_DAMPING_CLUSTER }),
    ] {
        let (chart, mut reduction, viscosity) = two_region_chart(corner, rest, aspect, 6, 100000.0);
        reduction.n_modes = 4;
        reduction.eigenvalue_window = (1.0, 1e8);
        let gate = CancelGate::new_clock_free();
        let image = chart.certified_viscous_radiators(&reduction, &viscosity, limits(), &gate).unwrap();
        let omega: Vec<f64> = image.radiators().iter().map(|r| r.omega).collect();
        let c = image.damping().matrix_per_s();
        let n = omega.len();
        let routed = image.certified_modal_damping(&budget).unwrap();
        let plain = fs_modal::certify_modal_damping(&omega, c, &budget).unwrap();
        let (w3, _) = chart.clone().with_kelvin_voigt_bending_loss().unwrap()
            .certified_radiators_with_force_ports(&reduction, &[], 0).unwrap();
        let mut omega_sq = vec![0.0; n * n];
        let mut identity = vec![0.0; n * n];
        for i in 0..n {
            omega_sq[n * i + i] = omega[i] * omega[i];
            identity[n * i + i] = 1.0;
        }
        let truth: Vec<C64> = fs_modal::quadratic_eigenvalues(&identity, c, &omega_sq, n).unwrap()
            .into_iter().filter(|p| p.im > 0.0).map(|p| refined_pole(&omega, c, p)).collect();
        let mut diagonal = 0;
        for (r, p) in routed.blocks.iter().zip(&plain.blocks) {
            if let fs_modal::ModalDampingBlock::Diagonal { mode, zeta, eigenvalue, radius } = r {
                diagonal += 1;
                assert_eq!(w3[*mode].omega.to_bits(), omega[*mode].to_bits(), "{label}: the same eigensolve");
                assert_eq!(zeta.to_bits(), w3[*mode].zeta.to_bits(), "{label}: mode {mode} zeta is W3's");
                let distance = truth.iter().map(|t| (*t - *eigenvalue).abs()).fold(f64::INFINITY, f64::min);
                assert!(distance <= *radius, "{label}: mode {mode} truth {distance:e} outside {radius:e}");
                let fs_modal::ModalDampingBlock::Diagonal { radius: base, .. } = p else {
                    panic!("{label}: routing changed the block structure");
                };
                rows.push(format!("{{\"chart\":\"{label}\",\"mode\":{mode},\"radius_over_decay\":{:e},\"unrouted_radius_over_decay\":{:e},\"truth_distance_over_radius\":{:e}}}",
                    radius / eigenvalue.re.abs(), base / eigenvalue.re.abs(), distance / radius));
            } else {
                assert_eq!(r, p, "{label}: full blocks are untouched");
            }
        }
        assert_eq!(diagonal, n, "{label}: every mode certified diagonal under this budget");
    }
    println!("{{\"suite\":\"irp-s0\",\"bead\":\"M0.2-W5\",\"test\":\"g1_certified_plate_zeta_is_the_w3_evaluator\",\"tier\":\"G1\",\"inputs\":{{\"fixture\":\"quadrant/rest isotropic plate: proportional and clustered\"}},\"provenance\":\"named fixture\",\"observed\":[{}],\"tolerance\":{{\"zeta\":\"bitwise to M0.2-W3\",\"enclosure\":\"refined dense truth inside the routed radius\"}},\"verdict\":\"pass\"}}",
        rows.join(","));
}
