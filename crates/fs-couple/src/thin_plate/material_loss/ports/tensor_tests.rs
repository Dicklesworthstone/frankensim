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
            relaxation_bending:None,polarization_detune:0.0,moving_end:false},0.0003,&state(1.0,None)).unwrap()
        .with_kelvin_voigt_bending_loss().unwrap())).to_vec();
    let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap();
    let cavity=Some(MaterialBodyCavity {cavity:HelmholtzModalCavity::new(HelmholtzCavity {
        volume_m3:0.004,neck_radius_m:0.01,neck_length_m:0.02},gas,CavityNeckLoss::Resistance(100.0)).unwrap(),
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
