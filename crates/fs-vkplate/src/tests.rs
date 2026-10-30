use super::*;

#[test]
fn shell_compatibility_retains_motion_below_the_rest_height_ulp() {
    // Exact dyadic rest heights: at the origin the curvatures are (1,1,0),
    // while the (1,1) stencil has mixed rest curvature -1/4. E*t = h = 1.
    let rest = |r: f64| if r < 0.5 { 0.0 } else if r < 1.1 { 0.5 } else { 1.0 };
    let mut body = Body::plate("small shell motion", &PlateSpec {
        radius: 3.5, hole: 0.0, grid: 1.0, edge: Edge::Free,
        thickness: &|_| 1.0 / 128.0, rest_height: &rest, young: 128.0,
        poisson: 0.3, density: 1000.0, tension: 0.0, bending: true,
        nonlinear: Nonlinear::VonKarman, t60: NO_LOSS, support: None, muffle: None,
    }).unwrap();
    assert_eq!(body.vk.as_ref().unwrap().eh, 1.0);
    let node = |x, y| body.xy.iter().position(|&p| p == [x, y]).unwrap();
    let origin = node(0.0, 0.0);
    let diagonal = node(1.0, 1.0);
    let core = |center| body.vk.as_ref().unwrap().nb.iter()
        .position(|s| s[0] == Some(center)).unwrap();
    let center_core = core(origin);
    let mixed_core = core(diagonal);
    let east = body.vk.as_ref().unwrap().nb[center_core][1].unwrap();
    let north_east = body.vk.as_ref().unwrap().nb[mixed_core][5].unwrap();
    let delta = 2.0_f64.powi(-55);
    let check_tangent = |body: &mut Body| {
        let potential = body.nonlinear_gradient();
        assert!(potential > 0.0);
        let work: f64 = body.w.iter().zip(&body.g).map(|(w, g)| w * g).sum();
        // Euler's identity for the quadratic rest-tangent energy. Higher
        // orders are smaller by O(delta); missing rest terms in g fail this.
        assert!((work - 2.0 * potential).abs() <= 1e-12 * (2.0 * potential),
            "the shell gradient must retain its linear rest tangent");
    };
    for displacement in [delta, -delta] {
        // The former nodal w+z0 evaluation loses this motion before it even
        // reaches the curvature/product subtraction. This is the negative.
        assert_eq!(0.5 + displacement, 0.5);
        body.w.fill(0.0);
        body.w[east] = displacement;
        check_tangent(&mut body);
        assert_eq!(body.vk.as_ref().unwrap().q[center_core], -displacement);
        body.w.fill(0.0);
        body.w[north_east] = displacement;
        check_tangent(&mut body);
        // w_xy = displacement/4, z_xy = -1/4. The quadratic term is below
        // half an ULP of the exact linear term -displacement/8.
        assert_eq!(body.vk.as_ref().unwrap().q[mixed_core], -displacement / 8.0);
    }
    body.w.fill(0.0);
    assert_eq!(body.nonlinear_gradient(), 0.0);
    assert!(body.vk.as_ref().unwrap().q.iter().all(|&q| q == 0.0));
}

#[test]
fn live_airy_source_system_matches_its_actual_solved_stress() {
    assert!(Body::point("no Airy system",1.0).airy_system().is_none());
    let mut live = live::Live::with_hihat_pedal(7.0,0.0).unwrap();
    live.strike(0.85*0.1778,0.0,2.5);
    let mut observed_nonzero = false;
    for _ in 0..8 {
        live.render(&mut [0.0;128],1.0);
        for b in &live.instrument().bodies {
            let Some((a,q,f)) = b.airy_system() else { continue; };
            assert_eq!(a.n,q.len()); assert_eq!(a.n,f.len());
            let mut af = vec![0.0;a.n]; a.mul_add(f,1.0,&mut af);
            let residual = af.iter().zip(q).map(|(a,b)|(a-b).powi(2)).sum::<f64>().sqrt();
            let qnorm = q.iter().map(|x|x*x).sum::<f64>().sqrt();
            let scale = qnorm + a.val.iter().map(|x|x*x).sum::<f64>().sqrt()*f.iter().map(|x|x*x).sum::<f64>().sqrt();
            assert!(residual.is_finite() && residual <= 1e-11*scale.max(1e-300), "Airy backward error {residual:e}/{scale:e}");
            observed_nonzero |= qnorm > 0.0;
        }
    }
    assert!(observed_nonzero,"exercise the actual struck source, not just its zero state");
}

#[test]
fn g0_affine_mount_energy_and_admission() {
    let mut body=Body::point("attachment",1.0);
    body.n=4;
    let angular=[7.0,2.0,9.0]; let damping=[0.4,0.1,0.8];
    for scale in [0.01,0.2,1.0] {
        body.xy=vec![[-scale,-2.0*scale],[scale,-2.0*scale],[-scale,2.0*scale],[scale,2.0*scale]];
        let nodes=vec![(0,2.0),(1,2.0),(2,2.0),(3,2.0)];
        let links=Link::affine_mount(0,&body,&nodes,Some((1,0)),[37.0,2.0],angular,damping).unwrap();
        let form=|w:&[f64],anchor:f64,damp:bool| links.iter().map(|l| {
            let wb=match l.b {Some((0,j))=>w[j],_=>anchor};
            (if damp {l.c}else{l.k})*(w[l.a.1]-wb).powi(2)
        }).sum::<f64>();
        // Independent orthogonal fields on a rectangle: translation, two tilts,
        // and the alternating non-affine saddle. The anchor translates too.
        let (t,tx,ty,q,anchor)=(0.003,0.2,-0.1,0.002,0.07);
        let w:Vec<_>=body.xy.iter().map(|[x,y]|anchor+t+tx*x+ty*y+q*x*y/(2.0*scale*scale)).collect();
        for (damp,ax,a) in [(false,37.0,angular),(true,2.0,damping)] {
            let expected=ax*(t*t+q*q)+a[0]*tx*tx+2.0*a[1]*tx*ty+a[2]*ty*ty;
            assert!((form(&w,anchor,damp)-expected).abs()<1e-10);
            for seed in 0..32 {
                let field:Vec<_>=(0..4).map(|i|((seed*4+i) as f64*1.731).sin()).collect();
                assert!(form(&field,0.3,damp)>=-1e-10,"assembled mount must be passive");
            }
        }
        let mut uneven=nodes.clone(); uneven[0].1=7.0;
        let center: [f64;2]=std::array::from_fn(|d| uneven.iter().map(|&(i,w)|w*body.xy[i][d]/13.0).sum());
        let links=Link::affine_mount(0,&body,&uneven,None,[37.0,2.0],angular,damping).unwrap();
        let w:Vec<_>=body.xy.iter().map(|[x,y]|t+tx*(x-center[0])+ty*(y-center[1])).collect();
        let work:f64=links.iter().map(|l| l.k*(w[l.a.1]-l.b.map_or(0.0,|(_,i)|w[i])).powi(2)).sum();
        assert!((work-(37.0*t*t+angular[0]*tx*tx+2.0*angular[1]*tx*ty+angular[2]*ty*ty)).abs()<1e-10);
    }
    let valid=vec![(0,1.0),(1,1.0),(2,1.0),(3,1.0)];
    for nodes in [vec![],vec![(0,1.0),(0,1.0),(2,1.0)],vec![(0,-1.0),(1,1.0),(2,1.0)],vec![(0,1.0),(1,f64::NAN),(2,1.0)],vec![(0,1.0),(1,1.0),(4,1.0)]] {
        assert!(Link::affine_mount(0,&body,&nodes,None,[1.0,1.0],angular,damping).is_err());
    }
    for bad in [[1.0,2.0,1.0],[-1.0,0.0,1.0],[1.0,f64::NAN,1.0]] {
        assert!(Link::affine_mount(0,&body,&valid,None,[1.0,1.0],bad,damping).is_err());
        assert!(Link::affine_mount(0,&body,&valid,None,[1.0,1.0],angular,bad).is_err());
    }
    assert!(Link::affine_mount(0,&body,&valid,Some((0,0)),[1.0,1.0],angular,damping).is_err());
    assert!(Link::affine_mount(0,&body,&valid,None,[-1.0,1.0],angular,damping).is_err());
    body.xy=vec![[0.0,0.0],[1.0,0.0],[2.0,0.0],[3.0,0.0]];
    assert!(Link::affine_mount(0,&body,&valid,None,[1.0,1.0],angular,damping).is_err());
}

#[test]
fn g3_streaming_filter_matches_offline_phase_and_clones() {
    for factor in [1,2,3,4,6,8] {
        let input:Vec<_>=(0..512*factor).map(|i| (0.31*i as f64).sin()+0.2*(1.13*i as f64).cos()).collect();
        let reference=decimate(&input,factor);
        let mut filter=StreamingDecimator::new(factor);
        let mut output=vec![];
        for (i,x) in input.iter().enumerate() {
            filter.push(*x);
            if i==73 {
                let mut copied=filter.clone();
                copied.push(0.125);
                let mut control=filter.clone(); control.push(0.125);
                assert_eq!(copied.value().to_bits(),control.value().to_bits());
            }
            if i%factor==0 { output.push(filter.value()); }
        }
        let delay=if factor==1 {0}else{16};
        for (i,x) in output.iter().enumerate().skip(delay) {
            assert_eq!(x.to_bits(),reference[i-delay].to_bits(),"factor {factor}, frame {i}");
        }
    }
}

#[test]
fn g3_live_internal_clock_matches_offline_and_arbitrary_output_blocks() {
    let mut a=crate::live::Live::new("snare-off",8.0).unwrap();
    assert!(a.ins.rate>=96_000.0);
    assert_eq!(a.filter_latency_frames(),16);
    a.strike(0.3*a.radius,0.0,2.0);
    let mut b=a.clone();
    let mut raw=a.ins.clone();
    let factor=(a.ins.rate/kit::OUT_RATE) as usize;
    let mut source=vec![0.0;512*factor]; let mut sample=[0.0];
    for x in &mut source { raw.step(&mut sample); *x=sample[0]; }
    let reference=decimate(&source,factor);
    let mut x=vec![0.0_f32;512];a.render(&mut x,1.0);
    assert!(x.iter().any(|x| *x!=0.0));
    for i in 16..512 { assert_eq!(x[i].to_bits(),(reference[i-16] as f32).to_bits()); }
    let mut y=vec![0.0_f32;512];
    b.render(&mut y[..193],1.0);
    let mut copied=b.clone();
    b.render(&mut y[193..194],1.0);
    b.render(&mut y[194..],1.0);
    let mut tail=vec![0.0;319];copied.render(&mut tail,1.0);
    assert_eq!(tail,y[193..]);
    assert_eq!(x,y);
    assert_eq!(a.ins.step_index,raw.step_index);
    assert_eq!(a.kinetic.to_bits(),b.kinetic.to_bits());
    assert_eq!(a.asleep(),b.asleep());
}


#[test]
fn an_axis_floor_reflects_bit_for_bit_as_the_floor_height_did() {
    // the snare and toms keep their bytes (ruling 8307): floor_z(h) images [x, y, 2 h - z] exactly
    for (h, p) in [(-0.65, [0.35, 0.25, 0.9]), (-0.75, [0.35, 0.25, 0.9]), (-0.5, [0.4, -0.2, 1.1]), (-0.55, [0.0, 0.0, 0.3])] {
        assert_eq!(floor_z(h).image(p).map(f64::to_bits), [p[0], p[1], 2.0 * h - p[2]].map(f64::to_bits));
    }
    let p = [0.0, 1.221, 0.1];
    assert_eq!(floor_y(-0.241).image(p).map(f64::to_bits), [p[0], 2.0 * -0.241 - p[1], p[2]].map(f64::to_bits));
}

#[test]
fn g3_microphone_source_height_matches_rigid_translation() {
    for facing in [1.0, -1.0] {
        let mut ins = Instrument::new(48_000.0);
        for raised in [false, true] {
            let mut b = Body::plate("source", &PlateSpec {
                radius: 0.06, hole: 0.0, grid: 0.015, edge: Edge::Free,
                thickness: &|_| 1e-3, rest_height: &|_| if raised { 0.02 } else { 0.0 },
                young: 1e11, poisson: 0.33, density: 8600.0, tension: 0.0,
                bending: true, nonlinear: Nonlinear::None, t60: NO_LOSS, support: None, muffle: None,
            }).unwrap();
            b.facing = facing;
            if !raised { b.origin[2] = facing * 0.02; }
            ins.add(b);
        }
        let mut mic = Mic::new([0.4, 0.2, 1.1]);
        mic.floor = Some(floor_z(-0.5));
        ins.mics.push(mic);
        ins.prepare();
        let taps = &ins.mics[0].taps;
        let a: Vec<_> = taps.iter().filter(|t| t.0 == 0).map(|t| (t.1,t.2,t.3.to_bits(),t.4.to_bits())).collect();
        let b: Vec<_> = taps.iter().filter(|t| t.0 == 1).map(|t| (t.1,t.2,t.3.to_bits(),t.4.to_bits())).collect();
        assert!(!a.is_empty());
        assert_eq!(a,b,"supplied source height must equal the same physical translation");
    }
}

#[test]
fn g0_hihat_contacts_preserve_opposing_geometry_and_requested_opening() {
    let law = Law::Power { k: 2e9, alpha: 1.5, mu: 0.08 };
    for scale in [1.0,4.0] {
        let closed = kit::hihat_grid_mu(0.0,0.0,&[],law,0.8,scale,0.0).unwrap();
        let open = kit::hihat_grid_mu(0.002,0.0,&[],law,0.8,scale,0.0).unwrap();
        let (top,bottom) = (&closed.instrument.bodies[0],&closed.instrument.bodies[1]);
        assert_eq!(top.xy,bottom.xy);
        assert!(top.rest_z.iter().any(|z| *z > 0.01));
        for (a,b) in top.rest_z.iter().zip(&bottom.rest_z) { assert_eq!(*a,-*b); }
        assert!((top.origin[2]-bottom.origin[2]-0.5*(0.85e-3+1.1e-3)).abs()<1e-15);
        let mut widest=0.0_f64;
        for (c,o) in closed.instrument.contacts.iter().zip(&open.instrument.contacts) {
            let (_,ib,_) = c.a[0]; let (_,it,_) = c.b[0];
            let [x,y] = top.xy[it];
            let waviness=7.5e-5*(1.0-(x*x-y*y)/(x*x+y*y));
            // Independent face clearance from physical mass, area and density.
            let tt=top.mass[it]/(kit::BRONZE.2*top.h*top.h);
            let tb=bottom.mass[ib]/(kit::BRONZE.2*bottom.h*bottom.h);
            let clearance=top.origin[2]+top.rest_z[it]-0.5*tt
                -(bottom.origin[2]+bottom.rest_z[ib]+0.5*tb);
            assert!((c.gap-clearance-waviness).abs()<1e-14);
            assert!((o.gap-c.gap-0.002).abs()<1e-14,"requested opening was capped");
            assert!(c.gap>=0.0);
            widest=widest.max(c.gap);
        }
        assert!(widest>0.003,"curved bow must not become a constant rim gap");
    }
    for gap in [-0.001,f64::NAN,f64::INFINITY] { assert!(kit::hihat(gap,0.0,&[]).is_err()); }
}

#[test]
fn g3_hihat_refinement_preserves_joint_load_and_rim_traction() {
    let law = Law::Power { k: 2e9, alpha: 1.5, mu: 0.08 };
    let compression=1e-5_f64; let speed=0.25;
    let pressure=2e9/5.5e-3_f64.powi(2)*compression.powf(1.5)*(1.0+0.08*speed);
    for scale in [0.9,1.0,1.5,4.0,8.0] {
        let ins=kit::hihat_grid_mu(0.0,0.0,&[],law,0.8,scale,0.0).unwrap().instrument;
        let rod=ins.bodies.iter().position(|b| b.name=="rod").unwrap();
        let clutch:Vec<_>=ins.links.iter().filter(|l| l.a.0==0&&l.b==Some((rod,0))).collect();
        assert!(clutch.len()>=4);
        assert!((clutch.iter().map(|l| l.k*1e-5).sum::<f64>()-80.0).abs()<1e-10);
        assert!((clutch.iter().map(|l| l.c*speed).sum::<f64>()-20.0).abs()<1e-10);
        // Generalized tilt torque must not grow with the clutch grid radius.
        let angular_k=2.0*PI*1.65e5*(0.5_f64*(0.0065+0.02)).powi(3);
        let cup:Vec<_>=ins.links.iter().filter(|l|l.a.0==1&&l.b.is_none()).collect();
        let area=PI*(0.02_f64.powi(2)-0.0065_f64.powi(2));
        assert!((cup.iter().map(|l|l.k).sum::<f64>()/(9.0e7*area)-1.0).abs()<1e-12);
        assert!((cup.iter().map(|l|l.c).sum::<f64>()/(2.7e3*area)-1.0).abs()<1e-12);
        for (body,rot_k) in [(0,angular_k),(1,0.5*angular_k)] {
            for axis in 0..2 {
                let mut torque=0.0; let mut net=0.0;
                for l in ins.links.iter().filter(|l|l.a.0==body) {
                    let xa=ins.bodies[body].xy[l.a.1][axis];
                    let xb=l.b.filter(|&(b,_)|b==body).map_or(0.0,|(_,i)|ins.bodies[body].xy[i][axis]);
                    let f=-l.k*0.01*(xa-xb);
                    net+=f; torque+=f*xa;
                    if l.b.is_some_and(|(b,_)|b==body) {net-=f; torque-=f*xb;}
                }
                assert!(net.abs()<1e-8);
                assert!((torque+rot_k*0.01).abs()<1e-9,"scale {scale}, body {body}, axis {axis}, torque {torque}");
            }
        }
        let force: f64=ins.contacts.iter().map(|c| {
            let Law::Power{k,alpha,mu}=c.law;
            k*compression.powf(alpha)*(1.0+mu*speed)
        }).sum();
        let area=ins.contacts.len() as f64*ins.bodies[0].h.powi(2);
        assert!((force/(pressure*area)-1.0).abs()<1e-12);
    }
}

#[test]
fn g0_shell_admission_includes_curvature_without_changing_live_state() {
    // Isolate the shell tangent from bending: a product preset may have a
    // bending-dominated highest mode even when its rest curvature matters.
    let mut b=Body::plate("shell admission", &PlateSpec {
        radius: 0.1, hole: 0.02, grid: 0.01, edge: Edge::Free,
        thickness: &|_| 0.001, rest_height: &|r| 0.2*r*r,
        young: 1e11, poisson: 0.3, density: 8000.0, tension: 0.0,
        bending: false, nonlinear: Nonlinear::VonKarman,
        t60: NO_LOSS, support: None, muffle: None,
    }).unwrap();
    let linear=b.max_rate2(); let shell=b.max_rate2_with_shell();
    assert_eq!(linear,0.0);
    assert!(shell.is_finite() && shell>0.0,"curved shell must have a tangent");
    let mut instrument=Instrument::new((shell/(4.0*0.5)).sqrt());
    instrument.add(b);
    assert!(instrument.stability().0<=0.6);
    let old_rate=((linear+shell)/(8.0*0.6)).sqrt();
    // Negative control: the former linear-only admission accepts this rate.
    assert!(linear/(4.0*old_rate.powi(2))<=0.6);
    instrument.rate=old_rate;
    instrument.bodies[0].f[0]+=0.5;
    let before=instrument.bodies[0].clone();
    assert!(instrument.stability().0>0.6);
    assert_eq!(instrument.bodies[0].f,before.f);
    assert_eq!(instrument.bodies[0].w,before.w);
    assert_eq!(instrument.bodies[0].wp,before.wp);
    assert_eq!(instrument.bodies[0].lw,before.lw);
}


fn crossings_hz(x: &[f64], rate: f64) -> f64 {
    let mean = x.iter().sum::<f64>() / x.len() as f64;
    let mut idx = vec![];
    for i in 1..x.len() {
        let (a, b) = (x[i - 1] - mean, x[i] - mean);
        if a < 0.0 && b >= 0.0 {
            idx.push(i as f64 - 1.0 + a / (a - b));
        }
    }
    assert!(idx.len() >= 3, "too few cycles");
    (idx.len() - 1) as f64 / ((idx[idx.len() - 1] - idx[0]) / rate)
}

const NO_LOSS: [(f64, f64); 2] = [(100.0, 1e12), (1000.0, 1e12)];

fn point_contact_plate() -> Body {
    Body::plate("contact", &PlateSpec {
        radius: 0.1, hole: 0.0, grid: 0.02, edge: Edge::Free,
        thickness: &|_| 1e-3, rest_height: &|_| 0.0, young: 1e11,
        poisson: 0.33, density: 8600.0, tension: 0.0, bending: true,
        nonlinear: Nonlinear::None, t60: NO_LOSS, support: None, muffle: None,
    }).unwrap()
}

#[test]
fn point_contact_retains_its_position_and_reconstructs_an_affine_surface() {
    let b = point_contact_plate();
    let (x, y) = (1.35 * b.h, -0.4 * b.h);
    let footprint = b.footprint(x, y, 0.0);
    assert_eq!(footprint.len(), 4);
    assert!(footprint.iter().all(|&(_, w)| w > 0.0 && w < 1.0));
    let sum: f64 = footprint.iter().map(|&(_, w)| w).sum();
    let px: f64 = footprint.iter().map(|&(i, w)| w * b.xy[i][0]).sum();
    let py: f64 = footprint.iter().map(|&(i, w)| w * b.xy[i][1]).sum();
    let surface: f64 = footprint.iter().map(|&(i, w)| w * (2.0 + 3.0 * b.xy[i][0] - 4.0 * b.xy[i][1])).sum();
    assert!((sum - 1.0).abs() < 1e-14);
    assert!((px - x).abs() < 1e-14 && (py - y).abs() < 1e-14);
    assert!((surface - (2.0 + 3.0 * x - 4.0 * y)).abs() < 1e-14);
}

#[test]
fn point_contact_moves_continuously_across_the_nearest_node_boundary() {
    let b = point_contact_plate();
    let centres: Vec<f64> = [0.49, 0.51].iter().map(|x| {
        b.footprint(x * b.h, 0.0, 0.0).iter().map(|&(i, w)| w * b.xy[i][0]).sum()
    }).collect();
    assert!((centres[1] - centres[0] - 0.02 * b.h).abs() < 1e-14);
    let node = b.node_at(b.h, -b.h).unwrap();
    assert_eq!(b.footprint(b.h, -b.h, 0.0), vec![(node, 1.0)]);
}

#[test]
fn interpolated_point_force_and_surface_velocity_do_the_same_work() {
    let b = point_contact_plate();
    let (x, y, force) = (-0.3 * b.h, 1.7 * b.h, 23.0);
    let footprint = b.footprint(x, y, 0.0);
    let nodal_work: f64 = footprint.iter().map(|&(i, w)| {
        force * w * (0.3 + 2.0 * b.xy[i][0] - b.xy[i][1])
    }).sum();
    let point_work = force * (0.3 + 2.0 * x - y);
    assert!((nodal_work - point_work).abs() < 1e-13);
}

#[test]
fn finite_contact_faces_and_incomplete_grid_cells_keep_their_existing_rules() {
    let b = point_contact_plate();
    let face = b.footprint(0.0, 0.0, 1.1 * b.h);
    assert_eq!(face.len(), 5);
    assert!(face.iter().all(|&(_, w)| w == 0.2));
    // The free-edge plate retains the union of complete curvature stencils;
    // its axial grid boundary is one node inside this geometric radius.
    let nearest = b.node_at(4.0 * b.h, 0.0).unwrap();
    assert_eq!(b.footprint(6.0 * b.h, 0.0, 0.0), vec![(nearest, 1.0)]);
}

#[test]
fn membrane_fundamental_matches_bessel_zero() {
    let (r, t, th) = (0.15, 2000.0, 2.54e-4);
    let thick = |_r: f64| th;
    let z = |_r: f64| 0.0;
    let mut b = Body::plate("m", &PlateSpec {
        radius: r, hole: 0.0, grid: 2.0e-3, edge: Edge::Clamped, thickness: &thick, rest_height: &z,
        young: 4e9, poisson: 0.38, density: 1390.0, tension: t, bending: false, nonlinear: Nonlinear::None,
        t60: NO_LOSS, support: None, muffle: None,
    }).unwrap();
    for i in 0..b.n {
        let p = b.xy[i];
        let rr = (p[0] * p[0] + p[1] * p[1]).sqrt() / r;
        // close to J0(2.405 r/R): mostly the (0,1) mode
        let v = 1e-4 * (1.0 - rr * rr) * (1.0 - 0.3 * rr * rr);
        b.w[i] = v;
        b.wp[i] = v;
    }
    b.sync();
    let mut ins = Instrument::new(96_000.0);
    let c = ins.add(b);
    let centre = ins.bodies[c].node_at(0.0, 0.0).unwrap();
    let mut out = [];
    let mut trace = vec![];
    for _ in 0..9600 {
        ins.step(&mut out);
        trace.push(ins.bodies[c].w[centre]);
    }
    let f = crossings_hz(&trace, 96_000.0);
    let exact = 2.404826 * (t / (1390.0 * th)).sqrt() / (2.0 * PI * r);
    assert!((f / exact - 1.0).abs() < 0.015, "f01 {f:.2} Hz vs {exact:.2} Hz");
}

#[test]
fn free_plate_lowest_elastic_mode_matches_disc_theory() {
    // free circular plate, nodal diameters 2, no nodal circle: lambda^2 = 5.253 (nu = 0.33)
    let (r, th, e, nu, rho) = (0.2, 1.0e-3, 1.0e11, 0.33, 8600.0);
    let thick = |_r: f64| th;
    let z = |_r: f64| 0.0;
    let mut b = Body::plate("p", &PlateSpec {
        radius: r, hole: 0.0, grid: 8.0e-3, edge: Edge::Free, thickness: &thick, rest_height: &z,
        young: e, poisson: nu, density: rho, tension: 0.0, bending: true, nonlinear: Nonlinear::None,
        t60: NO_LOSS, support: None, muffle: None,
    }).unwrap();
    for i in 0..b.n {
        let p = b.xy[i];
        let v = 1e-5 * (p[0] * p[0] - p[1] * p[1]) / (r * r); // r^2 cos(2 theta)
        b.w[i] = v;
        b.wp[i] = v;
    }
    let probe = b.node_at(0.15, 0.0).unwrap();
    b.sync();
    let mut ins = Instrument::new(96_000.0);
    let c = ins.add(b);
    let (ratio, _) = ins.stability();
    assert!(ratio < 1.0, "unstable: {ratio}");
    let mut out = [];
    let mut trace = vec![];
    for _ in 0..19200 {
        ins.step(&mut out);
        trace.push(ins.bodies[c].w[probe]);
    }
    let f = crossings_hz(&trace, 96_000.0);
    let d = e * th.powi(3) / (12.0 * (1.0 - nu * nu));
    let exact = 5.253 / (2.0 * PI * r * r) * (d / (rho * th)).sqrt();
    assert!((f / exact - 1.0).abs() < 0.04, "(2,0) {f:.2} Hz vs {exact:.2} Hz");
}

#[test]
fn shallow_shell_struck_hard_stays_bounded_and_spreads_energy_upward() {
    let spec = kit::crash();
    let mut p = kit::cymbal(&spec, &[(0.0, 8.0, 0.95, 0.008)]).unwrap();
    let out = kit::render(&mut p, 0.12, |_| {});
    let x = &out[0];
    assert!(x.iter().all(|v| v.is_finite()));
    let peak = x.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
    assert!(peak > 1e-3 && peak < 1e3, "peak {peak}");
}

#[test]
fn a_new_stroke_retires_the_sticks_earlier_contact() {
    // a hi-hat adds one stick contact per stroke; the stale one must stop pushing on the shell
    let mut p = kit::hihat(0.002, 0.0, &[(0.002, 2.0), (0.012, 2.0)]).unwrap();
    let first: Vec<usize> = (0..p.instrument.contacts.len()).filter(|&i| p.instrument.contacts[i].b.iter().any(|&(b, _, _)| p.instrument.bodies[b].name == "stick")).collect();
    assert_eq!(first.len(), 2, "one stick contact per stroke");
    let _ = kit::render(&mut p, 0.02, |_| {});
    let c = &p.instrument.contacts;
    assert!(!c[first[0]].active, "the first stroke's contact is still active after the second launch");
    assert!(c[first[1]].active && c[first[1]].peak_force > 0.0, "the second stroke never landed");
}

/// `image` (the export for the GPU renderer) only reads: exporting before the first step leaves every
/// output sample unchanged, and the exported Airy solve is symmetric positive definite, as A^-1 must be.
#[test]
fn exporting_an_instrument_changes_nothing_and_exports_a_spd_airy_inverse() {
    let run = |export: bool| {
        let mut p = live::Live::new("ride", 2.5).unwrap();
        p.strike(0.6 * p.radius, 0.0, 2.5);
        if export {
            let img = p.ins.image();
            assert_eq!(img.bodies.len(), p.ins.bodies.len());
            assert!(img.bodies.iter().all(|b| b.cf.len() == b.n && b.lap.ptr.len() == b.n + 1));
            assert!(img.mics.iter().all(|m| m.ring_len.is_power_of_two() && !m.taps.is_empty()));
        }
        let mut blk = vec![0.0f32; 960];
        p.render(&mut blk, 1.0);
        (blk.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), p)
    };
    let (plain, _) = run(false);
    let (exported, mut p) = run(true);
    assert_eq!(plain, exported, "exporting changed the output");
    let img = p.ins.image();
    let (bi, nc) = img.bodies.iter().enumerate().find_map(|(i, b)| b.vk.as_ref().map(|v| (i, v.f.len()))).unwrap();
    assert!(nc > 10, "no von Karman core found");
    let shell = &p.ins.bodies[bi];
    assert!(!shell.airy_solve(&mut vec![0.0; nc + 1]), "a wrong-length solve must refuse");
    let col = |j: usize| {
        let mut x = vec![0.0; nc];
        x[j] = 1.0;
        assert!(shell.airy_solve(&mut x));
        x
    };
    let (a, b) = (col(3), col(nc / 2));
    assert!(a[3] > 0.0 && b[nc / 2] > 0.0, "diagonal of A^-1 not positive");
    let rel = (a[nc / 2] - b[3]).abs() / a[nc / 2].abs().max(b[3].abs());
    assert!(rel < 1e-9, "A^-1 not symmetric: {rel:e}");
    // positive definite along a few rough vectors too (x' A^-1 x > 0)
    for seed in 1..4u64 {
        let x: Vec<f64> = (0..nc as u64).map(|i| (((i + 1) * 2_654_435_761 * seed) % 1000) as f64 / 500.0 - 1.0).collect();
        let mut y = x.clone();
        assert!(shell.airy_solve(&mut y));
        let q: f64 = x.iter().zip(&y).map(|(a, b)| a * b).sum();
        assert!(q > 0.0, "x' A^-1 x = {q:e} for seed {seed}");
    }
}

#[test]
fn g3_shell_reanchoring_is_continuous_at_zero_previous_potential() {
    let mut spec = kit::splash();
    spec.grid = 0.025;
    let shell = kit::cymbal(&spec, &[]).unwrap().instrument.bodies.remove(0);
    let mut targets = vec![];
    for fraction in [0.0, 1e-6] {
        let mut b = shell.clone();
        let q: Vec<_> = (0..b.n).map(|i| 1e-5 * (0.31 * i as f64).sin()).collect();
        // A prescribed checkpoint, independent of Body::sync's initial-state
        // policy: previous potential is exactly zero or approaches it.
        for i in 0..b.n { b.w[i] = fraction * q[i]; }
        b.r_prev = (2.0 * b.nonlinear_gradient()).sqrt();
        b.wp.clone_from(&b.w);
        b.w.clone_from(&q);
        b.sync();
        assert_eq!(b.r_prev == 0.0, fraction == 0.0);
        b.budget = b.nonlinear_gradient();
        b.linear_forces();
        b.predict(1.0 / 48_000.0);
        targets.push(b.psi);
    }
    assert!((targets[0] / targets[1] - 1.0).abs() < 1e-4,
        "zero previous potential caused a discontinuous correction: {targets:?}");
}

#[test]
fn g0_shell_loss_credit_uses_signed_total_viscous_work() {
    let mut spec = kit::splash();
    spec.grid = 0.025;
    let shell = kit::cymbal(&spec, &[]).unwrap().instrument.bodies.remove(0);
    let q: Vec<_> = (0..shell.n).map(|i| 1e-5 * (0.31 * i as f64).sin()).collect();
    let mass_norm: f64 = shell.mass.iter().zip(&q).map(|(m, q)| m * q * q).sum();
    // Independent free-edge flux: count each edge once, rather than calling
    // the production CSR damping or credit calculation in the oracle.
    let mut edge_norm = 0.0;
    for (&(i, j), &a) in &shell.grid_index {
        for (di, dj) in [(1, 0), (0, 1)] {
            if let Some(&b) = shell.grid_index.get(&(i + di, j + dj)) {
                edge_norm += 0.5 * (shell.mass[a] + shell.mass[b]) * (q[a] - q[b]).powi(2);
            }
        }
    }
    assert!(mass_norm > 0.0 && edge_norm > 0.0);
    let (dt, s1) = (1e-5, 0.02);
    for fraction in [0.5, 2.0] {
        for direction in [-1.0, 1.0] {
            let mut b = shell.clone();
            let s0 = 2.0 * fraction * s1 * edge_norm / (b.h * b.h * mass_norm);
            b.s0.fill(s0);
            b.s1 = s1;
            b.w.clone_from(&q);
            b.wp.fill(0.0);
            b.sync();
            b.linear_forces();
            b.predict(dt);
            assert!(b.sav);
            // Manufactured contact work reverses centered motion relative to
            // the backward velocity; test either sign of total viscous work.
            for i in 0..b.n { b.wn[i] = direction * q[i]; }
            let work0 = s0 * mass_norm / (2.0 * dt);
            let work1 = direction * s1 * edge_norm / (dt * b.h * b.h);
            let expected = (work0 + work1).max(0.0);
            b.finish(dt);
            assert!((b.budget - expected).abs() < 1e-11 * (work0 + work1.abs()),
                "s0 {work0} J, s1 {work1} J: credited {} J instead of net {expected} J", b.budget);
        }
    }
    let mut invalid = shell;
    invalid.wn[0] = f64::NAN;
    invalid.finish(dt);
    assert!(invalid.budget.is_nan(), "invalid viscous work became zero credit");
}

#[test]
fn g0_shell_loss_credit_refreshes_when_the_nonlinear_gradient_vanishes() {
    let mut spec = kit::splash();
    spec.grid = 0.025;
    let mut b = kit::cymbal(&spec, &[]).unwrap().instrument.bodies.remove(0);
    let dt = 1.0 / 48_000.0;
    b.s0.fill(10.0);
    b.s1 = 0.0;
    for i in 0..b.n {
        b.w[i] = 1e-5 * (0.31 * i as f64).sin();
        b.wp[i] = 2.0 * b.w[i];
    }
    b.sync();
    b.linear_forces();
    b.predict(dt);
    assert!(b.sav);
    // A contact correction reaches zero deformation after real viscous work.
    b.wn.fill(0.0);
    b.finish(dt);
    assert!(b.budget > 0.0);
    b.rotate();
    b.linear_forces();
    b.predict(dt);
    assert!(!b.sav, "undeformed shell must have zero nonlinear gradient");
    for speed in [0.0, 1e-4] {
        let mut step = b.clone();
        let mut expected = 0.0;
        for i in 0..step.n {
            let velocity = speed * (0.47 * i as f64).cos();
            step.wn[i] = step.wp[i] + 2.0 * dt * velocity;
            expected += 2.0 * dt * step.mass[i] * step.s0[i] * velocity * velocity;
        }
        step.finish(dt);
        assert!((step.budget - expected).abs() <= 1e-12 * expected.max(1e-30),
            "inactive SAV reused old credit: {} J, current work {expected} J", step.budget);
    }
}

/// The auxiliary variable must keep tracking sqrt(2V): a shallow shell carries its curvature
/// stiffness in V, so psi / sqrt(2V) scales that stiffness. Without the relaxation this ratio
/// reached 1.48 and then 0.52 within 2 s on this live-grid ride (and the 8-16 kHz bands grew back).
#[test]
fn shell_auxiliary_variable_tracks_its_potential() {
    let mut p = live::Live::new("ride", 2.5).unwrap();
    p.strike(0.6 * p.radius, 0.0, 2.5);
    let mut blk = vec![0.0f32; 4800];
    let mut worst = 0.0f64;
    for _ in 0..20 {
        p.render(&mut blk, 1.0);
        let (psi, r) = p.ins.bodies[0].sav_state();
        worst = worst.max((psi / r - 1.0).abs());
    }
    assert!(worst < 0.08, "psi / sqrt(2V) drifted {:.1} % from 1 within 2 s", 100.0 * worst);
}

#[test]
fn stencil_runs_sum_every_row_exactly_as_its_csr_row() {
    // the grid interior is computed as stencils; every row must come out bit for bit as the CSR row
    // (same products, same ascending-column order), signed zeros included
    let mut covered = 0;
    for (name, scale) in [("kick", 2.5), ("snare", 7.0), ("hihat-closed", 5.0), ("crash", 4.0), ("ride", 3.0)] {
        let p = live::Live::new(name, scale).unwrap();
        for b in &p.ins.bodies {
            for m in [&b.lap, &b.stiff] {
                covered += m.runs.iter().map(|r| r.len).sum::<usize>();
                let val = |i: usize, salt: usize| match (i * 7919 + salt) % 13 {
                    0 => 0.0,
                    1 => -0.0,
                    k => (k as f64 - 6.5) * 1.37e-3 / (1.0 + (i % 5) as f64),
                };
                let x: Vec<f64> = (0..b.n).map(|i| val(i, 3)).collect();
                let y0: Vec<f64> = (0..b.n).map(|i| val(i, 11)).collect();
                for s in [1.0, -1.0, 0.37] {
                    let (mut ya, mut yb, mut yc) = (y0.clone(), y0.clone(), y0.clone());
                    m.mul_add(&x, s, &mut ya);
                    // Keep the scalar CSR specification independent of both optimized paths.
                    for i in 0..m.n {
                        let mut acc = 0.0;
                        for p in m.ptr[i]..m.ptr[i + 1] {
                            acc += m.val[p] * x[m.col[p]];
                        }
                        yb[i] += s * acc;
                    }
                    m.rows_add(0, m.n, &x, s, &mut yc);
                    for i in 0..b.n {
                        assert_eq!(ya[i].to_bits(), yb[i].to_bits(), "{name}@{scale} {} row {i}", b.name);
                        assert_eq!(yc[i].to_bits(), yb[i].to_bits(), "{name}@{scale} {} scalar row {i}", b.name);
                    }
                }
            }
        }
    }
    assert!(covered > 1000, "the stencil runs cover only {covered} rows");
}

#[test]
fn centered_contact_preserves_elastic_impact_energy_and_dissipates_only_declared_loss() {
    // G0/G3: the old endpoint force retained only 5.89% of relative kinetic
    // energy at 48 kHz, even with mu=0. A zero-loss impact must be elastic.
    for alpha in [1.0, 1.5, 2.0, 2.5] {
        for rate in [48_000.0, 96_000.0, 192_000.0, 768_000.0] {
            for mu in [0.0, 0.05] {
                let (ma, reduced, speed) = (0.00724314, 0.0028, 2.5);
                let mb = 1.0/(1.0/reduced-1.0/ma);
                let dt = 1.0/rate;
                let kk = 1_118_132_165.064117*5e-5_f64.powf(1.5-alpha);
                let mut ins = Instrument::new(rate);
                let mut striker = Body::point("stick", ma);
                striker.w[0] = -2e-4;
                striker.wp[0] = striker.w[0]-speed*dt;
                ins.add(striker);
                ins.add(Body::point("free surface mass", mb));
                ins.contacts.push(Contact::new(vec![(0,0,1.0)],vec![(1,0,1.0)],0.0,
                    Law::Power { k:kk, alpha, mu }));
                let potential = |d: f64| kk*d.max(0.0).powf(alpha+1.0)/(alpha+1.0);
                let initial = 0.5*reduced*speed*speed;
                let (mut touched, mut separated) = (false,false);
                let mut last_energy = initial;
                for _ in 0..5000 {
                    ins.step(&mut []);
                    let (a,b) = (&ins.bodies[0],&ins.bodies[1]);
                    let (now,previous) = (a.w[0]-b.w[0],a.wp[0]-b.wp[0]);
                    let (va,vb) = ((a.w[0]-a.wp[0])/dt,(b.w[0]-b.wp[0])/dt);
                    let velocity = (now-previous)/dt;
                    let energy = 0.5*reduced*velocity*velocity+0.5*(potential(now)+potential(previous));
                    assert!((ma*va+mb*vb-ma*speed).abs()<1e-12,"contact must conserve momentum");
                    if mu == 0.0 {
                        assert!((energy/initial-1.0).abs()<1e-9,
                            "alpha {alpha}, rate {rate}: undeclared impact energy loss {energy}/{initial}");
                    } else {
                        assert!(energy<=last_energy+1e-11*initial,"declared loss added energy");
                    }
                    last_energy=energy;
                    touched |= ins.contacts[0].force>0.0;
                    // The secant still unloads while only the previous pose penetrates.
                    if touched && now<0.0 && previous<0.0 && velocity<0.0 {
                        if mu==0.0 { assert!((-velocity/speed-1.0).abs()<1e-9); }
                        else { assert!(-velocity<speed); }
                        separated=true;
                        break;
                    }
                }
                assert!(touched && separated,"the test must include a complete impact");
            }
        }
    }
}

#[test]
fn centered_contact_unloads_beyond_zero_penetration_and_the_old_bracket() {
    let law=Law::Power { k:1e9, alpha:1.5, mu:0.0 };
    for predicted in [-1e-4,1e-5] {
        let (previous,comp,dt)=(2e-4,1e-7,1e-5);
        let f=contact_force(predicted,comp,0.0,previous,law,true,dt);
        let next=predicted-comp*f;
        let expected=0.4e9*(next.max(0.0).powf(2.5)-previous.powf(2.5))/(next-previous);
        assert!(f>0.0 && next<0.0,"unloading must retain the stored contact energy");
        assert!((f/expected-1.0).abs()<1e-11);
        if predicted>0.0 { assert!(f>predicted/comp,"the old upper bracket cannot contain this root"); }
    }
    for alpha in [1.0,1.5,2.0,2.3,100.0] {
        let d=0.75_f64;
        let (value,derivative)=power_secant(d,d,alpha);
        assert!((value/d.powf(alpha)-1.0).abs()<1e-13);
        assert!((derivative/(0.5*alpha*d.powf(alpha-1.0))-1.0).abs()<1e-13);
    }
}

#[test]
fn queued_contact_solves_match_solving_one_at_a_time() {
    // the reference: every active contact in list order, each seeing the earlier ones' updates
    fn one_at_a_time(bodies: &mut [Body], contacts: &mut [Contact], k: f64) {
        for c in contacts.iter_mut().filter(|c| c.active) {
            let (mut d_pred, mut previous, mut comp) = (-c.gap, -c.gap, 0.0);
            for &(b, i, w) in &c.a {
                d_pred += w * bodies[b].wn[i];
                previous += w * bodies[b].wp[i];
                comp += w * w * bodies[b].compliance(i, k);
            }
            for &(b, i, w) in &c.b {
                d_pred -= w * bodies[b].wn[i];
                previous -= w * bodies[b].wp[i];
                comp += w * w * bodies[b].compliance(i, k);
            }
            let f = contact_force(d_pred, comp, c.force, previous, c.law, c.fixed_point_exit, k);
            apply_force(bodies, c, f, d_pred, comp, k);
        }
    }
    // a snare whose wire list also holds contacts that share nodes with others (sqrt-law copies, which
    // must close the batch, and power-law copies, which are solved alone), struck and ringing; the
    // copies leave the fixed-point exit off, so a batch mixes lanes with and without it
    let build = || {
        let mut p = live::Live::new("snare", 7.0).unwrap();
        p.strike(0.06, 0.02, 4.5);
        let n = p.ins.contacts.len();
        for (at, from, alpha) in [(3, 9, 1.5), (11, 4, 2.0), (n - 1, 20, 1.5), (7, 7, 2.5)] {
            let c = &p.ins.contacts[from];
            let Law::Power { k, mu, .. } = c.law;
            let copy = Contact::new(c.a.clone(), c.b.clone(), c.gap + 1e-7, Law::Power { k, alpha, mu });
            p.ins.contacts.insert(at, copy);
            // Inserting contacts moves the scheduled impact too. Exercise Instrument directly:
            // Live owns a private striker-slot index that this fixture must not invalidate.
            for s in &mut p.ins.strikes { if s.contact >= at { s.contact += 1; } }
        }
        p.ins
    };
    let (mut a, mut b) = (build(), build());
    let k = a.k();
    let mut out = [0.0f64; 1];
    let (mut solved, mut pressed) = (0, 0);
    for round in 0..120 {
        for _ in 0..64 { a.step(&mut out); b.step(&mut out); }
        // the same predicted positions in both: free flight plus a small push, some apart, some pressed
        for (ins, solve) in [(&mut a, true), (&mut b, false)] {
            for (bi, body) in ins.bodies.iter_mut().enumerate() {
                for i in 0..body.n {
                    let salt = ((i * 2654435761 + bi * 97 + round * 31) % 1009) as f64 / 1009.0 - 0.45;
                    // Two rounds use the ringing predictor; two impose known open/loaded
                    // poses so coverage does not depend on the live preset's internal clock.
                    body.wn[i] = if round % 4 < 2 { 0.0 } else { 2.0 * body.w[i] - body.wp[i] } + salt * 4e-7;
                    if body.name.starts_with("wire ") {
                        body.wn[i] += match round % 4 { 0 => 1e-4, 1 => -1e-4, _ => 0.0 };
                    }
                }
            }
            if round % 4 < 2 {
                // These imposed poses start at zero penetration; do not pair a changed pose
                // with an unrelated Hunt-Crossley penetration history.
                for contact in &mut ins.contacts { contact.d_now = 0.0; }
            }
            if solve {
                solve_contacts(&mut ins.bodies, &mut ins.contacts, &mut ins.cq, k);
            } else {
                one_at_a_time(&mut ins.bodies, &mut ins.contacts, k);
            }
        }
        for (ca, cb) in a.contacts.iter().zip(&b.contacts) {
            assert_eq!(ca.force.to_bits(), cb.force.to_bits(), "round {round}: force");
            assert_eq!(ca.d_now.to_bits(), cb.d_now.to_bits(), "round {round}: d_now");
            assert_eq!(ca.peak_force.to_bits(), cb.peak_force.to_bits(), "round {round}: peak");
            solved += 1;
            pressed += usize::from(ca.force > 0.0);
        }
        for (ba, bb) in a.bodies.iter().zip(&b.bodies) {
            assert!(ba.wn.iter().zip(&bb.wn).all(|(x, y)| x.to_bits() == y.to_bits()), "round {round}: {} moved differently", ba.name);
        }
    }
    assert!(pressed > solved / 4 && pressed < solved, "only {pressed} of {solved} solves pressed");
}

#[test]
fn empty_matrix_and_nan_penetration_do_nothing() {
    // an empty matrix adds nothing (a default Csr has no ptr[0] for the CSR rows to slice)
    let mut y = [0.25];
    Csr::default().mul_add(&[], 1.0, &mut []);
    Csr::default().mul_add(&[1.0], 2.0, &mut y);
    assert_eq!(y, [0.25]);
    // a NaN penetration exerts no force, solved alone or queued
    let law = Law::Power { k: 3.0e7, alpha: 1.5, mu: 0.35 };
    assert_eq!(contact_force(f64::NAN, 1e-6, 0.0, 0.0, law, true, 1e-5).to_bits(), 0.0f64.to_bits());
    let mut bodies = vec![Body::point("a", 0.01), Body::point("b", 0.01)];
    let mut contacts = vec![Contact::new(vec![(0, 0, 1.0)], vec![(1, 0, 1.0)], 0.0, law)];
    contacts[0].force = 1.0;
    bodies[0].wn[0] = f64::NAN;
    solve_contacts(&mut bodies, &mut contacts, &mut ContactQueue::default(), 1e-5);
    assert_eq!(contacts[0].force.to_bits(), 0.0f64.to_bits());
}

#[test]
fn centered_contact_loss_uses_actual_previous_shared_node_positions() {
    // G0: another contact can move a shared node after this contact stores
    // d_now. The two-step secant must read the actual previous positions.
    for alpha in [1.5, 2.0] { // queued and scalar paths
        for stale in [-1e-3, 1e-3] {
            let mut bodies = vec![Body::point("upper", 0.01), Body::point("lower", 0.02)];
            bodies[0].w[0] = 2e-5;
            bodies[0].wp[0] = 1e-5;
            bodies[0].wn[0] = 2e-5;
            let law = Law::Power { k: 3e7, alpha, mu: 0.35 };
            let mut contact = Contact::new(vec![(0, 0, 1.0)], vec![(1, 0, 1.0)], 5e-6, law);
            contact.d_now = stale;
            let previous = bodies[0].wp[0] - bodies[1].wp[0] - contact.gap;
            let dt = 1e-5;
            let mut contacts = vec![contact];
            solve_contacts(&mut bodies, &mut contacts, &mut ContactQueue::default(), dt);
            let compression = bodies[0].wn[0] - bodies[1].wn[0] - contacts[0].gap;
            // Evaluate the constitutive equation directly, independently of
            // the scalar Newton implementation used by both solver paths.
            let potential = |d: f64| 3e7*d.max(0.0).powf(alpha+1.0)/(alpha+1.0);
            let expected = (potential(compression)-potential(previous))/(compression-previous)
                * (1.0 + 0.35 * (compression - previous) / (2.0*dt)).max(0.0);
            let force = contacts[0].force;
            assert!(expected > 0.0 && force > 0.0);
            assert!((force / expected - 1.0).abs() < 1e-11,
                "alpha {alpha}, stale {stale}: force {force} does not satisfy centered law {expected}");
            assert!((0.01 * (bodies[0].wn[0] - 2e-5) + 0.02 * bodies[1].wn[0]).abs() < 1e-20);
        }
    }
}

#[test]
fn endpoint_sliding_friction_opposes_final_slip_for_anisotropic_compliance() {
    // G0: unequal masses and a skew tangent basis give non-scalar compliance.
    // Recover applied force from momentum changes, not the solver's output.
    for skew in [0.0, 0.5, -0.5] {
        let mut bodies = vec![Body::point("x", 1.0), Body::point("y", 0.25)];
        bodies[0].wn[0] = 1.0;
        bodies[1].wn[0] = 1.0;
        let mut contact = Contact::new(vec![], vec![], 0.0, Law::Power { k: 1.0, alpha: 1.0, mu: 0.0 });
        contact.force = 1.0;
        let mut frictions = vec![Friction {
            contact: 0, sx: vec![(0, 0, 1.0)], sy: vec![(0, 0, skew), (1, 0, 1.0)],
            mu: 0.5, force: 0.0, stuck: false,
        }];
        solve_frictions(&mut bodies, &[contact], &mut frictions, 1.0);
        let (x, y) = (bodies[0].wn[0], bodies[1].wn[0]);
        let fy = 0.25 * (y - 1.0);
        let fx = x - 1.0 - skew * fy;
        let (sx, sy) = (x, skew * x + y);
        assert!(!frictions[0].stuck);
        assert!(((fx * fx + fy * fy).sqrt() - 0.5).abs() < 1e-12);
        assert!((fx * sy - fy * sx).abs() < 1e-12, "force is not opposite slip: skew={skew}");
        assert!(fx * sx + fy * sy < 0.0);
        assert!(0.5 * x * x + 0.125 * y * y < 0.625, "kinetic energy increased");
    }
}

#[test]
fn endpoint_friction_recovers_prescribed_cone_solutions() {
    // G0/G3: prescribe f and lambda, then construct d = -(C + lambda I) f.
    // This checks the constitutive law across units and rotated anisotropy
    // without using another iterative solver as the oracle.
    for unit in [1e-12, 1e-6, 1e2] {
        for ratio in [1.0, 4.0, 100.0] {
            for angle in [0.0f64, 0.21, 0.77] {
                let (s, c) = angle.sin_cos();
                let cxx = unit * (c*c + ratio*s*s);
                let cyy = unit * (s*s + ratio*c*c);
                let cxy = unit * (1.0-ratio)*s*c;
                for cap in [1e-4, 1.0, 1e3] {
                    for multiplier in [0.25, 1.0] {
                        for lambda in [0.0, 1e-8*unit, unit, 1e6*unit] {
                            if multiplier < 1.0 && lambda != 0.0 { continue; }
                            let expected = [-0.6*cap*multiplier, 0.8*cap*multiplier];
                            let dx = -(cxx+lambda)*expected[0] - cxy*expected[1];
                            let dy = -cxy*expected[0] - (cyy+lambda)*expected[1];
                            let (fx, fy, stuck) = friction_force(dx, dy, cxx, cyy, cxy, cap);
                            assert!((fx-expected[0]).hypot(fy-expected[1]) <= 5e-10*cap,
                                "unit {unit}, ratio {ratio}, angle {angle}, cap {cap}, lambda {lambda}: ({fx}, {fy}) != {expected:?}");
                            assert!(fx.hypot(fy) <= cap*(1.0+2e-12));
                            if multiplier < 1.0 { assert!(stuck); }
                            let slip = [dx+cxx*fx+cxy*fy, dy+cxy*fx+cyy*fy];
                            let budget = 1e-9*cap*cap*(unit*ratio+lambda);
                            assert!((fx*slip[1]-fy*slip[0]).abs() <= budget);
                            assert!(fx*slip[0]+fy*slip[1] <= budget);
                        }
                    }
                }
            }
        }
    }
    assert_eq!(friction_force(0.0, 0.0, 1.0, 4.0, 0.5, 1.0), (0.0, 0.0, true));
    assert_eq!(friction_force(1.0, -2.0, 1.0, 4.0, 0.5, 0.0), (0.0, 0.0, false));
    // A finite prescribed sliding force remains finite when squaring the
    // unconstrained sticking force overflows its norm.
    let (cxx, cyy, cxy) = (1e-100, 4e-100, 0.3e-100);
    for lambda in [1e80, 1e100] {
        let (dx, dy) = (0.6*(cxx+lambda)-0.8*cxy, 0.6*cxy-0.8*(cyy+lambda));
        let (fx, fy, stuck) = friction_force(dx, dy, cxx, cyy, cxy, 1.0);
        assert!(!stuck && fx.is_finite() && fy.is_finite());
        assert!((fx+0.6).hypot(fy-0.8) < 5e-10);
        assert!(fx.hypot(fy) <= 1.0+2e-12);
    }
}

#[test]
fn rim_friction_sticks_inside_its_cone_and_contracts_backward_slip() {
    // the closed hat's rim frictions on a coarse grid, each solved alone on a random step
    let rim = Law::Power { k: 2.0e9, alpha: 1.5, mu: 0.08 };
    let mut ins = kit::hihat_grid_mu(0.0, 80.0, &[], rim, 0.8, 6.0, 0.4).unwrap().instrument;
    let k = ins.k();
    for b in ins.bodies.iter_mut() {
        b.coefs(k);
    }
    assert!(ins.frictions.len() >= 8, "only {} rim frictions", ins.frictions.len());
    let mut seed = 7u64;
    let mut rnd = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    let slip = |bodies: &[Body], f: &Friction| {
        let s = |t: &[(usize, usize, f64)]| t.iter().map(|&(b, i, w)| w * (bodies[b].wn[i] - bodies[b].w[i])).sum::<f64>();
        (s(&f.sx), s(&f.sy))
    };
    let (mut stuck, mut slid) = (0, 0);
    for (normal, expect_stick) in [(1.0e7, true), (1.0e-6, false)] {
        for j in 0..ins.frictions.len() {
            for b in ins.bodies.iter_mut() {
                for i in 0..b.n {
                    b.w[i] = 0.0;
                    b.wn[i] = 1.0e-6 * rnd();
                }
            }
            let c = ins.frictions[j].contact;
            ins.contacts[c].active = true;
            ins.contacts[c].force = normal;
            let (bx, by) = slip(&ins.bodies, &ins.frictions[j]);
            solve_frictions(&mut ins.bodies, &ins.contacts, &mut ins.frictions[j..j + 1], k);
            let f = &ins.frictions[j];
            let (ax, ay) = slip(&ins.bodies, f);
            let before = (bx * bx + by * by).sqrt();
            assert!(before > 0.0);
            if expect_stick {
                assert!(f.stuck, "friction {j} slid under {normal} N");
                assert!((ax * ax + ay * ay).sqrt() <= 1e-9 * before, "stuck friction {j} left slip {ax:e}, {ay:e} of {before:e}");
                stuck += 1;
            } else {
                assert!(!f.stuck, "friction {j} stuck under {normal} N");
                assert!((f.force - 0.4 * normal).abs() <= 1e-12 * normal, "sliding force {} != mu N", f.force);
                // Anisotropic compliance can rotate the remaining slip; its
                // projection on the incoming slip stays nonnegative.
                assert!(ax * bx + ay * by >= 0.0 && (ax * ax + ay * ay).sqrt() <= before, "friction {j} reversed or grew the slip");
                slid += 1;
            }
        }
    }
    assert!(stuck > 0 && slid > 0);
}

#[test]
fn g0_head_viscous_layer_increases_modal_damping() {
    // The layer also extends the retained shapes, so compare loads on the SAME
    // extended shapes: only the layer differs between them.
    let gas = standard_gas();
    let table = radiation::RadiationTable::shipped().unwrap();
    let piece = kit::drum(&kit::snare_wires_off(), &[]).unwrap();
    let b = &piece.instrument.bodies[0];
    let stiffness = |v: &[f64]| {
        let mut out = vec![0.0; v.len()];
        b.stiff.mul_add(v, 1.0, &mut out);
        b.lap.mul_add(v, b.tension, &mut out);
        out
    };
    let piston = if b.air_mass > 0.0 { b.air_mass / (b.area * b.area) } else { 0.0 };
    let areal = b.mass.iter().sum::<f64>() / b.area;
    let load = |viscous, extend| {
        radiation::ModalAirLoad::new_with_losses(&table, 0.1703, b.h, &b.xy, &b.mass, &stiffness, &gas, None, viscous, extend)
            .unwrap()
    };
    let table_only = load(None, None);
    let (plain, layered) = (load(None, Some(piston)), load(Some((areal, 140_000.0)), Some(piston)));
    assert_eq!(plain.cols.len(), layered.cols.len());
    assert!(plain.cols.len() > table_only.cols.len(), "the extension retained no shape beyond the table's");
    let n = plain.cols.len();
    let trace = |l: &radiation::ModalAirLoad| (0..n).map(|i| l.damping[i * n + i]).sum::<f64>();
    assert!(trace(&layered) > trace(&plain), "viscous layer did not add loss: {} <= {}", trace(&layered), trace(&plain));
    assert!(layered.rate_bound > plain.rate_bound);
    // Beyond the table a shape takes its evanescent air mass; on the table's own high
    // shapes (j_mn >= 15) that law agrees with the baked load (-9 % to -0.5 % on the kit
    // heads), which a wrong factor in it would not.
    let summary = plain.receipt.iter().find(|l| l.contains("against the table's own shapes with j_mn >= 15")).unwrap();
    let range: Vec<f64> = summary.rsplit(": ").next().unwrap().split(" % to ")
        .map(|v| v.trim_end_matches(" %").parse().unwrap()).collect();
    assert!(range.len() == 2 && range.iter().all(|v: &f64| v.abs() < 15.0), "{summary}");
    // The public attachment builds exactly the layered load.
    let mut attached = kit::drum(&kit::snare_wires_off(), &[]).unwrap();
    attached.instrument.bodies[0].attach_head_radiation_with_viscous(&table, 0.1703, &gas, None, 140_000.0).unwrap();
    assert_eq!(attached.instrument.bodies[0].head_radiation().unwrap().damping, layered.damping);
}

#[test]
fn g0_floor_tom_viscous_layer_refines_its_clock() {
    let gas = standard_gas();
    let point = materials::temperature_point("T", gas.temperature).unwrap()
        .with_quantity("moisture_content", fs_qty::QuantitySpec::dimensional(fs_qty::Dims::NONE), 0.12)
        .unwrap();
    let cards = materials::KitMaterials::from_card_packs(
        include_bytes!("../packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk"),
        &[],
        include_bytes!("../packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk"),
        &point,
        0.00730108512,
        0.005,
        0.05,
    ).and_then(materials::KitMaterials::with_stated_unknown_film_loss).unwrap();
    let piece = kit::drum_with_materials(&kit::tom(14.0), &[], gas, &cards).unwrap();
    assert!(piece.instrument.rate > kit::OUT_RATE, "viscous layer left the floor tom at 48 kHz");
    let (lam, who) = piece.instrument.admission_rate2();
    let ratio = lam / (4.0 * piece.instrument.rate.powi(2));
    assert!(ratio < 0.95, "floor tom was not admitted: {ratio} ({who})");
}
