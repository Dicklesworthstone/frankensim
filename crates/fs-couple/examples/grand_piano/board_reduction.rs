//! Cold bridge-driven Ritz projection of an already certified FE modal slice.
//! The low eigenpairs are copied exactly. Tail certificates belong to the
//! projected pencil; the original FE certificates remain in ReductionReport.

use super::{
    BoardGeometry, ReductionReport, cubic_triangle_shape, local_shape_displacement,
    nodal_displacement, motion::SourceBridgeFrame, ritz,
};
use fs_math::det;
use fs_modal::{ModePair, SliceReport};
use fs_plate::{PlateModel, ShellMesh};
use std::f64::consts::TAU;

pub(crate) struct ReducedBoard {
    pub modes: Vec<ModePair>,
    pub physical_damping: Vec<f64>,
    pub report: ReductionReport,
}

fn validate_source(free: usize, source: &SliceReport,
    apply_mass: &impl Fn(&[f64], &mut [f64])) -> Result<(), String> {
    if free == 0 || source.below_low != 0 || source.expected != source.modes.len()
        || source.below_high.checked_sub(source.below_low) != Some(source.expected)
        || source.modes.is_empty() || source.modes.len() > ritz::MAX_SOURCE_MODES {
        return Err("board reduction requires a complete certified source slice of 1..=512 modes".into());
    }
    let mut m_phi = vec![0.0; free];
    for (i, pair) in source.modes.iter().enumerate() {
        if pair.phi.len() != free || pair.phi.iter().any(|x| !x.is_finite())
            || !pair.lambda.is_finite() || pair.lambda <= 0.0
            || !pair.residual.is_finite() || pair.residual < 0.0
            || !pair.interval.0.is_finite() || pair.interval.0 <= 0.0
            || !pair.interval.1.is_finite() || pair.interval.1 < pair.interval.0
            || pair.lambda < pair.interval.0 || pair.lambda > pair.interval.1
            || (i != 0 && pair.lambda < source.modes[i - 1].lambda) {
            return Err("invalid or nonpositive certified source soundboard eigenpair".into());
        }
        apply_mass(&pair.phi, &mut m_phi);
        let norm: f64 = pair.phi.iter().zip(&m_phi).map(|(a, b)| a * b).sum();
        if !norm.is_finite() || (norm - 1.0).abs() > 1e-7 {
            return Err("source soundboard eigenvector is not mass normalized".into());
        }
        for other in &source.modes[..i] {
            let product: f64 = other.phi.iter().zip(&m_phi).map(|(a, b)| a * b).sum();
            if !product.is_finite() || product.abs() > 1e-7 {
                return Err("source soundboard modes are not mass orthogonal".into());
            }
        }
    }
    Ok(())
}

pub(super) fn project(
    geometry: &BoardGeometry,
    model: &PlateModel,
    source: &SliceReport,
    keys: &[u8],
    edge_cubic_mass: bool,
    options: &ritz::RitzOptions,
    frames: Option<&[SourceBridgeFrame]>,
) -> Result<ReducedBoard, String> {
    // Port construction indexes physical DOFs; the complete certificate and
    // mass-orthogonality checks follow in the shared modal projection.
    if source.modes.iter().any(|pair| pair.phi.len() != model.free) {
        return Err("invalid source soundboard eigenvector dimension".into());
    }
    let mesh = &geometry.chart.mesh;
    let mut ports = Vec::with_capacity(keys.len());
    for &key in keys {
        let site = geometry.bridge_sites.iter().find(|site| site.midi == key)
            .ok_or_else(|| format!("geometry has no bridge station for key {key}"))?;
        let tri = mesh.tris[site.triangle];
        let cubic = edge_cubic_mass.then(|| cubic_triangle_shape(mesh, tri, site.weights));
        ports.push(source.modes.iter().map(|pair| {
            if let Some(shape) = &cubic {
                local_shape_displacement(model, &pair.phi, tri, shape)
            } else {
                (0..3).map(|i| site.weights[i]
                    * nodal_displacement(model, &pair.phi, tri[i])).sum()
            }
        }).collect::<Vec<f64>>());
    }
    if let Some(frames)=frames {
        // Only source force projection needs this small geometry carrier;
        // never allocate a runtime MotionSurface for the 512-mode source.
        let shell=ShellMesh::new(mesh.nodes.iter().map(|&(x,y)|[x,y,0.]).collect(),
            mesh.tris.clone()).map_err(|e|e.to_string())?;
        append_secondary_ports(keys,frames,&shell,edge_cubic_mass,source,&mut ports,
            |phi,node| {
                let at=|c:usize|model.dof_map[3*node+c].map_or(0.,|i|phi[i]);
                [0.,0.,at(0),at(2),-at(1),0.]
            })?;
    }
    project_modal(model.free, geometry.damping_ratio, source, &ports, options,
        |x, y| model.m.spmv(x, y))
}

pub(crate) fn validate_frames(keys:&[u8],frames:&[SourceBridgeFrame])->Result<(),String> {
    if frames.len()!=keys.len() || frames.iter().zip(keys).any(|(frame,key)|frame.midi!=*key) {
        return Err("source bridge frames must match the complete admitted key order".into());
    }
    for frame in frames {
        if frame.primary.triangle!=frame.secondary.triangle
            || frame.primary.weights!=frame.secondary.weights
            || frame.primary.arm_m!=frame.secondary.arm_m {
            return Err("the two source bridge directions must share one physical site and arm".into());
        }
        let product:f64=frame.primary.direction.iter().zip(frame.secondary.direction)
            .map(|(a,b)|a*b).sum();
        if !product.is_finite() || product.abs()>1e-10 {
            return Err("source bridge directions must be orthogonal".into());
        }
    }
    Ok(())
}

/// Validate each supplied primary in the complete source basis, then append
/// its actual secondary force. The caller supplies its own nodal DOF map;
/// the site interpolation is shared with final runtime motion projection.
pub(crate) fn append_secondary_ports(
    keys:&[u8],frames:&[SourceBridgeFrame],mesh:&ShellMesh,edge_cubic:bool,
    source:&SliceReport,ports:&mut Vec<Vec<f64>>,
    nodal:impl Fn(&[f64],usize)->[f64;6],
)->Result<(),String> {
    validate_frames(keys,frames)?;
    if ports.len()!=keys.len() || ports.iter().any(|row|row.len()!=source.modes.len()) {
        return Err("source bridge coefficient count differs from the certified modal slice".into());
    }
    for (course,frame) in frames.iter().enumerate() {
        let primary=frame.primary.prepare(mesh,edge_cubic)?;
        let secondary=frame.secondary.prepare(mesh,edge_cubic)?;
        let mut lateral=Vec::with_capacity(source.modes.len());
        for (i,pair) in source.modes.iter().enumerate() {
            let (value,scale)=primary.project_nodal(primary.nodes().map(|node|nodal(&pair.phi,node)))?;
            let expected=ports[course][i];
            if !expected.is_finite()
                || (value-expected).abs()>1e-10*scale.max(expected.abs()).max(f64::MIN_POSITIVE) {
                return Err(format!("string polarization key {}, source mode {i}: supplied hammer-plane projection does not match the existing bridge site/basis",frame.midi));
            }
            let (value,_)=secondary.project_nodal(secondary.nodes().map(|node|nodal(&pair.phi,node)))?;
            lateral.push(value);
        }
        ports.push(lateral);
    }
    Ok(())
}

/// Project one already solved FE slice, independent of plate/shell DOF layout.
/// The supplied mass action belongs to that source pencil. Every source pair
/// is checked before rank selection, including pairs absent from all ports.
/// Callers retain the resulting full nodal vectors for every physical output.
pub(crate) fn project_modal(
    free: usize,
    damping_ratio: f64,
    source: &SliceReport,
    ports: &[Vec<f64>],
    options: &ritz::RitzOptions,
    apply_mass: impl Fn(&[f64], &mut [f64]),
) -> Result<ReducedBoard, String> {
    validate_source(free, source, &apply_mass)?;
    if !damping_ratio.is_finite() || damping_ratio < 0.0 {
        return Err("invalid source soundboard damping ratio".into());
    }
    let lambda: Vec<_> = source.modes.iter().map(|pair| pair.lambda).collect();
    let damping: Vec<_> = lambda.iter()
        .map(|value| 2.0 * damping_ratio * det::sqrt(*value)).collect();
    let basis = ritz::bridge_basis(&lambda, &damping, ports, options)?;
    let n = source.modes.len();
    let r = basis.columns.len();
    let low = options.keep_low_modes;
    let tail = r - low;
    let mut tail_k = vec![0.0; tail * tail];
    let mut identity = vec![0.0; tail * tail];
    for i in 0..tail {
        identity[i * tail + i] = 1.0;
        for j in 0..tail {
            tail_k[i * tail + j] = basis.stiffness[(low + i) * r + low + j];
        }
    }
    // Only the tail rotates. Even repeated protected low modes retain the
    // exact nodal eigenvectors, signs, eigenvalues and original certificates.
    let tail_pairs = fs_modal::eigh_gen_dense(&tail_k, &identity, tail)
        .map_err(|e| format!("reduced soundboard tail: {e}"))?;
    let mut modes = source.modes[..low].to_vec();
    let mut columns = basis.columns[..low].to_vec();
    for pair in tail_pairs {
        if !pair.lambda.is_finite() || pair.lambda <= 0.0
            || pair.phi.iter().any(|x| !x.is_finite())
            || !pair.interval.0.is_finite() || pair.interval.0 <= 0.0
            || !pair.interval.1.is_finite() || pair.interval.1 < pair.interval.0
            || !pair.residual.is_finite() || pair.residual < 0.0 {
            return Err("invalid or nonpositive projected soundboard eigenpair".into());
        }
        let mut column = vec![0.0; n];
        for (j, &weight) in pair.phi.iter().enumerate() {
            for (out, &entry) in column[low..].iter_mut().zip(&basis.columns[low + j][low..]) {
                *out += weight * entry;
            }
        }
        // Lower the final source-modal column to physical nodal DOFs once.
        // All bridge, volume, skin and motion projections use this same phi.
        let mut phi = vec![0.0; free];
        for (source_pair, &weight) in source.modes.iter().zip(&column) {
            if weight != 0.0 {
                for (out, &entry) in phi.iter_mut().zip(&source_pair.phi) {
                    *out += weight * entry;
                }
            }
        }
        modes.push(ModePair {
            lambda: pair.lambda, phi, residual: pair.residual, interval: pair.interval,
        });
        columns.push(column);
    }
    // The full congruence uses the final columns, including the tail rotation.
    // Computing from diagonal source C keeps symmetry exact and each diagonal
    // a sum of nonnegative terms. Off-diagonal physical loss is retained.
    let mut physical_damping = vec![0.0; r * r];
    for i in 0..r {
        for j in 0..=i {
            let value: f64 = (0..n)
                .map(|a| damping[a] * columns[i][a] * columns[j][a]).sum();
            if !value.is_finite() || (i == j && value < 0.0) {
                return Err("projected physical board damping overflow".into());
            }
            physical_damping[i * r + j] = value;
            physical_damping[j * r + i] = value;
        }
    }
    Ok(ReducedBoard {
        modes, physical_damping,
        report: ReductionReport {
            source_modes: n,
            protected_low_modes: low,
            source_frequency_intervals_hz: source.modes.iter().map(|pair|
                (det::sqrt(pair.interval.0) / TAU, det::sqrt(pair.interval.1) / TAU)).collect(),
            sample_hz: options.sample_hz.clone(),
            snapshot_count: basis.snapshot_count,
            max_relative_snapshot_error: basis.max_relative_snapshot_error,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn secondary_source_force_keeps_lateral_motion_and_bad_dropped_primary_refuses() {
        // G0: two uncoupled physical source coordinates at a known shell node.
        // The primary force cannot see lateral translation, but its supplied
        // orthogonal force must retain that mode in the response basis.
        let mesh=ShellMesh::new(vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],vec![[0,1,2]]).unwrap();
        let mut vertical=vec![0.;18];vertical[2]=1.;
        let mut lateral=vec![0.;18];lateral[0]=1.;
        let source=SliceReport {window:(0.,10.),below_low:0,below_high:2,expected:2,
            modes:vec![ModePair {lambda:1.,phi:vertical,residual:0.,interval:(1.,1.)},
                ModePair {lambda:9.,phi:lateral,residual:0.,interval:(9.,9.)}],
            stats:fs_modal::SliceStats {shift:0.,factorizations:0,lanczos_iters:0,
                restarts:0,factor_nnz_l:0,factor_peak_bytes:0,pivots_delayed:0}};
        let primary=super::super::motion::SourceBridgePort {triangle:0,weights:[1.,0.,0.],
            arm_m:[0.;3],direction:[0.,0.,1.]};
        let frame=SourceBridgeFrame {midi:69,primary,
            secondary:super::super::motion::SourceBridgePort {direction:[1.,0.,0.],..primary}};
        let original=vec![vec![1.,0.]];
        let mut ports=original.clone();
        let nodal=|phi:&[f64],node:usize|std::array::from_fn(|c|phi[6*node+c]);
        append_secondary_ports(&[69],&[frame],&mesh,false,&source,&mut ports,nodal).unwrap();
        assert_eq!(ports,vec![vec![1.,0.],vec![0.,1.]]);
        let options=ritz::RitzOptions {max_modes:2,keep_low_modes:0,sample_hz:vec![0.3]};
        let reduce=|ports:&[Vec<f64>]|project_modal(18,0.01,&source,ports,&options,
            |x,y|y.copy_from_slice(x)).unwrap();
        let one=reduce(&original);let both=reduce(&ports);
        assert_eq!(one.modes.len(),1);assert_eq!(both.modes.len(),2);
        assert!(one.modes.iter().all(|mode|mode.phi[0]==0.));
        assert!(both.modes.iter().any(|mode|mode.phi[0].abs()>0.9));
        assert!(both.report.max_relative_snapshot_error<1e-12);

        // Now the second source mode is rotation about y. A false x arm
        // leaves the retained vertical mode unchanged but changes this tail's
        // primary displacement; validating only the reduced basis would hide it.
        let mut rotated=source.clone();rotated.modes[1].phi.fill(0.);
        rotated.modes[1].phi[4]=1.;
        let mut wrong=frame;wrong.primary.arm_m[0]=0.02;wrong.secondary.arm_m[0]=0.02;
        let error=append_secondary_ports(&[69],&[wrong],&mesh,false,&rotated,
            &mut original.clone(),nodal).unwrap_err();
        assert!(error.contains("source mode 1"),"{error}");
        let error=append_secondary_ports(&[70],&[frame],&mesh,false,&source,
            &mut original.clone(),nodal).unwrap_err();
        assert!(error.contains("key order"),"{error}");
    }

    #[test]
    fn complete_source_validation_cannot_hide_bad_unselected_modes_above_runtime_rank() {
        // G0: the protected one-mode result still validates every one of the
        // 129 supplied source pairs, even with zero forcing on their tail.
        let n=129;
        let source=SliceReport {
            window:(0.,130.),below_low:0,below_high:n,expected:n,
            modes:(0..n).map(|i| {
                let mut phi=vec![0.;n];phi[i]=1.;let lambda=(i+1) as f64;
                ModePair {lambda,phi,residual:0.,interval:(lambda,lambda)}
            }).collect(),
            stats:fs_modal::SliceStats {shift:0.,factorizations:0,lanczos_iters:0,
                restarts:0,factor_nnz_l:0,factor_peak_bytes:0,pivots_delayed:0},
        };
        let options=ritz::RitzOptions {max_modes:1,keep_low_modes:1,sample_hz:vec![1.]};
        let mut port=vec![0.;n];port[0]=1.;let ports=vec![port];
        let run=|source:&SliceReport|project_modal(n,0.01,source,&ports,&options,
            |x,y|y.copy_from_slice(x));
        let reduced=run(&source).unwrap();
        assert_eq!(reduced.report.source_modes,n);assert_eq!(reduced.modes.len(),1);
        assert_eq!(reduced.modes[0].phi,source.modes[0].phi);
        assert_eq!(reduced.report.source_frequency_intervals_hz.len(),n);
        let mut wrong=source.clone();wrong.modes[n-1].phi[n-1]=2.;
        assert!(run(&wrong).is_err());
        wrong=source.clone();wrong.modes[n-1].phi=source.modes[0].phi.clone();
        assert!(run(&wrong).is_err());
        wrong=source.clone();wrong.modes[n-1].interval.1=128.;
        assert!(run(&wrong).is_err());
        wrong=source.clone();wrong.modes[n-1].residual=f64::NAN;
        assert!(run(&wrong).is_err());
        wrong=source.clone();wrong.expected-=1;
        assert!(run(&wrong).is_err());
        wrong=source.clone();wrong.below_low=1;
        assert!(run(&wrong).is_err());
    }
}
