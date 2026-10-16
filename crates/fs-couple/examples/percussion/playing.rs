//! Explicit playing inputs. They set mechanics, never pressure gain or pitch.
use super::Error;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stroke {
    pub speed_m_s: f64,
    // None preserves the original example's declared contact-location choice.
    pub position_m: Option<[f64;2]>,
}
impl Default for Stroke {
    fn default() -> Self { Self { speed_m_s: 0.8, position_m: None } }
}

/// Extract named physical inputs without changing existing positional syntax.
pub fn parse(args: Vec<String>) -> Result<(Vec<String>, Stroke), Error> {
    let mut positionals = Vec::new(); let mut stroke = Stroke::default();
    let mut seen_speed = false;
    let mut args = args.into_iter();
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--strike-speed-m-s" => {
                if seen_speed { return Err("strike speed may be specified only once".into()); }
                seen_speed = true;
                stroke.speed_m_s = args.next().ok_or("strike speed needs a value")?.parse()?;
                if !stroke.speed_m_s.is_finite() || !(0.0..=20.0).contains(&stroke.speed_m_s) {
                    return Err("strike speed must be finite in 0..=20 m/s; zero launches no stroke".into());
                }
            }
            "--strike-position-m" => {
                if stroke.position_m.is_some() { return Err("strike position may be specified only once".into()); }
                let x: f64 = args.next().ok_or("strike position needs x and y")?.parse()?;
                let y: f64 = args.next().ok_or("strike position needs x and y")?.parse()?;
                if !x.is_finite() || !y.is_finite() { return Err("strike position must be finite metres".into()); }
                stroke.position_m = Some([x,y]);
            }
            _ if arg.starts_with("--") => return Err(format!("unknown physical playing option: {arg}").into()),
            _ => positionals.push(arg),
        }
    }
    Ok((positionals,stroke))
}

/// Locate a vertical strike on the actual reference shell's XY chart. Unlike
/// the legacy example's nearest-facet choice, an explicit location is neither
/// snapped nor extrapolated. The existing point_port supplies the final shape.
pub fn shell_location(nodes: &[[f64;3]], triangles: &[[usize;3]], point: [f64;2])
    -> Result<(usize,[f64;3]),Error>
{
    if point.iter().any(|v| !v.is_finite()) { return Err("strike position must be finite".into()); }
    let mut hit:Option<(usize,[f64;3])>=None;
    for (index,t) in triangles.iter().enumerate() {
        let [a,b,c] = t.map(|i| nodes[i]);
        let det = (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
        if !det.is_finite() || det == 0.0 { return Err("strike needs a finite nondegenerate shell XY chart".into()); }
        let u = ((point[0]-a[0])*(c[1]-a[1])-(point[1]-a[1])*(c[0]-a[0]))/det;
        let v = ((b[0]-a[0])*(point[1]-a[1])-(b[1]-a[1])*(point[0]-a[0]))/det;
        let mut bary = [1.0-u-v,u,v];
        if bary.iter().all(|v| v.is_finite() && *v >= -32.0*f64::EPSILON) {
            for b in &mut bary { *b = b.max(0.0); }
            let sum = bary.iter().sum::<f64>(); for b in &mut bary { *b /= sum; }
            if let Some((prior,weights))=hit {
                // Adjacent facets legitimately meet at a shared edge/vertex.
                // An unrelated projected sheet must not win by file order.
                let previous=&triangles[prior];
                let on_common_feature=t.iter().enumerate().all(|(k,n)|
                    bary[k]<=64.0*f64::EPSILON || previous.contains(n))
                    && previous.iter().enumerate().all(|(k,n)|
                        weights[k]<=64.0*f64::EPSILON || t.contains(n));
                if !on_common_feature {
                    return Err("shell XY station is ambiguous across overlapping facets".into());
                }
            }else{hit=Some((index,bary));}
        }
    }
    hit.ok_or_else(||"strike position lies outside the physical shell surface (including its mounting hole)".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(s: &str) -> Vec<String> { s.split_whitespace().map(str::to_string).collect() }
    #[test]
    fn playing_options_preserve_legacy_arguments_and_are_explicit_si_values() {
        let old = args("drum-modal-mic 100 20 0.1 0.2 0.3");
        assert_eq!(parse(old.clone()).unwrap(),(old.clone(),Stroke::default()));
        let (out,stroke) = parse(args("--strike-speed-m-s 1.6 drum-modal-mic 100 --strike-position-m -0.06 0.02 20 0.1 0.2 0.3")).unwrap();
        assert_eq!(out,old); assert_eq!(stroke.speed_m_s,1.6); assert_eq!(stroke.position_m,Some([-0.06,0.02]));
        for bad in ["--strike-speed-m-s", "--strike-speed-m-s NaN", "--strike-speed-m-s -1",
            "--strike-speed-m-s 21", "--strike-speed-m-s 1 --strike-speed-m-s 2",
            "--strike-position-m 1", "--strike-position-m 1 inf", "--made-up 2"] {
            assert!(parse(args(bad)).is_err(),"{bad}");
        }
    }
    #[test]
    fn explicit_shell_position_is_interpolated_not_snapped_to_a_facet_center() {
        let nodes = [[0.0,0.0,0.0],[1.0,0.0,0.2],[0.0,1.0,0.1]];
        let (face,bary) = shell_location(&nodes,&[[0,1,2]],[0.2,0.3]).unwrap();
        assert_eq!(face,0);
        for (actual,expected) in bary.into_iter().zip([0.5,0.2,0.3]) { assert!((actual-expected).abs()<1e-14); }
        assert!(shell_location(&nodes,&[[0,1,2]],[0.8,0.8]).is_err());
    }
    #[test]
    fn shell_chart_accepts_shared_edges_but_refuses_overlapping_sheets() {
        let nodes=[[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
            [0.,0.,0.2],[1.,0.,0.2],[1.,1.,0.2]];
        let faces=[[0,1,2],[0,2,3]];
        for point in [[0.5,0.5],[0.,0.],[1.,1.]] {
            assert_eq!(shell_location(&nodes,&faces,point).unwrap().0,0);
        }
        for faces in [[[0,1,2],[4,5,6]],[[4,5,6],[0,1,2]]] {
            assert!(shell_location(&nodes,&faces,[0.6,0.2]).is_err());
            assert!(shell_location(&nodes,&faces,[0.5,0.5]).is_err());
        }
    }

}
