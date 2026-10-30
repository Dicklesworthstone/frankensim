//! Physical string and shell-bridge loads share one source-timed force program.
//! No modal participation vector or output gain is accepted by this grammar.
use super::*;

pub(super) fn parse(
    reader: &mut Reader<'_>,
    strings: usize,
    two_plane: impl Fn(usize) -> bool,
    rate: u32,
    samples: u64,
    block: usize,
) -> Result<Option<MaterialForceProgram>> {
    let last_string = strings
        .checked_sub(1)
        .ok_or("force input requires material strings")?;
    let mut row = reader.row("forces")?;
    let initialization = match row.word()? {
        "retain" => ForceInitialization::RetainState,
        "preload" => ForceInitialization::StaticPreload,
        _ => return Err("force initialization must be retain or preload".into()),
    };
    let ports = row.count(64)?;
    let event_count = row.count(MAX_EVENTS)?;
    row.end()?;
    let mut locations = Vec::with_capacity(ports);
    let mut initial_force_n = Vec::with_capacity(ports);
    for _ in 0..ports {
        let mut row = reader.row("port")?;
        let port = if row.words.clone().next() == Some("bridge") {
            row.word()?;
            // Both native shell profile versions admit at most eight original
            // bridge columns. Actual column existence remains source admission.
            MaterialForcePort::ShellBridge {
                bridge: row.count(7)?,
            }
        } else {
            let string = row.count(last_string)?;
            let station_fraction = row.scalar()?;
            let coupling = match row.word()? {
                "centreline" => StringForceCoupling::Centreline,
                "surface" => StringForceCoupling::Surface,
                "along" => {
                    let direction = TransverseDirection::try_new([row.scalar()?, row.scalar()?])?;
                    if !two_plane(string) && direction.components()[1] != 0.0 {
                        return Err("second-plane force requires polarizations two".into());
                    }
                    StringForceCoupling::CentrelineAlong(direction)
                }
                _ => return Err("force coupling must be centreline, surface or along DX DY".into()),
            };
            if !(0.0 < station_fraction && station_fraction < 1.0) {
                return Err("load station must be interior".into());
            }
            MaterialForcePort::String(StringForcePort {
                location: StringBowLocation {
                    string,
                    station_fraction,
                },
                coupling,
            })
        };
        locations.push(port);
        initial_force_n.push(row.scalar()?);
        row.end()?;
    }
    let mut events = Vec::with_capacity(event_count);
    for _ in 0..event_count {
        let mut row = reader.row("force")?;
        let sample: u64 = row.parse()?;
        let port = row.count(64)?;
        let force_n = row.scalar()?;
        row.end()?;
        if sample >= samples || port >= ports {
            return Err("force event exceeds its physical ports or observed window".into());
        }
        events.push(ModalForceEvent {
            sample,
            voice: 0,
            port,
            force_n,
        });
    }
    if ports == 0 && (event_count != 0 || initialization == ForceInitialization::StaticPreload) {
        return Err("preload and force events require an explicit load port".into());
    }
    Ok((ports != 0).then_some(MaterialForceProgram {
        ports: locations,
        initial_force_n,
        events,
        initialization,
        limits: ForceRenderConfig {
            sample_rate_hz: rate,
            max_block: block,
            max_events: MAX_EVENTS,
            max_controls: MAX_CONTROLS,
            max_projection_terms: 16_777_216,
        },
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(text: &str, two_plane: bool) -> Result<Option<MaterialForceProgram>> {
        let mut reader = Reader {
            lines: text.lines(),
            line: 0,
        };
        let value = parse(&mut reader, 2, |_| two_plane, 48_000, 16, 4)?;
        if reader.lines.next().is_some() {
            return Err("unexpected test record".into());
        }
        Ok(value)
    }

    #[test]
    fn mixed_string_and_bridge_loads_retain_physical_port_order_and_assignments() {
        let value = read("forces retain 3 3\nport 0 0.2 surface 0.002\nport bridge 1 -0.01\nport bridge 1 0.02\nforce 0 1 0.03\nforce 0 1 -0.04\nforce 9 2 0\n", false).unwrap().unwrap();
        assert!(
            matches!(value.ports[0], MaterialForcePort::String(StringForcePort {
            location: StringBowLocation { string: 0, station_fraction },
            coupling: StringForceCoupling::Surface,
        }) if station_fraction == 0.2)
        );
        assert_eq!(value.ports[1], MaterialForcePort::ShellBridge { bridge: 1 });
        assert_eq!(value.ports[2], value.ports[1]);
        assert_eq!(value.initial_force_n, [0.002, -0.01, 0.02]);
        assert_eq!(
            value
                .events
                .iter()
                .map(|e| (e.sample, e.voice, e.port, e.force_n))
                .collect::<Vec<_>>(),
            [(0, 0, 1, 0.03), (0, 0, 1, -0.04), (9, 0, 2, 0.0)]
        );
        assert_eq!(value.limits.sample_rate_hz, 48_000);
        assert_eq!(value.limits.max_block, 4);
    }

    #[test]
    fn original_string_syntax_and_direction_admission_remain_available() {
        let text = "forces preload 2 0\nport 0 0.3 centreline -0.01\nport 1 0.4 along 0 1 0.02\n";
        let value = read(text, true).unwrap().unwrap();
        assert_eq!(value.initialization, ForceInitialization::StaticPreload);
        assert!(
            matches!(value.ports[1], MaterialForcePort::String(StringForcePort {
            coupling: StringForceCoupling::CentrelineAlong(direction), ..
        }) if direction.components() == [0.0, 1.0])
        );
        assert!(read(text, false).is_err());
        assert!(read(&text.replace("along 0 1", "along 1 1"), true).is_err());
    }

    #[test]
    fn bridge_ports_require_finite_loads_valid_indices_and_complete_rows() {
        for port in [
            "bridge 8 0",
            "bridge -1 0",
            "bridge 0 NaN",
            "bridge 0 inf",
            "bridge 0",
            "bridge 0 1 extra",
            "2 0.3 surface 1",
            "0 1 surface 1",
        ] {
            assert!(
                read(&format!("forces retain 1 0\nport {port}\n"), false).is_err(),
                "{port}"
            );
        }
        // Source binding, not string-count arithmetic, owns actual bridge existence.
        let value = read("forces retain 1 0\nport bridge 7 -1\n", false)
            .unwrap()
            .unwrap();
        assert_eq!(value.ports[0], MaterialForcePort::ShellBridge { bridge: 7 });
    }

    #[test]
    fn empty_and_bounded_force_programs_keep_existing_admission() {
        assert!(read("forces retain 0 0\n", false).unwrap().is_none());
        for text in [
            "forces preload 0 0\n",
            "forces retain 65 0\n",
            "forces retain 0 16385\n",
            "forces retain 1 1\nport bridge 0 0\nforce 16 0 0\n",
            "forces retain 1 1\nport bridge 0 0\nforce 0 1 0\n",
            "forces retain 1 1\nport bridge 0 0\nforce 0 0 NaN\n",
        ] {
            assert!(read(text, false).is_err(), "{text}");
        }
    }
}
