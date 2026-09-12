use askama::Template;
use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;
use std::{env, fs};

const N_MAX: usize = 12;
const N_COEFF: usize = N_MAX * (N_MAX + 3) / 2; // 90, n = 1..=12

const WMM_FILES: [&str; 2] = ["ncei.noaa.gov/WMM2020.COF", "ncei.noaa.gov/WMM2025.COF"];
const WMM_DATA_PATH: &str = "src/wmm_data.rs";

type Result<T> = std::result::Result<T, Box<dyn Error>>;

struct WmmModel {
    epoch_year: i32,
    g: [f32; N_COEFF],
    h: [f32; N_COEFF],
    g_dot: [f32; N_COEFF],
    h_dot: [f32; N_COEFF],
}

#[derive(Template)]
#[template(path = "wmm_data.rs.j2")]
struct WmmDataTemplate {
    models: Vec<WmmModel>,
}

/// Flat index of WMM Gauss coefficient \(g_n^m\) / \(h_n^m\) (\(n \ge 1\)).
fn coeff_index(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m - 1
}

fn parse_next<'a, T: FromStr>(parts: &mut impl Iterator<Item = &'a str>, field: &str) -> Result<T>
where
    T::Err: Error + 'static,
{
    let token = parts.next().ok_or_else(|| format!("missing {field}"))?;
    token
        .parse()
        .map_err(|error| format!("unable to parse {field}: {error}").into())
}

fn parse_file(path: &Path) -> Result<WmmModel> {
    let wmm_file = fs::File::open(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let mut wmm_file_lines = BufReader::new(wmm_file).lines();

    let header_line = wmm_file_lines
        .next()
        .ok_or_else(|| format!("{} is empty", path.display()))??;
    let mut header_line_parts = header_line.split_whitespace();
    let epoch_year = parse_next::<f32>(&mut header_line_parts, "epoch year")? as i32;

    let mut g_map = HashMap::new();
    let mut h_map = HashMap::new();
    let mut g_dot_map = HashMap::new();
    let mut h_dot_map = HashMap::new();

    for line in wmm_file_lines {
        let content = line.map_err(|error| format!("{}: {error}", path.display()))?;
        if content.chars().all(|ch| ch == '9') {
            break;
        }

        let mut elements = content.split_whitespace();
        let n: usize = parse_next(&mut elements, "n index")?;
        let m: usize = parse_next(&mut elements, "m index")?;
        let g: f32 = parse_next(&mut elements, "g main field coefficient")?;
        let h: f32 = parse_next(&mut elements, "h main field coefficient")?;
        let g_dot: f32 = parse_next(&mut elements, "g secular variation coefficient")?;
        let h_dot: f32 = parse_next(&mut elements, "h secular variation coefficient")?;

        g_map.insert((n, m), g);
        h_map.insert((n, m), h);
        g_dot_map.insert((n, m), g_dot);
        h_dot_map.insert((n, m), h_dot);
    }

    let mut g = [0.0f32; N_COEFF];
    let mut h = [0.0f32; N_COEFF];
    let mut g_dot = [0.0f32; N_COEFF];
    let mut h_dot = [0.0f32; N_COEFF];

    for n in 1..=N_MAX {
        for m in 0..=n {
            g[coeff_index(n, m)] = *g_map
                .get(&(n, m))
                .ok_or_else(|| format!("no g coefficient for n={n} m={m}"))?;
            h[coeff_index(n, m)] = *h_map
                .get(&(n, m))
                .ok_or_else(|| format!("no h coefficient for n={n} m={m}"))?;
            g_dot[coeff_index(n, m)] = *g_dot_map
                .get(&(n, m))
                .ok_or_else(|| format!("no g secular variation for n={n} m={m}"))?;
            h_dot[coeff_index(n, m)] = *h_dot_map
                .get(&(n, m))
                .ok_or_else(|| format!("no h secular variation for n={n} m={m}"))?;
        }
    }

    Ok(WmmModel {
        epoch_year,
        g,
        h,
        g_dot,
        h_dot,
    })
}

fn main() -> Result<()> {
    println!("cargo:rerun-if-changed=ncei.noaa.gov");
    println!("cargo:rustc-check-cfg=cfg(has_ncei_test_file, values(any()))");
    if let Ok(entries) = fs::read_dir("ncei.noaa.gov") {
        for entry in entries.flatten() {
            let name = entry.file_name();
            let Some(name) = name.to_str() else {
                continue;
            };
            if name.ends_with("_TestValues.txt") && entry.path().is_file() {
                println!("cargo:rustc-cfg=has_ncei_test_file=\"{name}\"");
            }
        }
    }

    if env::var("GEN_WMM_SRC").is_ok_and(|value| value == "YES") {
        let models = WMM_FILES
            .into_iter()
            .map(|file| parse_file(Path::new(file)))
            .collect::<Result<Vec<_>>>()?;
        fs::write(
            Path::new(WMM_DATA_PATH),
            WmmDataTemplate { models }.render()?,
        )?;
    }
    Ok(())
}
