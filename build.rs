use rinja::Template;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::{env, fs};

const WMM_FILES: [&str; 2] = ["ncei.noaa.gov/WMM2020.COF", "ncei.noaa.gov/WMM2025.COF"];
const WMM_RENDER_FILE: &str = "src/wmm.rs";

struct WmmModel {
    model_version: i32,
    g_mfc: [f32; 90],
    h_mfc: [f32; 90],
    g_svc: [f32; 90],
    h_svc: [f32; 90],
}

#[derive(Template)]
#[template(path = "wmm.rs.j2")]
struct WMMTemplate {
    models: Vec<WmmModel>,
}

fn ix(n: usize, m: usize) -> usize {
    (1..=n).sum::<usize>() + m - 1
}

fn parse_file(path: &Path) -> WmmModel {
    let wmm_file = fs::File::open(path).expect("WMM file not found!");
    let mut wmm_file_lines = BufReader::new(wmm_file).lines();

    let header_line = wmm_file_lines
        .next()
        .expect("WMM file empty!")
        .expect("Read error!");
    let mut header_line_parts = header_line.split_whitespace();
    let model_version = header_line_parts
        .next()
        .expect("No year defined in WMM file!")
        .parse::<f32>()
        .expect("Invalid year format!") as i32;

    let mut g_mfc_map = HashMap::new();
    let mut h_mfc_map = HashMap::new();
    let mut g_svc_map = HashMap::new();
    let mut h_svc_map = HashMap::new();

    for line in wmm_file_lines {
        match line {
            Ok(content) => {
                if content.chars().all(|ch| ch == '9') {
                    break;
                }

                let mut elements = content.split_whitespace();
                let n: usize = elements
                    .next()
                    .expect("Lack of n index!")
                    .parse()
                    .expect("Unable to parse n index!");
                let m: usize = elements
                    .next()
                    .expect("Lack of m index!")
                    .parse()
                    .expect("Unable to parse m index!");
                let g_mfc: f32 = elements
                    .next()
                    .expect("Lack of g main field coefficient!")
                    .parse()
                    .expect("Unable to parse g main field coefficient!");
                let h_mfc: f32 = elements
                    .next()
                    .expect("Lack of h main field coefficient!")
                    .parse()
                    .expect("Unable to parse h main field coefficient!");
                let g_svc: f32 = elements
                    .next()
                    .expect("Lack of g secular variation coefficient!")
                    .parse()
                    .expect("Unable to parse g secular variation coefficient!");
                let h_svc: f32 = elements
                    .next()
                    .expect("Lack of h secular variation coefficient!")
                    .parse()
                    .expect("Unable to parse h secular variation coefficient!");

                g_mfc_map.insert((n, m), g_mfc);
                h_mfc_map.insert((n, m), h_mfc);
                g_svc_map.insert((n, m), g_svc);
                h_svc_map.insert((n, m), h_svc);
            }

            Err(error) => {
                panic!("Error reading lines: {}!", error)
            }
        }
    }

    let mut g_mfc = [0.0f32; 90];
    let mut h_mfc = [0.0f32; 90];
    let mut g_svc = [0.0f32; 90];
    let mut h_svc = [0.0f32; 90];

    for n in 1..=12 {
        for m in 0..=n {
            g_mfc[ix(n, m)] = *g_mfc_map
                .get(&(n, m))
                .expect("No value for g main field coefficient!");
            h_mfc[ix(n, m)] = *h_mfc_map
                .get(&(n, m))
                .expect("No value for h main field coefficient!");
            g_svc[ix(n, m)] = *g_svc_map
                .get(&(n, m))
                .expect("No value for g secular variation coefficient!");
            h_svc[ix(n, m)] = *h_svc_map
                .get(&(n, m))
                .expect("No value for h secular variation coefficient!");
        }
    }

    WmmModel {
        model_version,
        g_mfc,
        h_mfc,
        g_svc,
        h_svc,
    }
}

fn main() {
    if env::var("GEN_WMM_SRC").is_ok_and(|value| value == "YES") {
        let wmm = WMMTemplate {
            models: WMM_FILES.map(|file| parse_file(Path::new(file))).into(),
        };
        fs::write(Path::new(&WMM_RENDER_FILE), wmm.render().unwrap()).unwrap();
    }
}
