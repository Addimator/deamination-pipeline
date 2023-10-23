
use std::fs::File;
use std::io::{BufReader, Result, BufRead};
use std::path::PathBuf;


pub fn filter_mutations(mut vcf_positions: Vec<(String, u32)>, mutations: PathBuf) -> Result<Vec<(String, u32)>> {
    let mutations_file = File::open(mutations)?;
    let vcf_reader = BufReader::new(mutations_file);
    let mut to_remove = Vec::new();

    for line in vcf_reader.lines() {
        let line = line?;
        if !line.starts_with('#') {
            // Skip comment lines
            println!("LINE: {:?}", line); 
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 2 {
                let chrom = fields[0].to_string();
                if let Ok(pos) = fields[1].parse::<u32>() {
                    let reference = fields[3].to_string();
                    let alternate = fields[4].to_string();
                    if vcf_positions.contains(&(chrom.clone(), pos)) && reference == "C" {
                        to_remove.push((chrom, pos));
                    }
                }
            }
        }
    }
    // Use the `retain` method to keep only the elements not in `to_remove`
    vcf_positions.retain(|(chrom, pos)| !to_remove.contains(&(chrom.to_owned(), *pos)));

    Ok(vcf_positions)
}


pub fn filter_inconsistent(mut vcf_positions: Vec<(String, u32)>, consistent: PathBuf) -> Result<Vec<(String, u32)>> {
    let bedgraph_file = File::open(consistent).expect("Unable to open bedGraph file");
    let bedgraph_reader = BufReader::new(bedgraph_file);
    let mut vcf_positions_filtered = Vec::new();

    for bed_line in bedgraph_reader.lines() {
        let bed_line = bed_line.expect("Error reading bedGraph line");
        if bed_line.starts_with("track") {
            continue;
        }
        let bed_fields: Vec<&str> = bed_line.split('\t').collect();
    
        // Extrahiere relevante Informationen aus dem Bedgraph
        let chrom_orig = bed_fields[0].to_string();
        let chrom = if let Some(s) = chrom_orig.strip_prefix("chr") {
            s.to_string()
        } else {
            chrom_orig
        };
        let start = bed_fields[1].parse::<usize>().expect("Invalid position value") + 1;
        let end = bed_fields[2].parse::<usize>().expect("Invalid position value");
        for pos in start..end {
            if vcf_positions.contains(&(chrom.clone(), pos as u32)) {
                vcf_positions_filtered.push((chrom.clone(), pos as u32));
            }
        }
    }

    Ok(vcf_positions_filtered)
        
}