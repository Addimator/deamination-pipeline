
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


pub fn filter_inconsistent() {

}