use std::collections::{HashMap, BTreeMap};
use std::fs::File;
use std::io::{BufReader, Result, BufRead};
use std::collections::HashSet;
use std::path::PathBuf;


pub fn extract_vcf_positions(vcf_file_path: PathBuf) -> Result<HashSet<(String, u32)>> {
    let vcf_file = File::open(vcf_file_path)?;
    let vcf_reader = BufReader::new(vcf_file);

    let mut vcf_positions = HashSet::new();

    for line in vcf_reader.lines() {
        let line = line?;
        if !line.starts_with('#') {
            // Skip comment lines

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 2 {
                if let Ok(chrom) = fields[0].parse::<String>() {
                    if let Ok(position) = fields[1].parse::<u32>() {
                        vcf_positions.insert((chrom, position));
                    }
                }
            }
        }
    }

    Ok(vcf_positions)
}


pub fn count_bases_in_reads(sam_file_path: PathBuf, vcf_positions: &HashSet<(String, u32)>) -> Result<BTreeMap<(String, u32, char), HashMap<char, u32>>> {
    // Open the SAM file for reading
    let sam_file = File::open(sam_file_path)?;
    let sam_reader = BufReader::new(sam_file);

    // Initialize a HashMap to store the counts
    let mut position_counts: BTreeMap<(String, u32, char), HashMap<char, u32>> = BTreeMap::new();

    // Process the SAM file line by line
    for line in sam_reader.lines() {
        let line = line.unwrap();

        // Skip header lines
        if line.starts_with('@') {
            continue;
        }

        // Parse the SAM line
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue; // Malformed SAM line, skip
        }

        let flag = fields[1].parse::<u16>().unwrap();
        let chrom = fields[2].to_string();
        let position = fields[3].parse::<u32>().unwrap();
        let sequence = fields[9].as_bytes();

        let reverse_read = read_reverse_strand(flag, false);
        let direction = if reverse_read { 'r' } else { 'f' };

    
        for (base_index, base) in sequence.iter().enumerate() {
            let base_pos;
            if !reverse_read {
                base_pos = position + base_index as u32;
            }
            else {
                base_pos = position + base_index as u32 - 1;

            }
            if vcf_positions.contains(&(chrom.clone(), base_pos)) {
                let entry = position_counts
                    .entry((chrom.clone(), base_pos, direction))
                    .or_insert(HashMap::new());
                *entry.entry(*base as char).or_insert(0) += 1;
            }
        }
    }


    Ok(position_counts)
}


fn read_reverse_strand(flag:u16, paired: bool) -> bool {
    let read_reverse = 0b10000;
    let mate_reverse = 0b100000;
    let first_in_pair = 0b1000000;
    let second_in_pair = 0b10000000;
    if paired{
        if (flag & read_reverse) != 0 && (flag & first_in_pair) != 0 {
            return true
        }
        else if (flag & mate_reverse) != 0 && (flag & second_in_pair) != 0 {
            return true
        }
    }
    else {
        if (flag & read_reverse) != 0 {
            return true
        }
    }
    false
    // read.inner.core.flag == 163 || read.inner.core.flag == 83 || read.inner.core.flag == 16
    
}