use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead,  Write, BufWriter};
use std::collections::HashSet;
use std::env;


fn extract_vcf_positions(vcf_file_path: &str) -> io::Result<HashSet<u32>> {
    let vcf_file = File::open(vcf_file_path)?;
    let vcf_reader = io::BufReader::new(vcf_file);

    let mut vcf_positions = HashSet::new();

    for line in vcf_reader.lines() {
        let line = line?;
        if !line.starts_with('#') {
            // Skip comment lines

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 2 {
                if let Ok(position) = fields[1].parse::<u32>() {
                    vcf_positions.insert(position);
                }
            }
        }
    }

    Ok(vcf_positions)
}


fn count_bases_in_reads(sam_file_path: &str, vcf_positions: &HashSet<u32>) -> io::Result<HashMap<(String, u32), HashMap<char, u32>>> {
    // Open the SAM file for reading
    let sam_file = File::open(sam_file_path)?;
    let sam_reader = io::BufReader::new(sam_file);

    // Initialize a HashMap to store the counts
    let mut position_counts: HashMap<(String, u32), HashMap<char, u32>> = HashMap::new();

    // Process the SAM file line by line
    for line in sam_reader.lines() {
        let line = line?;
        println!("{}", line);

        // Skip header lines
        if line.starts_with('@') {
            continue;
        }

        // Parse the SAM line
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue; // Malformed SAM line, skip
        }

        let reference_name = fields[2].to_string();
        let position = fields[3].parse::<u32>().unwrap();
        let sequence = fields[9].as_bytes();

    
        for (buchstaben_index, buchstabe) in sequence.iter().enumerate() {
            let buchstaben_position = position + buchstaben_index as u32;
            if vcf_positions.contains(&buchstaben_position) {
                let entry = position_counts
                    .entry((reference_name.clone(), buchstaben_position))
                    .or_insert(HashMap::new());
                *entry.entry(*buchstabe as char).or_insert(0) += 1;
            }
        }
    }

    Ok(position_counts)
}

fn main() -> io::Result<()> {
    // Erhalte die Befehlszeilenargumente als Iterator
    let args: Vec<String> = env::args().collect();
    println!("Test");
    // Überprüfe, ob genügend Argumente übergeben wurden
    if args.len() < 4 {
        eprintln!("Usage: {} <sam_file_path> <vcf_file_path> <output_file_path>", args[0]);
        return Err(io::Error::new(io::ErrorKind::Other, "Not enough arguments"));
    }

    // Die Argumente sind nullbasiert, args[0] ist der Name des Programms selbst
    let sam_file_path = &args[1];
    let vcf_file_path = &args[2];
    let output_file_path = &args[3];

    // Jetzt können Sie die Argumente verwenden
    println!("SAM File: {}", sam_file_path);
    println!("VCF File: {}", vcf_file_path);
    println!("Output File: {}", output_file_path);

    let vcf_positions = extract_vcf_positions(vcf_file_path)?;
    let position_counts = count_bases_in_reads(sam_file_path, &vcf_positions)?;

    // Öffnen Sie die Ausgabedatei zum Schreiben
    let output_file = File::create(output_file_path)?;
    let mut writer = BufWriter::new(output_file);

    // Schreiben Sie die Ergebnisse in die Datei
    for ((reference_name, position), counts) in &position_counts {
        writeln!(
            &mut writer,
            "Reference: {}, Position: {}",
            reference_name, position
        )?;
        for (base, count) in counts {
            writeln!(&mut writer, "Base: {}, Count: {}", base, count)?;
        }
        writeln!(&mut writer)?;
    }

    Ok(())
}
